##!/usr/bin/env python
import sys
import os
import re
import sqlite3
import glob
import time
from datetime import datetime
import gzip
import pandas as pd
import numpy as np
from itertools import chain
import math
from multiprocessing import Pool


# Information about commmand line arguments required by python script.
def help_info():
    print(
        "\n\n--vepfile=filename or -f=filename \t\t Input file containing results",
        "of VEP from Ensembl analysis. Mandatory.\n\n",
    )
    print("--output=directory or -o=directory \t\t Execution output dir. Default ./output.\n\n")
    print("--root=rootname or -r=rootname \t\t Root name for the output folder. Default: no root name")
    print(
        "--int=jobID or -i=jobID \t\t\t Job ID code (when executing from sequencingAP)."
        " Default: Generated during execution.\n\n"
    )
    print("--databases=directory or -d=directory \t\t\t Absolute path to databases directory. Mandatory.\n\n")
    print("--forkp \t\t\t\t\t\t Number of forks to improve runtime. Default: 4.\n\n")
    print(
        "--filtering=True/False \t\t Boolean value indicating whether variants must be"
        " filtered based ontheir impact. In case this parameter is set to True, the"
        " 'vep_data_sorted.csv' output file willonly contain variants with HIGH or"
        " MODERATE impact. Default: True. \n\n"
    )
    print(
        "\ni.e. VEP_parser.pl -f=file.vcf -o=/home/user/z13-222 -r=analysis"
        " -i=20140213_000000-d=/home/epineiro/Programs/PCDA/databases --forkp=8"
        " --filtering=False\n\n"
    )
    sys.exit()


# Generate a jobid if not provided by the user.
def get_runid():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return root_name + timestamp + "_VEP"


# Print all the steps being performed by the python script in the log file.
def printl(step):
    try:
        logfile = open(f"{outexp}{jobid}.log", "a+")
    except IOError:
        print(f"Could not open file {logfile}")
    logfile.write(step + "\n")
    logfile.close()
    print(step)
    return


# Extract NCBI Gene ID, KEGG pathway description and ID from genes_ids,
# kegg_gene_pathway_db and pathw_desc_db databases, respectively.
def get_kegg_path_desc(kegg_paths):
    results = ""
    for path_id in kegg_paths:
        if path_id in pathw_desc_db.index:
            results += f'{pathw_desc_db.loc[path_id, "KEGG Pathway desc"]}|'

    keggpaths = "|".join(kegg_paths)
    results = results.rstrip("|")
    return results, keggpaths


def get_kegg_path(kegg_id):
    if kegg_id in kegg_gene_pathway_db.index:
        kegg_path_desc = f'{kegg_gene_pathway_db.loc[kegg_id, "KEGG Pathway ID"]}'.split("|")
        return get_kegg_path_desc(kegg_path_desc)
    else:
        return ("", "")


def get_kegg_id_sym(hgnc_symbol):
    if hgnc_symbol in genes_ids_db.index:
        kegg_pathways = get_kegg_path(int(genes_ids_db.loc[hgnc_symbol, "NCBI Gene ID"]))
        return (hgnc_symbol, kegg_pathways[0], kegg_pathways[1])
    else:
        return ("", "", "")


# VSCORE creation for each variant.
def create_vscore(q_data_line, annotations):
    linedata = "".join(annotations).rstrip("\n").split("\t")[:50]
    score = 0
    scored_columns = {}
    # fmt: off
    scored_columns_list = ["Chr", "Loc", "mut", "poly_score", "sift_score", "CADD_phred", 
                           "gene_role", "cosmic_id", "protein_position", "mut_cosmic_freq", 
                           "gene_cosmic_freq", "consequence", "impact", "GMAF", "gnomAD", 
                           "pfam", "interpro", "clinvar_clinical_significance", "INFO","gene_hgnc"]
    # fmt:on
    scored_columns = {value: index for index, value in enumerate(linedata) if value in scored_columns_list}

    # Decide the branch for the calculation according to the variation or gene definition.
    components = q_data_line[scored_columns["gene_role"]].split("; ")
    genetype = list(chain(*[re.findall(r"ONC|TSG", gr) for gr in components if re.findall(r"ONC|TSG", gr)]))
    genetype_string = ":".join(list(set(genetype)))
    branch = "UNCLASSIFIED"

    if genetype_string == "ONC":
        branch = "ONC"
    elif genetype_string == "TSG":
        branch = "TSG"

    if branch == "UNCLASSIFIED" and re.search(
        "stop_gain|stop_lost|frameshift_variant|splice_donor_variant|splice_acceptor_variant|splice_region_variant",
        q_data_line[scored_columns["consequence"]],
    ):
        branch = "TSG"

    # Add scores
    prediction_damaging = 0

    # Cosmic ID
    cosmic_id = q_data_line[scored_columns["cosmic_id"]]
    mut_cosmic_freq = q_data_line[scored_columns["mut_cosmic_freq"]]
    gene_cosmic_freq = q_data_line[scored_columns["gene_cosmic_freq"]]

    if cosmic_id.startswith("COSV"):
        if "PATHOGENIC" in cosmic_id:
            prediction_damaging += 1

        # Match mutation frequency
        match_mut_freq = re.match(r"(\d+) \/ (\d+)", str(mut_cosmic_freq))
        if match_mut_freq and (branch == "ONC" or branch == "UNCLASSIFIED"):
            num, denom = int(match_mut_freq.group(1)), int(match_mut_freq.group(2))
            ss = 0.125 / 2 if num >= 100 else (0.125 / 2) * ((math.log(num)) / (math.log(denom)))
            score += ss

        # Match gene frequency
        match_gene_freq = re.match(r"(\d+) \/ (\d+)", str(gene_cosmic_freq))
        if match_gene_freq:
            num, denom = int(match_gene_freq.group(1)), int(match_gene_freq.group(2))
            ss = 0.125 / 2 if num >= 100 else (0.125 / 2) * ((math.log(num)) / (math.log(denom)))
            score += ss

    # Prediction score
    if q_data_line[scored_columns["poly_score"]] and float(q_data_line[scored_columns["poly_score"]]) > 0.435:
        prediction_damaging += 1

    if q_data_line[scored_columns["sift_score"]] and float(q_data_line[scored_columns["sift_score"]]) <= 0.05:
        prediction_damaging += 1

    if q_data_line[scored_columns["CADD_phred"]] and float(q_data_line[scored_columns["CADD_phred"]]) > 20:
        prediction_damaging += 1

    if prediction_damaging >= 3:
        score += 0.125
    elif prediction_damaging == 2:
        score += 0.08
    elif prediction_damaging == 1:
        score += 0.04

    # Mutation type
    if q_data_line[scored_columns["impact"]] == "HIGH":
        score += 0.125

    # Frequencies
    gmaf = q_data_line[scored_columns["GMAF"]]
    gnomAD = q_data_line[scored_columns["gnomAD"]]
    if (gmaf and float(gmaf) < 1) or (gmaf == ""):
        score += 0.125 / 2
    if (gnomAD and float(gnomAD) < 1) or gnomAD == "":
        score += 0.125 / 2
    # fmt:off
    # Domains
    pfam = q_data_line[scored_columns["pfam"]]
    interpro = q_data_line[scored_columns["interpro"]]
    match_pfam = re.match(r"^(\w+)\.", pfam)
    if (match_pfam and match_pfam.group(1) in cancer_domain_db.index or 
        interpro == "Mutation previous last protein domain"):
        score += 0.125
    elif pfam or interpro:
        score += 0.125 / 2

    # Clinvar
    if "Pathogenic" in q_data_line[scored_columns["clinvar_clinical_significance"]]:
        if zygosity["_".join(q_data_line[scored_columns[key]] for key in ["Chr", "Loc", "mut"])]:
            score += 0.125
        else:
            if branch in ["ONC", "UNCLASSIFIED"]:
                score += 0.250
            else:
                score += 0.3125

    # Homozigous
    if (zygosity["_".join(q_data_line[scored_columns[key]] for key in ["Chr", "Loc", "mut"])] == "Homozygous"):
        if branch in ["ONC", "UNCLASSIFIED"]:
            score += 0.125
        else:
            score += 0.1875

    # GScore
    if q_data_line[scored_columns["gene_hgnc"]] in gscore_db.index:
        score += 0.125 * gscore_db.loc[q_data_line[scored_columns["gene_hgnc"]], "gscore"]

    q_data_line.extend(["{:.4f}".format(score), branch])

    return q_data_line


## VEP_parser_per_chr(chr_num).
def VEP_parser_per_chr(chr_num):
    printl(f"Processing chromosome {chr_num}...")

    input_chr = open(f"{outexp}/ensembl_vep_{chr_num}.csv")
    out_chr = open(f"{outexp}/vep_data_{chr_num}.csv", "w")
    outsort_chr = open(f"{outexp}/vep_data_sorted_{chr_num}.csv", "w")

    out_chr.write("".join(annotations))
    outsort_chr.write("".join(annotations_sort))

    count = 1
    last_gene = ["", "", ""]
    for line in input_chr:
        if not re.match(r"^CHROM", line):
            line = line.rstrip("\n")
            line = line.split("\t")

            line[0] = re.sub(r"chr", "", line[0])
            line[2] = "" if line[2] == "." else line[2]

            VCF_pos = line[1]
            VCF_ref = line[3]
            VCF_alt = line[4]

            if (len(line[3]) != len(line[4])) and ("," not in line[4]):
                if len(line[3]) > len(line[4]):
                    # There might be an internal deletion that would not be
                    # properly processed
                    if line[3].startswith(line[4]):
                        line[3] = re.sub(rf"^{line[4]}", "", line[3])
                        line[1] += len(line[4])
                        line[4] = "-"
                else:
                    # There might be an internal insertion that would not be
                    # properly processed (ej. TGCTCTACC/TATAGATCGGAAGCTCTACC)
                    if line[4].startswith(line[3]):
                        line[4] = re.sub(rf"^{line[3]}", "", line[4])
                        line[1] += len(line[3])
                        line[3] = "-"

            num_fields = 12
            add_fields = num_fields - len(line)

            # Extend the list with empty elements.
            for field in range(add_fields):
                line.append("")

            line[10] = line[3] + "/" + line[4]

            line[11] = line[0] + ":" + line[1]

            if line[9] != "":
                gentype = line[9].split(":")
                if gentype[0] == "1/1":
                    zygosity[f"{line[0]}_{line[1]}_{line[10]}"] = "Homozygous"
                else:
                    zygosity[f"{line[0]}_{line[1]}_{line[10]}"] = "Heterozygous"
            else:
                zygosity[f"{line[0]}_{line[1]}_{line[10]}"] = ""

            q_data = []

            if re.findall(r"CSQ=(.+)", line[7]):
                match = "".join(re.findall(r"CSQ=(.+)", line[7]))
                # Extract variant annotation per each transcript.
                transcripts = match.split(",")

                for trans in transcripts:
                    fields = trans.split("|")
                    # fmt:on
                    fields[pos["Consequence"]] = fields[pos["Consequence"]].replace("&", ",")
                    if fields[pos["Existing_variation"]]:
                        fields[pos["Existing_variation"]] = fields[pos["Existing_variation"]].replace("&", ",")

                    HGVSp = ""
                    if "HGVSp" in pos.keys() and fields[pos["HGVSp"]] and re.findall(r":(.+)", fields[pos["HGVSp"]]):
                        HGVSp = "".join(re.findall(r":(.+)", fields[pos["HGVSp"]]))

                    HGVSc = ""
                    if "HGVSc" in pos.keys() and fields[pos["HGVSc"]] and re.findall(r":(.+)", fields[pos["HGVSc"]]):
                        HGVSc = "".join(re.findall(r":(.+)", fields[pos["HGVSc"]]))

                    pol_cons, pol_sco = "", ""
                    if (
                        "PolyPhen" in pos.keys()
                        and fields[pos["PolyPhen"]]
                        and re.findall(r"(\w+)\((\d+\.*\d*)\)", fields[pos["PolyPhen"]])
                    ):
                        PolyPhen_match = re.findall(r"(\w+)\((\d+\.*\d*)\)", fields[pos["PolyPhen"]])
                        pol_cons = PolyPhen_match[0][0]
                        pol_sco = PolyPhen_match[0][1]

                    # elif re.match('stop_gained', fields[pos["Consequence"]]) or
                    #     re.match('frameshift_variant', fields[pos["Consequence"]]):
                    #     pol_cons="inferred"
                    #     pol_sco=1

                    sift_cons, sift_sco = "", ""
                    if (
                        "SIFT" in pos.keys()
                        and fields[pos["SIFT"]]
                        and re.findall(r"(\w+)\((\d+\.*\d*)\)", fields[pos["SIFT"]])
                    ):
                        SIFT_match = re.findall(r"(\w+)\((\d+\.*\d*)\)", fields[pos["SIFT"]])
                        sift_cons = SIFT_match[0][0]
                        sift_sco = SIFT_match[0][1]

                    # elif re.match('stop_gained', fields[pos["Consequence"]]) or
                    #     re.match('frameshift_variant', fields[pos["Consequence"]]):
                    #     sift_cons="inferred"
                    #     sift_sco=0

                    CADD_phred, CADD_raw = "", ""
                    if "CADD_PHRED" in pos.keys():
                        CADD_phred = fields[pos["CADD_PHRED"]]
                    if "CADD_RAW" in pos.keys():
                        CADD_raw = fields[pos["CADD_RAW"]]

                    cosmic_id, cosmic_fathmm, gene_freq, mut_freq, cosmic_total = ("", "", 0, 0, "")
                    if "HGVSc" in pos.keys() and fields[pos["HGVSc"]]:
                        Cosmic_key = f'{fields[pos["SYMBOL"]]}:{fields[pos["Feature"]]}:{HGVSc}'
                        cosmic_query = cosmic_db.execute(
                            "SELECT cosmic_id, FATHMM, Gene_freq, Mut_freq, Total FROM cosmic_table WHERE"
                            f' Cosmic_key="{Cosmic_key}";'
                        ).fetchall()
                        if cosmic_query:
                            (cosmic_id, cosmic_fathmm, gene_freq, mut_freq, cosmic_total) = cosmic_query[0]
                            mut_freq = f"{mut_freq} / {cosmic_total}"
                        if cosmic_fathmm != "":
                            cosmic_id += f":{cosmic_fathmm}"

                    if fields[pos["SYMBOL"]] in cosmic_gene_freq_db.index:
                        gene_freq = (
                            f'{cosmic_gene_freq_db.loc[fields[pos["SYMBOL"]]].values[0]} /'
                            f' {cosmic_gene_freq_db.loc[fields[pos["SYMBOL"]]].values[1]}'
                        )

                    kegg_data, kegg_ids = "", ""
                    if fields[pos["SYMBOL"]]:
                        if fields[pos["SYMBOL"]] != last_gene[0]:
                            # Get information about gene symbol: pathway_description, pathway_ids and entrez_gene_id
                            last_gene = get_kegg_id_sym((fields[pos["SYMBOL"]]).upper())
                        kegg_data = last_gene[1]
                        kegg_ids = last_gene[2]

                    var_type = fields[pos["VARIANT_CLASS"]]

                    GMAF = ""
                    if "GMAF" in pos.keys() and fields[pos["GMAF"]]:
                        GMAF_a = fields[pos["GMAF"]].split("&")
                        for gf in GMAF_a:
                            if gf != "":
                                GMAF = float(gf) * 100

                    clinvar_acc, clinvar_dis, clinvar_pat = "", "", ""
                    query = f"{line[0]}:{VCF_pos}:{VCF_ref}:{VCF_alt}"
                    if query in clinvar_db.index:
                        clinvar_acc = "; ".join(clinvar_db.loc[query, "Acc"])
                        clinvar_dis = "; ".join(clinvar_db.loc[query, "Trait"])
                        clinvar_pat = "; ".join(clinvar_db.loc[query, "Significance"])

                    prot_pos = fields[pos["Protein_position"]]

                    uniprot, pfam, interpro = "", "", ""

                    if fields[pos["SYMBOL"]] in uniprot_b_db.index:
                        ident = uniprot_b_db.loc[fields[pos["SYMBOL"]], "PROTEIN_ID"].split(";")[0]

                        prot_end = 0

                        if re.match(r"(\d+)\-(\d+)", str(prot_pos)):
                            match = re.match(r"(\d+)\-(\d+)", str(prot_pos))
                            prot_pos = int(match.group(1))
                            prot_end = int(match.group(2))

                        if re.match(r"(\d+)\-(\?)", str(prot_pos)):
                            match = re.match(r"(\d+)\-(\?)", str(prot_pos))
                            prot_pos = int(match.group(1))
                            prot_end = 0

                        if re.match(r"(\?)\-(\d+)", str(prot_pos)):
                            match = re.match(r"(\?)\-(\d+)", str(prot_pos))
                            prot_pos = 0
                            prot_end = int(match.group(2))
                        start = time.time()
                        if prot_pos != "":
                            prot_pos = int(prot_pos)
                            # Extract Pfam accesion number and domain description.
                            if ident in pfam_a_db.index:
                                pfam_domains = pfam_a_db.loc[[ident]]
                                for index, domain in pfam_domains.iterrows():
                                    if (prot_pos >= domain["START"] and prot_pos <= domain["END"]) or (
                                        prot_end >= domain["START"] and prot_end <= domain["END"]
                                    ):
                                        pfam = f'{domain["PFAM_ACC"]}: {domain["DOMAIN_DESCRIPT"]}'

                            # Extract Interpro domain ID and domain description.
                            if ident in interpro_db.index:
                                interpro_domains = interpro_db.loc[[ident]]
                                for index, domain in interpro_domains.iterrows():
                                    if (prot_pos >= domain["START"] and prot_pos <= domain["END"]) or (
                                        prot_end >= domain["START"] and prot_end <= domain["END"]
                                    ):
                                        interpro = f'{domain["DOMAIN_ID"]}: {domain["DOMAIN_DESCRIPT"]}'
                            # fmt:off
                            if (re.match(r"(stop_gained|frameshift_variant)",fields[pos["Consequence"]]) and interpro == ""
                                and ident in last_domain_db.index and prot_pos <= last_domain_db.loc[ident, "START"]):
                                interpro = "Mutation previous last protein domain"

                    gene_role = ""
                    if fields[pos["SYMBOL"]] in generole_db.index:
                        generole_dict = generole_db.loc[fields[pos["SYMBOL"]]].dropna(axis=0).to_dict()
                        if generole_dict:
                            for key, generole in generole_dict.items():
                                gene_role += f"{key}:{generole}; "
                            gene_role = gene_role.rstrip("; ")

                    gnomAD = ""
                    if "gnomADe_AF" in pos.keys() and fields[pos["gnomADe_AF"]]:
                        gnomAD_list = fields[pos["gnomADe_AF"]].split("&")
                        for ex in gnomAD_list:
                            if ex:
                                gnomAD = float(ex) * 100

                    gnomAD_NFE = ""
                    if "gnomADe_NFE" in pos.keys() and fields[pos["gnomADe_NFE"]]:
                        gnomAD_list = fields[pos["gnomADe_NFE"]].split("&")
                        for ex in gnomAD_list:
                            if ex:
                                gnomAD_NFE = float(ex) * 100

                    # fmt: off
                    q_data_line = [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8],
                                 line[9], line[10], line[11], fields[pos["Allele"]], fields[pos["Gene"]],
                                 fields[pos["Feature"]], fields[pos["Feature_type"]], fields[pos["Consequence"]],
                                 fields[pos["IMPACT"]], fields[pos["cDNA_position"]], fields[pos["CDS_position"]],
                                 fields[pos["Protein_position"]], fields[pos["Amino_acids"]], fields[pos["Codons"]],
                                 fields[pos["Existing_variation"]], "", fields[pos["APPRIS"]], pol_cons, pol_sco,
                                 sift_cons, sift_sco, CADD_phred, CADD_raw, fields[pos["SYMBOL"]], gene_role, cosmic_id,
                                 kegg_data, kegg_ids, clinvar_acc, clinvar_dis, clinvar_pat, var_type, HGVSc, HGVSp,
                                 GMAF, gnomAD, gnomAD_NFE, pfam, interpro, gene_freq, mut_freq]

                    q_data_line = create_vscore(q_data_line, annotations)
                    q_data_line = [str(annot) for annot in q_data_line]
                    q_data.append(q_data_line)
                    out_chr.write("\t".join(q_data_line) + "\n")

                count += 1

            # fmt:off
            q_data.sort(key=lambda x: (x[0], x[1], x[10], x[13], x[42], x[14]))
            for q_data_line in q_data:
                if filtering:
                    if re.match(r"HIGH|MODERATE", q_data_line[17]):
                        outsort_chr.write("\t".join(q_data_line[0:2] + [q_data_line[10]] + q_data_line[13:18] 
                                                    + q_data_line[25:37] + q_data_line[20:22] + q_data_line[37:]) + "\n")
                else:
                    outsort_chr.write("\t".join(q_data_line[0:2] + [q_data_line[10]] + q_data_line[13:18] 
                                                + q_data_line[25:37] + q_data_line[20:22] + q_data_line[37:])+ "\n")

    input_chr.close()
    out_chr.close()
    outsort_chr.close()
    printl(f"Chromosome {chr_num} processed!")
    return


# Assignation of default values for variables
outdir = "./output/"
outpath = os.getcwd()
root_name = ""

max_processes = 8
filtering = True
zygosity = {}

vep_results_file = None
dbdir = None
jobid = None

# Command line arguments handle
if len(sys.argv) == 1 or any(re.match(r"^((\-\-help)|(\-h))$", arg) for arg in sys.argv[1:]):
    help_info()

for arg in sys.argv[1:]:
    # Input file
    if re.match(r"^((\-\-vepfile=)|(\-f=))", arg):
        match1 = re.match(r"\-(\-vepfile|f)=(.*)", arg)
        if match1:
            vep_results_file = (
                match1.group(2)
                if match1.group(2)
                else sys.exit("\nEmpty argument. Please enter the parameter information.\n\neg. -f=file.vcf\n\n")
            )
            vep_results_file_tmp = glob.glob(vep_results_file)
            if len(vep_results_file_tmp) != 0:
                vep_results_file = vep_results_file_tmp[0]
            else:
                vep_results_file = None

    # Output file
    elif re.match(r"^((\-\-output=)|(\-o=))", arg):
        match2 = re.match(r"\-(\-output|o)=(.*)", arg)
        if match2:
            outdir = (
                match2.group(2)
                if match2.group(2)
                else sys.exit(
                    "\nEmpty argument. Please enter the parameter information.\n\neg. -o=/home/user/z13-222\n\n"
                )
            )
            outpath = ""

    # Root for the file
    elif re.match(r"^((\-\-root=)|(\-r=))", arg):
        match3 = re.match(r"\-(\-root|r)=(.*)", arg)
        if match3:
            root_name = (
                match3.group(2) + "_"
                if match3.group(2)
                else sys.exit("\nEmpty argument. Please enter the parameter information.\n\neg. -r=analysis\n\n")
            )

    # JobID
    elif re.match(r"^((\-\-int=)|(\-i=))", arg):
        match4 = re.match(r"\-(\-int|i)=(.*)", arg)
        if match4:
            jobid = (
                match4.group(2)
                if match4.group(2)
                else sys.exit("\nEmpty argument. Please enter the parameter information.\n\neg. -i=20140213_000000\n\n")
            )

    # databases path
    elif re.match(r"^((\-\-databases=)|(\-d=))", arg):
        match5 = re.match(r"\-(\-databases|d)=(.*)", arg)
        if match5:
            dbdir = (
                match5.group(2)
                if match5.group(2)
                else sys.exit(
                    "\nEmpty argument. Please enter the parameter information.\n\neg."
                    " -d=/home/epineiro/Programs/PCDA/databases\n\n"
                )
            )

    # Fork processes
    elif re.match(r"^(\-\-forkp=)", arg):
        match6 = re.match(r"\-(\-forkp)=(.*)", arg)
        if match6:
            max_processes = (
                match6.group(2)
                if match6.group(2)
                else sys.exit("\nEmpty argument. Please enter the parameter information.\n\neg. --forkp=8\n\n")
            )

    # Filter variants by impact
    elif re.match(r"^(\-\-filtering=)", arg):
        match7 = re.match(r"\-(\-filtering)=(.*)", arg)
        if match7:
            filtering = (
                eval(match7.group(2))
                if match7.group(2)
                else sys.exit("\nEmpty argument. Please enter the parameter information.\n\neg. --filtering=True")
            )

    else:
        sys.exit(f"\nArgument {arg} not valid.\n\n")


if vep_results_file == None:
    sys.exit(
        "\nVEP file not indicated. Please enter the absolute path to VEP sorted output file.\n\neg. -f=z13-222\n\n"
    )
if dbdir == None:
    sys.exit("\nPath to databases not indicated. Please, enter the databases path.\n\neg. -d=databases\n\n")

# Create folders
os.makedirs(dbdir, exist_ok=True)

# Start time counter
start = time.time()

# get new experiment id
if jobid:
    jobid = root_name + jobid + "_VEP"
else:
    jobid = get_runid()

outexp = outdir + "/" + jobid + "/"

# Create experiment folder in output folder
os.makedirs(outexp)

# fmt:off
## Load files into variables.
# Larger files has been previous converted into sqlite3 database format.
cosmic_conn = sqlite3.connect(f"{dbdir}/cosmic.db")
cosmic_db = cosmic_conn.cursor()

# Rest of files are read as pandas dataframes.
cosmic_gene_freq_db = pd.read_csv(f"{dbdir}/cosmic_gene_freq.tsv", index_col=0, sep="\t")
genes_ids_db = pd.read_csv(f"{dbdir}/custom.tsv", index_col=0, sep="\t").dropna()
kegg_gene_pathway_db = pd.read_csv(f"{dbdir}/gene_pathway.tsv", index_col=0, sep="\t")
pathw_desc_db = pd.read_csv(f"{dbdir}/pathway_desc.tsv", index_col=0, sep="\t")
pfam_a_db = pd.read_csv(f"{dbdir}/Pfam-A.full.tsv", index_col=4, sep="\t")
uniprot_b_db = pd.read_csv(f"{dbdir}/Uniprot.tsv", index_col=1, sep="\t")
gscore_db = pd.read_csv(f"{dbdir}/gscore_Ene_2023.tsv", index_col=0, sep="\t")
cancer_domain_db = pd.read_csv(f"{dbdir}/domains.tsv", index_col=4, sep="\t")
generole_db = pd.read_csv(f"{dbdir}/generole.tsv", index_col="gene", sep="\t")
interpro_db = pd.read_csv(f"{dbdir}/Interpro.tsv", sep="\t")
last_domain_db = interpro_db.sort_values("START", ascending=False).groupby("PROTEIN_ACC").first()
interpro_db = interpro_db.set_index("PROTEIN_ACC")
clinvar_db = pd.read_csv(f"{dbdir}/Clinvar.tsv", sep="\t")
clinvar_db["Clinvar_key"] = (clinvar_db["Chr"] + ":" + clinvar_db["VCF_pos"].astype(str)+ ":" 
                     + clinvar_db["VCF_ref"] + ":" + clinvar_db["VCF_alt"])
clinvar_db.set_index('Clinvar_key', inplace=True, drop=False)
# fmt:on


# Variable intialization
data = ""
# Open vcf file
try:
    with gzip.open(vep_results_file, "rt") as infile:
        rfile = infile.readlines()
        printl(f"\nProcessing file {vep_results_file}...\n")
except IOError:
    printl(f"Could not open file {vep_results_file}")
    sys.exit(1)

pos = {}
count = 0

# fmt: off
possible_vep_fields = [
    "Consequence", "IMPACT", "Existing_variation", "Feature", "PolyPhen","SIFT",
    "CADD_PHRED", "CADD_RAW", "SYMBOL", "Protein_position", "Amino_acids", "HGVSc",
    "HGVSp", "GMAF", "CDS_position", "Allele","Gene", "Feature_type", "cDNA_position",
    "Codons", "VARIANT_CLASS", "gnomADe_AF", "gnomADe_NFE_AF", "EXON","APPRIS"
    ]
# fmt: on

# Remove lines of metainformation (with '##'), while keeping the header line names
# and the data of each variant.
for n in range(len(rfile)):
    if re.match("^#[^#]", rfile[n]):
        rfile[n] = rfile[n].replace("#", "")
        data += rfile[n]

    elif re.match("^[^#]", rfile[n]):
        data += rfile[n]
        count += 1
    # Create a dictionary with Ensembl VEP annotations.
    elif re.match('^##INFO=<ID=CSQ.+Format: (.+)">', rfile[n]):
        match_fields = re.match('^##INFO=<ID=CSQ.+Format: (.+)">', rfile[n])
        vep_fields = match_fields.group(1).split("|")

        for field_index in range(len(vep_fields)):
            for field_name in possible_vep_fields:
                if vep_fields[field_index] == field_name:
                    pos[field_name] = field_index


# Save modifications to ensembl_vep.csv
outexp_path = re.sub(r"^\.", "", outexp)
outpath += outexp_path
outpathfile = outpath + "ensembl_vep.csv"

try:
    with open(outpathfile, "w") as sfile:
        sfile.write(data)
except IOError:
    printl("Could not save temp file\n")
    sys.exit(1)

printl("ensembl_vep.csv file created!")


# Split ensembl_vep.csv into multiple files based on the chromosome number.
with open(f"{outexp}/ensembl_vep.csv", "r") as input:
    ensemblhead = ""
    chr_num = ""
    filechr = ""
    chr_names = []
    for line in input:
        if re.match("^CHROM", line):
            ensemblhead = line
        else:
            line_splitted = line.split("\t")
            if line_splitted[0] != chr_num:
                chr_num = line_splitted[0]
                chr_names.append(chr_num)
                if filechr != "":
                    outputchr.close()
                filechr = f"{outexp}/ensembl_vep_{chr_num}.csv"
                outputchr = open(filechr, "w")
                outputchr.write(ensemblhead)
                outputchr.write(line)
            else:
                outputchr.write(line)
    outputchr.close()

printl(f"\n\nCreating annotations for {count} variants...")

# fmt: off
annotations = [
    "Chr\tLoc\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\tmut\tlocation\t",
    "allele\tgene\tfeature\tfeature_type\tconsequence\timpact\tcdna_position\t",
    "cds_position\tprotein_position\tamino_acids\tcodons\texisting_variation\t",
    "extra\tprincipal\tpoly_effect\tpoly_score\tsift_effect\tsift_score\tCADD_phred\t",
    "CADD_raw\tgene_hgnc\tgene_role\tcosmic_id\tKEGG_data\tKEGG_path_id\tclinvar_acc\t",
    "clinvar_disease\tclinvar_clinical_significance\tvariation_type\tHGVS_cDNA\t",
    "HGVS_protein\tGMAF\tgnomAD\tgnomAD_NFE\tpfam\tinterpro\tgene_cosmic_freq\t",
    "mut_cosmic_freq\tvscore\tbranch\n",
]

annotations_sort = [
    "Chr\tLoc\tmut\tgene\tfeature\tfeature_type\tconsequence\timpact\tprincipal\t",
    "poly_effect\tpoly_score\tsift_effect\tsift_score\tCADD_phred\tCADD_raw\t",
    "gene_hgnc\tgene_role\tcosmic_id\tkegg_data\tkegg_path_id\tprotein_position\t",
    "amino_acids\tclinvar_acc\tclinvar_disease\tclinvar_clinical_significance\t",
    "variation_type\tHGVS_cDNA\tHGVS_protein\tGMAF\tgnomAD\tgnomAD_NFE\tpfam\t",
    "interpro\tgene_cosmic_freq\tmut_cosmic_freq\tvscore\tbranch\n",
]


# VEP_parser_per_chr('chr1')
# fmt:off
with Pool(processes=max_processes) as p:
    p.map(VEP_parser_per_chr, chr_names)
cosmic_conn.close()

# Create a single file with all variants located in different chromosomes.
greatout = open(f"{outexp}/vep_data.csv", "w")
greatout.write("".join(annotations))
for chr in chr_names:
    input_chr = open(f"{outexp}/vep_data_{chr}.csv", "r")
    for line in input_chr:
        if not re.match(r"^Chr", line):
            greatout.write(line)
    input_chr.close()
    os.remove(f"{outexp}/vep_data_{chr}.csv")
    os.remove(f"{outexp}/ensembl_vep_{chr}.csv")
greatout.close()
printl("\nvep_data.csv file created!\n")


# Create a single sorted file with all variants located in different chromosomes.
greatout_sorted = open(f"{outexp}/vep_data_sorted.csv", "w")
greatout_sorted.write("".join(annotations_sort))

for chr in chr_names:
    input_chr = open(f"{outexp}/vep_data_sorted_{chr}.csv", "r")
    for line in input_chr:
        if not re.match(r"^Chr", line):
            greatout_sorted.write(line)
    input_chr.close()
    os.remove(f"{outexp}/vep_data_sorted_{chr}.csv")
greatout_sorted.close()
printl("vep_data_sorted.csv file created\n")


# Select the highest score for the principal isoform.
vepsorted_list = []
try:
    vep_sorted_file = open(f"{outpath}vep_data_sorted.csv")
except IOError:
    printl(f"Could not open file {outpath}vep_data_sorted.csv")
    sys.exit(1)

for line in vep_sorted_file:
    line = line.rstrip("\n").split("\t")
    vepsorted_list.append(line)
vep_sorted_file.close()

isoform = {}
appris_isoforms = ["P1", "P2", "P3", "P4", "P5", "A1|A2|^$"]
for variant in vepsorted_list:
    if re.match(r"HIGH|MODERATE", variant[7]):
        if variant[15] not in isoform.keys():
            isoform[variant[15]] = [0] * 6

        for index, isof in enumerate(appris_isoforms):
            if re.match(isof, variant[8]):
                isoform[variant[15]][index] = 1

for gene_name in isoform.keys():
    for index, isof in enumerate(appris_isoforms):
        if isinstance(isoform[gene_name], list) and isoform[gene_name][index] == 1:
            isoform[gene_name] = isof
            if isoform[gene_name] == "A1|A2|^$":
                isoform[gene_name] = "A1|A2|"

genes_affected = {}

for index in range(len(vepsorted_list)):
    if index == 0:
        head = "\t".join(vepsorted_list[index])
    else:
        if re.match(r"HIGH|MODERATE", vepsorted_list[index][7]):
            if re.match(f"^{isoform[vepsorted_list[index][15]]}$", vepsorted_list[index][8]):
                if vepsorted_list[index][15] in genes_affected.keys():
                    if vepsorted_list[index][35] >= genes_affected[vepsorted_list[index][15]][0]:
                        genes_affected[vepsorted_list[index][15]] = [vepsorted_list[index][35],vepsorted_list[index][36]]
                        genes_affected[vepsorted_list[index][15]].extend(vepsorted_list[index])
                else:
                    genes_affected[vepsorted_list[index][15]] = [vepsorted_list[index][35], vepsorted_list[index][36]]
                    genes_affected[vepsorted_list[index][15]].extend(vepsorted_list[index])

genes_affected_file = open(f"{outpath}genes_affected.csv", "w")
genes_affected_file.write(f"gene_hgnc\tmax(vscore)\tbranch\t{head}\n")
for gene_name in genes_affected.keys():
    genes_affected[gene_name] = "\t".join(genes_affected[gene_name])
    genes_affected_file.write(f"{gene_name}\t{genes_affected[gene_name]}\n")

genes_affected_file.close()
printl("\ngenes_affected.csv file created!\n\n")


end = time.time()
printl(f'\nTotal time: {"{:.2f}".format(end - start)} seconds\n')
