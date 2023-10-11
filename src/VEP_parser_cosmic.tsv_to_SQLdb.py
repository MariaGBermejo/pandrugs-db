import sys
import re
import pandas as pd
import sqlite3


# Information about commmand line arguments required by python script.
def help_info():
    print("\n\n--databases=directory or -d=directory \t\t\t Absolute path to databases directory. Mandatory.\n\n")

    print("\ni.e. VEP_parser.py -d=/home/epineiro/Programs/PCDA/databases\n\n")
    sys.exit()


# Command line arguments handle
if len(sys.argv) == 1 or any(re.match(r"^((\-\-help)|(\-h))$", arg) for arg in sys.argv[1:]):
    help_info()

if re.match(r"^((\-\-databases=)|(\-d=))", sys.argv[1]):
    match = re.match(r"\-(\-databases|d)=(.*)", sys.argv[1])
    if match:
        dbdir = (
            match.group(2)
            if match.group(2)
            else sys.exit(
                "\nEmpty argument. Please enter the parameter information.\n\neg."
                " -d=/home/epineiro/Programs/PCDA/databases\n\n"
            )
        )
else:
    sys.exit(f"\nArgument {sys.argv[1]} not valid.\n\n")

# Transform COSMIC.tsv file into a sqlite3 database.
cosmic_df = pd.read_csv(f"{dbdir}/COSMIC.tsv", index_col=0, sep="\t")
cosmic_db = sqlite3.connect(f"{dbdir}/cosmic.db")
cosmic_df.to_sql("cosmic_table", cosmic_db, if_exists="replace")
cosmic_db.close()
