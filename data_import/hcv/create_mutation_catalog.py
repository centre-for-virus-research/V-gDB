import pandas as pd
import sqlite3

conn = sqlite3.connect("/Volumes/My Passport/CVR/gdb/HCV/HCV_full_new_Usher.db")
conn.row_factory = sqlite3.Row  
cur = conn.cursor()

# # Drop old table if it exists
cur.execute("DROP TABLE IF EXISTS mutation_catalog;")

# Read TSV
df = pd.read_csv("/Volumes/My Passport/CVR/GLUE/generalized_mutation_catalog_with_extra_info.tsv", sep="\t")



# Create table and insert data
df.to_sql(
    "mutation_catalog",
    conn,
    if_exists="replace",
    index=False
)

cur.execute("DROP TABLE IF EXISTS drug_regimen;")

df = pd.read_csv("/Volumes/My Passport/CVR/GLUE/output.csv")
# Create table and insert data
df.to_sql(
    "drug_regimen",
    conn,
    if_exists="replace",
    index=False
)


cur.execute("DROP TABLE IF EXISTS result_regimen;")

df = pd.read_csv("/Volumes/My Passport/CVR/GLUE/mysql_exports/phdr_result_regimen.csv")
# Create table and insert data
df.to_sql(
    "result_regimen",
    conn,
    if_exists="replace",
    index=False
)

cur.execute("DROP TABLE IF EXISTS result_trial;")

# df = pd.read_csv("/Volumes/My Passport/CVR/GLUE/mysql_exports/phdr_result_trial.csv")
df = pd.read_csv("/Volumes/My Passport/CVR/GLUE/merged_result_trial.csv")

# id_merged,phdr_in_vivo_result_id,phdr_clinical_trial_id,id,trial_display_name,trial_nct_id
columns = [
    "id",
    "phdr_in_vivo_result_id",
    "phdr_clinical_trial_id",
    "trial_display_name",
    "trial_nct_id"
]

df = df[columns]
# Create table and insert data
df.to_sql(
    "result_trial",
    conn,
    if_exists="replace",
    index=False
)

df = pd.read_csv("/Volumes/My Passport/CVR/GLUE/mysql_exports/phdr_in_vivo_result.csv")

# id,phdr_in_vivo_result_id,phdr_clinical_trial_id
for row in df.itertuples(index=False):
    cur.execute("""
              INSERT INTO result_trial (
                id, 
                phdr_in_vivo_result_id,
                phdr_clinical_trial_id
              )
              VALUES (?, ?, ?)
          """, (
              row.id,
              row.phdr_resistance_finding_id,
              row.cohort_description
          ))

conn.commit()
conn.close()