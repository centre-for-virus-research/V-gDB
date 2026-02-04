import sqlite3

conn = sqlite3.connect("/Volumes/My Passport/CVR/gdb/Flu/flu-gDB_dec02.db")
conn.row_factory = sqlite3.Row  
cur = conn.cursor()

cur.execute("CREATE INDEX IF NOT EXISTS idx_meta_data_isolate ON meta_data(isolate);")
cur.execute("CREATE INDEX IF NOT EXISTS idx_meta_data_primary_accession ON meta_data (primary_accession)")
cur.execute("CREATE INDEX IF NOT EXISTS idx_meta_data_segment ON meta_data (segment)")
cur.execute("CREATE INDEX IF NOT EXISTS idx_meta_data_host ON meta_data (host)")
cur.execute("CREATE INDEX IF NOT EXISTS idx_meta_data_length ON meta_data (length)")
cur.execute("CREATE INDEX IF NOT EXISTS idx_meta_data_Parsed_strain ON meta_data (Parsed_strain)")
cur.execute("CREATE INDEX IF NOT EXISTS idx_meta_data_accession_type ON meta_data (accession_type)")
cur.execute("CREATE INDEX IF NOT EXISTS idx_meta_data_strain ON meta_data (strain)")

cur.execute("CREATE INDEX IF NOT EXISTS idx_isolates_strain ON isolates (strain)")
cur.execute("CREATE INDEX IF NOT EXISTS idx_sequences_header ON sequences (header)")
cur.execute("CREATE INDEX IF NOT EXISTS idx_sequence_alignment_sequence_id ON sequence_alignment (sequence_id)")
cur.execute("CREATE INDEX IF NOT EXISTS idx_sequence_alignment_alignment_name ON sequence_alignment (alignment_name)")


cur.execute("CREATE INDEX IF NOT EXISTS idx_features_accession ON features (accession)")
cur.execute("CREATE INDEX IF NOT EXISTS idx_features_reference_accession ON features (reference_accession)")

cur.execute("CREATE INDEX IF NOT EXISTS idx_cluster_members_primary_accession ON cluster_members (primary_accession)")

conn.commit()
conn.close()
