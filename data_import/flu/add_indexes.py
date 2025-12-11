import sqlite3

conn = sqlite3.connect("/Volumes/My Passport/CVR/gdb/Flu/flu-gDB_dec02.db")
conn.row_factory = sqlite3.Row  
cur = conn.cursor()

# # Drop old table if it exists
# cur.execute("CREATE INDEX idx_meta_data_isolate ON meta_data (isolate)")
# cur.execute("CREATE INDEX idx_meta_data_primary_accession ON meta_data (primary_accession)")
# cur.execute("CREATE INDEX idx_meta_data_segment ON meta_data (segment)")
# cur.execute("CREATE INDEX idx_meta_data_host ON meta_data (host)")
# cur.execute("CREATE INDEX idx_meta_data_length ON meta_data (length)")
cur.execute("CREATE INDEX idx_meta_data_Parsed_strain ON meta_data (Parsed_strain)")


# cur.execute("CREATE INDEX idx_isolates_strain ON isolates (strain)")
# cur.execute("CREATE INDEX idx_sequences_header ON sequences (header)")
# cur.execute("CREATE INDEX idx_sequence_alignment_sequence_id ON sequence_alignment (sequence_id)")
# cur.execute("CREATE INDEX idx_sequence_alignment_alignment_name ON sequence_alignment (alignment_name)")


# cur.execute("CREATE INDEX idx_features_accession ON features (accession)")
# cur.execute("CREATE INDEX idx_features_reference_accession ON features (reference_accession)")



conn.commit()
conn.close()
