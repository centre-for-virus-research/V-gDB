import csv

def get_codon_labeling(refStart, refEnd):

    codonStart = None
    codonEnd = None

    try:
        refStart = int(refStart)
        refEnd = int(refEnd)
    except (TypeError, ValueError):
        return [None, None]

    if (refStart - (refEnd + 1)) % 3 == 0:
        codonStart = 1
        codonEnd = (refEnd + 1 - refStart) // 3

    return [codonStart, codonEnd]


def dictfetchall(cursor):
    """Returns all rows from a cursor as a list of dictionaries."""
    columns = [col[0] for col in cursor.description]
    return [
        dict(zip(columns, row))
        for row in cursor.fetchall()
    ]


def fetch_one(cursor, query, params):
    print(query, params)
    cursor.execute(query, params)
    results = dictfetchall(cursor)

    if len(results) > 0:
        results = results[0]
    return results



def fetch_all(cursor, query, params):
    print(query, params)
    cursor.execute(query, params)
    return dictfetchall(cursor)


def build_csv_file(data, file_name):
    """
    Export metadata to a CSV file.

    Args:
        data (list): List of dictionaries containing metadata to write.
        file_name (str): Path to the output CSV file.
    """
    print("DATA", data)
    with open(file_name, "w", newline="") as f:
        writer = csv.writer(f)
        # Write header
        writer.writerow(data[0].keys())
        # Write rows
        for row in data:
            writer.writerow(row.values())


def build_fasta_file(data, filename):

    ofile = open(filename, "w")

    for row in data:
        ofile.write(row)
    ofile.close()

    return
    

def translateCodon(codon):
    """
    Translate a DNA codon into an amino acid using the standard genetic code.

    Args:
        codon (str): DNA triplet string.

    Returns:
        str: Single-letter amino acid or '-' if invalid.
    """
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
    }
    return table.get(codon, "-")
