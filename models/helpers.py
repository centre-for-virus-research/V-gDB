import csv

def get_codon_labeling(refStart, refEnd):
    
    condonStart = None
    codonEnd = None
    if ( refStart - (refEnd + 1) ) % 3 == 0:
        condonStart = 1
        codonEnd = round(((refEnd+1) - (refStart )) / 3,0)

    return [condonStart, codonEnd]


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
    print(results)
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

def build_clade_tree(clades):
    nodes = {}
    
    for clade in clades:
        major = clade['major_clade']
        minor = clade['minor_clade']
        if major not in nodes:
            nodes[major] = {
                'name':major,
                'text':major,
                'parent':None,
                'nodes':[]
            }
        if minor != None:
            nodes[major]['nodes'].append({
                'name':minor,
                'text':minor,
                'parent':major
            })
    
    tree = []
    for clade in nodes.values():
        if len(clade['nodes']) == 0:
            clade['nodes'] = None
    
    tree = list(nodes.values())

    return(tree)



def build_feature_tree(features):

    tree = []
    nodes = {}

    for feature in features:
        name = feature['name']
        if name == 'whole_genome,NULL':
            name = 'whole_genome'
        display_name = feature['description'] if feature['description'] is not None else name
        parent_name = feature['parent_name']
        nodes[name] = {
            'name': name,
            'text': display_name,
            'parent': parent_name,
            'nodes': []
        }
    for feature in features:
        name = feature['name']
        if name == 'whole_genome,NULL':
            name = 'whole_genome'
        parent_name = feature['parent_name']
        display_name = feature['description']

        if parent_name is None:
            tree.append(nodes[name])
        else:
            nodes[parent_name]['nodes'].append(nodes[name])

    def clean_tree(node):
        if len(node['nodes']) == 0:
            node['nodes'] = None
        else:
            for child in node['nodes']:
                clean_tree(child)

    for root in tree:
        clean_tree(root)
    

    return(tree)
