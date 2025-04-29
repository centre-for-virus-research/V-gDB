from django.db import connections
from models.helper_functions import dictfetchall  # Ensure this function is imported
from models.codon_labeling import get_kuiken2006_codon_labeling

class Mutations:
    def __init__(self, database, hosts=None, region=None, codons=None, reference_sequence=None):
        self.database = database
        self.hosts = hosts
        self.region = region
        self.codons = codons
        self.reference_sequence = reference_sequence


    def get_mutations(self):

        if not self.reference_sequence:
            self.reference_sequence = 'NC_001542'

        if not self.hosts:
            raise ValueError("No hosts chosen")
        
        if not self.region:
            raise ValueError("No region chosen")
        
        if not self.codons:
            raise ValueError("No codons chosen")
       
        alignments = self.__get_alignments()
        master_reference = self.__get_master_reference()

        mutations = self.__parse_mutations(alignments, master_reference)

        return mutations
    
    def get_mutation_regions_and_codons(self):

        if not self.reference_sequence:
            self.reference_sequence = 'NC_001542'

        master_reference = self.__get_master_reference()
        for feature in master_reference:
            codon_start, codon_end = get_kuiken2006_codon_labeling(feature["ref_start"], feature["ref_end"])
            feature["codon_start"] = codon_start
            feature["codon_end"] = codon_end

        return master_reference
    
    def get_sequence_mutations(self, primary_accession, test_sequence):

        if not self.reference_sequence:
            self.reference_sequence = 'NC_001542'

        with connections[self.database].cursor() as cursor:
            cursor.execute('SELECT * FROM sequence_alignment WHERE sequence_id = %s;', [primary_accession])
            alignments = dictfetchall(cursor)

        master_reference = self.__get_master_reference()
        mutations = self.__compare_sequences(master_reference, alignments[0]["alignment"], test_sequence)

        return mutations

    def __compare_sequences(self, master_reference, sequence, reference):
        mutations = {}
        mutations["codons"] = []
        mutations["nucleotides"] = []

        for i in range(len(master_reference)):
            ref_start = master_reference[i]["cds_start"]
            ref_end = master_reference[i]["cds_end"]

            sub_seq = sequence[ref_start:ref_end+1]
            sub_ref = reference[ref_start:ref_end+1]
            seq_codons = [sub_seq[i:i+3] for i in range(0, len(sub_seq), 3)]
            ref_codons = [sub_ref[i:i+3] for i in range(0, len(sub_ref), 3)]
            
            for codon in range(1, len(seq_codons)):
                seq_codon = seq_codons[codon]
                ref_codon = ref_codons[codon]
                if seq_codon != '---' and ref_codon != '---':
                    if seq_codon != ref_codon:
                        seq_protein = self.__translate(seq_codon)
                        ref_protein = self.__translate(ref_codon)
                        if (seq_protein != ref_protein):
                            mutations["codons"].append({"codon":codon, "seq_protein":seq_protein, "ref_protein":ref_protein})
            # print(ref_start)
            for nucleotide in range(ref_start, ref_end+1, 1):
                seq_nucleotide = sequence[nucleotide]
                ref_nucleotide = reference[nucleotide]

                if seq_nucleotide != ref_nucleotide:
                    if seq_nucleotide != '-':
                        mutations["nucleotides"].append({"nucleotide":nucleotide, "seq_nucleotide":seq_nucleotide, "ref_nucleotide":ref_nucleotide})

        return mutations
    
    # PRIVATE FUNCTIONS
    def __parse_mutations(self, alignments, master_reference):
        print(master_reference)
        ref_start = master_reference[0]["cds_start"]
        ref_end = master_reference[0]["cds_end"]
        
        for i in range(len(alignments)):
            alignments[i]["mutations"] = {}
            sub_seq = alignments[i]["alignment"][ref_start:ref_end+1]
            codons_new = [sub_seq[i:i+3] for i in range(0, len(sub_seq), 3)]
            if len(codons_new) > 0:
                for codon in self.codons:
                    selected_codons = codons_new[int(codon)]
                    sub_seq = ''.join(selected_codons)
                    if "-" not in sub_seq:
                        if len(sub_seq)%3 == 0:
                            protein = self.__translate(sub_seq)
                            alignments[i]["mutations"][str(codon)] = protein
                    else:
                        if sub_seq != "---":
                            alignments[i]["mutations"][str(codon)] = "X"

        return alignments


    def __get_alignments(self):
        formatted_hosts = ', '.join(['%s'] * len(self.hosts))  
        query = f'SELECT s.*, m.* FROM sequence_alignment s LEFT JOIN meta_data m on s.sequence_id = m.primary_accession WHERE m.host IN ({formatted_hosts});'

        with connections[self.database].cursor() as cursor:
            cursor.execute(query, self.hosts)
            alignments = dictfetchall(cursor)

        return alignments
    
    
    # THIS FUNCTION CAN BE REFACTORED MORE TO JUST GRAB EVERYTHING
    def __get_master_reference(self):
        print(self.reference_sequence)
        with connections[self.database].cursor() as cursor:
            if self.region:
                cursor.execute("SELECT cds_start, cds_end FROM features WHERE master_ref_accession=%s AND product=%s", [self.reference_sequence, self.region])
            else: 
                cursor.execute("SELECT * FROM features WHERE master_ref_accession=%s", [self.reference_sequence])
            master_reference = dictfetchall(cursor)

        return master_reference
    
    def __translate(self, codon): 
       
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
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
        } 
        protein = "" 
        protein = table.get(codon, "-")
        return protein 