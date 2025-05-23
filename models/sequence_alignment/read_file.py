from itertools import groupby
def fasta(fasta_name):
	read_file = open(fasta_name)
	faiter = (x[1] for x in groupby(read_file, lambda line: line[0] == ">"))
	for header in faiter:
		#header = header.next()[1:].strip()
		header = next(header)[1:].strip()
		seq = "".join(s.strip() for s in next(faiter))
		yield header, seq