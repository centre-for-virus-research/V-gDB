def get_kuiken2006_codon_labeling(refStart, refEnd):
    
    condonStart = None
    codonEnd = None
    if ( refStart - (refEnd + 1) ) % 3 == 0:
        condonStart = 1
        codonEnd = round(((refEnd+1) - (refStart )) / 3,0)

    return [condonStart, codonEnd]