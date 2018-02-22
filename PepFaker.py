
def yes():
    #ThePeptideList = ['PHM732251','''GINKLRKNDIDVGCALVMSKTNMNHTDEIYDFMFENKLPFNVIPLNKSGNA
    #VDNYQDLGLGPDEYYLPWSKLYDKWFYAEPKNYIFVSDFVRKTQAILAGRAADCIGMAQCGNTSFSVDPVGDLYPCASLSAQPD
    #MKYGNLQNNSILELMSSTRATIYRTREGTDSCKKCKWQHVCHGGCPARAYKYHDNDISHKDYYCPSLYKM''']
    ThePeptideList = ['1g6rP', 'SIYRYYGL']
    name = ThePeptideList[0]
    peptide = ThePeptideList[1]
    #number = 1
    haha=pepfaker(name,peptide)
    return(haha)



def pepfaker(name,peptide):

    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML


    # Define output path
    name = name[:]
    new_path = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/Blast'
    new_path = new_path + name + '.xml'

    # Perform BLASTp on peptide (Adjusted for small peptides)
    blast_handle = NCBIWWW.qblast('blastp', 'nr', peptide, matrix_name=('PAM30'), gapcosts='9 1', expect=200000, word_size=2, perc_ident=100)

    pID = []

    # Read and find peptides
    records = NCBIXML.parse(blast_handle)
    item = next(records)
    for alignment in item.alignments:
        for hsp in alignment.hsps:

            # Decide which source peptide to use
            if hsp.align_length == 8: # Alignment has no gaps
                if hsp.identities == 8: # Every aa match (doesn't account for gaps)
                    if alignment.length > 100: # Picks relevant source peptides (not too small)
                        print('****Alignment****')
                        print('ID:', alignment.hit_id)
                        print('length:', alignment.length)
                        print ('ident:', hsp.identities)

                        IDID = alignment.hit_id
                        IDID = IDID.split('|')
                        pID.append(IDID[1])


    # Save FULL BLASTp result (Not necessary - can be deleted)
    blast_handle.seek(0)
    blast_file = open(new_path, 'w')
    blast_file.write(blast_handle.read())
    blast_file.close()

    # Choose the right sequence if multiple fits
    if len(pID) > 1:
        pNumber = pID[0]
    elif len(pID) == 0:
        pNumber = False
    else:
        pNumber = pID[0]

    return(pNumber)






