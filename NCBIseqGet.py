from Bio.Blast import NCBIXML
from Bio import Entrez

Entrez.email = 'ohsaibot@mail.com'
from Bio import SeqIO


def yay():
    xmlfil = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/BlastPHM7322511.xml'
    blast_handle = open(xmlfil, "r")


    pID = []

    # Read and find peptides
    records = NCBIXML.parse(blast_handle)
    item = next(records)
    for alignment in item.alignments:
        for hsp in alignment.hsps:

            # Decide the e-val limit for being picked 3e-147
            if hsp.expect < 3e-100:
                #print('****Alignment****')
                #print('ID:', alignment.hit_id)
                IDID = alignment.hit_id
                IDID = IDID.split('|')
                pID.append(IDID[1])

    # Choose the right sequence if multiple fits
    if len(pID)!=1:
        pNumber = pID[0]
    else:
        pNumber = pID[0]

    return pNumber
    #seq=NCBIseqGet('1269676692')
    #print(seq)
proteinIDnumber = "1269676692"

def NCBIseqGet(proteinIDnumber):
    # Float where output is saved
    gottenID = []

    # The use of Entrez to find and get sequence for ID
    request = Entrez.epost("protein", id=proteinIDnumber)
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    handle = Entrez.efetch(db="protein", retmode="xml", webenv=webEnv, query_key=queryKey)
    for r in Entrez.parse(handle):
        # Grab the GI
        try:
            gi = int([x for x in r['GBSeq_other-seqids'] if "gi" in x][0].split("|")[1])
        except ValueError:
            gi = None
        gottenID.append(gi)
        gottenID.append(r["GBSeq_sequence"])

        return gottenID



