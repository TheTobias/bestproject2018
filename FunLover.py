from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.SubsMat import MatrixInfo
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio import Entrez
import xlrd
import os
import matplotlib.pyplot as plt
Entrez.email = 'ohsaibot@mail.com'
blosum = MatrixInfo.blosum62


class PeptideFun(object):

    def __init__(self, name='', peptide='', proID=''):

        self.n = name
        self.p = peptide
        self.id = proID

# Blasts a peptide and finds ID of the best match for source protein ######
# Only works for small peptide ############################
    @property
    def blastIDfinder(self):


        # Perform BLASTp on peptide (Adjusted for small peptides)
        blast_handle = NCBIWWW.qblast('blastp', 'nr', self.p, matrix_name='PAM30', gapcosts='9 1', expect=200000,
                                      word_size=2, perc_ident=100)

        # Saving xml file (only relevant if you work with the code and don't want to run the blast all the time)
        new_path = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/thisisblasttest.xml'
        with open(new_path, "w") as out_handle:
           out_handle.write(blast_handle.read())
        blast_handle.close()
        blast_handle = open(new_path)

        # Defining output
        pep2=PeptideFun('','')

        # Read and find peptides
        records = NCBIXML.parse(blast_handle)
        item = next(records)
        ident = 0
        mismatch = [0, 0]
        for alignment in item.alignments:
            for hsp in alignment.hsps:

                # Decide which source peptide to use
                if hsp.align_length == len(self.p):  # Alignment has no gaps

                    if alignment.length > 100:  # Picks relevant source peptides (not too small)

                        if hsp.identities == len(self.p):  # Every aa match (doesn't account for gaps)
                            IDID = alignment.hit_id
                            pep2.n = alignment.hit_def
                            IDID = IDID.split('|')
                            pep2.id = IDID[1]
                            mismatch = [0, 0]
                            return(pep2,mismatch)
                        elif hsp.identities > ident:
                            ident=hsp.identities
                            IDID = alignment.hit_id
                            pep2.n = alignment.hit_def
                            IDID = IDID.split('|')
                            mismatch[0] = hsp.query
                            mismatch[1] = hsp.sbjct
                            pep2.id = IDID[1]



        return(pep2, mismatch)


# Switches your peptide out with another in a given complex ###########
    def pepswitch(self, complex1):

        # Split complex file and peptide file and insert the peptide
        newcomplex = complex1.split()
        newcomplex[2] = self.n
        newcomplex[3] = self.p

        return newcomplex

# Trade ID number for protein sequence from NCBI #######################
    def NCBIseqGet(self):

        # Float where output is saved
        pep2 = PeptideFun('', '')

        if self.id == '':
            print('You need an ID in self.ID to use the NCBIseqGet function')
            return pep2

        # The use of Entrez to find and get sequence for ID
        request = Entrez.epost("protein", id=self.id)
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
            pep2.id = gi
            pep2.p = r["GBSeq_sequence"]

            return pep2

# Finds peptides with binding affinity and gives blosum62#############
#compared to original peptide we are looking for (uses netMHCpan4 data as input)######
    def MHCpanBLOSUM(self, worksheet):

        score = []
        vals = []
        fakepep = []
        maxi = worksheet.nrows

        for i in range(maxi):

            if worksheet.cell(i, 9).value == 1 and str(worksheet.cell(i, 1).value) != self.p:
                seq2 = str(worksheet.cell(i, 1).value)
                score.append(nogapscore(self.p, seq2, blosum))
                fakepep.append(seq2)

                if float(worksheet.cell(i, 7).value) > 1:
                    vals.append(float(worksheet.cell(i, 7).value / 10000))
                    i += 1

                else:
                    vals.append(float(worksheet.cell(i, 7).value))
                    i += 1

            else:
                i += 1


        return(score, fakepep)

    # Takes a netmhcpan4 output and finds peptides with closets rank% to self peptide
    def MHCpanRank(self, worksheet):

        maxi = worksheet.nrows
        rankP = 0

        # Finds the rank% of the self peptide
        for i in range(maxi):
            if str(worksheet.cell(i, 1).value) == self.p:
                rankP = float(worksheet.cell(i, 7).value)
                if float(worksheet.cell(i, 7).value) > 1:
                    rankP = (float(worksheet.cell(i, 7).value / 10000))

        # Predefine the outputs
        diffs1 = float(100)
        diffs2 = float(100)
        diffs3 = float(100)
        peps = ['', '', '']

        # Go through every peptides rank
        for i in range(maxi):
            # Make sure we don't use the self peptide and only rows with data on it
            if str(worksheet.cell(i, 1).value) != self.p and worksheet.cell(i, 9).value in (0,1):

                # Some data is corrupt and needs to be adjusted back (different excel versions)
                rank = float(worksheet.cell(i, 7).value)
                if float(worksheet.cell(i, 7).value) > 1:
                    rank = (float(worksheet.cell(i, 7).value / 10000))

                # Calculate the difference in rank between self peptide and a peptide
                diffs = abs(rank - rankP)

                # If the difference is less than diffs3 it is changed
                if diffs < diffs3:
                    diffs3 = diffs
                    peps[2] = str(worksheet.cell(i, 1).value)

                    # This code should be made into a function
                    # If the difference is less than diffs2 it is changed
                    if diffs3 < diffs2:
                        changer = diffs2
                        diffs2 = diffs3
                        diffs3 = changer
                        changer1 = peps[1]
                        peps[1] = peps[2]
                        peps[2] = changer1

                        # If the difference is less than diffs1 it is changed
                        if diffs2 < diffs1:
                            changer = diffs1
                            diffs1 = diffs2
                            diffs2 = changer
                            changer1 = peps[0]
                            peps[0] = peps[1]
                            peps[1] = changer1

            i += 1

        # Rounding of the diffs and inserted as easier to work with output
        diffs1 = round(diffs1, 4)
        diffs2 = round(diffs2, 4)
        diffs3 = round(diffs3, 4)
        thediffs = [diffs1, diffs2, diffs3]
        return (rankP, peps, thediffs)

# Extra functions ########################################################
# Score a pair
def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]

# Score two sequences with no gaps
def nogapscore(seq1,seq2,matrix):
    score = 0

    for i in range(len(seq1)):
        pair = (seq1[i], seq2[i])
        score+=score_match(pair, matrix)

    return score

# Choose best BLOSUM scores and order them
def bestBLOSUM(blosumscore,peptides,numberofpeps):

    bests = sorted(zip(blosumscore, peptides), reverse=True)[:numberofpeps]

    return bests

# Removes gaps and xxx's from crystallography and returns aa sequence
def aaCleaner(aaSequence):

    # Define the correct characters (the 20 normal natural aa's)
    AAs20 = 'ARNDCEQGHILKMFPSTWYV'

    # Make the sequence which has to be cleaned into a list
    aaSequence = list(aaSequence)

    # If the character isn't in AAs20 it is switched with nothing ('')
    myseq = Seq(AAs20, generic_protein)
    i = 0
    while i < len(aaSequence):
        if aaSequence[i] in myseq:
            i += 1
        else:
            aaSequence[i] = ''
            i += 1

    # Make the list into a string again
    cleanedseq = ''.join(aaSequence)

    return(cleanedseq)

# Changes a sequence inside a longer sequence with another sequence
def seqchanger(fullseq, partchange, partchangeto):

    # Make sure all sequences are lower case (as lower and upper isn't the same)
    partchange = partchange.lower()
    fullseq = fullseq.lower()
    partchangeto = partchangeto.lower()

    # Remove the sequence you want to change by splitting it at that spot
    newseq = fullseq.split(partchange)

    # Insert the new sequence in between the split
    newseq.insert(1, partchangeto)

    # Join the list back together to creat the new sequence
    newseq = ''.join(newseq)

    return (newseq)

# Shortens a sequence to a length of under 20000 while keeping the important sequence
def shortener20000(sequence2short, importantPep):

    # Make sure both sequences are lowercase
    sequence2short = sequence2short.lower()
    importantPep = importantPep.lower()

    # Split at the important peptide
    newseq = sequence2short.split(importantPep)

    # Pick the shorten part which has the important peptide
    if len(newseq[0]) > len(newseq[1]):
        newseq = sequence2short[-19999:]
    else:
        newseq = sequence2short[19999:]

    return(newseq)

# Get the first worksheet from a xlsx file
def getworksheet(xlsxfile):
    if os.path.isfile(xlsxfile) == True:
        excelworkbook = xlrd.open_workbook(xlsxfile)
        excelworksheet = excelworkbook.sheet_by_index(0)

        return (excelworksheet)
    return(0)

# Switch the found peptides out with the real one
def PepSwitch3(complex,peptides,complexpath):

# Make amount of switches equal to amount of peptides you can switch
    for i in range (len(peptides)):
        names = '>P'+'{}'.format(i+1)
        pep1 = PeptideFun(name=names, peptide=peptides[i])

# Use the pepswitch function to do the switch
        newcomplex = pep1.pepswitch(complex)

# Define output path for the new complex
        new_path = complexpath.rstrip('.fsa')
        new_path = new_path + 'fake' + str(i+1) + '.fsa'

# Write new complex with switched peptide
        ncf = open(new_path, 'w')
        ncf.write('\n'.join(newcomplex))
        ncf.close()

# Scatter plot with multiple y's per x's (x's are tekst)
def plotmepls1(x,y,xlabel,ylabel):
    # For naming the x axis
    labels = x

    # plotting the y's for each x
    for xe, ye in zip(x, y):
        plt.scatter([xe] * len(ye), ye)

    # Defining the number of ticks (depending on the number of x's)
    ticks = []
    for i in range(len(x)):
        ticks.append(i)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


    # Inserting the ticks and rotating them
    plt.xticks(ticks, labels, rotation='vertical')
    plt.axes().set_xticklabels(x)

    plt.show()