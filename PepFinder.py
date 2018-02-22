from FunLover import *


# Runs upto netMHCpan4 once
def yeho():
    complexworksheet = getworksheet('/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/database_tpm_final1.xlsx')

    ThePeptideList = complexworksheet.cell(70, 0).value
    ThePeptideList = ThePeptideList.split(',')

    #a=PepSwitch1(ThePeptideList)
    #print(a)
    ##peptide=(PepSwitch2(a[0], a[1], a[2], ThePeptideList))
    ##peptide1=(peptide[0][1],peptide[1][1],peptide[2][1])

    ##complexpath = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/1g6r.fsa'
    # Open and Save complex file
    ##cf = open(complexpath, 'r')
    ##complex = cf.read()
    ##cf.close()

    # Switch
    ##PepSwitch3(complex,peptide1,complexpath)

# Script making the new fake complexes from netMHCpan4 results
def yehu():

    # Load excel with all complexes
    complexworksheet = getworksheet('/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/database_tpm_final1.xlsx')

    # Define the number of iterations
    maxi = complexworksheet.nrows
    # maxi = 5

    # Used to keep track of the number of complexes
    j = 1

    # Define outputs for plots
    toplotx = []
    toploty = []

    # Go through each excel row
    for i in range(maxi):

        # Skip empty excel rows
        if complexworksheet.cell(i, 0).value != '':

            # Take out the complex number i from the excel sheet
            ThePeptideList = complexworksheet.cell(i, 0).value
            ThePeptideList = ThePeptideList.split(',')

            # Use PeptideFun to define the peptide as an object
            pep1 = PeptideFun(name=ThePeptideList[0], peptide=ThePeptideList[10])

            # Define output path so every complex has a different but relevant name
            directorypath = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/FakeComplexResults/'
            new_path = directorypath + str(ThePeptideList[0]) + 'real.fsa'

            # Make the real fasta from the excel
            name = '>' + ThePeptideList[0]
            ncf = open(new_path, 'w')
            ncf.write(name + ':M' + '\n' + ThePeptideList[8] + '\n' + name + ':P' + '\n' + ThePeptideList[10] + '\n')
            ncf.write(name + ':A' + '\n' + ThePeptideList[3] + '\n' + name + ':B' + '\n' + ThePeptideList[5] + '\n')
            ncf.close()

            # Take the netMHCpan4 results corresponding to the peptide in the complex chosen above
            df = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/netMHCpanResults/netmhcpan'
            dfpath = df + str(j) + '.xlsx'
            worksheet = getworksheet(dfpath)

            # If the netMHCpan4 worksheet exists keep on going, otherwise skip
            if worksheet != 0:

                # Use MHCpanRank to find 3 peptides with the closets Rank% to self peptide
                differences = pep1.MHCpanRank(worksheet)

                # Save the self peptide together with the differences of the 3 peptides for plotting
                toplotx.append(ThePeptideList[0])
                toploty.append(differences[2])

                # Take out the sequences of the peptides with a rank close to
                peptide1 = (differences[1][0], differences[1][1], differences[1][2])

                # Open the output of the complex with the new peptides in it
                cf = open(new_path, 'r')
                complex = cf.read()
                cf.close()

                # The actual switch of the peptide
                PepSwitch3(complex, peptide1, new_path)

            else:
                print(0)

            j += 1


    # Plotting of the differences for a easy overview of the rank% differences
    xlabel = 'TCR-pMHC complex (number)'
    ylabel = 'Rank% difference from self peptide rank%'
    plotmepls1(toplotx,toploty,xlabel,ylabel)

# Finds the source protein for all the peptides and saves it in a fasta file
def yeha():

    # Get all complexes
    complexworksheet = getworksheet('/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/database_tpm_final1.xlsx')

    # Define output path for the new complex
    new_path = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/exampletest3.fsa'
    ncf = open(new_path, 'w')

    maxi = complexworksheet.nrows
    #maxi = 5
    j=1
    for i in range(maxi):
        if (complexworksheet.cell(i, 0).value) != '':
            ThePeptideList = complexworksheet.cell(i, 0).value
            ThePeptideList = ThePeptideList.split(',')
            pep1 = PeptideFun(name=ThePeptideList[0], peptide=ThePeptideList[10])
            print(j)

            # Find the ID and the part that needs to be changed
            IDresult = pep1.blastIDfinder
            sourceP = IDresult[0]
            seqtochange = IDresult[1]

            # Uses the ID to get the source protein
            sourceP = sourceP.NCBIseqGet()

            # If part in source has to be changed to match original peptide
            if seqtochange[0] != 0:
                sourceP.p = seqchanger(sourceP.p, seqtochange[1], seqtochange[0])
                print('Changed sourceP')

            # If source protein is larger than 20000 it need to be shorten (netMHCpan4 limit is 20000)
            if len(sourceP.p) > 20000:
                sourceP.p = shortener20000(sourceP.p, pep1.p)
                print('Shorten sourceP')

            print(sourceP.p)

            # Define inputs (source protein, HLA type and length of cleaved peptides)
            HLAtype = ThePeptideList[7]
            pepLength = len(ThePeptideList[10])

            # Write new complex with switched peptide
            ncf.write(str(j) + '|' + HLAtype + '|' + str(pepLength) + '\n')
            ncf.write(sourceP.p + '\n' )
            j += 1



    ncf.close()


