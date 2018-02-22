from FunLover import *
import matplotlib.pyplot as plt


def haha():

    #get all complexes
    complexworksheet = getworksheet('/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/database_tpm_final1.xlsx')

    ThePeptideList = complexworksheet.cell(120, 0).value
    ThePeptideList = ThePeptideList.split(',')
    pep1 = PeptideFun(name=ThePeptideList[0], peptide=ThePeptideList[10])

    name = '>' + ThePeptideList[0]

    # Define output path so every complex has a different but relevant name
    directorypath = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/FakeComplexResults/'
    new_path = directorypath + str(ThePeptideList[0]) + 'real.fsa'

    # Make the real fasta from the excel
    ncf = open(new_path, 'w')
    ncf.write(name + ':M' + '\n' + ThePeptideList[8] + '\n' + name + ':P' + '\n' + ThePeptideList[10] + '\n')
    ncf.write(name + ':A' + '\n' + ThePeptideList[3] + '\n' + name + ':B' + '\n' + ThePeptideList[5] + '\n')
    ncf.close()

    # netMHCpan4 results
    df = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/netMHCpanResults/netmhcpan'
    dfpath = df + str(61) + '.xlsx'

    #j += 1

    #peptide = PepSwitch2(dfpath, ThePeptideList)

    worksheet=getworksheet(dfpath)

    if worksheet != 0:

        #RankP = MHCpanselfRank(pep1,worksheet)

        differences = pep1.MHCpanRank(worksheet)


        plotmepls(differences[0], differences[2])


        ##peptide1 = (differences[1][0], differences[1][1], differences[1][2])
        # Switch
        ##cf = open(new_path, 'r')
        ##complex = cf.read()
        ##cf.close()

        ##PepSwitch3(complex, peptide1, new_path)
    else:
        print(0)


def plotmepls(mainpoint, diffpoints):
    x = [(mainpoint)]
    y = [(diffpoints[0],diffpoints[1],diffpoints[2])]

    print(x)
    print(y)

    for xe, ye in zip(x, y):
        plt.scatter([xe] * len(ye), ye)
    #plt.plot(x, y, s=60, c='red', marker='^')


    plt.ylim(0, 10)
    plt.xticks([1])
    #plt.axes().set_xticklabels([1])


    plt.show()






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





haha()