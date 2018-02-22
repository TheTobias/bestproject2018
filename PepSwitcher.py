
complexpath = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/1g6r.fsa'
peptidepath = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/4abc.fsa'


def pepswitch(complexpath,peptidepath,number):


    #Define output path
    new_path = complexpath.rstrip('.fsa')
    new_path = new_path+'fake'+str(number)+'.fsa'

    #Open and Save complex file
    cf = open(complexpath,'r')
    complex = cf.read()
    cf.close()

    #Open and Save peptide file
    pf = open(peptidepath,'r')
    peptide = pf.read()
    pf.close()

    #Split complex file and peptide file and insert the peptide
    newcomplex=complex.split()
    peptide=peptide.split()
    newcomplex[2]=peptide[0]
    newcomplex[3]=peptide[1]

    #Write new complex with false peptide
    ncf = open(new_path,'w')
    ncf.write('\n'.join(newcomplex))
    ncf.close()

    return

