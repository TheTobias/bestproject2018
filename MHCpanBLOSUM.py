

import xlrd
from FunLover import *
from Bio.SubsMat import MatrixInfo
blosum = MatrixInfo.blosum62

def yolo():
    df = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/32124_NetMHCpan.xls'
    workbook = xlrd.open_workbook(df)
    worksheet = workbook.sheet_by_index(0)
    self = 'SIYRYYGL'
    a = MHCpanBLOSUM(self,worksheet)
    return a


def MHCpanBLOSUM(self, worksheet):

    score = []
    vals = []
    maxi = worksheet.nrows


    for i in range (maxi):

        if worksheet.cell(i, 9).value == 1 and str(worksheet.cell(i, 1).value) != self:
            seq2 = str(worksheet.cell(i, 1).value)
            score.append(nogapscore(self, seq2, blosum))

            if float(worksheet.cell(i, 7).value) > 1:
                vals.append(float(worksheet.cell(i, 7).value / 10000))
                i += 1

            else:
                vals.append(float(worksheet.cell(i, 7).value))
                i += 1

        else:
            i += 1

    return(score,vals)

print(yolo())

#Define output path
new_path = '/Users/TobiasHO/MasterPrograms/THEFOLDER/testoutput.fsa'

    #Write new complex with false peptide
ncf = open(new_path,'w')
ncf.write(str(yolo()))
ncf.close()
