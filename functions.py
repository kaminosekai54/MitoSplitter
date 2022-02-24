# Import
import os
import pandas as pd
from Bio import SeqIO
import re
from datetime import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

################################################################################

# function setup,
# this function will create all the folder need for the app to start if the don't 
# already exist
def setup():
    if not os.path.isdir("./logs"):
        os.makedirs("./logs")

    if not os.path.isdir("./ressources"):
        os.makedirs("./ressources")

    if not os.path.isdir("./ressources/raw_data"):
        os.makedirs("./ressources/raw_data")

    if not os.path.isdir("./ressources/results"):
        os.makedirs("./ressources/results")


# function getCSVFiles,
# this function return a list of all the csv files
# found in the indicated folder
# @param,
#@ path, the path of the folder yu want to parth
# @ extension, extension of the file we search for, here "csv by default"
def getCSVFiles(path, extension = ".csv"):
    fileList = os.listdir(path)
    return [ csvFile for csvFile in fileList if csvFile.endswith( extension) ]

# function getFASTAFiles,
# this function return a list of all the fasta files
# found in the indicated folder
# @param,
#@ path, the path of the folder yu want to parth
# @ extension, extension of the file we search for, here "fasta by default"
def getFASTAFiles(path, extension = ".fasta"):
    fileList = os.listdir(path)
    return [ fastaFile for fastaFile in fileList if fastaFile.endswith( extension) ]


# log function
def writeLog(output_path, logToWrite, firstTime = False):
    if not os.path.isfile(output_path + "log.txt"):
        f = open(output_path + "log.txt", "w")
    else:
        f = open(output_path + "log.txt", "a")
    now = datetime.now().strftime('%d-%m-%Y:%H:%M:%S')
    if firstTime:
        f.write(str(now ) + "\n" + logToWrite + "\n")
    else:
        f.write(logToWrite + "\n")
    f.close()


def getMythogenome(fastaFile):
    for record in SeqIO.parse(fastaFile, "fasta"):
        return (record.id, record.seq)

# function correctMinMaxInputError,
# This function will return only a number extracted from the string in param
# @param,
# @val, the string to parth
# @columnName, name of the column we want to check
# @mtName, name of the mythogenome we are treating
# @geneName, name of the column we are checking
def correctMinMaxInputError(val, columnName, mtName, geneName):
    if not re.match("^[0-9]*$", val):
        val = val.replace(",", "").replace(".", "")
        res = re.findall(r'\b\d+\b', val)[0]
        log = "A correction have been made  on " + mtName + " : " + geneName + " on the column : " + columnName + "\n" + val + " was corrected in : " + res
        print(log)
        writeLog("./logs/",log)
        return int(res)
    else:
        return int(val)

def writeRecords(listRecords, mtName, destinationPath = "./ressources/results/"):
    file = open(destinationPath + mtName + "_classique.fasta", "w")
    writer = SeqIO.FastaIO.FastaWriter(file)
    writer.write_file(listRecords)
    file.close()


    for record in listRecords:
        fileName = destinationPath  + record.name + ".fasta"
        if not os.path.isfile(fileName):
            file = open(fileName, "w")
        else:
            file = open(fileName, "a")
            
            
        writer = SeqIO.FastaIO.FastaWriter(file)
        writer.write_record(record)
        file.close()

        # writer = SeqIO.FastaIO.FastaWriter(file)
        # writer.write_file(rec)
        # file.close()



def extractSeq(csv, fasta):
    df = pd.read_csv(csv, sep=",")
    mythogenomeName, mythogenomeSeq = getMythogenome(fasta)
    Subseq = []
    records = []

    for i in range(len(df)):
        min = correctMinMaxInputError(str(df.Minimum[i]), "Minimum", mythogenomeName, df.Name[i])
        max = correctMinMaxInputError(str(df.Maximum[i]), "Maximum", mythogenomeName, df.Name[i])
        if min > 0 and max < len(mythogenomeSeq):
            Subseq.append(mythogenomeSeq[min : max])
            record = SeqRecord(mythogenomeSeq[min : max], id=mythogenomeName + "-" + df.Name[i], name=df.Name[i], description="")
            records.append(record)
        else:
            Subseq.append("Not Found")

    df = df.assign(Subsequence = Subseq)
    df = df.assign(FromMythogenome= mythogenomeName)
    writeRecords(records, mythogenomeName)
    
def run():
    fastaFiles = getFASTAFiles("./ressources/raw_data")
    csvFiles = getCSVFiles("./ressources/raw_data")
    couple = []
    for f in fastaFiles:
        for c in csvFiles:
            if f[:-5] in c:
                couple.append((c,f))
    for c, f in couple:
        extractSeq("./ressources/raw_data/" + c, "./ressources/raw_data/" + f)

