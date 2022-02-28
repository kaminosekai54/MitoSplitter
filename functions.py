# Import
from http.client import METHOD_NOT_ALLOWED
from setting import *
import os
import re
from datetime import datetime
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

################################################################################

# get the settings
settings = getSettings()


# function setup,
# this function will create all the folder need for the app to start if the don't 
# already exist
def setup():
    # checking if the logs folder exist and creat it if not
    if not os.path.isdir("./logs"):
        os.makedirs("./logs")

# checking if the ressources folder exist and creat it if not
    if not os.path.isdir("./ressources"):
        os.makedirs("./ressources")

# checking if the raw_data folder exist and creat it if not
    if not os.path.isdir("./ressources/raw_data"):
        os.makedirs("./ressources/raw_data")

# checking if the results folder exist and creat it if not
    if not os.path.isdir("./ressources/results"):
        os.makedirs("./ressources/results")

    # checking if the results/csv folder exist and creat it if not
    if not os.path.isdir("./ressources/results/csv"):
        os.makedirs("./ressources/results/csv")


# log function
def writeLog(logToWrite, output_path = settings["logPath"], firstTime = False):
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

# function getCSVFiles,
# this function return a list of all the csv files
# found in the indicated folder
# @param,
#@ path, the path of the folder you want to parth
# @ extension, extension of the file we search for, here "csv by default"
def getCSVFiles(path = settings["rawFilePath"], extension = ".csv"):
    fileList = os.listdir(path)
    return [ csvFile for csvFile in fileList if csvFile.endswith( extension) ]

# function getFASTAFiles,
# this function return a list of all the fasta files
# found in the indicated folder
# @param,
#@ path, the path of the folder you want to parth
# @ extension, extension of the file we search for, here "fasta by default"
def getFASTAFiles(path = settings["rawFilePath"], extension = ".fasta"):
    fileList = os.listdir(path)
    return [ fastaFile for fastaFile in fileList if fastaFile.endswith( extension) ]


#getMitogenome,
#this function return a tuple containing the id and seq of the mitogenome
# @param
#@fastafile, the fasta fle containing the mitogenome
def getMitogenome(fastaFile):
    for record in SeqIO.parse(fastaFile, "fasta"):
        return (record.id, record.seq)

# function correctMinMaxInputError,
# This function will return only a number extracted from the string in param
# @param,
# @val, the string to parth
# @columnName, name of the column we want to check
# @mtName, name of the mitogenome we are treating
# @geneName, name of the column we are checking
def correctMinMaxInputError(val, columnName, mtName, geneName):
    if not re.match("^[0-9]*$", val):
        val = val.replace(",", "").replace(".", "")
        res = re.findall(r'\b\d+\b', val)[0]
        log = "A correction have been made  on " + mtName + " : " + geneName + " on the column : " + columnName + "\n" + val + " was corrected in : " + res
        print(log)
        writeLog(log)
        return int(res)
    else:
        return int(val)

# function writeRecords,
# this function write a liste of record in the apropriate fasta file
# @param,
# @listRecords, the list of record to write 
# @mtName, the name of the mitogenome 
# @destinationPath, the path where to write the file 
def writeRecords(listRecords, mtName, destinationPath = settings["resultPath"]):
    # writing the classic fasta file
    file = open(destinationPath + mtName + "_classique.fasta", "w")
    writer = SeqIO.FastaIO.FastaWriter(file)
    writer.write_file(listRecords)
    file.close()

#  writing all different gene in the corresponding fasta file
    for record in listRecords:
        fileName = destinationPath  + record.name + ".fasta"
        if not os.path.isfile(fileName):
            file = open(fileName, "w")
        else:
            file = open(fileName, "a")
            
            
        writer = SeqIO.FastaIO.FastaWriter(file)
        writer.write_record(record)
        file.close()

# function extractSeq,
# this function will extract the subseq from the mitogenome
#  based on the anotation
# @param,
# @csv, the csv file
# @fasta, the fastafile
def extractSeq(csv, fasta):
    df = pd.read_csv(csv, sep=",")
    mitogenomeName, mitogenomeSeq = getMitogenome(fasta)
    Subseq = []
    records = []

    for i in range(len(df)):
        name = df[settings["nameColName"]][i]
        min = correctMinMaxInputError(str(df[settings["minColName"]][i]), settings["minColName"], mitogenomeName, name) -1
        max = correctMinMaxInputError(str(df[settings["maxColName"]][i]), settings["maxColName"], mitogenomeName, name)
        if min > 0 and max < len(mitogenomeSeq):
            Subseq.append(mitogenomeSeq[min : max])
            record = SeqRecord(mitogenomeSeq[min : max], id=mitogenomeName + "-" + name, name=name, description="")
            records.append(record)
        else:
            Subseq.append("Not Found")
            log = "The sequence : " + name + " from the mitogenome " + mitogenomeName + " from the file " + csv + "has not been found"
            print(log)
            writeLog(log)

    df = df.assign(Subsequence = Subseq)
    df = df.assign(FromMitogenome= mitogenomeName)
    df.to_csv(settings["csvResultPath"] + csv[csv.rfind("/")+1:].replace(".csv","") + "_with_sequence.csv", index= False)

# writting the fasta
    writeRecords(records, mitogenomeName)
    
# run,
# this function is the main function to run
def run():
    writeLog("Starting treatement", firstTime=True)
    setup()
    fastaFiles = getFASTAFiles(settings["rawFilePath"])
    csvFiles = getCSVFiles(settings["rawFilePath"])
    couple = []
    for f in fastaFiles:
        for c in csvFiles:
            if f[:-5] in c:
                couple.append((c,f))

    for c, f in couple:
        extractSeq(settings["rawFilePath"] + c, settings["rawFilePath"] + f)