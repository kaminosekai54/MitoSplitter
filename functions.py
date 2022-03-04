# Import
from setting import *
import os
import re
from datetime import datetime
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import GenBank


################################################################################

# get the settings and already computed data
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

        # checking if the results/genes fasta folder exist and creat it if not
    if not os.path.isdir("./ressources/results/genes fasta"):
        os.makedirs("./ressources/results/genes fasta")

        # checking if the results/classic fasta folder exist and creat it if not
    if not os.path.isdir("./ressources/results/classic fasta"):
        os.makedirs("./ressources/results/classic fasta")


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

# getGBFiles
# this function return a list of all the .gb files
# found in the indicated folder
# @param,
#@ path, the path of the folder you want to parth
# @ extension, extension of the file we search for, here ".gb" by default"
def getGBFiles(path = settings["rawFilePath"], extension = ".gb"):
    fileList = os.listdir(path)
    return [ gbFile for gbFile in fileList if gbFile.endswith( extension) ]

# getGeneDict,
# this function return a dictionary where
#  the keys are the genes  
# and the value are
# the mitogenome that have this gene
# @param
# @path, the path where to look for the gene fasta file
def getGeneDict(path = settings["genesFastaResultPath"]):
    geneDict = {}
    fileList = os.listdir(path)
    for file in fileList:
        geneName = file[file.rfind("/")+1:].replace(".fasta","")
        if not geneName in geneDict.keys():
            geneDict[geneName] = []

        for record in SeqIO.parse(path+file, "fasta"):
            geneDict[geneName].append((record.id, record.seq, len(record.seq)))

    return geneDict

def isGenomeInGeneDict(gene, valToCheck):
    for taxon, seq, length in geneDict[gene]:
        if taxon == valToCheck:
            return True

    return False

def getLengthInGeneDict(gene, mitogenomeName):
    for taxon, seq, length in geneDict[gene]:
        if taxon == mitogenomeName:
            return length

    return pd.NA

#getMitogenome,
#this function return a tuple containing the id and seq of the mitogenome
# @param
#@fastafile, the fasta fle containing the mitogenome
def getMitogenome(fastaFile):
    for record in SeqIO.parse(fastaFile, "fasta"):
        mitogenomeName = fastaFile[fastaFile.rfind("/")+1:].replace(".fasta","").replace(" ", "")
        if "_" in mitogenomeName:
            mitogenomeName= mitogenomeName[0:mitogenomeName.find("_")]
        return (mitogenomeName, record.seq, record.id)

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
# @destinationPath, the path where to write the classic file
# @destinationPath2, the path where to write the gene fasta file

def writeRecords(listRecords, mtName, destinationPath = settings["classicFastaResultPath"], destinationPath2 = settings["genesFastaResultPath"]):
    # writing the classic fasta file
    with open(destinationPath + mtName + "_classique.fasta", "w") as file:
        writer = SeqIO.FastaIO.FastaWriter(file)
        writer.write_file(listRecords)


#  writing all different gene in the corresponding fasta file
    for record in listRecords:
        fileName = destinationPath2  + record.name + ".fasta"
        if not record.name in geneDict.keys():
            geneDict[record.name] = []
        if not isGenomeInGeneDict(record.name, record.id) :
            if not os.path.isfile(fileName):
                file = open(fileName, "w")
            else:
                file = open(fileName, "a")
        
            writer = SeqIO.FastaIO.FastaWriter(file)
            writer.write_record(record)
            geneDict[record.name].append((record.id, record.seq, len(record.seq)))

            file.close()
#checkSuperposition
# this function check if their is a overlaping betwing two sequence position
#@param
# @start, the start position of the sequence
# @prevEnd, the previous end position of the previous sequence
# @mitogenomeName, name of the mitogenomeName, only used ofr log
# @prevName, the name of the previous sequence, only used for logs
# @name, the name of the current sequence, only used for logs
def checkSuperposition(start, prevEnd, mitogenomeName, prevName, name):
    superposedPosition = []
    dif = start - prevEnd
    for i in range(dif):
        superposedPosition.append(prevEnd+i)
    if dif == 0:
        superposedPosition.append(prevEnd)
        
    log = "Superposition found  for " + mitogenomeName + " between " + prevName + " and " + name + " on the following position : \n" + str(superposedPosition)
    writeLog(log)
    print(log)
    return (prevName, name, superposedPosition)

# function extractSeqFromCSV,
# this function will extract the subseq from the mitogenome fasta file paired with a csv anotation file
#  based on the anotation
# @param,
# @csv, the csv file
# @fasta, the fastafile
def extractSeqFromCSV(csv, fasta):
    df = pd.read_csv(csv, sep=",")
    df = df[df[settings["typeColName"]] != "gene"]
    df[settings["nameColName"]] = df[settings["nameColName"]].str.replace("\s(.*)", "")

    mitogenomeName, mitogenomeSeq, mitogenomeId = getMitogenome(fasta)
    Subseq = []
    superposedSeq = []
    records = []
    lastMax = -1
    lastIndex = -1
    nbTreated = 0

    for i in df.index:
        name = df[settings["nameColName"]][i]
        min = correctMinMaxInputError(str(df[settings["minColName"]][i]), settings["minColName"], mitogenomeName, name) -1
        max = correctMinMaxInputError(str(df[settings["maxColName"]][i]), settings["maxColName"], mitogenomeName, name)
        if min > 0 and max < len(mitogenomeSeq):
            if lastMax == -1 : lastMax = max-1

            # treatement for superposed position
            if min >= lastMax and nbTreated >0:
                superposedSeq.append(checkSuperposition(min, lastMax, mitogenomeName, df[settings["nameColName"]][lastIndex], name))

                #  common treatement
            Subseq.append(mitogenomeSeq[min : max])
            record = SeqRecord(mitogenomeSeq[min : max], id=mitogenomeName, name=name, description="")
            records.append(record)
            lastMax = max
            lastIndex=i
            nbTreated+=1
        else:
            Subseq.append(pd.NA)
            log = "The sequence : " + name + " from the mitogenome " + mitogenomeName + " from the file " + csv + "has not been found"
            print(log)
            writeLog(log)

    df = df.assign(Subsequence = Subseq)
    df.to_csv(settings["csvResultPath"] + csv[csv.rfind("/")+1:].replace(".csv","") + "_with_sequence.csv", index= False)

# writting the fasta
    writeRecords(records, mitogenomeName)
    return (mitogenomeName, pd.NA)

#getMitogenomeFromGBFile
# this function return the mitogenome name, accessionID, and seq from a .gb file
# @param
# @gbFile, the path to the genbankFile
def getMitogenomeFromGBFile(gbFile):
    with open(gbFile) as gb:
        record = GenBank.read(gb)
        return (record.organism.replace(" ", ""), record.sequence, record.accession[0])


#extractSeqFromGBFile
# this function will extract the mitogenome from a genbank (.gb) file
# it will also rename it if need it
# @param
# @gbFile, the genbank file to deal with
def extractSeqFromGBFile(gbFile):
    mitogenomeName, mitogenomeSeq, accessionID = getMitogenomeFromGBFile(gbFile)
    listGene = []
    listRecords = []

    with open(gbFile) as gb:
        for record in SeqIO.parse(gb, "genbank"):
            name = ""
            prevType = ""
            prevEnd = -1
            nbTreated =0
            for gene in record.features:
                if gene.type == "tRNA":
                    name = gene.qualifiers["product"][0]
                elif gene.type == "CDS":
                    name = gene.qualifiers["gene"][0]
                elif gene.type == "rRNA":
                    name = gene.qualifiers["product"][0][:gene.qualifiers["product"][0].find(" ")]

                
                if gene.type != "gene" and gene.type != "source" and gene.type != "misc_feature":
                    start = correctMinMaxInputError(str(gene.location.start),"Minimum",mitogenomeName,name)
                    end= correctMinMaxInputError(str(gene.location.end), "Maximum",mitogenomeName,name)
                    seq = mitogenomeSeq [start:end]
                    record = SeqRecord(Seq(seq), id=mitogenomeName, name=name, description="")
                    listRecords.append(record)
                    listGene.append(name)

                    if nbTreated > 0 and prevType !=  "" and prevType != "source" and gene.type != "misc_feature":
                        checkSuperposition(start, prevEnd, mitogenomeName, listGene[-2], name)

                    nbTreated+=1
                    prevEnd=gene.location.end
                    prevType= gene.type
                
    writeRecords(listRecords,mitogenomeName)

    if gbFile[gbFile.rfind("/")+1:-3] != mitogenomeName:
        os.rename(gbFile, gbFile[:gbFile.rfind("/")+1] + mitogenomeName+ ".gb")

    return (mitogenomeName, accessionID )

def generatePresenceSummary(mitogenomeDict, csvPath= settings["csvResultPath"]):
    data={"Taxon":mitogenomeDict.keys()}
    for genome in mitogenomeDict.keys():
        for gene in geneDict.keys():
            if not gene in data.keys():
                data[gene] = []
            
            data[gene].append(int(isGenomeInGeneDict(gene, genome)))

    df = pd.DataFrame.from_dict(data)
    df.to_csv(csvPath+ "summary_presence.csv", index=False)

def generateLengthSummary(mitogenomeDict, csvPath= settings["csvResultPath"]):
    data={"Taxon":mitogenomeDict.keys()}
    for genome in mitogenomeDict.keys():
        for gene in geneDict.keys():
            if not gene in data.keys():
                data[gene] = []
            
            data[gene].append(getLengthInGeneDict(gene, genome))

    df = pd.DataFrame.from_dict(data)
    df.to_csv(csvPath+ "summary_length.csv", index=False)



def generateAccessionIDSummary(mitogenomeDict, csvPath= settings["csvResultPath"]):
    data={"Taxon":mitogenomeDict.keys(), "AccessionID":mitogenomeDict.values()}
    df = pd.DataFrame.from_dict(data)
    df.to_csv(csvPath+ "summary_accessionID.csv", index=False)

# run,
# this function is the main function to run
def run():
    
    writeLog("Starting treatement", firstTime=True)
    setup()
    mitogenomeDict = {}
    fastaFiles = getFASTAFiles()
    csvFiles = getCSVFiles()
    gbFiles = getGBFiles()
    couple = []
    for f in fastaFiles:
        for c in csvFiles:
            if f[:-5] in c:
                couple.append((c,f))
    
    for c, f in couple:
        mitogenomeName, accessionID = extractSeqFromCSV(settings["rawFilePath"] + c, settings["rawFilePath"] + f)
        if mitogenomeName not in mitogenomeDict.keys(): mitogenomeDict[mitogenomeName] = accessionID
    for file in gbFiles:
        mitogenomeName, accessionID = extractSeqFromGBFile(settings["rawFilePath"] +file)
        if mitogenomeName not in mitogenomeDict.keys(): mitogenomeDict[mitogenomeName] = accessionID
    
    generatePresenceSummary(mitogenomeDict)
    generateLengthSummary(mitogenomeDict)
    generateAccessionIDSummary(mitogenomeDict)

################################################################################
# initialisation of global variable
geneDict = getGeneDict()