# Import
from setting import *
import sys, os, platform, subprocess, re
from datetime import datetime
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import GenBank
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO


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

# checking if the raw_data folder exist and creat it if not
    if not os.path.isdir("./raw_data"):
        os.makedirs("./raw_data")

# checking if the results folder exist and creat it if not
    if not os.path.isdir("./results"):
        os.makedirs("./results")

    # checking if the results/csv folder exist and creat it if not
    if not os.path.isdir("./results/csv"):
        os.makedirs("./results/csv")

        # checking if the results/genes fasta folder exist and creat it if not
    if not os.path.isdir("./results/genes-fasta"):
        os.makedirs("./results/genes-fasta")

        # checking if the results/classic fasta folder exist and creat it if not
    if not os.path.isdir("./results/classic-fasta"):
        os.makedirs("./results/classic-fasta")

    if not os.path.isdir("./results/alignement"):
        os.makedirs("./results/alignement")

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
            geneDict[geneName].append((record.id, record.seq, len(record.seq), ""))

    return geneDict

def isGenomeInSuperpositionDict(gene, valToCheck):
    for taxon, length in superpositionDict[gene]:
        if taxon == valToCheck: return True

    return False

def getLengthInSuperpositionDict(gene, mitogenomeName):
    for taxon, length  in superpositionDict[gene]:
        if taxon == mitogenomeName: return length

    return pd.NA



def isGenomeInGeneDict(gene, valToCheck):
    for taxon, seq, length, id  in geneDict[gene]:
        if taxon == valToCheck:
            return True

    return False

def getAccesIDInGeneDict(gene, mitogenomeName):
    for taxon, seq, length, id in geneDict[gene]:
        if taxon == mitogenomeName:
            if id !="" : return id

    return pd.NA

def getLengthInGeneDict(gene, mitogenomeName):
    for taxon, seq, length, id  in geneDict[gene]:
        if taxon == mitogenomeName:
            return length

    return pd.NA

#getMitogenome,
#this function return a tuple containing the id and seq of the mitogenome
# @param
#@fastafile, the fasta fle containing the mitogenome
def getMitogenome(fastaFile):
    for record in SeqIO.parse(fastaFile, "fasta"):
        mitogenomeName = fastaFile[fastaFile.rfind("/")+1:].replace(".fasta","")
        if "_" in mitogenomeName:
            mitogenomeName= mitogenomeName[0:mitogenomeName.find("_")]

        if " " in mitogenomeName:
            mitogenomeName = mitogenomeName.replace(mitogenomeName[0], str.upper(mitogenomeNmae[0]), 1).replace(mitogenomeName[1:], str.lower(mitogenomeName[1:]), 1).replace(" ", "-")
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
    log = "Writing classic fasta file for " + mtName
    print(log)
    writeLog(log)
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
        
            accessId=record.description
            record.description = ""
            writer = SeqIO.FastaIO.FastaWriter(file)
            writer.write_record(record)
            geneDict[record.name].append((record.id, record.seq, len(record.seq), accessId))

        else:
            for i in range(len(geneDict[record.name])):
                mt, seq, length, id = geneDict[record.name][i]
                if mt == record.id:
                    geneDict[record.name][i] = (record.id, record.seq, len(record.seq), record.description)

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
        
    # log = "Superposition found  for " + mitogenomeName + " between " + prevName + " and " + name + " on the following position : \n" + str(superposedPosition)
    log = str(len(superposedPosition)) + " Superposition found  for " + mitogenomeName + " between " + prevName + " and " + name 
    writeLog(log)
    print(log)

    superpositionName = prevName +"/"+ name
    if not superpositionName in superpositionDict.keys():
        superpositionDict[superpositionName] = []
    if not isGenomeInSuperpositionDict(superpositionName, mitogenomeName) : superpositionDict[superpositionName].append((mitogenomeName, len(superposedPosition)))
    return (prevName, name, superposedPosition)

# function extractSeqFromCSV,
# this function will extract the subseq from the mitogenome fasta file paired with a csv anotation file
#  based on the anotation
# @param,
# @csv, the csv file
# @fasta, the fastafile
def extractSeqFromCSV(csv, fasta):
    log = "Starting extraction for " + csv
    print(log)
    writeLog(log)
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
        mitogenomeName = record.organism
        if " " in mitogenomeName:
            mitogenomeName = mitogenomeName.replace(mitogenomeName[0], str.upper(mitogenomeName[0]), 1).replace(mitogenomeName[1:], str.lower(mitogenomeName[1:]), 1).replace(" ", "-")

        return (mitogenomeName, record.sequence, record.accession[0])


#extractSeqFromGBFile
# this function will extract the mitogenome from a genbank (.gb) file
# it will also rename it if need it
# @param
# @gbFile, the genbank file to deal with
def extractSeqFromGBFile(gbFile):
    log = "Starting extraction for " + gbFile
    print(log)
    writeLog(log)
    # mitogenomeName, mitogenomeSeq, accessionID = getMitogenomeFromGBFile(gbFile)
    listGene = []
    listRecords = []
    listAccession = []
    needRename=False
    mitogenomeName = ""

    with open(gbFile) as gb:
        for record in SeqIO.parse(gb, "genbank"):
            mitogenomeName = record.features[0].qualifiers["organism"][0]
            if " " in mitogenomeName:
                mitogenomeName = mitogenomeName.replace(mitogenomeName[0], str.upper(mitogenomeName[0]), 1).replace(mitogenomeName[1:], str.lower(mitogenomeName[1:]), 1).replace(" ", "-")

            mitogenomeSeq = record.seq
            accessionID = record.id
            needRename=True

            if not accessionID in listAccession: listAccession.append(accessionID)
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
                    name = str.upper(name)
                    record = SeqRecord(Seq(seq), id=mitogenomeName, name=name, description= accessionID)
                    listRecords.append(record)
                    listGene.append(name)

                    if nbTreated > 0 and prevType !=  "" and prevType != "source" and gene.type != "misc_feature":
                        checkSuperposition(start, prevEnd, mitogenomeName, listGene[-2], name)

                    nbTreated+=1
                    prevEnd=gene.location.end
                    prevType= gene.type
                
    writeRecords(listRecords,mitogenomeName)

    # check for rename
    if gbFile[gbFile.rfind("/")+1:-3] != mitogenomeName and needRename:
        newName = gbFile[:gbFile.rfind("/")+1] + mitogenomeName+ ".gb"
        if os.path.isfile(newName):
            i = 2
            while os.path.isfile(newName):
                newName = gbFile[:gbFile.rfind("/")+1] + mitogenomeName+ "_" + str(i) + ".gb"
                i+=1

        os.rename(gbFile, newName)

    return (mitogenomeName, listAccession)

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
    data={"Taxon":mitogenomeDict.keys()}
    for genome in mitogenomeDict.keys():
        for gene in geneDict.keys():
            if not gene in data.keys():
                data[gene] = []
            
            data[gene].append(getAccesIDInGeneDict(gene, genome))
    df = pd.DataFrame.from_dict(data)
    df.to_csv(csvPath+ "summary_accessionID.csv", index=False)


def generateSuperpositionSummary(mitogenomeDict, csvPath= settings["csvResultPath"]):
    data={"Taxon":mitogenomeDict.keys()}
    for genome in mitogenomeDict.keys():
        for gene in superpositionDict.keys():
            if not gene in data.keys():
                data[gene] = []
            
            data[gene].append(getLengthInSuperpositionDict(gene, genome))

    df = pd.DataFrame.from_dict(data)
    df.to_csv(csvPath+ "summary_Superposition.csv", index=False)


################################################################################
# Sequence alignement
def aligneSequence(fasta, outputLocation = settings["sequenceAlignementResultPath"], muscleLocation = settings["musclePath"] ):
    osName = platform.system()
    tmpFile = outputLocation + "tmp.aln"
    outputFile = outputLocation + fasta[fasta.rfind("/")+1:-5]+ "_align.phy"
    muscleEXE = ""
    if osName == "Linux":
        muscleEXE= muscleLocation  + "muscle5.1.linux_intel64"
    elif osName == "Windows":
        muscleEXE= muscleLocation+ "muscle5.1.win64.exe"
    elif osName == "Darwin":
        muscleEXE=muscleLocation+"muscle5.1.macos_intel64"

    muscle_cline= muscleEXE + " -align " + fasta + " -output " + outputFile
    child= subprocess.Popen(str(muscle_cline), stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32"))
    child.wait()
    # print(child.stdout)
    # print(child.stdin)
    # print(child.stderr)
    # print(muscle_cline)


    # aligns= []
    # with open(outputFile) as align_handle:
        # aligns= AlignIO.pa(align_handle,"fasta")

    # outfile=open(outputFile , "w")
    # AlignIO.write([align], outfile, 'phylip')
    # outfile.close()
    # os.remove(tmpFile)




    
# run,
# this function is the main function to run
def run():
    setup()
    writeLog("Starting treatement", firstTime=True)
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
        if mitogenomeName not in mitogenomeDict.keys(): 
            mitogenomeDict[mitogenomeName] = [accessionID]
        else:
            for acces in accessionID:
                if not acces in mitogenomeDict[mitogenomeName]: mitogenomeDict[mitogenomeName].append(acces)

    generatePresenceSummary(mitogenomeDict)
    generateLengthSummary(mitogenomeDict)
    generateAccessionIDSummary(mitogenomeDict)
    generateSuperpositionSummary(mitogenomeDict)
    print("Finish")
    writeLog("Finish")

    for fasta in getFASTAFiles(path=settings ["genesFastaResultPath"]):
        aligneSequence(settings ["genesFastaResultPath"] + fasta)

################################################################################
# initialisation of global variable and environement
setup()
geneDict = getGeneDict()
superpositionDict={}