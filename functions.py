# Import
from audioop import mul
from setting import *
import sys, os, platform, subprocess, re
from datetime import datetime
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import GenBank
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator


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
            if id !="" and id !="-1" and not "csv_" in id: return id

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
        multipleFile = -1
        mitogenomeName = fastaFile[fastaFile.rfind("/")+1:].replace(".fasta","")
        if "_" in mitogenomeName:
            multipleFile= "csv" + mitogenomeName[mitogenomeName.find("_"):]
            mitogenomeName= mitogenomeName[0:mitogenomeName.find("_")]

        if " " in mitogenomeName:
            mitogenomeName = mitogenomeName.replace(mitogenomeName[0], str.upper(mitogenomeNmae[0]), 1).replace(mitogenomeName[1:], str.lower(mitogenomeName[1:]), 1).replace(" ", "-")
        return (mitogenomeName, record.seq, record.id, multipleFile)

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
                    if id == "":
                        geneDict[record.name][i] = (record.id, seq, len(seq), record.description)

                    elif id == "-1":
                        geneDict[record.name][i] = (record.id, record.seq, len(record.seq), record.description)
                    elif "csv_" in id:
                        x = int(id[id.find("_")+1:])
                        y = int(record.description[record.description.find("_")+1:])
                        if y == x+1:
                            newSeq = seq + record.seq
                            modifySeqInFasta(fileName,mt, newSeq)
                            geneDict[record.name][i] = (record.id, newSeq, len(newSeq), id)
                            geneDict[record.name].append((record.id, newSeq, len(newSeq), record.description))



            file.close()


def modifySeqInFasta(fileName,mitogenomeName, newSeq):
    listRec = []
    for record in SeqIO.parse(fileName, "fasta"):
        if record.id == mitogenomeName:
            record.seq = newSeq

        listRec.append(record)

    with open(fileName, "w") as file:
        writer = SeqIO.FastaIO.FastaWriter(file)
        writer.write_file(listRec)


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
    df = df[df[settings["typeColName"]] != "source"]
    df[settings["nameColName"]] = df[settings["nameColName"]].str.replace("\s(.*)", "")
    df[settings["nameColName"]] = df[settings["nameColName"]].str.replace(")", "-")
    df[settings["nameColName"]] = df[settings["nameColName"]].str.replace("(", "-")
    df[settings["nameColName"]] = df[settings["nameColName"]].str.upper()

    mitogenomeName, mitogenomeSeq, mitogenomeId, fileNumber = getMitogenome(fasta)
    Subseq = []
    superposedSeq = []
    records = []
    lastMax = -1
    lastIndex = -1
    nbTreated = 0
    prevName =""
    nextName = ""
    listName = []

    for i in df.index:
        name = df[settings["nameColName"]][i]
        min = correctMinMaxInputError(str(df[settings["minColName"]][i]), settings["minColName"], mitogenomeName, name) -1
        max = correctMinMaxInputError(str(df[settings["maxColName"]][i]), settings["maxColName"], mitogenomeName, name)

        if i+1 < len(df.index) and i+1 in df.index :
            nextName = df[settings["nameColName"]][i+1]

        if min >= 0 and max <= len(mitogenomeSeq):
            if lastMax == -1 : lastMax = max-1

            # treatement for superposed position
            if min >= lastMax and nbTreated >0:
                superposedSeq.append(checkSuperposition(min, lastMax, mitogenomeName, df[settings["nameColName"]][lastIndex], name))

                #  common treatement
            if name == "TRNA-LEU":
                if prevName == "COX1" and nextName == "COX2": name = name +"1"
                elif prevName == "COX1" : name = name +"1"
                elif nextName == "COX2": name = name +"1"
                elif (prevName == "ND1" or prevName == "NAD1")and nextName == "16S" : name = name +"2"
                elif (prevName == "ND1" or prevName == "NAD1"): name = name +"2"
                elif nextName == "16S" : name = name +"2"

            if name == "TRNA-SER":
                if prevName == "TRNA-ASN" and nextName == "TRNA-GLU" : name = name +"1"
                elif prevName == "TRNA-ASN" : name = name +"1"
                elif nextName == "TRNA-GLU" : name = name +"1"
                elif prevName == "CYTB" and (nextName == "ND1" or nextName == "NAD1") : name = name +"2"
                elif prevName == "CYTB": name = name +"2"
                elif (nextName == "ND1" or nextName == "NAD1") : name = name +"2"
            Subseq.append(mitogenomeSeq[min : max])
            record = SeqRecord(mitogenomeSeq[min : max], id=mitogenomeName, name=name, description=str(fileNumber))
            records.append(record)
            listName.append(name)
            lastMax = max
            lastIndex=i
            nbTreated+=1
            prevName= name
        else:
            Subseq.append(pd.NA)
            log = "The sequence : " + name + " from the mitogenome " + mitogenomeName + " from the file " + csv + "has not been found"
            print(log)
            writeLog(log)

    df = df.assign(Subsequence = Subseq)
    df.to_csv(settings["csvResultPath"] + csv[csv.rfind("/")+1:].replace(".csv","") + "_with_sequence.csv", index= False)

    leu1Found = "TRNA-LEU1" in listName
    leu2Found = "TRNA-LEU2" in listName
    ser1Found= "TRNA-SER1" in listName
    ser2Found= "TRNA-SER2" in listName
    for record in records :
        if record.name == "TRNA-LEU" and leu1Found and not leu2Found: record.name = record.name +"2"
        if record.name == "TRNA-LEU" and not leu1Found and leu2Found: record.name = record.name +"1"
        if record.name == "TRNA-SER" and not ser1Found and ser2Found : record.name = record.name +"1"
        if record.name == "TRNA-SER" and ser1Found and not ser2Found : record.name = record.name +"2"



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
            listName = []

            if not accessionID in listAccession: listAccession.append(accessionID)
            name = ""
            prevType = ""
            prevEnd = -1
            nbTreated =0
            prevName = ""
            nextName = ""
            nextIndex =1
            for gene in record.features:
                if gene.type == "tRNA":
                    name = gene.qualifiers["product"][0]
                elif gene.type == "CDS":
                    name = gene.qualifiers["gene"][0]
                elif gene.type == "rRNA":
                    name = gene.qualifiers["product"][0][:gene.qualifiers["product"][0].find(" ")]

                if nextIndex < len(record.features):
                    if record.features[nextIndex].type == "tRNA":
                        nextName= str.upper(record.features[nextIndex].qualifiers["product"][0])
                    elif record.features[nextIndex].type == "CDS":
                        nextName= str.upper(record.features[nextIndex].qualifiers["gene"][0])
                    elif record.features[nextIndex].type == "rRNA":
                        nextName= str.upper(record.features[nextIndex].qualifiers["product"][0][:record.features[nextIndex].qualifiers["product"][0].find(" ")])
                if gene.type != "gene" and gene.type != "source" and gene.type != "misc_feature":
                    start = correctMinMaxInputError(str(gene.location.start),"Minimum",mitogenomeName,name)
                    end= correctMinMaxInputError(str(gene.location.end), "Maximum",mitogenomeName,name)
                    strand = gene.location.strand
                    seq = Seq(mitogenomeSeq [start:end])
                    # if strand == -1:
                        # print("need revers")
                        # seq = seq.reverse_complement()
                    name = str.upper(name)

                    if name == "TRNA-LEU":
                        if prevName == "COX1" and nextName == "COX2": name = name +"1"
                        elif prevName == "COX1" : name = name +"1"
                        elif nextName == "COX2": name = name +"1"
                        elif (prevName == "ND1" or prevName == "NAD1")and nextName == "16S" : name = name +"2"
                        elif (prevName == "ND1" or prevName == "NAD1"): name = name +"2"
                        elif nextName == "16S" : name = name +"2"

                    if name == "TRNA-SER":
                        if prevName == "TRNA-ASN" and nextName == "TRNA-GLU" : name = name +"1"
                        elif prevName == "TRNA-ASN" : name = name +"1"
                        elif nextName == "TRNA-GLU" : name = name +"1"
                        elif prevName == "CYTB" and (nextName == "ND1" or nextName == "NAD1") : name = name +"2"
                        elif prevName == "CYTB": name = name +"2"
                        elif (nextName == "ND1" or nextName == "NAD1") : name = name +"2"

                    record = SeqRecord(seq, id=mitogenomeName, name=name, description= accessionID)
                    listRecords.append(record)
                    listGene.append(name)

                    if nbTreated > 0 and prevType !=  "" and prevType != "source" and gene.type != "misc_feature":
                        checkSuperposition(start, prevEnd, mitogenomeName, listGene[-2], name)

                    nbTreated+=1
                    prevEnd=gene.location.end
                    prevType= gene.type
                    prevName= name
                    nextIndex+=1
                    listName.append(name)

    
    leu1Found = "TRNA-LEU1" in listName
    leu2Found = "TRNA-LEU2" in listName
    ser1Found= "TRNA-SER1" in listName
    ser2Found= "TRNA-SER2" in listName
    for record in listRecords:
        if record.name == "TRNA-LEU" and leu1Found and not leu2Found: record.name = record.name +"2"
        if record.name == "TRNA-LEU" and not leu1Found and leu2Found: record.name = record.name +"1"
        if record.name == "TRNA-SER" and not ser1Found and ser2Found : record.name = record.name +"1"
        if record.name == "TRNA-SER" and ser1Found and not ser2Found : record.name = record.name +"2"

    writeRecords(listRecords, mitogenomeName)

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
def aligneSequence(fasta, useMuscle =settings["useMuscle"], useMafft=settings["useMafft"],  outputLocation = settings["sequenceAlignementResultPath"], muscleLocation = settings["musclePath"], mafftLocation = settings["mafftPath"]  ):
    osName = platform.system()
    outputFile =""
    if useMuscle:
        log = "Processing Muscle alignement for " + fasta
        print(log)
        outputFile = outputLocation + fasta[fasta.rfind("/")+1:-6]+ "_muscle_align.phy"
        tmpFile= outputLocation + fasta[fasta.rfind("/")+1:-6]+ "_muscle_align.fasta"

        muscleEXE = ""
        if osName == "Linux":
            muscleEXE= muscleLocation  + "muscle5.1.linux_intel64"
            subprocess.Popen("chmod +x " + muscleEXE, stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32")).wait()
        elif osName == "Windows":
            muscleEXE= muscleLocation+ "muscle5.1.win64.exe"
        elif osName == "Darwin":
            muscleEXE=muscleLocation+"muscle5.1.macos_intel64"
            subprocess.Popen("chmod +x " + muscleEXE, stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32")).wait()

        muscle_cline= muscleEXE + " -align " + os.path.normpath(fasta) + " -output " + os.path.normpath(tmpFile)
        child= subprocess.Popen(str(muscle_cline), stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32"))
        child.wait()
        AlignIO.convert(tmpFile, "fasta", outputFile, "phylip-relaxed")
        os.remove(tmpFile)


    # output try with mafft
    if useMafft:
        log = "Processing Mafft alignement for " + fasta
        print(log)
        outputFile = outputLocation + fasta[fasta.rfind("/")+1:-6]+ "_mafft_align.phy"
        tmpFile= outputLocation + fasta[fasta.rfind("/")+1:-6] + "_mafft_align.fasta"
        mafftEXE =""

        if osName == "Linux":
            mafftEXE= mafftLocation  + "muscle5.1.linux_intel64"
            subprocess.Popen("chmod +x " + mafftEXE, stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32")).wait()
        elif osName == "Windows":
            mafftEXE= mafftLocation + "mafft-win/mafft.bat"
            mafft_cline= "cmd.exe /C " + os.path.abspath(mafftEXE) + " --auto --out "+ os.path.normpath(tmpFile) + " " + os.path.normpath(fasta)
        elif osName == "Darwin":
            mafftEXE= mafftLocation +"mafft-mac/mafft.bat"
            subprocess.Popen("chmod +x " + mafftEXE, stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32")).wait()
            mafft_cline= os.path.abspath(mafftEXE) + " --auto --out "+ os.path.normpath(tmpFile) + " " + os.path.normpath(fasta)

        print("Commande : \n" + mafft_cline)
        child= subprocess.Popen(str(mafft_cline), stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32"))
        child.wait()
        AlignIO.convert(tmpFile, "fasta", outputFile, "phylip-relaxed")
        os.remove(tmpFile)
        
    return outputFile


def writeSingleRecordList(tmpFasta, listRec):
    with open(tmpFasta, "w") as file:
        writer = SeqIO.FastaIO.FastaWriter(file)
        writer.write_file(listRec)


def getReversedRecordList(fastaFile, taxonName):
    listRec = []
    for record in SeqIO.parse(fastaFile, "fasta"):
                if record.id == taxonName:
                    rec = SeqRecord(record.seq.reverse_complement(), id=record.id, name=record.name, description="")
                else:
                    rec = record
                listRec.append(rec)

    return listRec

def getPDistMatrix(alignementFile):
    aln = list(AlignIO.parse(alignementFile, "phylip-relaxed"))
    calculator = DistanceCalculator('identity')
    return calculator.get_distance(aln[0])


def checkAlignement(alignementFile, alignementPath= settings ["sequenceAlignementResultPath"], pathToFasta = settings["genesFastaResultPath"]):
    print("alignement check initialised")
    dm = getPDistMatrix(alignementFile)
    fastaFile = pathToFasta + alignementFile[alignementFile.rfind("/")+1:alignementFile.find("_")] + ".fasta"
    tmpFasta = fastaFile[:-6] + "_tmp.fasta"
    listRec= []

    for i in range(len(dm.matrix)):
        taxonName= dm.names[i]
        prevPDist = dm.matrix[i][0]
        if dm.matrix[i][0] > 0.5:
            # print(dm.matrix[i][0])
            listRec = getReversedRecordList(fastaFile, taxonName)
            writeSingleRecordList(tmpFasta, listRec)
            tmpAlignement = aligneSequence(tmpFasta, False, True)
            dm2 = getPDistMatrix(tmpAlignement)
            for j in range(len(dm2.matrix)):
                if dm2.names[j] == taxonName:
                    if dm2.matrix[j][0] < prevPDist:
                        os.remove(fastaFile)
                        os.remove(alignementFile)
                        os.rename(tmpFasta, fastaFile)
                        os.rename(tmpAlignement, alignementFile)
                        return checkAlignement(alignementFile)
                    else:
                        os.remove(tmpFasta)
                        os.remove(tmpAlignement)


    print("alignementCheck finished")
    return



    
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
        alignedFile = aligneSequence(settings ["genesFastaResultPath"] + fasta, useMuscle=settings ["useMuscle"], useMafft = settings ["useMafft"] )
        checkAlignement(alignedFile)

################################################################################
# initialisation of global variable and environement
setup()
geneDict = getGeneDict()
superpositionDict={}