# Import
from setting import *
import sys, os, platform, subprocess, re, time
from datetime import datetime
from colorama import Fore, Back, Style, init
import numpy as np
import pandas as pd
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import GenBank
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment


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

    # checking if the results/alignement folder exist and creat it if not
    if not os.path.isdir("./results/alignement"):
        os.makedirs("./results/alignement")

    # checking if the results/tree folder exist and creat it if not
    if not os.path.isdir("./results/tree"):
        os.makedirs("./results/tree")


# utility function to print in color
def prRed(skk): print("\033[91m {}\033[00m" .format(skk))
def prGreen(skk): print("\033[92m {}\033[00m" .format(skk))
def prYellow(skk): print("\033[93m {}\033[00m" .format(skk))
def prLightPurple(skk): print("\033[94m {}\033[00m" .format(skk))
def prPurple(skk): print("\033[95m {}\033[00m" .format(skk))
def prCyan(skk): print("\033[96m {}\033[00m" .format(skk))
def prLightGray(skk): print("\033[97m {}\033[00m" .format(skk))
def prBlack(skk): print("\033[98m {}\033[00m" .format(skk))

# log function
def writeLog(logToWrite, output_path = settings["logPath"], firstTime = False):
    if not os.path.isfile(output_path + "log.txt"):
        f = open(output_path + "log.txt", "w")
    else:
        if firstTime: f = open(output_path + "log.txt", "w")
        else: f = open(output_path + "log.txt", "a")
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
# getCSVDelimiter,
#  This function return the delimiter of a csv File
# @param
# @csvFile, the csv file
def getCSVDelimiter(csvFile):
    with open(csvFile, 'r') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.readline())
        return str(dialect.delimiter)


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



# function isGenomeInSuperpositionDict
# This function return a bollean value, checking if the genome existe in the dictionary
# @param
# @gene, the gene corresponding to the key
# @mitogenomeName, name of the taxon we want to check the presence
def isGenomeInSuperpositionDict(gene, mitogenomeName):
    for taxon, length in superpositionDict[gene]:
        if taxon == mitogenomeName: return True

    return False

# function getLengthInSuperpositionDict
# This function return the length of a sequence, for a given gene and taxonName
# Or NA, is it's not found
# @param
# @gene, the gene corresponding to the key
# @mitogenomeName, name of the taxon we want to get the length
def getLengthInSuperpositionDict(gene, mitogenomeName):
    for taxon, length  in superpositionDict[gene]:
        if taxon == mitogenomeName: return length

    return pd.NA



# function isGenomeInGeneDict
# This function return a bollean value, checking if the genome existe in the dictionary
# @param
# @gene, the gene corresponding to the key
# @taxonName, name of the taxon we want to check the presence
def isGenomeInGeneDict(gene, taxonName):
    for taxon, seq, length, id  in geneDict[gene]:
        if taxon == taxonName: return True

    return False

# function getAccesIDInGeneDict
# This function return the accessID of a sequence, for a given gene and taxonName
# Or NA, is it's not found
# @param
# @gene, the gene corresponding to the key
# @mitogenomeName, name of the taxon we want to get the accessID
def getAccesIDInGeneDict(gene, mitogenomeName):
    for taxon, seq, length, id in geneDict[gene]:
        if taxon == mitogenomeName:
            if id !="" and id !="-1" and not "csv_" in id: return id

    return pd.NA

# function getLengthInGeneDict
# This function return the length of a sequence, for a given gene and taxonName
# Or NA, is it's not found
# @param
# @gene, the gene corresponding to the key
# @mitogenomeName, name of the taxon we want to get the length
def getLengthInGeneDict(gene, mitogenomeName):
    for taxon, seq, length, id  in geneDict[gene]:
        if taxon == mitogenomeName:
            return length

    return pd.NA




#getMitogenome,
#this function return a tuple containing the id and seq of the mitogenome for a fasta file coupled with a csv file
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
# @param
# @val, the string to parth
# @columnName, name of the column we want to check
# @mtName, name of the mitogenome we are treating
# @geneName, name of the column we are checking
def correctMinMaxInputError(val, columnName, mtName, geneName):
    if not re.match("^[0-9]*$", val):
        val = val.replace(",", "").replace(".", "")
        res = re.findall(r'\b\d+\b', val)[0]
        log = "A correction have been made  on " + mtName + " : " + geneName + " on the column : " + columnName + "\n" + val + " was corrected in : " + res
        # print(log)
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
        seqToModify={}
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
                        x = int(id[id.rfind("_")+1:])
                        y = int(record.description[record.description.rfind("_")+1:])
                        # print(x)
                        # print(y)
                        if y == x+1:
                            newSeq = seq + record.seq
                            seqToModify[mt]= newSeq
                            geneDict[record.name][i] = (record.id, newSeq, len(newSeq), id[0: id.rfind("_")+1]+str(y))
                            # geneDict[record.name].append((record.id, newSeq, len(newSeq), record.description))
                            # print(geneDict[record.name])



            file.close()
            modifySeqInFasta(fileName, seqToModify)

# function modifySeqInFasta
# This function will modify a sequence for a given mitogenome in a given file and 
#  rewrite the file, overwriting the previous one
# @param
# @fileName, name of the file to read and write 
# @mitogenomeName, name of the mitogenome, equal to the seqID, used to find the correct seq to modify
# @newSeq, the new sequence to write
def modifySeqInFasta(fileName, seqToModify):
    listRec = []
    for record in SeqIO.parse(fileName, "fasta"):
        if record.id in seqToModify.keys():
            record.seq = seqToModify[record.id]

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
    # print(log)

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
def extractSeqFromCSV(csv, fasta, destinationPath = settings["classicFastaResultPath"], destinationPath2 = settings["genesFastaResultPath"], csvPath= settings["csvResultPath"]):
    log = "Starting extraction for " + csv
    print(log)
    writeLog(log)
    separator = getCSVDelimiter(csv)
    df = pd.read_csv(csv, sep=separator)
    df = df[df[settings["typeColName"]] != "gene"]
    df = df[df[settings["typeColName"]] != "source"]
    df[settings["nameColName"]] = df[settings["nameColName"]].str.replace("\s(.*)", "",  regex=True)
    df[settings["nameColName"]] = df[settings["nameColName"]].str.replace(r")", "-", regex=False)
    df[settings["nameColName"]] = df[settings["nameColName"]].str.replace(r"(", "-", regex=False)
    df[settings["nameColName"]] = df[settings["nameColName"]].str.upper()
    for realName, aliasList in settings["geneAlias"].items():
            for alias in aliasList: 
                if alias in df[settings["nameColName"]] :     df[settings["nameColName"]] = df[settings["nameColName"]].str.replace(alias, realName)

    colToRemove = []

    for gene in df.Name:
        if not gene in settings["geneToDetect"] : colToRemove.append(gene)

    for gene in colToRemove:
        df = df[df[settings["nameColName"]] != gene]


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
    geneCounter= {}
    for gene in settings["geneToDetect"]:
        geneCounter[gene] = 0

    for i in df.index:
        name = df[settings["nameColName"]][i]
        min = correctMinMaxInputError(str(df[settings["minColName"]][i]), settings["minColName"], mitogenomeName, name) -1
        max = correctMinMaxInputError(str(df[settings["maxColName"]][i]), settings["maxColName"], mitogenomeName, name)

        if name in settings["geneToDetect"]:

            if i+1 in df.index :
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
                if name in geneCounter.keys():
                    geneCounter[name] +=1
                    if geneCounter[name]  > 1 :
                        if any(char.isdigit() for char in name):name = name + "-" + str(geneCounter[name]) 
                        else:name = name + str(geneCounter[name])
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

    # print(len(df))
    # print(len(Subseq))
    # for seq in df.Name:
        # if seq not in listName: print(seq)
    # print(df)
    # print(Subseq)
    df = df.assign(Subsequence = Subseq)
    df.to_csv(csvPath + csv[csv.rfind("/")+1:].replace(".csv","") + "_with_sequence.csv", index= False)

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
    writeRecords(records, mitogenomeName, destinationPath , destinationPath2 )
    return (mitogenomeName, pd.NA)


# function extractSeqFromSingleFasta
# This function will get the sequence from single fasta (Not paired with a csv)
# and add them to the global gene-fasta file to witch they correspond
# @param
# @fasta, the path to the fasta file to parth
# @destinationPath, path to the destion where to write the new fasta or to find it if it exist
def extractSeqFromSingleFasta(fasta, destinationPath = settings["classicFastaResultPath"], destinationPath2 = settings["genesFastaResultPath"]):
    log = "Starting extraction from " + fasta
    print(log)
    writeLog(log)
    name = str.upper(fasta[fasta.rfind("/")+1:].replace(".fasta", "").replace(" ", "-"))
    listMitogenome = []
    if "_" in name : name = name[:name.find("_")]
    name = getAlias(name)
    if name in settings["geneToDetect"]:
        if not name in geneDict .keys() :  geneDict[name] = []
        

        # for classic fasta
        if not os.path.isfile(destinationPath + name + ".fasta"):
            SeqIO.write(SeqIO.parse(fasta, "fasta"), destinationPath + name + ".fasta", "fasta")

        # For gene-fasta
        if not os.path.isfile(destinationPath2 + name + ".fasta"):
            # print("file don't exist")
            SeqIO.write(SeqIO.parse(fasta, "fasta"), destinationPath2 + name + ".fasta", "fasta")
            
        
            for record in SeqIO.parse(destinationPath2 + name + ".fasta", "fasta"):
                listMitogenome.append(record.id)
                if not isGenomeInGeneDict(name, record.id) : geneDict[name].append((record.id, record.seq, len(record.seq), "-1"))
                else:
                    log ="WARNING : unexpectable error happened : \n the gene : " + name + " seams to already exist for the taxon : " + record.id + " although the file for this gene just been created \n Please check if their isn't an isue in the name of your taxon or weerd things \n The error come from the file : " + fasta
                    prRed(log)
                    writeLog(log)

        # if the file already exist
        else:
            # print("file already exist")
            file = open(destinationPath2+ name + ".fasta", "a")
            # print("file ouvert")
            # print(list(records))
            for record  in SeqIO.parse(fasta, "fasta"):
                # print("ok")
                listMitogenome.append(record.id)
                if not isGenomeInGeneDict(name, record.id) : 
                    # print("c'est ok")
                    geneDict[name].append((record.id, record.seq, len(record.seq), "-1"))
                    writer = SeqIO.FastaIO.FastaWriter(file)
                    writer.write_record(record)
                    # print("should write the seq" + record.id)

                else:
                    log = "WARNING : Unexpected error : the gene " + name + " seams to already exist for the mitogenome : " + record.id + "\n  It's will be ignored in the file : " + fasta + "\n please check what is going on"
                    prRed(log)
                    writeLog(log)

            file.close()

    return (listMitogenome, "-1")

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
def extractSeqFromGBFile(gbFile, destinationPath = settings["classicFastaResultPath"], destinationPath2 = settings["genesFastaResultPath"]):
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
            geneCounter= {}
            for gene in settings["geneToDetect"]:
                geneCounter[gene] = 0


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

                if str.upper(gene.type) == "CDS" or str.upper(gene.type) == "TRNA" or str.upper(gene.type) == "RRNA": 
                    start = correctMinMaxInputError(str(gene.location.start),"Minimum",mitogenomeName,name)
                    end= correctMinMaxInputError(str(gene.location.end), "Maximum",mitogenomeName,name)
                    strand = gene.location.strand
                    seq = Seq(mitogenomeSeq [start:end])
                    name = str.upper(name)


                    #  name auto rename
                    name = getAlias(name)
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

                    if name in settings["geneToDetect"]:
                        if name in geneCounter.keys():
                            geneCounter[name] +=1
                            if geneCounter[name]  > 1 : 
                                if any(char.isdigit() for char in name):name = name + "-" + str(geneCounter[name]) 
                                else:name = name + str(geneCounter[name])
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

    writeRecords(listRecords, mitogenomeName, destinationPath, destinationPath2 )

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


# function getAlias,
# This function return the correct name if an alias is encontered
# @param,
# @name, the alias to search
def getAlias(name):
    for realName, aliasList in settings["geneAlias"].items():
        if name in aliasList: return str.upper(realName)

    return str.upper(name)


################################################################################
# generating csv summaris for a global view on the extracted data


# function generatePresenceSummary
# This function will creat a csv indicating with 1 or 0 the presence of a gene in a taxon
# @param
# @mitogenomeDict, a dictionary containing all the mitogenome name found during the data extraction
# used to creat the first column (Taxon) of the csv
# @csvPath, path where to write the csv
def generatePresenceSummary(mitogenomeDict, csvPath= settings["csvResultPath"]):
    data={"Taxon":mitogenomeDict.keys()}
    for genome in mitogenomeDict.keys():
        for gene in geneDict.keys():
            if not gene in data.keys():
                data[gene] = []
            
            data[gene].append(int(isGenomeInGeneDict(gene, genome)))

    df = pd.DataFrame.from_dict(data)
    df.to_csv(csvPath+ "summary_presence.csv", index=False)



# function generateLengthSummary
# This function will creat a csv containing the length of each sequence found for the gene in the taxon
# @param
# @mitogenomeDict, a dictionary containing all the mitogenome name found during the data extraction
# used to creat the first column (Taxon) of the csv
# @csvPath, path where to write the csv
def generateLengthSummary(mitogenomeDict, csvPath= settings["csvResultPath"]):
    data={"Taxon":mitogenomeDict.keys()}
    for genome in mitogenomeDict.keys():
        for gene in geneDict.keys():
            if not gene in data.keys():
                data[gene] = []
            
            data[gene].append(getLengthInGeneDict(gene, genome))

    df = pd.DataFrame.from_dict(data)
    df.to_csv(csvPath+ "summary_length.csv", index=False)

# function generateAccessionIDSummary
# this function will creat a csv that give the accessionID for the genes if it's exist
# @param
# @mitogenomeDict, a dictionary containing all the mitogenome name found during the data extraction
# used to creat the first column (Taxon) of the csv
# @csvPath, path where to write the csv
def generateAccessionIDSummary(mitogenomeDict, csvPath= settings["csvResultPath"]):
    data={"Taxon":mitogenomeDict.keys()}
    for genome in mitogenomeDict.keys():
        for gene in geneDict.keys():
            if not gene in data.keys():
                data[gene] = []
            
            data[gene].append(getAccesIDInGeneDict(gene, genome))
    df = pd.DataFrame.from_dict(data)
    df.to_csv(csvPath+ "summary_accessionID.csv", index=False)



# function generateSuperpositionSummary
# this function will creat a csv summarising all the superposition found, indicating the length of it
# @param
# @mitogenomeDict, a dictionary containing all the mitogenome name found during the data extraction
# used to creat the first column (Taxon) of the csv
# @csvPath, path where to write the csv
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


# function getMyPDistMatrix,
# this function return a dictionnary corresponding to 
# the simplified pdistance matrix for thhe alignement
# @param
# @alignementFile, the, the alignement file
def getMyPDistMatrix(alignementFile):
    matrix = {}
    firstSeq = ""
    listRec = []
    for record in AlignIO.read(alignementFile, "phylip-relaxed"):
        matrix[record.id] = []

        distList = []
        listRec.append(str(record.seq))
        # for seq in listRec:
            # distList.append(getMyPDist(str(record.seq), seq)/len(seq))
        distList.append(getMyPDist(str(record.seq), listRec[0])/len(listRec[0]))
        matrix[record.id] = distList

    return matrix

# function getMyPdist,
# this function compute the number of difference between two strings
#  ignoring the difference when a ? is encontered
def getMyPDist(s1, s2):
    nbDif = 0
    for i in range(len(s1)):
        if s1[i] != s2[i] and s1[i] != "?" and s2[i] != "?": nbDif +=1 

    return nbDif


# function aligneSequenceWithMuscleV3
# this function will creat a sequence alignement file (.phy) using the muscleV3 software
# @param
# @fasta, path to the fasta file on witch to process the alignement
# @outputLocation path where to write the alignement file
# @muscleLocation, path to the muscle executable file
def aligneSequenceWithMuscleV3(fasta, outputLocation = settings["sequenceAlignementResultPath"], muscleLocation = settings["musclePath"]):
    osName = platform.system()
    outputFile =""
    log = "Processing MuscleV3 alignement for " + fasta
    print(log)
    outputFile = outputLocation + fasta[fasta.rfind("/")+1:-6]+ "_muscle_align.phy"
    tmpFile= outputLocation + fasta[fasta.rfind("/")+1:-6]+ "_muscle_align.fasta"

    muscleEXE = ""
    if osName == "Linux":
        muscleEXE= muscleLocation  + "muscle3.8.31_i86linux64"
        subprocess.Popen("chmod +x " + muscleEXE, stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32")).wait()
    elif osName == "Windows":
        muscleEXE= muscleLocation+ "muscle3.8.31_i86win32.exe"
    elif osName == "Darwin":
        muscleEXE=muscleLocation+"muscle3.8.31_i86darwin64"
        subprocess.Popen("chmod +x " + muscleEXE, stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32")).wait()

    muscle_cline= muscleEXE + " -in " + os.path.abspath(fasta) + " -out " + os.path.abspath(tmpFile)
    if len(settings["muscleParameterValue"].keys()) >0:
        for k,v in settings["muscleParameterValue"].items():
            muscle_cline += " " + str(k) + " " + str(v)

    if settings["debugLog"] : print("Commande : \n" + muscle_cline)
    # child= subprocess.Popen(str(muscle_cline), stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32"))
    try:
        child= subprocess.check_output(str(muscle_cline), text=True)
    # child.wait()
    # rc= child.returncode
    # stdout = child.stdout
    # stderr = child.stderr
    # if rc == 0:
        correctGape(tmpFile)
        AlignIO.convert(tmpFile, "fasta", outputFile, "phylip-relaxed")
        os.remove(tmpFile)
        return outputFile
    # else:
    except subprocess.CalledProcessError as error: 

        prRed("Error while trying to align")
        log = "The executed command was : \n" + str(muscle_cline) + "\n and the error code return was : " + str(error.returncode) + "\n"
        log += "Please check the log file for more infos, or try to execute the above command directly \n"
        prRed(log)
        log += "Output for error : \n" + str(error.output)
        writeLog(log)
        return outputFile


# function aligneSequenceWithMuscle
# this function will creat a sequence alignement file (.phy) using the muscle software
# @param
# @fasta, path to the fasta file on witch to process the alignement
# @outputLocation path where to write the alignement file
# @muscleLocation, path to the muscle executable file
def aligneSequenceWithMuscle(fasta, outputLocation = settings["sequenceAlignementResultPath"], muscleLocation = settings["musclePath"]):
    if settings["useMuscleV3"] : return aligneSequenceWithMuscleV3(fasta, outputLocation =outputLocation, muscleLocation =muscleLocation  )
    osName = platform.system()
    outputFile =""
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

    muscle_cline= muscleEXE + " -align " + os.path.abspath(fasta) + " -output " + os.path.abspath(tmpFile)
    if settings["debugLog"] : print("Commande : \n" + muscle_cline)
    # child= subprocess.Popen(str(muscle_cline), stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32"))
    try:
        child= subprocess.check_output(str(muscle_cline), text=True)
    # child.wait()
    # rc= child.returncode
    # stdout = child.stdout
    # stderr = child.stderr
    # if rc == 0:
        correctGape(tmpFile)
        AlignIO.convert(tmpFile, "fasta", outputFile, "phylip-relaxed")
        os.remove(tmpFile)
        return outputFile
    # else:
    except subprocess.CalledProcessError as error:
        prRed("Error while trying to align")
        log = "The executed command was : \n" + str(muscle_cline) + "\n and the error code return was : " + str(error.returncode) + "\n"
        log += "Please check the log file for more infos, or try to execute the above command directly \n"
        prRed(log)
        log += "Output for error : \n" + str(error.output)
        writeLog(log)
        return outputFile



# function aligneSequenceWithMafft
# this function will creat a sequence alignement file (.phy) using the mafft software
# @param
# @fasta, path to 
# @outputLocation path where to write the alignement file
# @mafftLocation, path to the mafft executable file
def aligneSequenceWithMafft(fasta, outputLocation = settings["sequenceAlignementResultPath"], mafftLocation = settings["mafftPath"]  ):
    osName = platform.system()
    outputFile =""
    log = "Processing Mafft alignement for " + fasta
    print(log)
    outputFile = outputLocation + fasta[fasta.rfind("/")+1:-6]+ "_mafft_align.phy"
    tmpFile= outputLocation + fasta[fasta.rfind("/")+1:-6] + "_mafft_align.fasta"
    mafftEXE =""

    if osName == "Linux":
        mafftEXE= mafftLocation  + "muscle5.1.linux_intel64"
        # subprocess.Popen("chmod +x " + mafftEXE, stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32")).wait()
        mafft_cline= "/usr/bin/mafft --auto --out "+ os.path.abspath(tmpFile) + " " + os.path.abspath(fasta)
    elif osName == "Windows":
        mafftEXE= mafftLocation + "mafft-win/mafft.bat"
        mafft_cline= os.path.abspath(mafftEXE) + " --auto --out "+ os.path.abspath(tmpFile) + " " + os.path.abspath(fasta)
    elif osName == "Darwin":
        mafftEXE= mafftLocation +"mafft-mac/mafft.bat"
        subprocess.Popen("chmod +x " + mafftEXE, stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32")).wait()
        mafft_cline= os.path.abspath(mafftEXE) + " --auto --out "+ os.path.abspath(tmpFile) + " " + os.path.abspath(fasta)

    if settings["debugLog"] : print("Commande : \n" + mafft_cline)
    # child= subprocess.Popen(str(mafft_cline), stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32"))
    try:
        child= subprocess.check_output(str(mafft_cline), text=True)
    # child.wait()
    # rc= child.returncode
    # stdout = child.stdout
    # stderr = child.stderr
    # if rc == 0 and not "fatal:" in str(stderr.readline()) and os.path.isfile(tmpFile):
        correctGape(tmpFile)
        AlignIO.convert(tmpFile, "fasta", outputFile, "phylip-relaxed")
        os.remove(tmpFile)
        return outputFile
    # else:
    except subprocess.CalledProcessError as error:
        prRed("Error while trying to align")
        log = "The executed command was : \n" + str(muscle_cline) + "\n and the error code return was : " + str(error.returncode) + "\n"
        log += "Please check the log file for more infos, or try to execute the above command directly \n"
        prRed(log)
        log += "Output for error : \n" + str(error.output)
        writeLog(log)
        return outputFile

# function getAlignementDict
# this function will parth the alignement file to creat a dictionnary containing the alignement record
# @param
# @alignementType, "mafft" or "muscle", type of the alignement you want to work with
# @path, path where to find the alignement files 
def getAlignementDict(alignementType, path = settings["sequenceAlignementResultPath"]):
    alignementDict = {}
    fileList = [ file for file in os.listdir(path) if alignementType in file and not "Matrix" in file]
    for file in fileList:
        geneName = file[file.rfind("/")+1:file.find("_")].replace(".phy","").replace("_", "")
        if not geneName in alignementDict.keys(): alignementDict[geneName] = []

        aln =AlignIO.read(path+file, "phylip-relaxed")
        for record in aln:
            alignementDict[geneName].append((record.id, record.seq, len(record.seq), record))

    return alignementDict


# function getAlignedMitogenomeDict    
# this function will parse the alignement file to creat a dictionnary containing the mitogenme name foundable in the records
# @param
# @alignementType, "mafft" or "muscle", type of the alignement you want to work with
# @path, path where to find the alignement files 
def getAlignedMitogenomeDict    (alignementType, path = settings["sequenceAlignementResultPath"]):
    alignedMitogenomeDict    = {}
    fileList = [ file for file in os.listdir(path) if alignementType in file and not "Matrix" in file]
    for file in fileList:
        aln =AlignIO.read(path+file, "phylip-relaxed")
        for record in aln:
            if not record.id in alignedMitogenomeDict    .keys() : alignedMitogenomeDict    [record.id] = 1

    return alignedMitogenomeDict    

# function correctGape
# this function will replace the gape on the begining and 
# end of the align sequence by "?", for unknown data, for a better accuracy for the trees
# @param
# @alignementFile, the alignement file on witch to process the correction, (should be fasta)
def correctGape(alignementFile):
    nucléotyde = ["A", "T", "C", "G", "U", "N"]

    msa= AlignIO.read(alignementFile, "fasta")
    for record in msa:
        seq =str(record.seq).upper()
        if seq.startswith("-"):
            nucIndex = []
            for nuc in nucléotyde:
                if seq.find(nuc) != -1 : nucIndex.append(seq.find(nuc))

            nucIndex.sort()
            seq= seq[:nucIndex[0]].replace("-", "?") + seq[nucIndex[0]:]


        # correction of end gape
        if seq.endswith("-"):
            nucIndex = []
            for nuc in nucléotyde:
                if seq.rfind(nuc) != -1 : nucIndex.append(seq.rfind(nuc))

            nucIndex.sort()
            seq= seq[:nucIndex[-1]]  + seq[nucIndex[-1]:].replace("-", "?")


        record.seq= Seq(seq)

    AlignIO.write(msa, alignementFile,  "fasta")    



# function getseqInAlignementDict
# this function return the seq in the alignementDict if it exist
#  or a sequence of "?" of the correct length if it's not found
# @param
# @alignementDict, dictionarry to search in
# @gene, name of the gene, corresponding to the key
# @mitogenomeName, name of the mitogenome we want to get the sequence
def getseqInAlignementDict(alignementDict, gene, mitogenomeName):
    len = 0
    for taxon, seq, length, desc   in alignementDict[gene]:
        len = length
        if taxon == mitogenomeName: return seq

    return Seq("?"*len)

# function writeConcatenatedMatrix
# this function will output an alignement file containing all the alignement data
# to build a global concatenated matrix
# @param
# @alignementDict, a dictionnary containing the alignement data
# @mitogenomeDict, a dictionnary containing all the taxon name found during the data extraction
# @alignementType, "mafft" or "muscle", type of the alignement you want to work with
# @destinationPath path where to write the file
def writeConcatenatedMatrix(alignementDict, mitogenomeDict, alignementType, destinationPath = settings["sequenceAlignementResultPath"]):
    data ={}
    anotationDict = {}
    for genome in mitogenomeDict.keys():
        lengthCounter= 0
        if not genome in data.keys(): data[genome] = ""
        for gene in alignementDict.keys():
            alignement = str(getseqInAlignementDict(alignementDict, gene, genome))
            if not gene in anotationDict .keys() :
                anotationDict [gene] =[ str(lengthCounter+1) + "-" + str(lengthCounter+ len(alignement))]
            data[genome] += alignement
            lengthCounter= len(data[genome])

            
    msa =[]
    for id, seq in data.items():
        msa.append(SeqRecord(Seq(seq), id = id))

    msa = MultipleSeqAlignment(msa, annotations=anotationDict)
    file = AlignIO.write(msa, destinationPath+ "globalConcatenatedMatrix_" + alignementType + ".phy", "phylip-relaxed")
    AlignIO.convert(destinationPath+ "globalConcatenatedMatrix_" + alignementType + ".phy", "phylip-relaxed", destinationPath+ "globalConcatenatedMatrix_" + alignementType + ".nex", "nexus", molecule_type="DNA")
    infoLocation = "[BEGIN MrBayes;" + "\n" + "	outgroup OUT1;" + "\n"
    for gene, info in anotationDict.items():
        infoLocation += "	CHARSET " + str(gene) + " = " + info[0] + ";" + "\n"

    infoLocation+= "	partition markers = 7: " + str(anotationDict.keys()).replace("[", "").replace("]", "") + "\n" + "end;]"
    fileContent =""
    with open(destinationPath+ "globalConcatenatedMatrix_" + alignementType + ".nex", "r") as fileObj:
        fileContent = str(fileObj.read())
    with open(destinationPath+ "globalConcatenatedMatrix_" + alignementType + ".nex", "w") as fileObj:
        fileObj.write(fileContent.replace("end;", infoLocation))




# function writeSingleRecordList
# this function write a list of record
# @param
# @tmpFasta, path to the file to write
# @listRec, list of record to write
def writeSingleRecordList(tmpFasta, listRec):
    with open(tmpFasta, "w") as file:
        writer = SeqIO.FastaIO.FastaWriter(file)
        writer.write_file(listRec)


# function getReversedRecordList
# this function return a list of record, with reverse record for the one mentionned in parameter
# @param
# fastaFile, name of the fasta file to parth
# @taxonName, name of the taxon (seqId) to reverse
def getReversedRecordList(fastaFile, taxonToReverse):
    listRec = []
    for record in SeqIO.parse(fastaFile, "fasta"):
                if record.id in taxonToReverse:
                    rec = SeqRecord(record.seq.reverse_complement(), id=record.id, name=record.name, description="")
                else:
                    rec = record
                listRec.append(rec)

    return listRec

# function getPDistMatrix
# this function return a p-distance matrix for the alignement file in parameter
# @param
# @alignementFile, alignement file to compute the p-distance matrix
def getPDistMatrix(alignementFile):
    aln = list(AlignIO.parse(alignementFile, "phylip-relaxed"))
    calculator = DistanceCalculator('identity')
    return calculator.get_distance(aln[0])


# function checkMafftAlignement
# this function check the alignement for the mafft program
# it will check if the p-distance is > to a predefind value, and try to reverse the sequence, and realign the sequence to 
# check if it increase the alignement quality
# @param
# @alignementFile, alignement file to check
# @alignementPath, path to the alignement file are located, to read and write 
# @pathToFasta, path where to read / write fasta 
def checkMafftAlignement(alignementFile, alignementPath= settings ["sequenceAlignementResultPath"], pathToFasta = settings["genesFastaResultPath"], mafftLocation = settings["mafftPath"]):
    if not  os.path.isfile(alignementFile):
        log = "Aborting alignement check for " + alignementFile + "it seams that the basic alignement file has not been found, certainly due to an error during it"
        prRed(log)
        writeLog(log)
        return
    print("alignement check initialised for mafft ")
    dm = getMyPDistMatrix(alignementFile)
    fastaFile = pathToFasta + alignementFile[alignementFile.rfind("/")+1:alignementFile.find("_")] + ".fasta"
    tmpFasta = fastaFile[:-6] + "_tmp.fasta"
    tmpAlignement = ""
    listRec= []
    prevValues = {}
    taxonToRevers = []
    taxonToKeep = []
    gene = fastaFile[fastaFile.rfind("/")+1: fastaFile.rfind(".")]
    threshold = settings ["pDistThreshold"]["default"]
    if gene in settings["pDistThreshold"]: threshold = settings["pDistThreshold"][gene]

    for taxonName, pDistList in dm.items():
        taxonPDist = sum (pDistList) / len(pDistList)
        prevValues[taxonName]= taxonPDist
        if taxonPDist > threshold: taxonToRevers.append(taxonName)

    if len(taxonToRevers) > 0:
        if settings["debugLog"]: prCyan("Taxon that will be reverse complement for testing : "+ str(taxonToRevers))
        listRec = getReversedRecordList(fastaFile, taxonToRevers)
        writeSingleRecordList(tmpFasta, listRec)
        tmpAlignement = aligneSequenceWithMafft(tmpFasta, outputLocation=alignementPath, mafftLocation=mafftLocation)
        dm2 = getMyPDistMatrix(tmpAlignement)
        for taxonName, pDistList in dm.items():
            taxonPDist = sum (pDistList) / len(pDistList)
            if prevValues[taxonName] > taxonPDist  and taxonName in taxonToRevers:
                # print("pdist for " + taxonName+ " was " + str(prevValues[taxonName] )+ " and after reverse, it became "+ str(taxonPDist ))
                taxonToKeep.append(taxonName)

    if len(taxonToKeep) > 0:
        if settings["debugLog"]: prCyan("taxon where the reverse complement seams more efficiant : " + str(taxonToKeep))
        listRec = getReversedRecordList(fastaFile, taxonToKeep)
        writeSingleRecordList(fastaFile, listRec)
        finalAlignement = aligneSequenceWithMafft(fastaFile, outputLocation=alignementPath, mafftLocation=mafftLocation)
        dm3 = getMyPDistMatrix(finalAlignement )
        # print("values of pdist after reversing and keeping the best sequence")
        for taxonName, pDistList in dm.items():
            taxonPDist = sum (pDistList) / len(pDistList)
            # if taxonName in taxonToKeep : print(taxonName + " : " + str(taxonPDist))
            if taxonPDist > threshold:
                prRed("Warning : Even after reverse complement the taxon " + taxonName + " still have a p-distance > to the  defind threshold (" + str(threshold) + ") for the alignement file " + alignementFile)

    if os.path.isfile(tmpFasta): os.remove(tmpFasta)
    if os.path.isfile(tmpAlignement): os.remove(tmpAlignement)
    print("alignementCheck finished")

# function checkMuscleAlignement
# this function check the alignement for the mafft program
# it will check if the p-distance is > to a predefind value, reverse the concerne sequence, and re align the sequence to 
# check if it increase the alignement quality
# @param
# @alignementFile, alignement file to check
# @alignementPath, path to the alignement file are located, to read and write 
# @pathToFasta, path where to read / write fasta 
def checkMuscleAlignement(alignementFile, alignementPath= settings ["sequenceAlignementResultPath"], pathToFasta = settings["genesFastaResultPath"], muscleLocation = settings["musclePath"]):
    if not  os.path.isfile(alignementFile):
        log = "Aborting alignement check for " + alignementFile + "it seams that the basic alignement file has not been found, certainly due to an error during it"
        prRed(log)
        writeLog(log)
        return

    print("alignement check initialised for muscle  ")
    dm = getMyPDistMatrix(alignementFile)
    fastaFile = pathToFasta + alignementFile[alignementFile.rfind("/")+1:alignementFile.find("_")] + ".fasta"
    tmpFasta = fastaFile[:-6] + "_tmp.fasta"
    tmpAlignement = ""
    listRec= []
    prevValues = {}
    taxonToRevers = []
    taxonToKeep = []
    gene = fastaFile[fastaFile.rfind("/")+1: fastaFile.rfind(".")]
    threshold = settings ["pDistThreshold"]["default"]
    if gene in settings["pDistThreshold"]: threshold = settings["pDistThreshold"][gene]


    for taxonName, pDistList in dm.items():
        taxonPDist = sum (pDistList) / len(pDistList)
        prevValues[taxonName]= taxonPDist
        if taxonPDist >= threshold: taxonToRevers.append(taxonName)

    if len(taxonToRevers) > 0:
        print(taxonToRevers)
        listRec = getReversedRecordList(fastaFile, taxonToRevers)
        writeSingleRecordList(tmpFasta, listRec)
        tmpAlignement = aligneSequenceWithMuscle(tmpFasta, outputLocation=alignementPath, muscleLocation=muscleLocation)
        dm2 = getMyPDistMatrix(tmpAlignement)
        for taxonName, pDistList in dm.items():
            taxonPDist = sum (pDistList) / len(pDistList)
            if prevValues[taxonName] > taxonPDist  and taxonName in taxonToRevers:
                print("pdist for " + taxonName+ " was " + str(prevValues[taxonName] )+ " and after reverse, it became "+ str(taxonPDist ))
                taxonToKeep.append(taxonName)

    if len(taxonToKeep) > 0:
        print(taxonToKeep)
        listRec = getReversedRecordList(fastaFile, taxonToKeep)
        writeSingleRecordList(fastaFile, listRec)
        finalAlignement = aligneSequenceWithMuscle(fastaFile, outputLocation=alignementPath, muscleLocation=muscleLocation)
        dm3 = getMyPDistMatrix(finalAlignement )
        print("values of pdist after reversing and keeping the best sequence")
        for taxonName, pDistList in dm.items():
            taxonPDist = sum (pDistList) / len(pDistList)
            if taxonName in taxonToKeep : print(taxonName + " : " + str(taxonPDist))
            if taxonPDist > threshold:
                prRed("Warning : Even after reverse complement the taxon " + taxonName + " still have a p-distance > to the  defind threshold (" + str(threshold) + ") for the alignement file " + alignementFile)

    if os.path.isfile(tmpFasta): os.remove(tmpFasta)
    if os.path.isfile(tmpAlignement): os.remove(tmpAlignement)
    print("alignementCheck finished")

def generatePDistSummary(alignMitogenomeDict, alignementDict, type, pathToFile = settings["sequenceAlignementResultPath"], csvPath= settings["csvResultPath"]):
    data={"Taxon":alignMitogenomeDict.keys()}
    matrixDict= {}
    
    for gene in alignementDict.keys():
        if not gene in data.keys(): data[gene] = []
        if not gene in matrixDict.keys(): matrixDict[gene] = getMyPDistMatrix(pathToFile + gene + "_" +type + "_align.phy")

    for gene, matrix in matrixDict.items():
        for genome in alignMitogenomeDict.keys():
            if genome in matrix.keys() : data[gene].append(sum(matrix[genome]) / len(matrix[genome]))
            else: data[gene].append(pd.NA)

    df = pd.DataFrame.from_dict(data)
    df.to_csv(csvPath+ "summary_pDist.csv", index=False)
    





################################################################################

# function generateDistanceTree
# this function will creat distance tree 
# for cds gene, and rRNA, and one for the global concatenation alignement matrix
# @param
# @alignementType, the type of alignement, "mafft" or "muscle"
# @alignementLocation , where the alignement file are located, by default the value is the one set in the settings.py file
# @outputLocation, where to output the tree files, by default the value is the one set in the settings.py file
# tree section
def generateDistanceTree(alignementType, alignementLocation = settings["sequenceAlignementResultPath"], outputLocation = settings["treeResultPath"]):
    fileList = [ file for file in os.listdir(alignementLocation) if alignementType in file and file.endswith(".phy")]

    for alignementFile in fileList:
        if not alignementFile.startswith("TRNA"):
            gene = alignementFile[:alignementFile.find("_")]
            aln = AlignIO.read(alignementLocation + alignementFile, 'phylip-relaxed')

            print("Generating distance tree for : " + alignementFile)
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(aln)
            constructor = DistanceTreeConstructor(calculator, 'nj')
            tree =constructor.build_tree(aln)
            Phylo.write(tree, outputLocation+ gene +"_" + alignementType+ "_nj_dist_tree.tree", "newick")





# run,
# this function is the main function to run
def run():
    setup()
    writeLog("Starting treatement", firstTime=True)
    for dir in os.listdir("./results"):
        for file in os.listdir("results/"+dir):
            os.remove("results/"+dir+"/"+file)

    mitogenomeDict = {}
    tExtraction = time.time()
    fastaFiles = getFASTAFiles()
    csvFiles = getCSVFiles()
    gbFiles = getGBFiles()
    couple = []
    singleFasta = []
    for f in fastaFiles:
        for c in csvFiles:
            if f[:-5] in c:
                couple.append((c,f))
                singleFasta.append(f)
    
    singleFasta= [x for x in fastaFiles if x not in singleFasta]
    couple.sort()
    singleFasta.sort()
    gbFiles.sort()
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

    for file in singleFasta:
        mitogenomes, accessionID = extractSeqFromSingleFasta(settings["rawFilePath"]+file)
        for mitogenomeName  in mitogenomes:
            if mitogenomeName not in mitogenomeDict.keys(): mitogenomeDict[mitogenomeName] = accessionID
    tExtraction = time.time() - tExtraction
    tSummaryGeneration = time.time()

# generation of csv summary
    generatePresenceSummary(mitogenomeDict)
    generateLengthSummary(mitogenomeDict)
    generateAccessionIDSummary(mitogenomeDict)
    generateSuperpositionSummary(mitogenomeDict)
    tSummaryGeneration= time.time() - tSummaryGeneration

    tMuscleAlignement = time.time()
    tMafftAlignement = time.time()
    # Alignement of file
    if settings["useMuscle"]:
        for fasta in getFASTAFiles(path=settings ["genesFastaResultPath"]):
            alignedFile = aligneSequenceWithMuscle(settings ["genesFastaResultPath"] + fasta)
            if settings["checkAlignement"] : checkMuscleAlignement(alignedFile)
        tMuscleAlignement = time.time() - tMuscleAlignement

    if settings["useMafft"]:
        for fasta in getFASTAFiles(path=settings ["genesFastaResultPath"]):
            alignedFile = aligneSequenceWithMafft(settings ["genesFastaResultPath"] + fasta)
            if settings["checkAlignement"] : checkMafftAlignement(alignedFile)
        tMafftAlignement = time.time() - tMafftAlignement


    tMuscleTree = time.time()
    tMafftTree = time.time()
    if   settings["useMuscle"]:
        alignementDict= getAlignementDict("muscle")
        generatePDistSummary(getAlignedMitogenomeDict("muscle"), alignementDict, "muscle")
        writeConcatenatedMatrix(alignementDict, mitogenomeDict, "muscle")
        generateDistanceTree("muscle")
        tMuscleTree= time.time() - tMuscleTree

    if settings["useMafft"]:
        alignementDict= getAlignementDict("mafft")
        generatePDistSummary(getAlignedMitogenomeDict("mafft"), alignementDict, "mafft")
        writeConcatenatedMatrix(alignementDict, mitogenomeDict, "mafft")
        generateDistanceTree("mafft")
        tMafftTree = time.time() - tMafftTree


    prYellow("Data extraction was made in : "+ str(tExtraction) + " Secondes")

    if settings["useMuscle"]:
        prYellow("Muscle Alignement and their check was made in : " + str(tMuscleAlignement) + " Secondes")
        prYellow("Distance tree based on muscle alignement was made in : " + str(tMuscleTree) + " Secondes")

    if settings["useMafft"]:
        prYellow("Mafft Alignement and their check was made in : " + str(tMafftAlignement) + " Secondes")
        prYellow("Distance tree based on mafft alignement was made in : " + str(tMafftTree) + " Secondes")
    
    prGreen("Finish")
    writeLog("Finish")



################################################################################
# initialisation of global variable and environement
# setup()
geneDict = {}
superpositionDict={}
init() # to have color in the terminal





################################################################################
# debug function
# put what ever you want inside to debug
def debug():
    # for record in SeqIO.parse(settings ["genesFastaResultPath"] + "12S.fasta", "fasta"):
    # record = list(SeqIO.parse(settings ["genesFastaResultPath"] + "12S.fasta", "fasta"))
    # print(len(record))
        # if "BL" in record.id :print(len(record.seq))
        # print(record.id + ":"+str(len(record.seq)))checkMu
        # checkMuscleAlignement(settings ["sequenceAlignementResultPath"] + "12S_muscle_align.phy")

    # alignedFile = aligneSequenceWithMafft(settings ["genesFastaResultPath"] + "18S.fasta")
    # alignedFile = aligneSequenceWithMafft(settings ["genesFastaResultPath"] + "17S.fasta")
    # alignedFile = aligneSequenceWithMuscle(settings ["genesFastaResultPath"] + "18S.fasta")
    # alignedFile = aligneSequenceWithMuscle(settings ["genesFastaResultPath"] + "17S.fasta")
    # taxonOfIntrest = ["BL619", "BL862"]
    # listPdistMean = []
    # listRec = getReversedRecordList(settings ["genesFastaResultPath"] + "12S.fasta", taxonOfIntrest )
    # writeSingleRecordList(settings ["genesFastaResultPath"] + "12S_tmp.fasta", listRec)
    # aligneSequenceWithMuscle(settings ["genesFastaResultPath"] + "12S_tmp.fasta")
    # correctGape(settings ["sequenceAlignementResultPath"] + "12S_muscle_align.fasta")
    # AlignIO.convert(settings ["sequenceAlignementResultPath"] + "12S_muscle_align.fasta", "fasta", settings ["sequenceAlignementResultPath"] + "12S_muscle_align.phy", "phylip-relaxed")
    df=pd.read_csv(settings["csvResultPath"]+ "summary_presence.csv")
    data= {}
    geneList = []
    geneCount=[]
    print(len(df))
    for gene in settings["geneToDetect"]:
        if gene in df:
            if not gene in geneList:
                geneList.append(gene)
                geneCount.append(df[gene].sum())

            # if not gene in data.keys(): data[gene] = []
            # data[gene].append(df[gene].sum())

    print(data)
    # sumdf = pd.DataFrame.from_dict(data)
    sumdf = pd.DataFrame.from_dict({"gene_name":geneList, "species_count":geneCount})
    sumdf.to_csv("tmp.csv", index= False, sep=";")
    print(sumdf)

    # alignementDict= getAlignementDict("muscle")
    # generatePDistSummary(getAlignedMitogenomeDict("muscle"), alignementDict, "muscle")
    # dm = getMyPDistMatrix(settings ["sequenceAlignementResultPath"] + "12S_muscle_align.phy")
    # dm2 = getMyPDistMatrix(settings ["sequenceAlignementResultPath"] + "12S_tmp_muscle_align.phy")

    # print("first matrix")
    # for taxon, val in dm.items():
        # listPdistMean.append(sum(val) / len(val))
        # if taxon in taxonOfIntrest: 
            # print(taxon + " : " + str(val))
            # print(taxon + " : " + str(sum(val) / len(val)))

    # print("mean of pdist : " + str(sum(listPdistMean) / len(listPdistMean)))
    # listPdistMean= []
    # print("matrix after reverse")
    # for taxon, val in dm2.items():
        # listPdistMean.append(sum(val) / len(val))
        # if taxon in taxonOfIntrest: 
            # print(taxon + " : " + str(val))
            # print(taxon + " : " + str(sum(val) / len(val)))

    # print("mean of pdist : " + str(sum(listPdistMean) / len(listPdistMean)))

    

debug()