# import
import os, sys
# to be able to imort parent functions
path = os.path.abspath("../")
if not path in sys.path:
    sys.path.append(path)

from functions import *
from setting import *


################################################################################
# global variable
settings = getSettings()


################################################################################
# main function
def main():
    mitogenomeDict = {}
    geneDict = getGeneDict()
    tExtraction = time.time()
    fastaFiles = getFASTAFiles(path ="." +  settings["rawFilePath"])
    csvFiles = getCSVFiles(path ="." +  settings["rawFilePath"])
    gbFiles = getGBFiles(path ="." +  settings["rawFilePath"])
    couple = []
    singleFasta = []
    for f in fastaFiles:
        for c in csvFiles:
            if f[:-5] in c:
                couple.append((c,f))
                singleFasta.append(f)
    
    singleFasta= [x for x in fastaFiles if x not in singleFasta]
    for c, f in couple:
        mitogenomeName, accessionID = extractSeqFromCSV("." + settings["rawFilePath"] + c, "." + settings["rawFilePath"] + f, destinationPath = "." + settings["classicFastaResultPath"], destinationPath2 = "." + settings["genesFastaResultPath"], csvPath= "."+ settings["csvResultPath"])
        if mitogenomeName not in mitogenomeDict.keys(): mitogenomeDict[mitogenomeName] = accessionID
    
    for file in gbFiles:
        mitogenomeName, accessionID = extractSeqFromGBFile("." + settings["rawFilePath"] +file, destinationPath = "." +  settings["classicFastaResultPath"], destinationPath2 = "." + settings["genesFastaResultPath"])
        if mitogenomeName not in mitogenomeDict.keys(): 
            mitogenomeDict[mitogenomeName] = [accessionID]
        else:
            for acces in accessionID:
                if not acces in mitogenomeDict[mitogenomeName]: mitogenomeDict[mitogenomeName].append(acces)

    for file in singleFasta:
        mitogenomes, accessionID = extractSeqFromSingleFasta("." + settings["rawFilePath"]+file, destinationPath = "." +  settings["classicFastaResultPath"], destinationPath2 = "." + settings["genesFastaResultPath"])
        for mitogenomeName  in mitogenomes:
            if mitogenomeName not in mitogenomeDict.keys(): mitogenomeDict[mitogenomeName] = accessionID
    tExtraction = time.time() - tExtraction
    tSummaryGeneration = time.time()

    # generation of csv summary
    generatePresenceSummary(mitogenomeDict, csvPath= "."+ settings["csvResultPath"])
    generateLengthSummary(mitogenomeDict, csvPath= "."+ settings["csvResultPath"])
    generateAccessionIDSummary(mitogenomeDict, csvPath= "."+ settings["csvResultPath"])
    generateSuperpositionSummary(mitogenomeDict,csvPath= "."+ settings["csvResultPath"])
    tSummaryGeneration= time.time() - tSummaryGeneration

    prYellow("Data extraction was made in : "+ str(tExtraction) + " Secondes")

    if os.path.isdir("./raw_data") : os.removedirs("./raw_data")
    if os.path.isdir("./results") :
        for dir in os.listdir("./results"):
            os.removedirs("./results/"+dir)

    prGreen("Finish")


if __name__ == '__main__':
    main()