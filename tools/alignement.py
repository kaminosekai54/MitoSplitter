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
    tMuscleAlignement = time.time()
    tMafftAlignement = time.time()
    # Alignement of file
    if settings["useMuscle"]:
        for fasta in getFASTAFiles(path="." + settings ["genesFastaResultPath"]):
            alignedFile = aligneSequenceWithMuscle("."+ settings ["genesFastaResultPath"] + fasta, outputLocation ="." + settings["sequenceAlignementResultPath"], muscleLocation = "." + settings["musclePath"])
            if settings["checkAlignement"] : checkMuscleAlignement(alignedFile, alignementPath="."+ settings ["sequenceAlignementResultPath"], pathToFasta ="."+settings["genesFastaResultPath"], muscleLocation = "." + settings["musclePath"])
        tMuscleAlignement = time.time() - tMuscleAlignement

    if settings["useMafft"]:
        for fasta in getFASTAFiles(path="."+settings ["genesFastaResultPath"]):
            alignedFile = aligneSequenceWithMafft("."+ settings ["genesFastaResultPath"] + fasta, outputLocation ="." + settings["sequenceAlignementResultPath"], mafftLocation = "." + settings["mafftPath"])
            if settings["checkAlignement"] : checkMafftAlignement(alignedFile, alignementPath="."+ settings ["sequenceAlignementResultPath"], pathToFasta ="."+settings["genesFastaResultPath"], mafftLocation = "." + settings["mafftPath"])
        tMafftAlignement = time.time() - tMafftAlignement

    if settings["useMuscle"]:
        prYellow("Muscle Alignement and their check was made in : " + str(tMuscleAlignement) + " Secondes")

    if settings["useMafft"]:
        prYellow("Mafft Alignement and their check was made in : " + str(tMafftAlignement) + " Secondes")

    if os.path.isdir("./raw_data") : os.removedirs("./raw_data")
    if os.path.isdir("./results") :
        for dir in os.listdir("./results"):
            os.removedirs("./results/"+dir)
            
    prGreen("Finish")

if __name__ == '__main__':
    main()