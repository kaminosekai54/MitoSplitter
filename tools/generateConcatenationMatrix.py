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

    tMuscleTree = time.time()
    tMafftTree = time.time()
    if   settings["useMuscle"]:
        mitogenomeDict= getAlignedMitogenomeDict   ("muscle", path ="."+ settings["sequenceAlignementResultPath"])
        alignementDict= getAlignementDict("muscle", path ="."+ settings["sequenceAlignementResultPath"])
        generatePDistSummary(mitogenomeDict, alignementDict, "muscle", pathToFile="."+ settings["sequenceAlignementResultPath"], csvPath="."+ settings["csvResultPath"])
        writeConcatenatedMatrix(alignementDict, mitogenomeDict, "muscle", destinationPath = "." + settings["sequenceAlignementResultPath"])
        generateDistanceTree("muscle", alignementLocation ="." + settings["sequenceAlignementResultPath"], outputLocation = "." + settings["treeResultPath"])
        tMuscleTree= time.time() - tMuscleTree

    if settings["useMafft"]:
        mitogenomeDict= getAlignedMitogenomeDict    ("mafft", path ="."+ settings["sequenceAlignementResultPath"])
        alignementDict= getAlignementDict("mafft", path ="."+ settings["sequenceAlignementResultPath"])
        generatePDistSummary(mitogenomeDict, alignementDict, "mafft", pathToFile="."+ settings["sequenceAlignementResultPath"], csvPath="."+ settings["csvResultPath"])
        writeConcatenatedMatrix(alignementDict, mitogenomeDict, "mafft", destinationPath = "." + settings["sequenceAlignementResultPath"])
        generateDistanceTree("mafft",alignementLocation ="." + settings["sequenceAlignementResultPath"], outputLocation = "." + settings["treeResultPath"])
        tMafftTree = time.time() - tMafftTree

    if settings["useMuscle"]:
        prYellow("Distance tree and concatenation matrix based on muscle alignement was made in : " + str(tMuscleTree) + " Secondes")

    if settings["useMafft"]:
        prYellow("Distance tree and concatenation matrix based on mafft alignement was made in : " + str(tMafftTree) + " Secondes")
        
    if os.path.isdir("./raw_data") : os.removedirs("./raw_data")
    if os.path.isdir("./results") :
        for dir in os.listdir("./results"):
            os.removedirs("./results/"+dir)    
    prGreen("Finish")


if __name__ == '__main__':
    main()