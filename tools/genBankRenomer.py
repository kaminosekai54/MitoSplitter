# @author, Alexis Culpin 
# This is a simple stand alone script
# to rename .gb file according to the organism they are treating


################################################################################
# import
import os
from Bio import GenBank
from Bio import SeqIO

################################################################################

# global variable:
sourcePath = "./" # the path where to search the file


################################################################################
#  functions

# getGBFiles
# this function return a list of all the .gb files
# found in the indicated folder
# @param,
#@ path, the path of the folder you want to parth
# @ extension, extension of the file we search for, here ".gb" by default"
def getGBFiles(path, extension = ".gb"):
    fileList = os.listdir(path)
    return [ gbFile for gbFile in fileList if gbFile.endswith( extension) ]


# renameGBFile
# this function will rename a .gb file according to its 
# organism
# @param
# @gbFile, the path to the genbankFile
def renameGBFile(gbFile):
    mitogenomeName = ""
    newName = ""
    with open(gbFile) as gb:
        for record in SeqIO.parse(gb, "genbank"):
            if record.features[0].qualifiers["organism"]:
                mitogenomeName = record.features[0].qualifiers["organism"][0]
                if " " in mitogenomeName:
                    mitogenomeName = mitogenomeName.replace(mitogenomeName[0], str.upper(mitogenomeName[0]), 1).replace(mitogenomeName[1:], str.lower(mitogenomeName[1:]), 1).replace(" ", "-")

                break

    if mitogenomeName!= "" and   (gbFile[gbFile.rfind("/")+1:-3] != mitogenomeName):
        newName = gbFile[:gbFile.rfind("/")+1] + mitogenomeName+ ".gb"
        if os.path.isfile(newName):
            i = 2
            while os.path.isfile(newName):
                newName = gbFile[:gbFile.rfind("/")+1] + mitogenomeName+ "_" + str(i) + ".gb"
                i+=1

        
    if newName !="" :
        os.rename(gbFile, newName)

# main,
# this is the main function to run
def main():
    gbFiles = getGBFiles(sourcePath)
    for file in gbFiles:
        renameGBFile(sourcePath+ file)


if __name__ == '__main__':
    main()