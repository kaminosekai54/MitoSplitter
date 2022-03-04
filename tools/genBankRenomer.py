# @author, Alexis Culpin 
# This is a simple stand alone script
# to rename .gb file according to the organism they are treating


################################################################################
# import
import os
from Bio import GenBank

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
    with open(gbFile) as gb:
        record = GenBank.read(gb)
        mitogenomeName = str.upper(record.organism.replace(" ", "-"))

    if gbFile[gbFile.rfind("/")+1:-3] != mitogenomeName:
        os.rename(gbFile, gbFile[:gbFile.rfind("/")+1] + mitogenomeName+ ".gb")

# main,
# this is the main function to run
def main():
    gbFiles = getGBFiles(sourcePath)
    for file in gbFiles:
        renameGBFile(sourcePath+ file)


if __name__ == '__main__':
    main()