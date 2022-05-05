# importing of package
import sys,os, re
from Bio import Entrez
import pandas as pd
# to be able to imort parent functions
path = os.path.abspath("../")
if not path in sys.path:
    sys.path.append(path)

from functions import getCSVDelimiter


################################################################################
# definding the email used for quaries
Entrez.email = "alexis.culpin@cri-paris.org"


################################################################################
# function
def getIDListFromCSV(file = "./list_accessionID.csv"):
    df = ""
    if os.path.isfile(file):
        separator= getCSVDelimiter(file)
        df=pd.read_csv(file, sep=separator)
    else:
        print("File not found : " + file + "\n Please make sure your file is named correctly, have the correct extension, and   is at the correct position")
        return -1

    if "accessionID" in df:
        return df["accessionID"].tolist()
    else:
        print("Column not found : accessionID \n Please make sure it's name correctly in your file (Case sensitive)")
        return -1

def downloadFiles(idList, db = "nuccore", outputPath = "../raw_data/"):
    for id in idList:
        print("download file correponding to accessionID :"+ id)
        handle = Entrez.efetch(db=db, id=id, rettype="gb", retmode="text")
        record =handle.read()
        fileName = "sequence.gb"
        index = 0
        while os.path.isfile(outputPath+ fileName):
            index +=1
            fileName = "sequence_" + str(index) + ".gb"
        with open(outputPath + fileName, "w") as file:
            file.write(record)
        handle.close()


def main():
    idList = getIDListFromCSV()
    print(idList)
    if idList != -1:
        downloadFiles(idList=idList)
        if os.path.isdir("./raw_data") : os.removedirs("./raw_data")
        if os.path.isdir("./results") :
            for dir in os.listdir("./results"):
                os.removedirs("./results/"+dir)
        print("Download finish")
    else:
        if os.path.isdir("./raw_data") : os.removedirs("./raw_data")
        if os.path.isdir("./results") :
            for dir in os.listdir("./results"):
                os.removedirs("./results/"+dir)
        return

if __name__ == '__main__':
    main()