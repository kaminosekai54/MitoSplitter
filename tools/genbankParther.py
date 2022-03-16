# importing of package
import os, re
from xmlrpc.server import SimpleXMLRPCRequestHandler
from Bio import Entrez



################################################################################
# initialisation of global variable
# definding the email used for quaries
Entrez.email = "alexis.culpin@cri-paris.org"
# getting the last updated infos:
info = Entrez.einfo()

def search():
    dataBaseList = Entrez.read(info)["DbList"]
    print(dataBaseList)
    print("Please type the number corresponding to the database you want to search in :")
    for i in range(len(dataBaseList)):
        print(str(i) + " : "+ dataBaseList[i])

    databaseID = input("Number of database id")
    while not re.match("^[0-9]*$", databaseID) or not databaseID or int(databaseID) >= len(dataBaseList):
        print("error, please enter a valid number")
        databaseID = input("Number of database :")

    database = dataBaseList[int(databaseID)]

    print("please enter the key words you want to search")
    termes = input("Enter the terme you want to search for, just like you would do in genbank search bar")
    handle = Entrez.esearch(db=database, term=termes)
    record = Entrez.read(handle)
    handle.close()
    # print(record)
    return (database, record['IdList'])

def downloadFiles(db, idList, outputPath):
    for id in idList:
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





db, idList = search()
downloadFiles(db, idList, "../raw_data/")