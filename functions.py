# Import
import os  
import pandas as pd



################################################################################

# function setup,
# this function will create all the folder need for the app to start if the don't 
# already exist
def setup():
    if not os.path.isdir("./ressources"):
        os.makedirs("./ressources")

    if not os.path.isdir("./ressources/raw_data"):
        os.makedirs("./ressources/raw_data")

    if not os.path.isdir("./ressources/results"):
        os.makedirs("./ressources/results")


# function getCSVFiles,
# this function return a list of all the csv files
# found in the indicated folder
# @param,
#@ path, the path of the folder yu want to parth
# @ extension, extension of the file we search for, here "csv by default"
def getCSVFiles(path, extension = ".csv"):
    fileList = os.listdir(path)
    return [ csvFile for csvFile in fileList if csvFile.endswith( extension) ]

# function getFASTAFiles,
# this function return a list of all the fasta files
# found in the indicated folder
# @param,
#@ path, the path of the folder yu want to parth
# @ extension, extension of the file we search for, here "fasta by default"
def getFASTAFiles(path, extension = ".fasta"):
    fileList = os.listdir(path)
    return [ fastaFile for fastaFile in fileList if fastaFile.endswith( extension) ]




# log function
def writeLog(output_path, logToWrite):
    if not os.path.isfile(output_path + "log.txt"):
        f = open(output_path + "log.txt", "w")
    else:
        f = open(output_path + "log.txt", "a")
    date = datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
    f.write(str(date) + "\n" + logToWrite)
    f.close()

fastaFiles = getFASTAFiles("./ressources/raw_data")
csvFiles = getCSVFiles("./ressources/raw_data")
df =[]
print(csvFiles)
for csvFile in csvFiles:
    df=pd.read_csv("./ressources/raw_data/" + csvFile,sep=",")
    print("new")
    print(df.columns)
    print(df)

# print(df)