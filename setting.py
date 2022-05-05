settings = {
"rawFilePath": r"./raw_data/",
"resultPath" : r"./results/",
"csvResultPath" : r"./results/csv/",
"classicFastaResultPath" : r"./results/classic-fasta/",
"genesFastaResultPath": r"./results/genes-fasta/",
"sequenceAlignementResultPath": r"./results/alignement/",
"treeResultPath": r"./results/tree/",
"musclePath": r"./tools/muscle/",
"mafftPath": r"./tools/mafft/",
"logPath": "./logs/",
"minColName" : "Minimum",
"maxColName" : "Maximum",
"nameColName" : "Name",
"typeColName" : "Type",
"useMuscle" : True,
"useMafft" : False,
"checkAlignement" : False,
"debugLog" : True,
"geneToDetect":[
    "12S", "16S", "18S", "28S", 
    "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "H3",
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", 
    "TRNA-ALA", "TRNA-ARG", "TRNA-ASN", "TRNA-ASP", 
    "TRNA-CYS", "TRNA-GLN", "TRNA-GLU", "TRNA-GLY", "TRNA-HIS", 
    "TRNA-ILE", "TRNA-LEU1", "TRNA-LEU2", "TRNA-LEU", "TRNA-LYS", 
    "TRNA-MET", "TRNA-PHE", "TRNA-PRO", "TRNA-SER1", "TRNA-SER2", "TRNA-SER", 
    "TRNA-THR", "TRNA-TRP", "TRNA-TYR", "TRNA-VAL", 
    ],
}

def getSettings():
    return settings