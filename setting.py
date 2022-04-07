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
"useMafft" : True,
}

def getSettings():
    return settings