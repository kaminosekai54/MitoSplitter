# data-compiler

################################################################################

Installation
To run this program you will need the following package :

-- pandas
-- biopython
-- re

First make sure you have pip installed, or install it  by following the instruction here :
https://pip.pypa.io/en/stable/installation/
To install them via 
 pip open a command prompt ant type :
 pip install pandas
 pip install biopython
 pip install re
 
 
 If you are using conda  :
 conda install pandas
 conda install biopython
 conda install re
 
 
 
 ################################################################################
Utilisation :

To use this program, simply run the main.py file by double clicking on it or run 
py main.py
in a command prompt (At the same location of the file).



To use it correctly, you will need atleest two file
An anotation file (.csv) and your mitogenome file (.fasta) that contain only the full sequence of your mitogenome.

The anotation file and the mitogenome file will be paired together, so please make sure that they have the exact same name (Except the extension) as follows :
my_file.csv (Anotation file)
my_file.fasta (mitogenome file)


By default you need to put your list of paired file in the ressources/raw_data/ folder,
But feel free to edit the setting.py file to change the default location of your file.

Simply run the main.py file and let the magic do its work.