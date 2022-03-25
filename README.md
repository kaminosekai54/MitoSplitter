# data-compiler

################################################################################

Installation
To run this program you will need the following package :

-- pandas
-- biopython
-- re

First make sure you have pip installed, or install it  by following the instruction here :
https://pip.pypa.io/en/stable/installation/
To install them via pip
 pip open a command prompt ant type :
 pip install pandas
 pip install biopython
 pip install re
 
 
 If you are using conda  :
 conda install pandas
 conda install biopython
 conda install re
 
 
 This program also use  muscle 5.1 	ns mafft 7.5. They are included in the tools folder, but feel free to réinstall them if you enconter any issue.
 
On macOS and linux , you might have to allow execution with the following command :
chmod +x tools/mafft/mafft-mac/mafft.bat # change for the apropriate folder on linux
chmod +x tools/muscle/muscle5.1.macos_intel64 
or
chmod +x tools/muscle/muscle5.1.linux_intel64

Please don't change the name of the executable file, they are used in the script to find the command line program.
If you have those software already installed, feel free to change their location in the settings.py file
 
 ################################################################################
Utilisation :

To use this program, simply run the main.py file by double clicking on it or run 
py main.py
in a command prompt (At the same location of the file).



To use it correctly, you will need atleest two file or one if it's a .gb (genbank) file
An anotation file (.csv) and your mitogenome file (.fasta) that contain only the full sequence of your mitogenome.

The anotation file and the mitogenome file will be paired together, so please make sure that they have the exact same name (Except the extension) and their name refer to the orgnasim they are sequencing as follows :
organism.csv (Anotation file)
organism.fasta (mitogenome file)

if you have a mitogenome sequenced in two part, please name it normaly, but followed by _number_indicating_the_part:
organism_1.csv
organism_1.fasta
organism_2.csv
organism_2.fasta

The name of the file will be used to name the sequence in the gene fasta file


if you are using .gb file, no worries, they will be automaticly rename correctly, you have also a stand alone rename tool in "tools" folder.

By default you need to put your list of paired file in the ressources/raw_data/ folder,
But feel free to edit the setting.py file to change the default location of your file.

Simply run the main.py file and let the magic do its work.
To get file from genbank, you can find in the tools folder a utility script "genbankDownloader.py" that will download a liste of file based on accessionID stored in a csv file.
The csv file must be called "list_accessionID.csv" with atleest a column named "accessionID" that contain your accessionID.
The csv file should be at the same place as the genbankDownloader.
Simply execute the file by double clicking on it or by the command line
python genbankDownloader.py
It will automaticly download the file in the row_data folder.

WARNING :  on macOS you might need to install a default ssl certificate with the command :
cd /Applications/Python 3.6/
./Install Certificates.command
Change the python version for the one you are using.
