"""Impor the main class"""
from organism import Organism

"""How create an Organism"""
covid1 = Organism() #Empty organism

#This doesn't work as an empty organism
#gripa = Organism(' ','',' ') 

#Full organism
covid10 = Organism('NC_12345.0', 12, 'AGSUJDJGGAS', 'Description of an organism') 

#Can't be half full
#covid10 = Organism(None, 'AGUHSHDS', 'Description') 

"""Create organism from file FASTA"""

# 1) Create an empty Organism
covid2 = Organism()

# 2) Use createByFasta
covid2.createByFASTA('./samples/sequence2.fasta') 
# The path is relative to the directory in which it is run
# the program or can be an absolute path

# 3) Ready for instructions
covid2.printOrganism()

"""Export an organism to FASTA file"""
covid2.exportToFASTA()

"""Create organism from file JSON"""

# 1) Create an empty Organism
covid3 = Organism()

# 2) Use createByJSON
covid3.createByJSON('./samples/sequence1.json', 'organism', 
'nc', 'sequence', 'description')
# We need the name of the collection, the name
# of the NC attribute, the name of the sequence attribute and
# the name of the description attribute in the json file.

# 3) Ready for instructions
covid3.printOrganism()

"""Export an organism to JSON file"""
covid3.exportToJSON("hola.json")
# The route of the output file it's optional for every export method

"""Create organism from file CLUSTAL"""

# 1) Create an empty Organism
covid4 = Organism()

# 2) Use createByCLUSTAL
covid4.createByCLUSTAL('./samples/sequence3.aln')

# 3) Ready for instructions
covid4.printOrganism()

"""Export an organism to CLUSTAL file"""
covid4.exportToCLUSTAL()
