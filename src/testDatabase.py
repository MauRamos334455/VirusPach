from organism import Organism
import os
import re

# no --> new_organisms
# oo --> old_organisms 
# ro --> result_organisms

# 1) Create the route of the directory with all the new files
no_path = './'

# 2) Create the route of the directory with all the old files
# Also we can use the actual directory, where the program is executed
oo_path = os.getcwd()

# 3) Create arrays for the files
no_content = os.listdir(no_path)
oo_content = os.listdir(oo_path)

# 4) Create an array of all the new organisms
no = []
extension_pattern = re.compile(r'\.[a-z]*')
for file in no_content:
    extension = re.search(extension_pattern, file)
    filename = no_path+os.sep+file
    if (extension == '.fasta'):
        no = Organism.createManyByFASTA(filename, no)
    elif (extension == '.aln'):
        no = Organism.createManyByCLUSTAL(filename, no)
    elif (extension == '.json'):
        no = Organism.createManyByJSON(filename, no)
    else:
        continue

# 5) Create an array of all the old organisms
oo = []
extension_pattern = re.compile(r'\.[a-z]*')
for file in oo_content:
    extension = re.search(extension_pattern, file)
    filename = oo_path+os.sep+file
    if (extension == '.fasta'):
        oo = Organism.createManyByFASTA(filename, no)
    elif (extension == '.aln'):
        oo = Organism.createManyByCLUSTAL(filename, no)
    elif (extension == '.json'):
        oo = Organism.createManyByJSON(filename, no)
    else:
        continue

# 6) Create an array for the results of the comparison
ro = []
for i in oo:
    ro.append(i)

for i in no:
    for x in ro:
        if x == i:
            continue
        else:
            ro.append(i)

# 7) Watch the results
for x in ro:
    x.printOrganism()
