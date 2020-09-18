from organism import Organism
import os
import re

# no --> new_organisms
# oo --> old_organisms 
# ro --> result_organisms

# 1) Create the route of the directory with all the new files
no_path = './samples/new'

# 2) Create the route of the directory with all the old files
# Also we can use the actual directory, where the program is executed
oo_path = './samples/old'

# 3) Create arrays for the files
no_content = os.listdir(no_path)
oo_content = os.listdir(oo_path)

# 4) Create an array of all the new organisms
no = []
extension_pattern = re.compile(r'\.[a-z]*')
for file in no_content:
    extension = re.search(extension_pattern, file)
    extension = extension.group(0)
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
    extension = extension.group(0)
    filename = oo_path+os.sep+file
    if (extension == '.fasta'):
        oo = Organism.createManyByFASTA(filename, oo)
    elif (extension == '.aln'):
        oo = Organism.createManyByCLUSTAL(filename, oo)
    elif (extension == '.json'):
        oo = Organism.createManyByJSON(filename, oo)
    else:
        continue

# 6) Create an array for the results of the comparison
ro = []
print("-----OLD ORGANISMS----")
for i in oo:
    i.printId()
    if i in ro:
        for x in ro:
            if i == x:
                if not x.version:
                    ro.remove(x)
                    ro.append(i)
                elif not i.version:
                    continue
                elif i.version > x.version:
                    ro.remove(x)
                    ro.append(i)
                else:
                    continue
    else:
        ro.append(i)

print("-----NEW ORGANISMS----")
for i in no:
    i.printId()
    if i in ro:
        for x in ro:
            if i == x:
                if not x.version:
                    ro.remove(x)
                    ro.append(i)
                elif not i.version:
                    continue
                elif i.version > x.version:
                    ro.remove(x)
                    ro.append(i)
                else:
                    continue
    else:
        ro.append(i)


# 7) Watch the results
print("---RESULT ORGANISMS----")
for i in ro:
    i.printId()

