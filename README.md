# VirusPach
Application that provides the necessary tools to manage genetic sequences. Provides within its content a tool to compare databases based on FASTA, CLUSTAL and JSON files.
## _Version_
0.0.1 (Development)
## _Dependencies_
* You need install the BioPython library, more information:
    * _https://biopython.org/wiki/SeqIO_
* You need the json library, more information:
    * _https://docs.python.org/3/library/json.html_
## Instructions
In the **src** directory you can find the complete code of this project:

* The file _organism.py_ contains a class that provides an object
called Organism, this objects can uses to manipulate complete 
genetic sequences, from a collection of methods that read files FASTA, CLUSTAL and JSON.

* The file _testClass.py_ contains a program that, with the class
**Organism** compare a collection of files to updated them.

* The file _testQuery.py_ contains an example of how to do a
query on a web page (like the web page of the NCBI) for download
a collection of FASTA or CLUSTAL files.

* The file _testDatabase.py_ contains an example of how to do
a complete update from a collection of files recently downloaded
and a database that needs to be change, creating new organisms and replacing fields of the existing ones

**The test files must be run from src directory, for example:**
_/src$ python testClass.py_