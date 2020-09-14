# VirusPach
Application that provides the necessary tools to manage genetic sequences. Provides within its content a tool to compare databases based on FASTA, CLUSTAL and JSON files.
## _Version_
0.0.1 (Development)
## _Dependencies_
* You need install the BioPyhon library, more information:
    *_https://biopython.org/wiki/SeqIO_
* You need the json library, more information:
    *_https://docs.python.org/3/library/json.html_
## Instructions
The file _organism.py_ contains a class that provides an object
called Organism, this objects can uses to manipulate complete 
genetic sequences, from a collection of methods that read files FASTA, CLUSTAL and JSON.

The file _test.py_ contains a program that, with the class
**Organism** compare a collection of files to updated them.