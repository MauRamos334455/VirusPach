import json
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Organism:
    """ 
    A class used to represented a sequence organism details.

    ...

    Attributes
    ----------
    id: str
        The identifier number of the organism
    version: int
        The version of the document from which it was obtained
    sequence: str
        The genome sequence of the organism
    description: str
        The name and description of the organism

    Methods
    -------
    createByJSON(filename, collection='organisms', nc='nc', 
    sequence='sequence', description='description')
        Create a Organism object with the information
        contained in the file given (.JSON)
    createByFASTA(filename)
        Create a Organism object with the information
        contained in the file given (.fasta)
    createByCLUSTAL(filename)
        Create a Organism object with the information
        contained in the file given (.fasta)
    createManyByFASTA(filename, organisms=[])
        Create a Organism object with the information
        contained in the file given (.fasta)
    createManyByCLUSTAL(filename, organisms=[])
        Create a Organism object with the information
        contained in the file given (.fasta)
    createManyByJSON(filename, organisms=[], collection='organisms', nc='nc', 
    sequence='sequence', description='description')
        Create a Organism object with the information
        contained in the file given (.fasta)
    exportToFASTA(route)
        Create a file with extension .fasta (FASTA file)
        with the information of the organism
    exportToJSON(route)
        Create a file with extension .json
        with the information of the organism
    printOrganism()
        Prints on screen all attributes of the self organism
    printID()
        Prints on screen the identification number of the self
        organism

    """
    def __init__(self, id=None, version=None, sequence=None, description=None):
        """
        Parameters
        ----------
        id : str
            The identifier number
        version : int
            The number version of the document 
            from which it was obtained
        sequence : str
            The complete gnome
        description : str, optional
            The name or the description of the organism
        """
        if (id or version or sequence or description):    
            if not isinstance(id, str):
                raise TypeError("The identifier must be a integer")
            elif not isinstance(version, int):
                raise TypeError("The version must be a string")
            elif not isinstance(sequence, str):
                raise TypeError("The sequence must be a string")
            elif not isinstance(description, str):
                raise TypeError("The description must be a string")
        self.id = id
        self.version = version
        self.sequence = sequence
        self.description = description
    
    def __eq__(self, anorganism):
        """
        From another organism (anorganism) compares
        the identifier number, if is equal the result is true.

        Parameters
        ----------
        anorganism : Organism
            The organism to compare
        """
        if isinstance(anorganism, Organism) == False:
            raise TypeError("Just can compare organisms with organisms")
        if (self.id == anorganism.id):
            return True
        else:
            return False

    def printOrganism(self):
        """
        Prints all information of the Organism on screen
        """
        print("Description: "+self.description)
        print("ID: "+self.id)
        print("Version: "+ str(self.version))
        #print("Sequence: "+self.sequence)
    
    def printId(self):
        """
        Prints just the id of the Organism on screen
        """
        if self.version:
            print("ID: "+self.id+'.'+str(self.version))
        else:
            print("ID: "+self.id)

    def createByJSON(self, filename, collection='organisms', id='id', 
    sequence='sequence', description='description'):
        with open(filename) as handle:
            data = json.load(handle)
            for p in data[collection]:
                aux = p[id].split('.')
                self.id = aux[0]
                try:
                    self.version = int(aux[1])
                except IndexError:
                    self.version = None
                self.sequence = p[sequence]
                self.description = p[description]
                return

    def createByFASTA(self, filename):
        id_pattern = re.compile(r'(?<=>)[A-Z_]*[0-9]*\.?[0-9]*')
        with open(filename, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                id = re.search(id_pattern, record.id)
                if id:
                    id = id.group(0)
                    id = id.split('.')
                    self.id = id[0]
                    try:
                        self.version = int(id[1])
                    except IndexError:
                        self.version = None
                else:
                    print("Identifier not found")
                    return
                self.description = record.description
                self.sequence = str(record.seq)
                return

    def createByCLUSTAL(self, filename):
        id_pattern = re.compile(r'(?<=>)[A-Z_]*[0-9]*\.?[0-9]*')
        with open(filename, "rU") as handle:
            for record in SeqIO.parse(handle, "clustal"):
                id = re.search(id_pattern, record.id)
                if id:
                    id = id.group(0)
                    id = id.split('.')
                    self.id = id[0]
                    try:
                        self.version = int(id[1])
                    except IndexError:
                        self.version = None
                else:
                    print("Identifier not found")
                    return
                self.description = record.description
                self.sequence = str(record.seq)
                return
    
    def exportToCLUSTAL(self, route=None):
        if not (route):
            route = self.id + ".aln"
        if self.version:
            version = '.'+str(self.version)
        else:
           version = ''
        record = SeqRecord(
            Seq(self.sequence),
            id = self.id+version,
            name = self.description,
            description = self.description
        )
        sequences = []
        sequences.append(record)
        with open(route, "a") as output_handle:
            SeqIO.write(sequences, output_handle, "clustal")
        return

    def exportToFASTA(self, route=None):
        if not (route):
            route = self.id + ".fasta"
        if self.version:
            version = '.'+str(self.version)
        else:
           version = ''
        record = SeqRecord(
            Seq(self.sequence),
            id = self.id+version,
            description = self.description
        )
        sequences = []
        sequences.append(record)
        with open(route, "a") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")
            output_handle.write("\n")

    def exportToJSON(self, route=None):
        if not (route):
            route = self.id + ".json"
        if self.version:
            version = '.'+str(self.version)
        else:
           version = ''
        data = {}
        data['organisms'] = []
        data['organisms'].append({
            'id': self.id+version,
            'sequence': self.sequence,
            'description': self.description
        })
        with open(route, 'a') as outfile:
            json.dump(data, outfile)
    
    @classmethod
    def createManyByFASTA(cls, filename, organisms=[]):
        id_pattern = re.compile(r'\A[A-Z_]*[0-9]*\.?[0-9]*')
        with open(filename, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                x = Organism()
                id = re.search(id_pattern, record.id)
                if id:
                    id = id.group(0)
                    id = id.split('.')
                    x.id = id[0]
                    try:
                        x.version = int(id[1])
                    except IndexError:
                        x.version = None
                else:
                    print("ID not found")
                    continue
                x.description = record.description
                x.sequence = str(record.seq)
                organisms.append(x)
        return organisms

    @classmethod
    def createManyByCLUSTAL(cls, filename, organisms=[]):
        id_pattern = re.compile(r'(?<=>)[A-Z_]*[0-9]*\.?[0-9]*')
        with open(filename, "rU") as handle:
            for record in SeqIO.parse(handle, "clustal"):
                x = Organism()
                id = re.search(id_pattern, record.id)
                if id:
                    id = id.group(0)
                    id = id.split('.')
                    x.id = id[0]
                    try:
                        x.version = int(id[1])
                    except IndexError:
                        x.version = None
                else:
                    print("ID not found")
                    continue
                x.description = record.description
                x.sequence = str(record.seq)
                organisms.append(x)
        return organisms

    @classmethod
    def createManyByJSON(cls, filename, organisms=[], collection='organisms', 
    id='id', sequence='sequence', description='description'):
        with open(filename) as handle:
            data = json.load(handle)
            for p in data[collection]:
                x = Organism()
                aux = p[id].split('.')
                x.id = aux[0]
                try:
                    x.version = int(aux[1])
                except IndexError:
                    x.version = None
                x.sequence = p[sequence]
                x.description = p[description]
                organisms.append(x)
        return organisms
