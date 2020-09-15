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
    nc: int
        The identifier number of the organism
    sequence: str
        The genome sequence of the organism
    description: str
        The name and description of the organism

    Methods
    -------
    createByJSON(file)
        Create a Organism object with the information
        contained in the file given (.JSON)
    createByFASTA(file)
        Create a Organism object with the information
        contained in the file given (.fasta)
    exportToFASTA()
        Create a file with extension .fasta (FASTA file)
        with the information of the organism
    exportToJSON()
        Create a file with extension .json
        with the information of the organism

    """
    def __init__(self, nc=None, sequence=None, description=None):
        """
        Parameters
        ----------
        nc : int
            The identifier number
        sequence : str
            The complete gnome
        description : str, optional
            The name or the description of the organism
        """
        if (nc or sequence or description):    
            if not isinstance(nc, float):
                raise TypeError("The NC must be a float number")
            elif not isinstance(sequence, str):
                raise TypeError("The sequence must be a string")
            elif not isinstance(description, str):
                raise TypeError("The description must be a string")
        self.nc = nc
        self.sequence = sequence
        self.description = description
    
    def __eq__(self, anorganism):
        """
        From another organism (anorganism) compares
        the nc number, if is equal the result is true.

        Parameters
        ----------
        anorganism : Organism object
            The organism to compare
        """
        if isinstance(anorganism, Organism) == False:
            raise TypeError("Just can compare organisms with organisms")
        if self.nc == anorganism.nc:
            return True
        else:
            return False

    def printOrganism(self):
        print("Description: "+self.description)
        print("NC: "+str(self.nc))
        print("Sequence: "+self.sequence)

    def createByJSON(self, namefile, collection, nc, sequence, description):
        with open(namefile) as handle:
            data = json.load(handle)
            for p in data[collection]:
                self.nc = float(p[nc])
                self.sequence = p[sequence]
                self.description = p[description]
                return

    def createByFASTA(self, namefile):
        nc_pattern = re.compile(r'(?<=NC_)[0-9]*\.[0-9]*')
        with open(namefile, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                nc = re.search(nc_pattern, record.id)
                if nc:
                    self.nc = float(nc.group(0))
                else:
                    print("NC not found")
                    return
                self.description = record.description
                self.sequence = str(record.seq)
                return
    
    def createByCLUSTAL(self, namefile):
        nc_pattern = re.compile(r'(?<=NC_)[0-9]*\.[0-9]*')
        with open(namefile, "rU") as handle:
            for record in SeqIO.parse(handle, 'clustal'):
                nc = re.search(nc_pattern, record.id)
                if nc:
                    self.nc = float(nc.group(0))
                else:
                    print("NC not found")
                    return
                self.description = record.description
                self.sequence = str(record.seq)
                return
    
    def exportToCLUSTAL(self, route=None):
        if not (route):
            route="NC_"+str(self.nc)+".aln"
        record = SeqRecord(
            Seq(self.sequence),
            id="NC_"+str(self.nc),
            name=self.description,
            description=self.description
        )
        sequences = []
        sequences.append(record)
        with open(route, "w") as output_handle:
            SeqIO.write(sequences, output_handle, "clustal")
        return

    def exportToFASTA(self, route=None):
        if not (route):
            route="NC_"+str(self.nc)+".fasta"
        record = SeqRecord(
            Seq(self.sequence),
            id="NC_"+str(self.nc),
            name=self.description,
            description=self.description
        )
        sequences = []
        sequences.append(record)
        with open(route, "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")

    def exportToJSON(self, route=None):
        if not (route):
            route = "NC_"+str(self.nc)+".json"
        data = {}
        data['organisms'] = []
        data['organisms'].append({
            'nc': str(self.nc),
            'sequence': self.sequence,
            'description': self.description
        })
        with open(route, 'w') as outfile:
            json.dump(data, outfile)