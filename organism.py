import json

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
    createByJson(file)
        Create a Organism object with the information
        contained in the file given (.JSON)
    CreateByFasta(file)
        Create a Organism object with the information
        contained in the file given (.FNA)
    exportToFasta()
        Create a file with extension .fna (FASTA file)
        with the information of the organism
    exportToJson()
        Create a file with extension .json
        with the information of the organism

    """
    def __init__(self, nc, sequence, description='None'):
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
        if not isinstance(nc, int):
           raise TypeError("The NC must be a number")
        elif not isinstance(sequence, str):
           raise TypeError("The sequence must be a string")
        elif not isinstance(description, str):
            raise TypeError("The description must be a string")
        else:
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

    def createByJSON(self, file, collection, nc, sequence, description):
        data = json.load(file)
        for p in data[collection]:
            organism = Organism(int(p[nc]), p[sequence], p[description])
            return organism

    def exportToJson(self):
        data = {}
        data['organisms'] = []
        data['organisms'].append({
            'nc': str(self.nc),
            'sequence': self.sequence,
            'description': self.description
        })
        with open(str(self.nc)+'.json', 'w') as outfile:
            json.dump(data, outfile)