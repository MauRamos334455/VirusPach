from organism import Organism

covid10 = Organism(123, 'hola', 'hola')
covid20 = Organism(12, 'esto no', 'no lo se')

with open('123.json') as file:
    covid10.createByJSON(file, 'organisms', 'nc', 'sequence', 'description')
covid10.printOrganism()