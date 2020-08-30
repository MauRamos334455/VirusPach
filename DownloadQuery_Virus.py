#Programa que descarga todas la información de números GI obtenidos de una consulta a Nucleotide

#Libreías necesarias
import urllib.request
import urllib.parse
import re

#Hacer consulta a la DB de nucleotide
query='("Viruses"[Organism]+OR+viruses[All+Fields])+AND+(viruses[filter]+AND+refseq[filter])'
base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
url = base+'esearch.fcgi?db=nucleotide&term='+query+'&usehistory=y'
f = urllib.request.urlopen(url)
search = f.read().decode('utf-8')

#Obtener todos los GI numbers con ayuda de expresiones regulares
#Obtención de NCID
patron = re.compile(r'<WebEnv>\S*</WebEnv>', re.A)
flag = patron.search(search)
ncid=flag.group()
ncid=ncid.split('<WebEnv>')
ncid=ncid[1].split('</WebEnv>')
ncid=ncid[0]
print("ID de la busqueda: "+ ncid)

#Obtención del total de resultados
patron = re.compile(r'<Count>\b\d+\b</Count>', re.A)
flag = patron.search(search)
count=flag.group()
count=count.split('<Count>')
count=count[1].split('</Count>')
count=count[0]
print("Resultados: "+ count)

#Obtención de la llave de la consulta
patron = re.compile(r'<QueryKey>\b\d+\b</QueryKey>', re.A)
flag = patron.search(search)
key = flag.group()
key = key.split('<QueryKey>')
key = key[1].split('</QueryKey>')
key = key[0]
print("Llave de busqueda: "+ key)

#Consultar y descargar la información de cada GI number
retmax=500
for retstart in range (0, int(count), retmax):
    efetch_url = base+'efetch.fcgi?db=nucleotide&WebEnv='+ncid
    efetch_url += '&query_key='+key+'&retstart='+str(retstart) 
    efetch_url += '&retmax='+str(retmax)+'&rettype=fasta&retmode=text'
    urllib.request.urlretrieve(efetch_url, './virus'+str(retstart)+'.fna')

