import re

texto = "algo|NC_009874.12|:c0998"
patron_nc= re.compile(r'NC_[0-9]*\.?[0-9]*')
resultado = re.search(patron_nc, texto)
print(resultado.group(0))

"""
texto = "archivo.fasta"
patron_nc = re.compile(r'\.[a-z]*')
resultado = re.search(patron_nc, texto)
if resultado:
    print(resultado)
else:
    print("no hubo resultado")

"Hola"+"hola2"""
