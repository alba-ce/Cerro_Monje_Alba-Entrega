import os

from Bio import SeqIO



def parser(gb, output) :

	""" gb: genbank a parsear
	    output: archivo en el que se imprime multifasta

	    Parsea todas las proteinas de un archivo genbank de las que se tenga su traduccion
	    Las imprime en un multifasta que indica: proteinid product[organism] """

	# Lo que escribamos se append al archivo output si ya existe, si no se crea un archivo nuevo
	if os.path.exists(output) :
		file = open(output,'a+')
	else :
		file = open(output, 'w')
	
	
	with open(gb, 'r') as input_handle :
		for record in SeqIO.parse(input_handle, 'genbank') :
			pass
		for feature in record.features : 
			if feature.type == 'CDS' :
				try : # Si no tiene la traduccion pasamos 
					traduccion = feature.qualifiers['translation'][0]
					try : # por si no tiene un nombre para el producto
						product = feature.qualifiers['product'][0]
					except : 
						product = ' '

					file.write('>' + feature.qualifiers['protein_id'][0] + ' ' + product + ' [' + record.annotations['organism'] + ']\n')
					file.write(traduccion+'\n')
				except : 
					pass
	file.close()




def genbankParser(output, gb = None) :

	""" output: archivo en el que se imprime multifasta
	    gb: genbank que se desea parsear

	    Llama a parser() con el archivo gb si el usuario lo ha indicado
	    Si no, llama a parser() con todos los genbanks que haya en el directorio """	

	# Hacer directorio data y results si no esta hecho 
	if not os.path.exists('data') :
		os.mkdir('data')

	if not os.path.exists('results') :
		os.mkdir('results')

	# Si usuario no explicita un genbank, se usan todos los genbanks que haya en el directorio
	if gb == None :
		for file in os.listdir(os.getcwd()) :
			if file.endswith('.gbff') :
				print(file) # para facilitar al usuario el progreso del programa
				parser(file, output)
				os.rename(file, 'data/'+file) # mueve genbank a directorio data
	else :
		parser(gb, output)
		os.rename(gb, 'data/'+gb)

	os.rename(output, 'results/'+output)





