import os

from Entrega import genbank, blastMuscle, prosite

""" Inputs: 
	Tantos genbanks como quiera el usuario y un archivo fasta con proteinas query. 
    		Es necesario que esten en el mismo directorio que el script principal. 
	Parametros del blast (pident, cov y evalue)

    Genera un directorio data al que se mueven los inputs.
    Genera un directorio results que contiene un multifasta con todas las proteinas de todos los genbanks. Contiene ademas
	Un directorio blast, con archivos multifasta con todas las proteinas con las que hace hit cada query
	Un directorio alineamientos, con un alineamiento de Muscle por cada query
	Un directorio arboles, con un arbol Neighbour joining por cada query
	Un directorio prosite, con un archivo de texto que contiene informacion acerca de los dominios de cada query
"""



if __name__ == '__main__' :
	# 0. Datos
	# genbanks
	while True :
		print('Introduzca el nombre del genbank que desea parsear')
		gb = raw_input('Si no introduce ninguno, se parsearan todos los genbanks del directorio\n')
		if gb == '' or os.path.exists(gb) :
			break
		else :
			print('No existe ningun archivo con ese nombre')

	multifasta = raw_input('Introduzca el nombre para el archivo multifasta que se generara: ')

	# query de blast
	while True :
		query = raw_input('Introduzca el nombre del archivo con las queries para hacer blast: ')
		if os.path.exists(query) :
			break
		else :
			print('No existe ningun archivo von ese nombre')

	# parametros del blast. Si no se introduce nada, se pone el valor por defecto. Si no se introduce un numero, se vuelve a preguntar
	print('Parametros del blast:')
	while True :	
		pident = raw_input('\tLimite de identidad (0 por defecto): ')
		if pident == '' :
			pident = 0
		try :
			pident = float(pident)
			break
		except :
			print('Por favor, introduzca un numero.')

	while True :
		cover = raw_input('\tLimite de coverage (0 por defecto): ')
		if cover == '' :
			cover = 0
		try :
			cover = float(cover)
			break
		except :
			print('Por favor, introduzca un numero.')

	while True :
		eval = raw_input('\tLimite de evalue (1e-06 por defecto): ')
		if eval == '' : 
			eval = 1e-06
		try :
			eval = float(eval)
			break
		except :
			print('Por favor, introduzca un numero')



	#1. Multifasta a partir del genbank
	print('Creando un multifasta con todas las proteinas de los genbanks seleccionados...')
	if gb == '' : # si no se ha especificado uno, hacemos una lista con todos los genbanks para comprobar que hay alguno en el directorio
		list_gb = [] # lista con todos los archivos genbank
		for file in os.listdir('.') :
			if file.endswith('.gbff') :
				list_gb.append(file)
		if len(list_gb) == 0 : # si no hay ninguno se sale del programa para que el usuario meta gb en el directorio
			print('No hay ningun genbank en el directorio.')
			exit()	
		genbank.genbankParser(multifasta) # usa todos los genbanks del directorio 
	else :
		genbank.genbankParser(multifasta, gb = gb) # usa un genbank concreto
	print('Hecho')
	print(' ')



	# 2. Blast
	# tsv con todos los hits de todas las queries (id_query, id_subject, pident, cov, evalue)
	# devuelve ademas una lista (id) que contiene varias listas, una por cada query. Cada sublista contiene los identificadores
	# 	de las protenias con las que hace hit esa query.
	print('Haciendo blast...')
	id = blastMuscle.blast(query, 'results/'+multifasta, id = pident, cov = cover, evalue = eval) 

	# un multifasta para cada query con las secuencias con las que hace hit esa query
	blastMuscle.blast2Fasta(id, 'data/'+query, 'results/'+multifasta)
	print('Hecho')
	print(' ')



	# 3. Muscle y Prosite
	# alineamiento + busqueda de dominio: un alineamiento/resumen de dominios para cada query
	# se hacen a la vez porque funciones tienen el mismo input
	print('Alineando las secuencias...')
	print('Buscando dominios conservados...')
	for file in os.listdir('results/blast') :
		if file.endswith('.fasta') :
			query = file.replace('blast_result_', '').replace('.fasta', '') # nombre de la query que vamos a alinear ahora
			print(query)
			blastMuscle.align('results/blast/'+file, query)
			prosite.findDomains('results/blast/'+file, query)
	print(' ')

	# arboles filogeneticos: uno para cada query
	print('Creando arboles filogeneticos...')
	for file in os.listdir('results/alineamientos') :
		query = file.replace('alineamiento_', '').replace('.fasta','') # nombre de la query de la que vamos a hacer el arbol
		print(query)
		blastMuscle.makeTreeNJ('results/alineamientos/'+file, query)
	print(' ')
