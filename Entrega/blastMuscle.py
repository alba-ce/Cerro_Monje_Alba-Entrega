from subprocess import Popen, PIPE
import sys
import os


def blast(fasta, multifasta, id = 0, cov = 0, evalue = 1e-06) :

	""" Fasta: query
	    Multifasta: subject

	    Devuelve el resultado de hacer blast de la query contra el subject filtrado con los parametros id, cov y evalue
	    en un archivo blast_result.tsv 
	    Devuelve una lista de listas. Cada sublista esta compuesta por todos los sseqid con los que hace hit una misma query. 
		El primer elemento de esa lista es el qseqid """

	busqueda = ['blastp', '-query', fasta, '-subject', multifasta, '-outfmt', "6 qseqid sseqid pident qcovs evalue"]
	sseqid = [] # Lista con todos los ID, para luego hacer un fasta con todas las secuencias que hacen blast
	queries = [] # Lista con todos los ID de las queries, para luego hacer un fasta con cada query

	proceso = Popen(busqueda, stdout = PIPE, stderr = PIPE)
	blast = proceso.stdout.read()

	proceso.stdout.close()
	
	if not os.path.exists('data') :
		os.mkdir('data')

	if not os.path.exists('results/blast') :
		os.mkdir('results/blast')

	os.rename(fasta, 'data/'+fasta)


	prots = blast.split('\n') # Cada linea del blast es un componente de la lista prots

	file = open('results/blast/blast_result.tsv', 'w')
	file.write('>QueryID\tSubjectID\tIdent\tCoverage\tEvalue\n')
	for prot in prots :
		if prot != '' :
			campos = prot.split('\t')
			if float(campos[2]) > float(id) and float(campos[3]) > float(cov) and float(campos[4]) < float(evalue) : # filtrado	
				file.write(prot+'\n')
				# para hacer la lista
				if campos[0] in queries :
					for i in range(len(sseqid)) :
						if campos[0] == sseqid[i][0] :
							sseqid[i].append(campos[1]) # si ya tenemos una lista para esa query, metemos sseqid en la que corresponda
							break
				else :
					sseqid.append([campos[0], campos[1]]) # si no, hacemos una nueva lista
					queries.append(campos[0])

	file.close()
	return sseqid


def blast2Fasta(sseqid, query, multifasta) :

	""" sseqid: lista con todos los id de las proteinas con las que hace blast
	    query: contiene las secuencias query
	    multifasta: tiene las secuencias subject

	    Convierte resultado del blast a fasta """	
	
	if not os.path.exists('results/blast') :
		os.mkdir('results/blast')

	for i in range(len(sseqid)) :
		m = open(multifasta, 'r')	
		f = open('results/blast/blast_result_' + sseqid[i][0] + '.fasta','w') # creamos un fasta para cada query
		siguiente = False
		for line in m : # para cada linea del multifasta con todas las prots
			if line.startswith('>') :
				campos = line.split(' ')
				if campos[0].replace('>','') in sseqid[i] : # si el identificador corresponde con alguno de la lista
					f.write(line) # guardamos la cabecera y la siguiente linea (la proteina)
					siguiente = True
				else :
					siguiente = False

			elif siguiente == True :
				f.write(line)
				siguiente = False
		m.close()		
	
		siguiente = False
		q = open(query, 'r')
		for line in q : # para guardar la query tambien
			if line.startswith('>') and line.replace('>','').replace('\n', '') == sseqid[i][0] :
				f.write(line)
				siguiente = True
			elif siguiente == True :
				f.write(line)
				siguiente = False
		q.close()
		f.close()

	



def align(file, output = '') :
	""" File: fasta con las secuencias a alinear
	    output: nombre de la query que vamos a alinear """

	if not os.path.exists('results/alineamientos') :
		os.mkdir('results/alineamientos')

	output_file = str('results/alineamientos/alineamiento_' + output + '.fasta') # hacemos un nuevo archivo para cada query
	alin = [ 'muscle', '-in', file, '-out', output_file ]

	proceso = Popen(alin, stdout = PIPE, stderr = PIPE)
	muscle = proceso.stdout.read()
	error = proceso.stderr.read()

	proceso.stderr.close()
	proceso.stdout.close()



def makeTreeNJ(alin, output = '') :
	""" Alin: alineamiento
	    Output: nombre de la query """

	if not os.path.exists('results/arboles') :
		os.mkdir('results/arboles')

	output_file = str('results/arboles/NJtree_' + output + '.phy') # Un archivo distinto para cada query
	tree = [ 'muscle', '-maketree', '-in', alin, '-out', output_file, '-cluster', 'neighborjoining' ]

	proceso = Popen(tree, stdout = PIPE, stderr = PIPE)
	muscle = proceso.stdout.read()
	error = proceso.stderr.read()

	proceso.stdout.close()
	proceso.stderr.close()


