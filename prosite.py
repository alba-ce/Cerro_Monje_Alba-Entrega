import os

from Bio.ExPASy import Prosite, Prodoc
import re



def repl(re) :
	""" Cambia expresiones regulares de prosite a las del modulo re de python"""

	mal = ['.', '-', '{', '}', '(', ')', '<', '>', 'x', '>]'] # La de prosite
	bien = ['', '', '[^', ']', '{', '}', '^', '$', '[GAVLIMPFYWSCTNQDEKRH]', ']?$'] # La buena

	for i in range(len(mal)) :
		re = re.replace(mal[i], bien[i])

	return re




def findDomains(multifasta, output = '') :
	""" Multifasta: archivo con todas las proteinas en las que se van a buscar dominios
	    Output: nombre de la query """

	file = open(multifasta, 'r')
	
	if not os.path.exists('results/prosite') :
		os.mkdir('results/prosite')
	output_file = str('results/prosite/dominios_' + output + '.txt') # un archivo para cada query
	result = open(output_file, 'w')
	accession_bruto = [] # todos los numeros de acceison de dominios encontrados en el multifasta
	accession = [] # lo mismo pero eliminando repeticiones

	for line in file : 
		if line.startswith('>') :
			result.write('*************************************************************************************************************'+'\n')
			result.write(line.replace('>', '') + '\n') # titulo: nombre de la proteina
		else :
			handle = open('prosite.dat', 'r')
			records = Prosite.parse(handle)
			for record in records :
				patron = repl(record.pattern) # traduccion patron 
				if len(patron) !=0 and re.search(patron, line) : # si existe el patron y se encuentra
					result.write('Patron: ' + record.pattern + '\nName: ' + record.name + '\nAccession: ' + record.accession + '\nDescription: ' + record.description + '\n\n')
					accession_bruto.append(record.accession) # guardamos info y el numero de accesion (necesario para buscar en prodoc)
	
	for a in accession_bruto : # para eliminar los repetidos
		if a not in accession :
			accession.append(a)

	result.write('\n\n\n\nINFORMACION DE LOS DOMINIOS\n\n')
	handle = open('prosite.doc', 'r')
	records = Prodoc.parse(handle)
	for record in records :
		if len(record.prosite_refs) != 0 and record.prosite_refs[0][0] in accession : # el numero de accesion de prosite esta en accesion
			result.write(record.text + '\n\nAccession prodoc: ' + record.accession + '\nAccession prosite: ' + record.prosite_refs[0][0] + '\n')
			result.write('**************************************************************************************\n\n\n')



