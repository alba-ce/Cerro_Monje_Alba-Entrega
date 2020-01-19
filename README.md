# Entrega final programación

Alba Cerro Monje

Desarrollado en Python 2.7

El paquete 'Entrega' contiene los módulos que se explican a continuación. El programa main.py es un ejemplo de uso de todas 
las funciones que contiene.

GENBANKPARSER.PY

Requiere los módulos os y Biopython.

Parsea tantos archivos genbanks como desee el usuario, produciendo un archivo multifasta con todas las proteínas de todos 
los genbanks.

Contiene la función genbankParser(output, gb = None)
- output: el nombre del archivo multifasta que se creará
- gb: el nombre del genbank que se desea parsear. Si no se da este argumento, se parsean todos los genbanks del directorio.
Esta función crea además, si no existen:
- un directorio data, donde introduce los genbanks parseado
- un directorio results, donde introduce el multifasta producido

BLASTMUSCLE.PY

Requiere blast, muscle y los módulos subprocess, sys y os. 

Contiene las siguientes funciones:

* blast(fasta, multifasta, id = 0, cov = 0, evalue = 1e-06) 
  - fasta: query del blast
  - multifasta: subject
  - id: porcentaje de identidad
  - cov: coverage de la query
 Esta funcion realiza una busqueda blast de la query sobre el subject filtrando con los parámetros proporcionados.
 
 Devuelve el resultado en un tsv en formato 6: qseqid sseqid pident qcovs evalue.
 
 Crea además una lista con todos los id de las proteínas con las que hace hit la query. En el caso de uan búsqueda multiquery, 
 devuelve una lista de listas, una para cada query.
 
 Crea además, si no existen:
 - un directorio data, donde introduce los genbanks parseados
 - un directorio results/blast, donde introduce el tsv resultado del blast
 
* blast2Fasta(sseqid, query, multifasta) 
   - sseqid: lista con los id de las proteínas que han hecho hit
   - query: fasta con secuencia(s) query
   - multifasta: fasta con las secuencias subject
  Crea un archivo fasta para cada query con las secuencias de la query y con las que ha hecho hit
 
  Crea, si no existe, un directorio results/blast, donde introduce dichos fasta
  
* align (file, output = '') 
  - file: fasta con las secuencias a alinear
  - output: nombre de la query, para identificar el archivo producido
  Crea un alineamiento con el metodo Muscle
  Crea, si no existe, un directorio results/alineamientos, en el que introduce el alienamiento

* makeTreeNJ(alin, output = '') 
  - alin: alineamiento
  - output: nombre de la query, para identificar el archivo producido
  Crea con Muscle un arbol neighbour joining en formato phy.
 
  Crea, si no existe, un directorio results/arboles, en el que introduce el arbol.

PROSITE.PY

Requiere los módulos os, Biopython y re. También requiere los archivos prosite.dat y prosite.doc en el directorio del programa principal.

Contiene las siguientes funciones: 

* repl(re) 
  - re: patrón de Prosite
  Traduce un patrón de Prosite a una expresión regular que pueda entender el módulo re.

* findDomains(multifasta, output = '') 
  - multifasta: fasta con todas las proteínas en las qeu se van a buscar dominios
  - output: nombre de la query, para identificar el archivo producido
  Produce un archivo de texto con la siguiente información para cada proteína:
  - nombre de la proteína
  - patrones encontrados
  - nombre de los patrones
  - numero de accesion del patron
  - descripcion del patron
  Al final del todo, se incluye una breve descripcion bibliografica de todos los dominios encontrados.
  
 
