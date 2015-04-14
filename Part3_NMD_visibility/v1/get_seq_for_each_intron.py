import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
import tempfile # On importe tempfile pour créer des fichiers temporaires


def extract_coord(file_name):
	coord_by_intron = {}
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête "Intron	Coordinates	type	txRegion	PTCstatus"

	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		coord_by_intron[content[0]]=content[1]
	file_in.close()
	return(coord_by_intron)

# On défnit une fonction qui va extraire les coordonnées et les enregistrer au format BED
def BED_creation(file_name,coord_by_intron): 
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() 

	file_out = tempfile.NamedTemporaryFile(mode="w+",delete=False,suffix='.bedPos')
	file_out_name = file_out.name
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		if content[1] in coord_by_intron:
			coords = coord_by_intron[content[1]]
			regex = re.compile("^(chr[^\+-]+)[\+-][0-9]+:([0-9]+)_([0-9]+):[0-9]+")
			result = regex.findall(coords)
			Bed_coord = result[0][0]+"\t"+str(int(result[0][1])-1)+"\t"+result[0][2]+"\t"+content[0]+"\n"
			# Écriture du fichier BED
			file_out.write(Bed_coord)
	file_in.close() 
	file_out.close()
	return(file_out_name)


opts, args = getopt.getopt(sys.argv[1:],'',['gen_ref=','output_fasta=','intron_list=','intron_braunch='])
for elmts in opts:
	if elmts[0] == '--gen_ref':
		gen_ref = elmts[1]
	elif elmts[0] == '--output_fasta':
		output_fasta = elmts[1]
	elif elmts[0] == '--intron_list':
		intron_file = elmts[1]
	elif elmts[0] == '--intron_braunch':
		intron_braunch = elmts[1]

coord_by_intron = extract_coord(intron_braunch)
print ("Extraction des introns et écriture au format BED")
file_out_name = BED_creation(intron_file,coord_by_intron)
print("Écriture du fichier fasta avec twobitToFa")
command = "twobitToFa -bedPos -bed="+file_out_name+" "+gen_ref+" "+output_fasta
####Execution de la commande en ligne de commande
os.system(command)
os.unlink(file_out_name)