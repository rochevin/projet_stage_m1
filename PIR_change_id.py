####################################################################
#Auteur : ROCHER Vincent										   #
#																   #
#But : Modifie l'identifiant dans fichier PIR					   #
####################################################################

import re #On importe re pour expression régulière
import sys,getopt



def extract_id(file_name): 
	intron_id={}
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	file_in.readline() 
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		find_id = re.compile('^[^:]+:(.+)$')
		resultat = find_id.search(content[0])
		key_id = resultat.group(1)
		intron_id[key_id]=content[0] #Le dictionnaire contient en clé l'id au même format que le fichier PIR et en valeur l'id au format voulu
	file_in.close() 
	return(intron_id)

def change_name_id(file_name_in,file_name_out,intron_id):
	file_in=open(file_name_in,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() 
	file_out=open(file_name_out,"w")
	file_out.write(name)
	for line in file_in:
		content=line.split("\t") #On split le contenu dans une liste
		if content[0] in intron_id:
			content[0]=intron_id[content[0]]
			sep = "\t"
			new_line = sep.join(content)
			file_out.write(new_line)



####Début du programme ####
####Définition des fichiers passés en arguments####
opts, args = getopt.getopt(sys.argv[1:],'i:o:',['index='])
for elmts in opts:
	if elmts[0] == '-i':
		input_file = elmts[1] # fichier d'entré  contenant les PIR pour chaque introns avec un format id WASH7P:NR_024540:10

	elif elmts[0] == '-o':
		output_file = elmts[1] # fichier de sortie qui contiendra le nouveau fichier avec les ids corrigés

	elif elmts[0] == '--index':
		index = elmts[1] # fichier qui contient le bon id au format EIE000001:WASH7P:NR_024540:10

####Lancement de la fonction qui va récupérer tous les ids d'introns####
print('Extraction de chaque id ...')
intron_id = extract_id(index)
print("Ok !")
####Lancement de la fonction change_name_id qui va remplacer l'ancien id qui ne convient pas par le nouveau####
print('Remplacement des ids par les nouveaux')
change_name_id(input_file,output_file,intron_id)
print('Ok !')
message = "Résultats disponibles dans le nouveau fichier "+output_file
print(message)