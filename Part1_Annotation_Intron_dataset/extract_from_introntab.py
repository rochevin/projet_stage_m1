####################################################################
#Auteur : ROCHER Vincent										   #
#																   #
#But : Extraire les coordonnées des introns de chaque genes		   #
#	   avec -20 +20 à chaque extrémité et enregistrer au format BED#
####################################################################


import re #On importe re pour expression régulière
import sys

def extract_coord(file_name,file_out_name): # On défnit une fonction qui va extraire les coordonnées et les enregistrer au format BED
	retrieve_name = {} # Dictionnaire qui va contenir les coordonnées en clé et l'id en valeur
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête "Intron	Coordinates	type	txRegion	PTCstatus"
	file_out=open(file_out_name,"w") # On créer le fichier qi va contenir les données
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		intron_name=content[0] #On récupère le nom de l'intron
		#Expression régulière qui récupère le nom du chromosome, 
		#les coordonnées de la dernière base de l'exon au site d'épissage 5'
		#les coordonnées de la première base de l'exon au site d'épissage 3'
		regex = re.compile('^(chr[\w]+)([-\+])[0-9]+:([0-9]+)_([0-9]+):[0-9]+')
		#On enregistre les résultats dans une liste (deux dimensions)
		result = regex.findall(content[1])
		#On écris sous la forme "chr	debut-20	fin+20	ID"
		file_out.write(result[0][0]) #chromosome
		file_out.write("\t") #tabulation
		debut= int(result[0][2])-19 #debut - 20
		file_out.write(str(debut)) #On ecris sans oublier de convertir en chaine de caractère
		file_out.write("\t") #tabulation
		fin = int(result[0][3])+19 #fin + 20
		file_out.write(str(fin)) #On ecris en str
		file_out.write("\t") #tabulation
		file_out.write(intron_name) #On écris l'id de l'intron
		file_out.write("\n")

		key = result[0][0]+":"+str(debut)+"-"+str(fin) # Conversion des coordonnées au même format que le nom de la seq
		value = intron_name+result[0][1]
		
		if key in retrieve_name:
			retrieve_name[key].append(value)

		else:
			retrieve_name[key]=[] # On precise à python que la valeur est une liste
			retrieve_name[key].append(value) #Sachant que le nom du fichier fasta sont les coordonnées de la sequence
		
		#on enregistre pour pouvoir retrouver l'id correspondant
		#et on enregistre le signe à la fin pour savoir si reverse ou forward
	#On ferme les deux fichiers et on renvoit le dictionnaire pour retrouver l'id
	file_in.close() 
	file_out.close()
	return(retrieve_name)

