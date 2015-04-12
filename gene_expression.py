####################################################################
#Auteur : ROCHER Vincent										   #
#																   #
#But : Programme principal										   #
#	   Fait la jointure entre le taux d'expression				   #
#	   Et notre jeu d'intron									   #
####################################################################

import sys,getopt #Pour pouvoir faire passer des arguments en ligne de commande




def extract_tx(file_name,ens_output,creation): 
	ENS_id={}
	if creation== 1: #Si creation vaut 1, on veut un fichier de sortie avec tous les ID gene ensembl
		file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
		file_out=open(ens_output,"w")  # On créer le fichier de sortie qui contiendra la liste des genes ENS
		header = file_in.readline() 
		list_head = header.split("\t")
		list_head = list_head[1:] #Sans le premier élément
		for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
			content=line.split("\t") #On split le contenu dans une liste
			file_out.write(content[0]+"\n") # On écrit l'identfiant Ensembl du gene dans le fichier
			ID = content[0].replace('\n', '') # On enregistre l'identifiant Ensembl comme clé sans le saut de ligne
			ENS_id[ID]= content[1:] # On enregistre dans un dictionnaire chaque objet avec comme clé l'identifiant Ensembl pour le retrouver
		file_in.close() 

	else: # Sinon on se contente d'extraire les taux d'expression pour chaque gene
		file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
		header = file_in.readline() 
		list_head = header.split("\t")
		list_head = list_head[1:] #Sans le premier élément
		for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
			content=line.split("\t") #On split le contenu dans une liste
			ID = content[0].replace('\n', '') # On enregistre l'identifiant Ensembl comme clé sans le saut de ligne
			ENS_id[content[0]]= content[1:] # On enregistre dans un dictionnaire chaque objet avec comme clé l'identifiant Ensembl pour le retrouver
		file_in.close() 
	return(ENS_id,list_head)




def paste_file(file_name,file_out_name,ENS_id,list_head): # On défnit une fonction qui va extraire les coordonnées et les enregistrer au format BED

	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête
	file_out=open(file_out_name,"w") # On créer le fichier qui va contenir les données
	new_header = name[:-1]+"\t"+"\t".join(list_head)
	file_out.write(new_header)
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		new_line = line[:-1] #sans le saut de ligne
		ID = content[1].replace('\n', '')
		if ID in ENS_id:
			new_line += "\t"+"\t".join(ENS_id[ID])
			file_out.write(new_line)
		else:
			new_line+="\tNA" * len(list_head) # On ecrit NA autant de fois qu'il y a d'éléments dans la liste
			new_line+="\n"
			file_out.write(new_line)
	file_in.close()
	file_out.close()




####Début du programme ####
###Valeur par defaut creation qui détermine si on veut un fichier de sorti qui recense tous les identifiants ENS
creation = 0
ens_output = None
####Définition des fichiers passés en arguments####
opts, args = getopt.getopt(sys.argv[1:],'i:o:',['index=','ens_list='])
for elmts in opts:
	if elmts[0] == '-i':
		input_file = elmts[1] # fichier d'entré  contenant les données, initialement HumanExpression.cRPKM_perGene_UB141209.tab

	elif elmts[0] == '-o':
		output_file = elmts[1] # fichier de sortie qui contiendra la jointure entre l'index et les données

	elif elmts[0] == '--index':
		index = elmts[1] # fichier qui contient la référence entre l'id ENS du gene et l'id de chaque intron

	elif elmts[0] == '--ens_list':
		creation = 1
		ens_output = elmts[1]



####Lancement de la fonction extract_tx qui va récupérer le taux d'expression de chaque gene dans tous les tissus
print('Extraction des données d\'expression de chaque gene dans tous les tissus')
ENS_id,list_head = extract_tx(input_file,ens_output,creation)
print("Ok !")
####Lancement de la fonction paste_file qui va joindre ce taux d'expression à chaque identifiant d'intron contenu dans "l'index"
message = 'Jointure entre les deux fichiers '+input_file+' et '+index+' dans '+output_file
print(message)
paste_file(index,output_file,ENS_id,list_head)
