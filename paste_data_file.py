####################################################################
#Auteur : ROCHER Vincent										   #
#																   #
#But : Ajoute les données de GC rate + seq d'intéret			   #
#	   Au fichier principal contenant id + coordonnées génomique   #
#	   (intronsHuman.tab)										   #
####################################################################


import re #On importe re pour expression régulière


#Nouveaux fichier qui porte le nom de l'ancien avec "result_" devant
def paste_data(file_name1,file_out_name,dataset): #On créer une fonction qui va servir a coller les informations
	#On ouvre les fichiers
	position={} # Dictionnaire qui contiendra la position de chaque intron pour un id donné

	for key in dataset:
		regex = re.compile('^[^:]+:(.+):[0-9]{,2}')
		result = regex.findall(key)
		if result[0] in position:
			position[result[0]].append(key)

		else:
			position[result[0]]=[]
			position[result[0]].append(key)


	file_in=open(file_name1,"r") #En lecture
	file_out=open(file_out_name,"w") #En écriture pour créer le fichier
	for line in file_in: #On parcours le fichier ligne par ligne
		if re.match("Intron", line): #Si la ligne contient intron, c'est que c'est l'en tête
			#Format : "Intron	Coordinates	type	txRegion	PTCstatus	(-20/+30)	(-30/+20)	GCrate"
			newline = line[0:-1]+"\t(-20/+30)\t(-30/+20)\tIntron_length\tIntron_position\tGCrate\n" #On complète donc l'en tête avec les nouveaux titres
			file_out.write(newline) #On écris dans le fichier
		else: #Sinon c'est qu'on est dans une ligne informative
			content=line.split("\t") #On split le contenu par tabulation
			if content[0] in dataset: #Si l'id est contenu dans dataset (moins le signe), c'est qu'on a des données
				result = regex.findall(content[0]) # On reprend la même expression régulière pour déterminer le nombre d'introns
				number_intron=str(len(position[result[0]])) # On enregistre ce nombre dans une variable
				find_number = re.compile(':([0-9]{,2})$')
				resultat = find_number.search(content[0])
				num_intron = resultat.group(1)
				intron_position = num_intron+"/"+number_intron
				newcontent = line[0:-1]+"\t"+dataset[content[0]][1]+"\t"+dataset[content[0]][2]+"\t"+str(dataset[content[0]][3])+"\t"+intron_position+"\t"+str(dataset[content[0]][4])+"\n"
				#On concatene la ligne sans le saut de ligne avec les nouvelles informations
				file_out.write(newcontent) #On écris la nouvelle ligne dans le fichier
			else: #Si l'id n'est pas dans dataset, c'est qu'on a pas d'information dessus, on écris alors "NA"
				NAcontent = line[0:-1]+"\tNA\tNA\tNA\tNA\tNA\n"
				file_out.write(newcontent)
	#On ferme les deux fichiers
	file_in.close()
	file_out.close()