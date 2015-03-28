####################################################################
#Auteur : ROCHER Vincent										   #
#																   #
#But : Programme principal										   #
#	   Execute tous les autres									   #
####################################################################


from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq
import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
import tempfile # On importe tempfile pour créer des fichiers temporaires
# On importe chaque fonction issue des autres scripts


class IntronInfoB(object):
	"Classe contenant les données de base contenus dans TableS6_IntronsHuman.tab"
	def __init__(self,Intron_id,Intron_coord,Intron_type,Intron_txRegion,Intron_PTCstatus,nb_nt_exon):
		nb_nt_exon = int(nb_nt_exon)-1 # On veut connaitre le nombre de bp à prendre en compte dans les exons flanquants -1 (car ce sont les coordonnées d'exon qui sont données)
		self.id = Intron_id # Contient l'identifiant de l'intron
		self.coords = Intron_coord # Contient les coordonnées de l'intron

		#Expression régulière qui récupère le nom du chromosome, 
		#les coordonnées de la dernière base de l'exon au site d'épissage 5'
		#les coordonnées de la première base de l'exon au site d'épissage 3'
		regex = re.compile('^(chr[\w]+)([-\+])[0-9]+:([0-9]+)_([0-9]+):[0-9]+')
		#On enregistre les résultats dans une liste (deux dimensions)
		result = regex.findall(self.coords)
		# Enregistre dans variables les coordonnées à prendre en compte pour retrouver la séquence fasta (globalement coordonnées introns +/- nb_nt_exon bp des exons flanquants)
		self.chr, self.signe, self.debut, self.fin = [result[0][0],result[0][1],int(result[0][2])-nb_nt_exon,int(result[0][3])+nb_nt_exon] 



		self.type = Intron_type # Contient le type de l'intron
		self.txregion = Intron_txRegion # Contient txRegion
		self.PTC = Intron_PTCstatus # Content le status PTC

	def formating_coord(self): # Formate les coordonnées pour retrouver la séquence fasta (le nom de la séquence est sous ce format)
		return self.chr+":"+str(self.debut)+"-"+str(self.fin)
	def BED_coord(self): # Formate les coordonnées au format BED
		return self.chr+"\t"+str(self.debut)+"\t"+str(self.fin)+"\t"+self.id

class IntronInfoAnnotation(IntronInfoB):
		"Classe contenant les informations obtenues après parsing pour chaque Intron"

		def __init__(self,Intron_id,Intron_coord,Intron_type,Intron_txRegion,Intron_PTCstatus,nb_nt_exon,seq):
			IntronInfoB.__init__(self,Intron_id,Intron_coord,Intron_type,Intron_txRegion,Intron_PTCstatus,nb_nt_exon)
			self.seq = seq # Contient la séquence d'intéret, intron + nb_nt_exon bases des exons flanquants (déjà reverse complement si signe == -)
		# Calcul du taux de GC à partir de la séquence fasta
		def GCrate(self,total_len,nb_GC):
			if self.seq != "NA":
				# Détermination des bp à prendre en compte pour faire le calcul de GC : les bords de la séquence minimale de l'intron moins les nucleotides d'intéret
				intron_seq = self.seq[total_len:-total_len]
				A = intron_seq.count('A')
				C = intron_seq.count('C')
				G = intron_seq.count('G')
				T = intron_seq.count('T')
				# On détermine une valeur minimale de nucleotides pour calculer le taux de GC
				if ((A+T+C+G)>=nb_GC):
					#Calcul du taux de GC
					GCrate = ((G+C)/(A+T+G+C))*100
					GCrate = round(GCrate,2)
				else:
					GCrate = "NA"
				return GCrate
			else:
				return "NA"
		# Détermine nos séquences d'intéret, par exemple (-20/+30) et (-30/+20)
		def seq_begin(self,total_len):
			if self.seq != "NA":
				return self.seq[0:total_len]
			else:
				return "NA"
		def seq_end(self,total_len):
			if self.seq != "NA":
				taille = len(self.seq)
				return self.seq[taille - total_len :taille]
			else:
				return "NA"
		# Détermine la taille de l'intron pour affichage dans fichier final
		def intron_length(self,nb_nt_exon):
			if self.seq != "NA":
				return len(self.seq[nb_nt_exon:-nb_nt_exon])
			else:
				return "NA"
		# Détermine le numéro de l'intron pour connaitre sa position par rapport aux autres
		def intron_num(self):
			find_number = re.compile(':([0-9]{,3})$')
			resultat = find_number.search(self.id)
			return resultat.group(1)


# On défnit une fonction qui va extraire les coordonnées et les enregistrer au format BED
def extract_coord(file_name,nb_nt_exon = 20,nb_nt_intron = 30,nb_GC = 40,len_seq = 50): 
	retrieve_name = {} # dictionnaire qui va contenir les objets issues de IntronInfoB
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête "Intron	Coordinates	type	txRegion	PTCstatus"

	file_out = tempfile.NamedTemporaryFile(mode="w+",delete=False,suffix='.bedPos')
	file_out_name = file_out.name
	#file_out=open(file_out_name,"w") # On créer le fichier BED qui va contenir les données
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		intron_name = content[0] #On récupère le nom de l'intron
		intron_name = IntronInfoB(content[0],content[1],content[2],content[3],content[4],nb_nt_exon)
		
		retrieve_name[intron_name.id]=intron_name
		# Écriture du fichier BED
		file_out.write(intron_name.BED_coord())
		file_out.write("\n")
	file_in.close() 
	file_out.close()
	return(retrieve_name,file_out_name)



#On définit une fonction pour récupérer les informations : GC rate (-nb_nt_exon/+nb_nt_intron) (-nb_nt_intron/+nb_nt_exon)
def extract_fasta_info(filename,retrieve_name,nb_nt_exon = 20,nb_nt_intron = 30,nb_GC = 40,len_seq = 50): 
	
	seq_info = {} # On définit un dictionnaire qui contient la séquence en valeur et ses coordonnées en clé
	dataset={} # On définit un dictionnaire qui contiendra les données
	
	handle = open(filename, "rU") #On ouvre le fichier en mode lecture

	for record in SeqIO.parse(handle, "fasta") : #On parse le fichier avec seqIO
		seq_id = record.id #On enregistre le nom de la sequence, ici les coordonnées
		seq = record.seq #On enregistre la sequence fasta dans une variable
		seq= seq.upper() #On met en majuscule au cas ou la séquence ne l'est pas
		seq_info[seq_id]=seq # On enregistre la séquence dans le dictionnaire avec ses coordonnées comme clé
		
	handle.close()
	# On parcours notre dictionnaire contenant les identifiants de l'intron et ses annotations
	# Va enregistrer la séquence correspondant à l'intron à partir de ses coordonnées et le stocker dans un nouvel objet de type IntronInfoAnnotation
	for intron_id in retrieve_name.values():
		coords = intron_id.formating_coord()
		seq = seq_info[coords]
		if len(seq)>=len_seq:
			if intron_id.signe == "-":
				seq = seq.reverse_complement()
			intron_id = IntronInfoAnnotation(intron_id.id,intron_id.coords,intron_id.type,intron_id.txregion,intron_id.PTC,nb_nt_exon,str(seq))
		else:
			intron_id = IntronInfoAnnotation(intron_id.id,intron_id.coords,intron_id.type,intron_id.txregion,intron_id.PTC,nb_nt_exon,"NA")
		dataset[intron_id.id]=intron_id
	return(dataset)



# Fonction qui va joindre les données de bases à nos annotations dans un seul fichier
def paste_data(file_name1,file_out_name,dataset,nb_nt_exon = 20,nb_nt_intron = 30,nb_GC = 40,len_seq = 50):

	total_len = int(nb_nt_exon)+int(nb_nt_intron) # Détermine les bords de la séquence minimale de l'intron moins les nucleotides d'intéret

	#On ouvre les fichiers
	position={} # Dictionnaire qui contiendra la position de chaque intron pour un id donné
	#On parcours les identifiants pour déterminer le nombre d'introns présent dans notre jeu de données dans chaque gene
	regex = re.compile('^[^:]+:(.+):[0-9]{,3}')
	for key in dataset:
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
			#Format : "Intron	Coordinates	type	txRegion	PTCstatus	(-20/+30)	(-30/+20)	Intron_length	Intron_position 	Total_Intron 	GCrate"
			newline = line[0:-1]+"\t(-"+str(nb_nt_exon)+"/+"+str(nb_nt_intron)+")\t(-"+str(nb_nt_intron)+"/+"+str(nb_nt_exon)+")\tIntron_length\tIntron_position\tTotal_Intron\tGCrate\n" #On complète donc l'en tête avec les nouveaux titres
			file_out.write(newline) #On écris dans le fichier
		else: #Sinon c'est qu'on est dans une ligne informative
			content=line.split("\t") #On split le contenu par tabulation
			if content[0] in dataset: #Si l'id est contenu dans dataset (moins le signe), c'est qu'on a des données
				result = regex.findall(content[0]) # On fait l'expression régulière pour déterminer le nombre d'introns
				total_intron = str(len(position[result[0]])) # On enregistre ce nombre dans une variable
				intron_position = dataset[content[0]].intron_num()
				newcontent = line[0:-1]+"\t"+dataset[content[0]].seq_begin(total_len)+"\t"+dataset[content[0]].seq_end(total_len)+"\t"+str(dataset[content[0]].intron_length(nb_nt_exon))+"\t"+intron_position+"\t"+total_intron+"\t"+str(dataset[content[0]].GCrate(total_len,nb_GC))+"\n"
				file_out.write(newcontent)
			else: #Si l'id n'est pas dans dataset, c'est qu'on a pas d'information dessus, on écris alors "NA"
				NAcontent = line[0:-1]+"\tNA\tNA\tNA\tNA\tNA\tNA\n"
				file_out.write(NAcontent)
	#On ferme les deux fichiers
	file_in.close()
	file_out.close()

######DEBUT DU PROGRAMME######

###Parsing des éléments en argument ####
opts, args = getopt.getopt(sys.argv[1:],'i:g:o:',['exon_len=','intron_len=','len_GC=','len_seq='])
for elmts in opts:
	if elmts[0] == '-i':
		input_file = elmts[1] # nom du fichier d'origine (input)

	elif elmts[0] == '-o':
		output_file = elmts[1] # nom du fichier de sortie

	elif elmts[0] == '-g':
		gen_ref = elmts[1] # nom du fichier génome de référence (binaire)

	elif elmts[0] == '--exon_len':
		nb_nt_exon = int(elmts[1]) # détermine le nombre de nucleotides à prendre en compte dans les exons flanquants

	elif elmts[0] == '--intron_len':
		nb_nt_intron = int(elmts[1]) # détermine le nombre de nucleotides à prendre en compte aux extrémités de l'intron

	elif elmts[0] == '--len_GC':
		nb_GC = int(elmts[1]) # détermine la longueur minimale de l'intron -(nb_nt_exon+nb_nt_intron) pour calculer le taux de GC

	elif elmts[0] == '--len_seq':
		len_seq = int(elmts[1]) # détermine la taille minimale de la séquence pour annoter


####Lancement des fonctions une par une####

####Lancement de la fonction extract_coord
####Permet de créer un fichier temporaire au format BED pour twobittofa et un dictionnaire d'objet contenant les informations issues de la publication
print("1 - Création du fichier temporaire BED")
retrieve_name,file_out_name = extract_coord(input_file,nb_nt_exon,nb_nt_intron,nb_GC,len_seq)
####On récupère le dictionnaire et le nom du fichier temporaire
print('OK !')
####Création d'une commande terminal pour lancer twobitToFa en précisant le fichier BED temporaire et le genome de reference
print("2 - Lancement de twobitToFa -> création du fichier fasta (output.fa)")
command = "twobitToFa -bedPos -bed="+file_out_name+" "+gen_ref+" output.fa"
####Execution de la commande en ligne de commande
os.system(command)
print('OK !')
####On suprimme le fichier temporaire .bedPos généré
os.unlink(file_out_name)

####Lancement de la fonction extract_fasta_info
####On récupère le jeu de donné annoté avec la séquence en plus, les objets contenus dans le dictionnaire permettent d'agir dessus
print("3 - Récupération des séquences")
dataset =extract_fasta_info("output.fa",retrieve_name,nb_nt_exon,nb_nt_intron,nb_GC,len_seq)
print("Ok !")
####On suprimme le fichier fasta généré
os.unlink("output.fa")

message = "4 - Jointure du fichier "+input_file+" avec les données obtenues à partir de la séquence dans : "+output_file
print(message)
####Lancement de la fonction paste_data, celle-ci va générer un nouveau fichier qui contiendra les données de bases ainsi que nos données en jointure au format tabulé
paste_data(input_file,output_file,dataset,nb_nt_exon,nb_nt_intron,nb_GC,len_seq) 
print("Ok !")
message = "Résultats disponibles dans "+output_file
print(message)