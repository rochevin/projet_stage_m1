####################################################################
#Auteur : ROCHER Vincent										   #
#																   #
#But : Programme principal										   #
#	   Associe les données de Braunschweig + sequence intron	   #
#	   à l'id de chaque gene Ensembl							   #
####################################################################

from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq
import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
import tempfile # On importe tempfile pour créer des fichiers temporaires

class IntronInfoAnnotate(object):
	"Classe contenant les données des introns annotées générées via script précédent"
	def __init__(self,Intron_id,Intron_coord,Intron_type,Intron_txRegion,Intron_PTCstatus,seq_annotation_begin,seq_annotation_end,intron_len,intron_position,total_intron,GCrate):
		self.id = Intron_id # Contient l'identifiant de l'intron
		self.coords = Intron_coord # Contient les coordonnées de l'intron
		#Expression régulière qui récupère le nom du chromosome, 
		#les coordonnées de la dernière base de l'exon au site d'épissage 5'
		#les coordonnées de la première base de l'exon au site d'épissage 3'
		regex = re.compile('^(chr[\w]+)([-\+])[0-9]+:([0-9]+)_([0-9]+):[0-9]+')
		#On enregistre les résultats dans une liste (deux dimensions)
		result = regex.findall(self.coords)
		self.chr, self.signe, self.debut, self.fin = [result[0][0],result[0][1],int(result[0][2])+1,int(result[0][3])-1]


		self.type = Intron_type # Contient le type de l'intron
		self.txregion = Intron_txRegion # Contient txRegion
		self.PTC = Intron_PTCstatus # Contient le status PTC
		self.begin_annotation = seq_annotation_begin # Contient la séquence d'intéret de début
		self.end_annotation = seq_annotation_end # Contient la séquence d'intéret de fin
		self.len = intron_len # Contient la taille totale de l'intron
		self.pos = intron_position # Contient la position de l'intron par rapport aux autres dans le transcrit
		self.total = total_intron # Nombre total d'intron dans le transcrit
		self.GCrate = GCrate # Taux de GC de l'intron sans les séquences d'intéret

	def formating_coord(self):
		return self.chr+":"+str(self.debut)+"-"+str(self.fin)

	def flank_exon_left(self):
		regex = re.compile('^chr[\w]+[-\+]([0-9]+):([0-9]+)_[0-9]+:[0-9]+')
		result = regex.findall(self.coords)
		return result[0]
	def flank_exon_right(self):
		regex = re.compile('^chr[\w]+[-\+][0-9]+:[0-9]+_([0-9]+):([0-9]+)')
		result = regex.findall(self.coords)
		return result[0]


class IntronInfoComplete(IntronInfoAnnotate):
		"Classe contenant les informations des introns de la classe parente + seq + ID Ensembl"

		def __init__(self,Intron_id,Intron_coord,Intron_type,Intron_txRegion,Intron_PTCstatus,seq_annotation_begin,seq_annotation_end,intron_len,intron_position,total_intron,GCrate,seq,ensembl):
			IntronInfoAnnotate.__init__(self,Intron_id,Intron_coord,Intron_type,Intron_txRegion,Intron_PTCstatus,seq_annotation_begin,seq_annotation_end,intron_len,intron_position,total_intron,GCrate)
			self.seq = seq
			self.ensembl = ensembl

def extract_ensembl_id(file_name):
	IDlist = {} # Dictionnaire contenant en clé l'id de l'intron et en valeur l'id Ensembl
	list_Ensembl_ids = {} # Dictionnaire contenant tous les ids Ensembl en clé pour éviter la redondance
	file_in=open(file_name,"r")
	header = file_in.readline() # On enregistre l'en-tête
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t")
		IDlist[content[0]] = content[1].replace('\n', '')
		list_Ensembl_ids[content[1].replace('\n', '')] =""
	file_in.close()
	return(IDlist,list_Ensembl_ids)

def extract_coord(file_name): # On définit une fonction qui va extraire les annotations et enregistrer les coordonnées au format BED
	retrieve_name = {} # dictionnaire qui va contenir les objets issues de IntronInfoB
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête "Intron	Coordinates	type	txRegion	PTCstatus"

	file_out = tempfile.NamedTemporaryFile(mode="w+",delete=False,suffix='.bedPos')
	file_out_name = file_out.name
	#file_out=open(file_out_name,"w") # On créer le fichier BED qui va contenir les données
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		intron_name = content[0] #On récupère le nom de l'intron
		intron_name = IntronInfoAnnotate(content[0],content[1],content[2],content[3],content[4],content[5],content[6],content[7],content[8],content[9],content[10].replace('\n', ''))
		
		retrieve_name[intron_name.id]=intron_name
		# Écriture du fichier BED
		file_out.write(intron_name.chr)
		file_out.write("\t")
		file_out.write(str(intron_name.debut))
		file_out.write("\t")
		file_out.write(str(intron_name.fin))
		file_out.write("\t")
		file_out.write(intron_name.id)
		file_out.write("\n")
	file_in.close() 
	file_out.close()
	return(retrieve_name,file_out_name)

def extract_fasta_info(filename,retrieve_name,IDlist): #On définit une fonction pour récupérer la séquence fasta de chacun de nos introns
	
	seq_info = {} # On définit un dictionnaire qui contient la séquence en valeur et ses coordonnées en clé
	dataset={} # On définit un dictionnaire qui contiendra les données clé : id_intron, valeur = objet_intron
	
	handle = open(filename, "rU") #On ouvre le fichier en mode lecture

	for record in SeqIO.parse(handle, "fasta") : #On parse le fichier avec seqIO
		seq_id = record.id #On enregistre le nom de la sequence, ici les coordonnées
		seq = record.seq #On enregistre la sequence fasta dans une variable
		seq= seq.upper() #On met en majuscule au cas ou la séquence ne l'est pas
		seq_info[seq_id]=seq # On enregistre la séquence dans le dictionnaire avec ses coordonnées comme clé
		
	handle.close()

	for intron_id in retrieve_name.values():
		coords = intron_id.formating_coord() # On récupère les coordonées de façon formatés pour retrouver la sequence
		seq = seq_info[coords] # On va dans le dictionnaire qui contient les séquences et on récupère la séquence correspondant aux coordonnées de l'intron
		Ensembl_id = IDlist[intron_id.id] # On récupère l'id ensembl du gene de l'intron à partir de l'id de l'intron

		if len(seq)>=50:
			if intron_id.signe == "-":
				seq = seq.reverse_complement()
			intron_id = IntronInfoComplete(intron_id.id,intron_id.coords,intron_id.type,intron_id.txregion,intron_id.PTC,intron_id.begin_annotation,intron_id.end_annotation,intron_id.len,intron_id.pos,intron_id.total,intron_id.GCrate,str(seq),Ensembl_id)
		else:
			intron_id = IntronInfoComplete(intron_id.id,intron_id.coords,intron_id.type,intron_id.txregion,intron_id.PTC,intron_id.begin_annotation,intron_id.end_annotation,intron_id.len,intron_id.pos,intron_id.total,intron_id.GCrate,"NA",Ensembl_id)
		dataset[intron_id.id]=intron_id
	return(dataset)



IDlist,list_Ensembl_ids = extract_ensembl_id(sys.argv[1]) # Récupère l'id Gene Ensembl correspondant au gène de chaque intron
retrieve_name,file_out_name = extract_coord(sys.argv[2]) # Récupère les informations annotés issu de notre précédent parsing (new_main.py) dans un dico + fichier temporaire BED avec coordonnées de chaque intron annotés
command = "twobitToFa -bedPos -bed="+file_out_name+" "+sys.argv[3]+" output.fa" # Commande qui lance TwobitToFa pour récupérer la séquence de chaque coordonnée intronique
os.system(command) # Execute la commande précédente
os.unlink(file_out_name) # Suprimme le fichier BED temporaire
dataset = extract_fasta_info("output.fa",retrieve_name,IDlist) # Rajoute aux éléments annotés la séquence de l'intron ainsi que l'id Ensembl correspondant au gene dont fait parti l'intron
os.unlink("output.fa") # Suprimme le fichier fasta contenant les séquences des introns
print(list_Ensembl_ids)
