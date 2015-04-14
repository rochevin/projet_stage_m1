# Programme qui va récupérer la liste des exons codants et des introns pour chaque transcrit canonique ensembl
# A partir de ces annotations, va organiser chaque transcrit canonique de façon à ce qu'il contienne tous ses exons/introns
# A partir de cette liste, va construire le CDS de chaque transcrit, puis introduire un à un chaque intron 
# Pour chaque intron, on déterminera si celui-ci est NMD visible, et si il décale la phase, on indiquera donc la phase de l'intron


import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq
#Définition des classes

class ExonInfo(object):
	"On définit un classe qui contiendra les annotations de chaque exon"
	def __init__(self,exon_id,trans_id,gene_id,coords,seq):
		self.id = exon_id
		self.coords = coords
		self.gene_id = gene_id
		self.trans_id = trans_id
		self.seq = seq
		self.get_type = "exon"

		regex = re.compile("(chr[^\+-]+)([\+-])([0-9]+):([0-9]+)")
		result = regex.findall(self.coords)
		self.chr = result[0][0]
		self.strand = result[0][1]
		self.start = result[0][2]
		self.end = result[0][3]
	def get_seq(self):
		if self.strand == "-":
			return self.seq.reverse_complement()
		else:
			return self.seq


class IntronInfo(object):
	"On définit un classe qui contiendra les annotations de chaque intron"
	def __init__(self,intron_id,trans_id,gene_id,coords,seq):
		self.id = intron_id
		self.coords = coords
		self.gene_id = gene_id
		self.trans_id = trans_id
		self.seq = seq
		self.get_type = "intron"

		regex = re.compile("(chr[^\+-]+)([\+-])([0-9]+):([0-9]+)")
		result = regex.findall(self.coords)
		self.chr = result[0][0]
		self.strand = result[0][1]
		self.start = result[0][2]
		self.end = result[0][3]
	def get_seq(self):
		if self.strand == "-":
			return self.seq.reverse_complement()
		else:
			return self.seq

# Définition des fonctions 

#Fonction qui va récupérer les séquences fasta pour chaque exon/intron et les enregistrer dans un dictionnaire
def get_seq_fasta(file_name):
	seq_info = {}
	handle = open(file_name, "rU") #On ouvre le fichier en mode lecture

	for record in SeqIO.parse(handle, "fasta") : #On parse le fichier avec seqIO
		seq_id = record.id #On enregistre le nom de la sequence, ici les coordonnées
		seq = record.seq #On enregistre la sequence fasta dans une variable
		seq= seq.upper() #On met en majuscule au cas ou la séquence ne l'est pas
		seq_info[seq_id]=seq # On enregistre la séquence dans le dictionnaire 
	handle.close()
	return(seq_info)

# Fonction qui va enregistrer chaque exon dans un objet, puis l'ajouter dans un dictionnaire pour chaque transcrit
def get_all_cds_for_transcript(file_name,seq_info):

	CDS_by_transcript = {}
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête

	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") # On split la ligne dans une liste
		# On récupère chaque élément de la liste dans une variable
		exon_id = content[0]
		trans_id = content[1]
		gene_id = content[2]
		coords = content[3].replace('\n', '')
		# On récupère la séquence au format BioSeq via le dictionnaire
		exon_seq = seq_info[exon_id]
		# On enregistre le tout dans l'objet
		exon = ExonInfo(exon_id,trans_id,gene_id,coords,exon_seq)
		# Puis on ajoute l'objet dans un dictionnaire, annoté par transcrit
		if exon.trans_id in CDS_by_transcript:
			CDS_by_transcript[exon.trans_id].append((exon.coords,exon))
		else:
			CDS_by_transcript[exon.trans_id]=[]
			CDS_by_transcript[exon.trans_id].append((exon.coords,exon))
	file_in.close()
	return(CDS_by_transcript)

# Fonction qui va enregistrer chaque intron dans un objet, puis l'ajouter dans un dictionnaire pour chaque transcrit
def get_all_intron_for_transcript(file_name,seq_info):

	intron_by_transcript = {}
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête

	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t")
		intron_id = content[0]
		trans_id = content[1]
		gene_id = content[2]
		coords = content[3].replace('\n', '')
		intron_seq = seq_info[intron_id]
		intron = IntronInfo(intron_id,trans_id,gene_id,coords,intron_seq)
		if intron.trans_id in intron_by_transcript:
			intron_by_transcript[intron.trans_id].append((intron.coords,intron))
		else:
			intron_by_transcript[intron.trans_id]=[]
			intron_by_transcript[intron.trans_id].append((intron.coords,intron))
	file_in.close()
	return(intron_by_transcript)

# Fonction qui va associer pour chaque transcrit tous ses exons et introns, et les trier par coordonnées
def association_between_intron_and_exon(intron_by_transcript,CDS_by_transcript):
	transcript_complete = {}
	for trans_id, list_exon in CDS_by_transcript.items():
		if trans_id in intron_by_transcript:
			introns = intron_by_transcript[trans_id]
			full_list = []
			full_list.extend(list_exon)
			full_list.extend(introns)
			# Fonction de tri, on a des coordonnées celon la conventien : chr[+/-]debut:fin
			# On va donc prendre en compte uniquement la fin, donc on split la chaine pour obetenir fin, qu'on converti en entier, et qu'on tri dans l'ordre
			full_list.sort(key=lambda x: int(x[0].split(':')[1]))
			transcript_complete[trans_id]=full_list
	return(transcript_complete)

##########################################################
################FONCTIONS NMD VISIBILITY##################
##########################################################
# Fonction qui va construire les éléments du CDS à partir du premier exon codant jusqu'au dernier (sans les introns UTR)
def CDS_construction(transcript):
	CDS_elmt = [] # Liste des éléments du CDS -> sans les introns avant et après les premiers et dernier exons
	# Enregistre les coordonnées de chaque bord des exons +1 et -1 pour déterminer les introns présents dans le CDS
	# Globalement, le CDS va de la coordonnées la plus petite à la coordonnées la plus grande dans cette liste
	exons_coords = [int(elmt[1].start)-1 for elmt in transcript if elmt[1].get_type == "exon"]
	exons_coords.extend([int(elmt[1].end)+1 for elmt in transcript if elmt[1].get_type == "exon"])
	# On construit une liste contenant uniquement les introns afin de pouvoir les trier
	introns = [(elmt[1].start,elmt[1].end,elmt[1]) for elmt in transcript if elmt[1].get_type == "intron"]
	# Si un intron possède deux coordonnées présentes dans la liste exons_coords, alors c'est qu'il fait parti du CDS
	CDS_elmt = [intron for intron in introns if (int(intron[0]) in exons_coords and int(intron[1]) in exons_coords)]
	CDS_elmt.extend([(elmt[1].start,elmt[1].end,elmt[1]) for elmt in transcript if elmt[1].get_type == "exon"])
	CDS_elmt.sort(key=lambda x: int(x[1]))
	
	return(CDS_elmt)

# Fonction qui va lister les introns ne faisant pas parti du CDS
def intron_no_CDS(transcript,CDS_elmt):
	# On met dans une liste tous les introns qui font partis du CDS
	introns = [elmt[2].id for elmt in CDS_elmt if elmt[2].get_type == "intron"]
	# Puis on parcours chaque intron contenu dans le transcrit et on renvoit ceux qui ne sont pas dans le CDS
	return [(elmt[1].start,elmt[1].end,elmt[1]) for elmt in transcript if (elmt[1].get_type == "intron" and elmt[1].id not in introns)]

# Fonction qui va construire la séquence pour une liste d'identifiants
def seq_CDS(list_id):
	if list_id[0][2].strand == "-":
		return [str(elmt[2].seq) for elmt in reversed(list_id)]
	else:
		return [str(elmt[2].seq) for elmt in list_id]

# Fonction qui va déterminer combien de codon stop il y'a et la position comme control
def stop_control(CDS,transcript_id):
	complete_CDS = ""
	for elmt in CDS : complete_CDS+=elmt
	position = 0
	stop_control = 0
	while position<len(complete_CDS):
		triplet = complete_CDS[position:position+3]
		if triplet.endswith(("TGA","TAG","TAA")) : stop_control+=1
		position+=3
	if stop_control > 1 : 
		return 1
	else : return 0

# Fonction qui va enlever tous les introns pour construire le CDS
def CDS(CDS_elmt):
	return [(elmt[2].start,elmt[2].end,elmt[2]) for elmt in CDS_elmt if elmt[2].get_type == "exon"]

# Fonction qui va introduire un intron dans le CDS
def insert_intron(CDS_list,intron):
	CDS_intron = CDS_list
	CDS_intron.append(intron)
	CDS_intron.sort(key=lambda x: int(x[1]))
	return(CDS_intron)

# Fonction principale qui va parcourir chaque transcrit et déterminer la NMD visibilitée des introns du transcrit
def main_for_NMD(transcript_complete):
	# Dictionnaire qui va contenir l'information pour chacuns de nos introns si il est ou non NMD visible
	intron_annotation = {}
	comptor = 0
	# On parcours chaque transcrit de notre liste
	for transcript_id, transcript_content in transcript_complete.items():
		# On envoi le contenu du transcrit (intron et exons codants) à notre fonction
		# Va renvoyer la même liste sans les introns non contenus dans le CDS
		CDS_elmt = CDS_construction(transcript_content)
		# On récupère la liste des introns qui ne sont pas dans le CDS
		introns_in_UTR = intron_no_CDS(transcript_content,CDS_elmt)
		# On créer une liste contenant uniquement les exons -> CDS
		CDS_list = CDS(CDS_elmt)
		# On construit la séquence du CDS (sous forme de liste)
		# Afin de vérifier si le CDS possède un codon stop prématuré
		exon_CDS = seq_CDS(CDS_list)
		comptor += stop_control(exon_CDS,transcript_id)
		# Ajout de nos introns du CDS un a un :
		intron_list = [elmt[2].id for elmt in CDS_elmt if elmt[2].get_type == "intron"]
		# Pour chacun de nos introns contenus dans le CDS : 
		for intron in intron_list:
			# On insert l'intron dans le CDS : 
			CDS_intron = insert_intron(CDS_list,intron)
			print(CDS_intron)
			# Puis on fait la séquence 
			seq_CDS_intron = seq_CDS(CDS_intron)

	print(comptor)

# Interface avec l'utilisateur :
opts, args = getopt.getopt(sys.argv[1:],'',['liste_exon=','liste_intron=','fasta=',])
for elmts in opts:
	if elmts[0] == '--liste_exon':
		liste_exon = elmts[1] # Fichier contenant les coordonnées pour chaque exon : liste_exon.tab

	elif elmts[0] == '--liste_intron':
		liste_intron = elmts[1] # Fichier contenant les coordonnées pour chaque intron : list_intron.tab

	elif elmts[0] == '--fasta':
		fasta_file = elmts[1] # Fichier contenant les séquences au format fasta des exons et des introns : seq.fa

# Lancement des fonctions :
# Récupération des séquences fasta
seq_info = get_seq_fasta(fasta_file)
# Récupération des exons 
CDS_by_transcript = get_all_cds_for_transcript(liste_exon,seq_info)
# Récupération des introns
intron_by_transcript = get_all_intron_for_transcript(liste_intron,seq_info)
# Association des introns et des exons dans un transcrit
transcript_complete = association_between_intron_and_exon(intron_by_transcript,CDS_by_transcript)