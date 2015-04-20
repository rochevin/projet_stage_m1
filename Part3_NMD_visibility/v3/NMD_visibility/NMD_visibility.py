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

class TranscriptInfo(object):
	"Classe qui contiendra les annotations de chaque transcrits"
	def __init__(self,transcript_id,list_object):
		# Expression régulière qui va récupérer l'id, le debut et la fin du CDS dans la séquence du transcrit
		regex = re.compile('([\w]+)\(([0-9]+):([0-9]+)\)')
		result = regex.findall(transcript_id)
		self.id = result[0][0] # ID du transcrit
		self.transcript_content = list_object # Liste d'objets contenus dans le transcrits (introns + exons)
		self.exons = [(exon[1].start,exon[1].end,exon[1]) for exon in list_object if exon[1].get_type == "exon"] # Liste des exons sous forme d'objet
		self.introns = [(intron[1].start,intron[1].end,intron[1]) for intron in list_object if intron[1].get_type == "intron"] # Liste des introns sous forme d'objet
		self.gene_id = self.transcript_content[0][1].gene_id # On récupère l'ID du gene
		self.strand = self.strand() # On appelle la fonction qui va récupérer le brin
		# Codon start :
		self.cds_start = int(result[0][1])-1
		# Codon stop :
		self.cds_stop = int(result[0][2])
		# Séquences du transcrit :
		# Séquence complète :
		self.seq = self.transcript_seq()
		# Séquence du CDS :
		self.CDS_seq = self.CDS_seq()
		# Calcul de position pour tous les éléments du transcrit :
		# Exons :
		self.exon_pos = pos_in_transcript(self.exons)
		# Introns : 
		self.intron_pos = [intron_start(self.exons,one_intron) for one_intron in self.introns]
		# PTC :
		self.PTC = self.PTC_pos()
		# Type du transcrit :
		self.type = self.CDS_type()
	# On récupère le brin du transcrit à partir de l'annotation des éléments du transcrit
	def strand(self):
		return self.transcript_content[0][1].strand

	# On récupère la séquence du transcrit, c'est à dire la séquence de tous les exons
	def transcript_seq(self):
		self.exons.sort(key=lambda x: int(x[1]))
		return "".join([str(elmt[2].seq) for elmt in reversed(self.exons)]) if self.exons[0][2].strand == "-" else "".join([str(elmt[2].seq) for elmt in self.exons])

	# On récupère le CDS du transcrit à partir de la séquence entière, en bornant avec le codon initiateur et terminateur
	def CDS_seq(self):
		return self.transcript_seq()[self.cds_start:self.cds_stop]

	# On récupère sous forme de liste la position de tous les PTC du CDS, renvoit none si aucun PTC
	def PTC_pos(self):
		# On définit un tuple contenant tous les codons stop possible
		stop_codon = ("TGA","TAG","TAA")
		# On vérifie qu'il ne contient pas un codon stop prématuré
		codons = [self.CDS_seq[i:i+3] for i in range(0,len(self.CDS_seq),3)]
		index_codons = [ i*3+self.cds_start+3 for i in range(0,len(codons),1) if codons[i] in stop_codon][:-1]
		return index_codons if len(index_codons) > 0 else None
		
	# On récupère le type du transcrit, celui-ci peut être OK, partiel (sans start/stop ni multiple de 3), ou PTC (CaD avec un codon stop prematuré)
	def CDS_type(self):
		# On définit un tuple contenant tous les codons stop possible
		stop_codon = ("TGA","TAG","TAA")
		# On vérifie qu'il commence par un start et fini par un stop :
		if not (self.CDS_seq.endswith(stop_codon) or self.CDS_seq.startswith("ATG")) : 
			return "Partial"
		# On vérifie qu'il est bien multiple de 3 :
		if not len(self.CDS_seq)%3 == 0 : 
			return "Partial"
		# On vérifie si il contient un PTC
		if self.PTC != None:
			return "PTC"
		# Sinon on envoi l'annotation OK
		else:
			return "OK"
	# On récupère la position d'un exon dans la liste des exons
	def exon_position(self,one_exon):
		exon_in_list = self.exons.index(one_exon)
		return self.exon_pos[exon_in_list]


class ExonInfo(object):
	"Classe qui contiendra les annotations de chaque exon"
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
	"Classe qui contiendra les annotations de chaque intron"
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


class NMDInfo(object):
	"Classe qui contiendra les annotations de chaque intron pour la NMD visibilité"
	def __init__(self,intron_id,trans_id,gene_id,CDS_status,dist_next_intron,dist_next_stop,dist_CDS_stop,NMD_status,intron_position,intron_phase,intron_start):
		self.id = intron_id
		self.trans_id = trans_id
		self.gene_id = gene_id
		self.CDS_status = CDS_status
		self.dist_next_intron = dist_next_intron
		self.dist_next_stop = dist_next_stop
		self.dist_CDS_stop = dist_CDS_stop
		self.NMD_status = NMD_status
		self.position = intron_position
		self.phase = intron_phase
		self.start = intron_start
	def format_print(self):
		return self.id+"\t"+self.trans_id+"\t"+self.gene_id+"\t"+self.CDS_status+"\t"+str(self.dist_next_intron)+"\t"+str(self.dist_next_stop)+"\t"+str(self.dist_CDS_stop)+"\t"+self.NMD_status+"\t"+str(self.position)+"\t"+str(self.phase)+"\n"
# Fonctions qui servent à l'objet :

def intron_start(exon_list,one_intron):
		CDS_intron = []
		CDS_intron.extend(exon_list)
		CDS_intron.append(one_intron)
		CDS_intron.sort(key=lambda x: int(x[1]))
		index_intron = CDS_intron.index(one_intron)
		intron_pos = pos_in_transcript(CDS_intron)[index_intron]
		return(intron_pos)

def pos_in_transcript(elmt_list):
	start_list = [c for c in cumsum(elmt_list)]
	end_list = [start_list[i]+len(elmt_list[i][2].seq)-1 for i in range(0,len(elmt_list),1)]
	return [(start_list[i],end_list[i],elmt_list[i][2]) for i in range(0,len(elmt_list),1)]

def cumsum(liste):
	s = 0
	yield s
	for elmt in liste[:-1]:
		s += len(elmt[2].seq)
		yield s

# Fonctions pour déroulement du programme
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
# Puis va construire un objet transcrit qui contientra toutes nos données
def association_between_intron_and_exon(intron_by_transcript,CDS_by_transcript):
	transcript_complete = {}
	for trans_id, list_exon in CDS_by_transcript.items():
		real_trans_id = trans_id.split("(")[0]
		if trans_id in intron_by_transcript:
			introns = intron_by_transcript[trans_id]
			full_list = []
			full_list.extend(list_exon)
			full_list.extend(introns)
			# Fonction de tri, on a des coordonnées celon la conventien : chr[+/-]debut:fin
			# On va donc prendre en compte uniquement la fin, donc on split la chaine pour obetenir fin, qu'on converti en entier, et qu'on tri dans l'ordre
			full_list.sort(key=lambda x: int(x[0].split(':')[1]))
			# On créer l'objet transcrit, qui contiendra les introns et les exons du transcrit
			transcript_object = TranscriptInfo(trans_id,full_list)
			transcript_complete[real_trans_id]=transcript_object
	return(transcript_complete)

##########################################################
################FONCTIONS NMD VISIBILITY##################
##########################################################
# Fonction qui va construire la séquence du CDS avec l'intron retenu
def seq_with_intron_retention(seq,start,intron_seq,intron_pos):
	return seq[start:intron_pos]+intron_seq+seq[intron_pos:]

def seq_with_intron_retention_in_5UTR(seq,intron_seq,intron_pos):
	seq_with_intron = seq[:intron_pos]+intron_seq+seq[intron_pos:]
	new_start_pos = seq_with_intron.find('ATG')
	return seq_with_intron[new_start_pos:intron_pos]+intron_seq+seq_with_intron[intron_pos:]
		 
# Fonction qui va déterminer la phase de l'intron dans le CDS
def get_phase(seq,start,intron_pos):
	return len(seq[start:intron_pos])%3

# Fonction qui va déterminer si l'intron est NMD visible lorsqu'il est retenu
def NMD_visibility(seq_with_intron,cds_stop_with_intron):
	# On définit un tuple contenant tous les codons stop possible
	stop_codon =("TGA","TAG","TAA")
	# On créer une liste de tous les codons de la séquence
	codons = [seq_with_intron[i:i+3] for i in range(0,len(seq_with_intron),3)]
	index_codons = [ i*3 for i in range(0,len(codons),1) if codons[i] in stop_codon]
	stop_number = len(index_codons)
	# if not cds_stop_with_intron in index_codons :
	# 	return "NMD"
	if stop_number>1 :
		return "NMD"
	else :
		return "NA"



# Fonction qui va calculer le stop le plus proche de l'intron, si il n'y en a pas (donc intron dans région 3' UTR), la fonction retourne -1
def next_stop(PTC_list,intron_start,CDS_stop):
	list_for_fonction = []
	list_for_fonction.extend(PTC_list)
	list_for_fonction.append(CDS_stop)
	return min([PTC_pos for PTC_pos in list_for_fonction if PTC_pos > intron_start]) if max(list_for_fonction) > intron_start else -1

# Fonction qui va calculer la distance entre l'intron et le prochain, après retention de sa séquence
def dist_next_intron(one_intron,intron_list):
	if one_intron == intron_list[-1]:
		return -1
	else:
		next_intron = intron_list[intron_list.index(one_intron)+1]
		return int(next_intron[0])+len(one_intron[2].seq)-1

# Fonction qui va permettre de récupérer les annotations pour chaque introns du transcrit
def one_intron_retention(transcript_complete):
	NMD_dic = {}
	intron_NMD = 0
	no_NMD = 0
	five = 0
	three = 0
	total_intron = 0
	for trans_id,trans_object in transcript_complete.items():
		# On récupère la séquence du transcrit
		seq_for_transcript = trans_object.seq
		# On récupère les coordonnées transcrit du codon start et stop canonique 
		start = trans_object.cds_start
		stop = trans_object.cds_stop
		# Si on veut la position du codon stop dans le CDS, on doit soustraire la position du start du CDS, et enlever 3 pour obtenir le début du codon
		stop_position_in_CDS = stop-start-3
		# On récupère tous les introns du transcrit, ainsi que leur position sur le transcrit
		introns_for_transcript = trans_object.intron_pos
		# Pour chaque intron du transcrit :
		for intron in introns_for_transcript:
			total_intron +=1
			# On récupère le rang de l'intron dans le transcrit 
			position = introns_for_transcript.index(intron)+1
			# On récupère sa position dans le transcrit
			intron_start = intron[0]
			intron_object = intron[2]
			# On calcule la phase de l'intron
			phase = get_phase(seq_for_transcript,start,intron_start)
			# On calcule la distance en BP du stop canonique dans le CDS, après retention de l'intron dans le CDS
			dist_CDS_stop = stop_position_in_CDS+len(intron_object.seq)-1

			# Détermination de l'annotation de l'intron : CDS_OK, CDS_PTC (5' ou 3'), CDS_Partial

			# Si l'intron est contenu dans le CDS
			if intron_start > start and intron_start <= stop:
				# On simule la retention de l'intron dans le CDS
				seq_with_intron = seq_with_intron_retention(seq_for_transcript,start,str(intron_object.seq),intron_start)
				# On détermine si l'intron est NMD visible, c'est à dire si le transcrit va être détecter par le système NMD et dégradé lors de la retention de l'intron
				NMD_status = NMD_visibility(seq_with_intron,dist_CDS_stop)
				# Si il est NMD, on ajoute +1 au compteur
				if NMD_status == "NMD":
					intron_NMD+=1 
				else:
					print(intron_object.id)
					no_NMD +=1
				# Si le transcrit possède un codon stop prématuré dans la séquence, on annote l'intron comme étant CDS_PTC
				if trans_object.type == "PTC":
					CDS_status = "CDS_PTC_5'" if intron_start < min(trans_object.PTC) else "CDS_PTC_3'"
				# Si le status du transcrit est partiel, c'est à dire qu'il ne remplit pas une des conditions suivante :
				#	-Pas de codon start
				#	-Pas de codon stop
				#	-Pas multiple de 3
				# Alors on annote l'intron comme étant CDS_Partial
				elif trans_object.type =="Partial":
					CDS_status = "CDS_Partial"
				# Si le transcrit est OK, c'est à dire qu'il produit une protéine viable, on annote l'intron CDS_OK
				elif trans_object.type == "OK":
					CDS_status = "CDS_OK"
			# Si l'intron est dans le 5' UTR, on annote l'intron faisant parti de cette région et non détéctable par le système NMD
			elif intron_start <= start:
				# if len(intron_object.seq)%3 != 0:
				# 	seq_with_intron_UTR = seq_with_intron_retention_in_5UTR(seq_for_transcript,str(intron_object.seq),intron_start)
				# 	NMD_status = NMD_visibility(seq_with_intron_UTR,dist_CDS_stop)
				# 	print(NMD_status)
				# else:
				NMD_status = "NA" 
				CDS_status = "5'UTR"
				five +=1
			# Ou dans la région 3' UTR
			elif intron_start > stop:
				NMD_status = "NA"
				CDS_status = "3'UTR"
				three+=1
			# Dans le cas ou aucune condition n'est réunie :
			else:
				print(trans_id)
				print("start intron =",intron_start)
				print("start CDS=",start)
				print("stop CDS=",stop)
			# Si le transcrit possède un ou plusieurs codon stop prématuré
			if trans_object.PTC != None :
				# On détermine le codon stop le plus proche du transcrit
				intron_next_stop = next_stop(trans_object.PTC,intron_start,stop)
				# Puis on calcule sa distance en BP dans le CDS avec retention de l'intron
				dist_next_stop = intron_next_stop+len(intron_object.seq)-1
			# Sinon le stop le plus proche est le canonique, alors la distance du codon stop le plus proche sera égale à celle du canonique
			else :
				dist_next_stop = dist_CDS_stop
			
			# On calcule la distance entre notre intron et le prochain en cas de retention de notre intron
			dist_next = dist_next_intron(intron,introns_for_transcript)

			# On enregistre toutes nos annotations dans un objet NMD, qu'on enregistre dans un dictionnaire
			NMD_info = NMDInfo(intron[2].id,trans_object.id,trans_object.gene_id,CDS_status,dist_next,dist_next_stop,dist_CDS_stop,NMD_status,position,phase,intron_start)
			NMD_dic[NMD_info.id] = NMD_info

	print ("Introns NMD visibles :",intron_NMD,"sur",total_intron)
	print("No NMD :",no_NMD)
	print ("5'UTR :",five)
	print("3'UTR :",three)
	return NMD_dic

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
# Pour chaque transcrit, on annote chaque intron du transcrit pour vérifier si il est NMD visible
NMD_dic = one_intron_retention(transcript_complete)

count = 0
for key,value in CDS_PTC_5.items():
	if value.NMD_status != "NMD":
		count+=1
print(count)