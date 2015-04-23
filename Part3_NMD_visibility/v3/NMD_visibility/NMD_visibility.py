# Programme qui va récupérer la liste des exons codants et des introns pour chaque transcrit canonique ensembl
# A partir de ces annotations, va organiser chaque transcrit canonique de façon à ce qu'il contienne tous ses exons/introns
# A partir de cette liste, va construire le CDS de chaque transcrit, puis introduire un à un chaque intron 
# Pour chaque intron, on déterminera si celui-ci est NMD visible, et si il décale la phase, on indiquera donc la phase de l'intron


import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq
from math import *
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
		self.exon_pos = pos_in_transcript(self.exons,self.strand)
		# Introns : 
		self.intron_pos = [intron_start(self.exons,one_intron,self.strand) for one_intron in self.introns]
		# PTC :
		self.PTC = self.PTC_pos()
		# Type du transcrit :
		self.type = self.CDS_type()
		# Bornes 5' UTR : (la fin de la borne est défini par le début du codon start)
		self.five_UTR_region = (0,self.cds_start)
		# Bornes 3' UTR :
		self.three_UTR_region = (self.cds_stop,len(self.seq[self.cds_stop:])-1)
		# Fenetres de position pour les régions 5' et 3' UTR
		self.five_UTR_windows = [(c,c-29) for c in cumsum_for_UTR(self.five_UTR_region[1],floor(self.five_UTR_region[1]/30),-30) if c-29 >=0]
		self.three_UTR_windows = [(c,c+29) for c in cumsum_for_UTR(self.three_UTR_region[0],floor(self.three_UTR_region[1]/30),30) if c+29 <= self.three_UTR_region[1] ]
		# Introns dans les régions UTRs :
		self.intron_in_five_UTR = [elmt[0] for elmt in self.intron_pos if elmt[0] <= self.five_UTR_region[1]]
		self.intron_in_three_UTR = [elmt[0] for elmt in self.intron_pos if elmt[0] >= self.three_UTR_region[0]]
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


class PTCInfo(object):
	"Classe qui contiendra les annotations de chaque intron pour la NMD visibilité"
	def __init__(self,intron_id,transcript,gene_id,CDS_status,dist_next_intron,dist_next_stop,dist_CDS_stop,PTC_status,intron_position,intron_phase,intron_start):
		self.id = intron_id
		self.transcript = transcript
		self.gene_id = gene_id
		self.CDS_status = CDS_status
		self.dist_next_intron = dist_next_intron
		self.dist_next_stop = dist_next_stop
		self.dist_CDS_stop = dist_CDS_stop
		self.PTC_status = PTC_status
		self.position = intron_position
		self.phase = intron_phase
		self.start = intron_start
	def format_print(self):
		return self.id+"\t"+self.transcript.id+"\t"+self.gene_id+"\t"+self.CDS_status+"\t"+str(self.dist_next_intron)+"\t"+str(self.dist_next_stop)+"\t"+str(self.dist_CDS_stop)+"\t"+self.PTC_status+"\t"+str(self.position)+"\t"+str(self.phase)+"\n"


class TranscriptAnnotate(object):
	"Classe qui contiendra les annotations du transcrit"
	def __init__(trans_id,gene_id,status,intron_number,exon_number,intron_bf_stop):
		self.id = trans_id
		self.gene_id = gene_id
		self.status = status
		self.intron_number = intron_number
		self.exon_number = exon_number
		self.intron_bf_stop = intron_bf_stop

	def format_print(self):
		return self.id+"\t"+self.gene_id+"\t"+self.status+"\t"+self.intron_number+"\t"+self.exon_number+"\t"+self.intron_bf_stop+"\n"


# Fonctions qui servent à l'objet :

def intron_start(exon_list,one_intron,strand):
		CDS_intron = []
		CDS_intron.extend(exon_list)
		CDS_intron.append(one_intron)
		CDS_intron.sort(key=lambda x: int(x[1]))
		CDS_upon_intron_pos = pos_in_transcript(CDS_intron,strand)
		intron_pos = [elmt for elmt in CDS_upon_intron_pos if elmt[2].get_type == "intron"][0]
		return(intron_pos)

def pos_in_transcript(list_object,strand):
	elmt_list = list_object
	if strand == "-":
		elmt_list.reverse()
	start_list = [c for c in cumsum(elmt_list)]
	end_list = [start_list[i]+len(elmt_list[i][2].seq)-1 for i in range(0,len(elmt_list),1)]
	return [(start_list[i],end_list[i],elmt_list[i][2]) for i in range(0,len(elmt_list),1)]

def cumsum(liste):
	s = 0
	yield s
	for elmt in liste[:-1]:
		s += len(elmt[2].seq)
		yield s

def cumsum_for_UTR(s,max,window):
	yield s
	for i in range(max-1):
		s+=window
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
def get_all_exon_for_transcript(file_name,seq_info):

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
# Fonction qui va annoter chaque transcrit de notre jeu de données
# ANNOTATIONS TRANSCRITS
# -Gene
# -Status de son CDS
# -Nombres d'introns 
# -Nombres d'exons
# -Nombre d'introns avant le stop, suceptible d'avoir un impact sur le CDS
def CDS_annotation(transcript_complete):
	CDS_dic = {}
	for trans_id,trans_object in transcript_complete.items():
		# Identifiant du transcrit :
		transcript_id = trans_id
		# Identifiant du gène du transcrit :
		gene_id = trans_object.gene_id
		# Status du CDS du transcrit :
		CDS_status = trans_object.type
		# Nombre d'intron :
		intron_number = len(trans_object.introns)
		# Nombre d'exons : 
		exon_number = len(trans_object.exons)
		# Nombre d'introns avant le stop :
		intron_before_stop = len([elmt for elmt in trans_object.intron_pos if elmt[0]<= trans_object.cds_stop])
		# On enregistre les annotations dans un objet
		Transcript = TranscriptAnnotate(trans_id=transcript_id,gene_id=gene_id,status=CDS_status,intron_number=intron_number,exon_number=exon_number,intron_bf_stop=intron_before_stop)
		CDS_dic[Transcript.id]=Transcript


# Fonction principale, pour chaque transcrit, va annoter tous ses introns, et donner les informations suivantes :
# ANNOTATIONS INTRONS
# -Transcrit
# -Gene
# -Numéro de l'intron
# -Status du CDS
# -Phase de l'intron
# -Si PTC en cas de rétention
# -Distance en bp entre le stop du CDS et la position de l'intron en cas de rétention
# -Distance en bp du prochain intron en cas de rétention
# -Distance en bp du prochain stop en cas de rétention
def PTC_annotation(transcript_complete):
	PTC_dic = {}
	# On initialise les compteurs à 0
	total_intron = intron_NMD = no_NMD = five = three = 0
	for trans_id,trans_object in transcript_complete.items():
		#### ANNOTATIONS SUR LE TRANSCRIT
		# Identifiant du transcrit :
		transcript_id = trans_id
		# Identifiant du gène du transcrit :
		gene_id = trans_object.gene_id
		# Séquence du transcrit
		seq_for_transcript = trans_object.seq
		# On récupère les coordonnées transcrit du codon start et stop canonique 
		cds_start = trans_object.cds_start
		cds_stop = trans_object.cds_stop
		# On récupère tous les introns du transcrit, ainsi que leur position sur le transcrit
		introns_for_transcript = trans_object.intron_pos
		#### ANNOTATIONS SUR LES INTRONS
		# Pour chaque intron du transcrit :
		for intron in introns_for_transcript:
			# Compteur total :
			total_intron+=1
			# Identifiant de l'intron :
			intron_id = intron[2].id
			# On récupère le rang de l'intron dans le transcrit 
			intron_rank = introns_for_transcript.index(intron)+1
			# On récupère sa position dans le transcrit
			intron_start = intron[0]
			# On récupère sa séquence :
			intron_object = intron[2]
			intron_seq = intron_object.seq


			# Détermination de l'annotation de l'intron : CDS_OK, CDS_PTC (5' ou 3'), CDS_Partial
			# Si l'intron est contenu dans le CDS
			if intron_start > start and intron_start <= stop:
				# Si le transcrit est annté comme OK, c'est que son CDS a réussi les tests de control qualité (voir fonction objet CDS_type)
				if trans_object.type == "OK":
					CDS_status = "CDS_OK"
				# Si le transcrit possède un codon stop prématuré dans la séquence, on annote l'intron comme étant CDS_PTC
				elif trans_object.type == "PTC":
					CDS_status = "CDS_PTC_5'" if intron_start < min(trans_object.PTC) else "CDS_PTC_3'"
				# Si le status du transcrit est partiel, c'est à dire qu'il ne remplit pas une des conditions suivante :
				#	-Pas de codon start
				#	-Pas de codon stop
				#	-Pas multiple de 3
				# Alors on annote l'intron comme étant CDS_Partial
				elif trans_object.type =="Partial":
					CDS_status = "CDS_Partial"
			# Si l'intron est dans le 5' UTR, on annote l'intron faisant parti de cette région et non détéctable par le système NMD
			elif intron_start <= start:
				CDS_status = "5'UTR"
				five+=1
			# Ou dans la région 3' UTR
			elif intron_start > stop:
				CDS_status = "3'UTR"
				three+=1

			# Phase de l'intron :
			intron_phase = get_phase(status=CDS_status,seq=seq_for_transcript,start=cds_start,intron_pos=intron_start)

			# PTC status :
			PTC_status = annotation_PTC_status(status=CDS_status,seq=seq_for_transcript,start=cds_start,stop=cds_stop,intron_pos=intron_start,intron_seq=intron_seq)
			if PTC_status != "NA":
				intron_NMD+=1
			else:
				no_NMD+=1
			# Calcul de la distance entre le stop du cds et le début de l'intron en cas de rétention
			dist_intron_CDS_stop = dist_CDS_stop(status=CDS_status,stop=cds_stop,intron_pos=intron_start,len_intron=len(intron_seq)-1)

			# Calcul de la distance entre l'intron et l'intron suivant en cas de retention d'intron
			dist_intron_next_intron = dist_next_intron(intron,introns_for_transcript)

			# Calcul de la distance du stop le plus proche : 
			dist_next_stop = next_stop(seq=seq_for_transcript,intron_seq=intron_object.seq,intron_start=intron_start)

			# On enregistre toutes nos annotations dans un objet NMD, qu'on enregistre dans un dictionnaire
			PTC_info = NMDInfo(intron_id,transcript_id,gene_id,CDS_status,dist_intron_next_intron,dist_next_stop,dist_intron_CDS_stop,PTC_status,intron_rank,intron_phase,intron_start)
			PTC_dic[PTC_info.id] = PTC_info
	print ("Introns NMD visibles :",intron_NMD,"sur",total_intron)
	print("No NMD :",no_NMD)
	print ("5'UTR :",five)
	print("3'UTR :",three)
	return PTC_dic

# Fonction qui va déterminer la phase de l'intron dans le CDS
def get_phase(status,start,seq,intron_pos):
	if status =="5'UTR" or status == "3'UTR":
		return status
	else:
		return len(seq[start:intron_pos])%3

# Fonction qui va déterminer si l'intron est NMD visible lorsqu'il est retenu
def annotation_PTC_status(status,seq,start,stop,intron_pos,intron_seq):
	if status != "CDS_OK":
		return "NA"
	else:
		seq_with_intron = seq_with_intron_retention(seq,start,stop-3,intron_seq,intron_pos)
		stop_position = stop_position_in_seq(seq_with_intron)
		if stop_position != None and len(stop_position)>=1 :
			return "CDS_PTC_upon_insertion"
		else :
			return "CDS_noPTC_upon_insertion"

# Fonction qui va construire la séquence du CDS avec l'intron retenu
def seq_with_intron_retention(seq,start,stop,intron_seq,intron_pos):
	return seq[start:intron_pos]+intron_seq+seq[intron_pos:stop]

# Fonction qui va retourner la position dans une séquence de tous les codons stop
# Si on part du start, la position du codon stop dans le transcrit sera position_trouvée+position_start_transcrit+3
def stop_position_in_seq(any_seq):
	# On définit un tuple contenant tous les codons stop possible
	stop_codon =("TGA","TAG","TAA")
	# On créer une liste de tous les codons de la séquence
	codons = [any_seq[i:i+3] for i in range(0,len(any_seq),3)]
	index_codons = [ i*3 for i in range(0,len(codons),1) if codons[i] in stop_codon]
	if index_codons != None and len(index_codons) >= 1:
		return (index_codons)
	else:
		return None

# Fonction de calcul de la distance entre le stop canonique et le début de l'intron en cas de retention
def dist_CDS_stop(status,stop,intron_pos,len_intron):
	# Si l'intron est dans la région 3' UTR, la retention de l'intron n'a pas d'impact sur la position du stop
	if status == "3'UTR" :
		return stop-intron_pos
	# Sinon, on doit ajouter la taille de l'intron à la position du stop, car celle-ci sera décalée d'autant de bp que contient la séquence de l'intron
	else:
		return (stop+len_intron)-intron_pos

# Fonction qui va calculer la distance entre l'intron et l'intron suivant en cas de retention
def dist_next_intron(one_intron,intron_list):
	if one_intron == intron_list[-1]:
		return -1
	else:
		next_intron = intron_list[intron_list.index(one_intron)+1]
		return int(next_intron[0])-int(one_intron[0])+len(one_intron[2].seq)-1

# Fonction qui va calculer le stop le plus proche de l'intron, si il n'y en a pas (donc intron dans région 3' UTR), la fonction retourne -1
def next_stop(seq,intron_seq,intron_start):
	seq_intron_and_next_transcript = intron_seq+seq[intron_start:]
	stop_position_list = stop_position_in_seq(seq_intron_and_next_transcript)
	if stop_position_list != None:
		dist_between_stop_and_intron_start = min([PTC_pos for PTC_pos in stop_position_list if PTC_pos > intron_start]) if max(stop_position_list) > intron_start else -1
		return dist_between_stop_and_intron_start
	else:
		return -1

##########################################################
##################FONCTIONS BORNES UTR####################
##########################################################
def intron_density_in_UTR(transcript_complete):
	# On définit deux dictionnaires qui contiendrons les fenêtres des régions 5' et 3'
	five_UTR_density = {}
	three_UTR_density = {}
	# Pour chaque transcrit de notre jeu de données
	for trans_id,trans_object in transcript_complete.items():
		# UTR window :
		five_windows = trans_object.five_UTR_windows
		three_windows = trans_object.three_UTR_windows
		# Introns in UTR :
		five_introns = trans_object.intron_in_five_UTR
		five_introns.sort(reverse = True)
		three_introns = trans_object.intron_in_three_UTR
		three_introns.sort()

		# Parcours des fenetres et assignation des introns
		# Pour le 5' UTR
		for window in five_windows:
			# Si aucun intron dans le 5', on associe la valeur 0 à notre fenêtre
			if len(five_introns) == 0:
				five_UTR_density[window]=0
			# Sinon, on parcours les introns du 5' UTR
			else:
				for intron in five_introns:
					# Si l'intron est contenu dans la fenêtre
					if intron>=window[1] and intron<=window[0]:
						# Alors on l'ajoute +1 à la fenêtre
						if window in five_UTR_density:
							five_UTR_density[window]+=1
						else:
							five_UTR_density[window]=1
					# Sinon
					else :
						# On dit que la fenêtre ne contient aucun intron si celle si n'existe pas dans le dictionnaire
						if not window in five_UTR_density:
							five_UTR_density[window]=0
		# Pour le 3' UTR
		for window in three_windows:
			if len(three_introns) == 0:
				three_UTR_density[window]=0
			else:
				for intron in three_introns:
					if intron>=window[0] and intron<=window[1]:
						if window in three_UTR_density:
							three_UTR_density[window]+=1
						else:
							three_UTR_density[window]=1
					else :
						if not window in three_UTR_density:
							three_UTR_density[window]=0
	return(five_UTR_density,three_UTR_density)

def intron_density_in_UTR_test(transcript_complete):
	# On définit deux dictionnaires qui contiendrons les fenêtres des régions 5' et 3'
	five_UTR_density = {}
	three_UTR_density = {}
	# Pour chaque transcrit de notre jeu de données
	for trans_id,trans_object in transcript_complete.items():
		# UTR window :
		five_windows = trans_object.five_UTR_windows
		three_windows = trans_object.three_UTR_windows
		# On construit le dictionnaire en initialisant les valeurs à 0 :
		for elmt in five_windows:
			if not elmt in five_UTR_density:
				five_UTR_density[elmt]=0
		for elmt in three_windows:
			if not elmt in three_UTR_density:
				three_UTR_density[elmt]=0
		# Introns in UTR :
		five_introns = trans_object.intron_in_five_UTR
		three_introns = trans_object.intron_in_three_UTR
		# Assignation du nombre d'introns par fenêtres :
		# Pour le 5' UTR
		five_UTR_density = calcul_density(UTR_density=five_UTR_density,introns=five_introns,windows=five_windows,reverse=True)
		# Pour le 3' UTR
		three_UTR_density = calcul_density(UTR_density=three_UTR_density,introns=three_introns,windows=three_windows,reverse=False)
	return(five_UTR_density,three_UTR_density)

# Fonction qui calcule, pour deux listes données, le nombre d'éléments dans la deuxième liste présents entre les deux bornes de chaque élément de la première
def calcul_density(UTR_density,introns,windows,reverse):
	# En fonction des régions UTR, l'élement le plus grand des limites de chaque fenêtre peut être le premier ou le second
	# Dans le cas du 5' UTR, on part du start vers le début du transcrit, on a donc l'élément le plus grand à gauche et le plus petit à droite
	# Dans le cas du 3' UTR, on part du stop vers la fin, on a donc l'élement le plus petit à gauche et le plus grandà droite
	min_window_limit = 1 if reverse == True else 0
	max_window_limit = 0 if reverse == True else 1
	# On effectue un tri des deux listes dans l'ordre croissant
	windows.sort()
	introns.sort()
	# On définit nos compteurs à 0
	i = 0
	j = 0
	# Tant que le compteur i est inférieur à la taille de la liste des fenêtres, on continue
	while i<len(windows):
		# Si le compteur j est inférieur à la taille de la liste des introns, on peut vérifier l'élement de la liste des introns à la position j
		if j<len(introns):
			# Si la position de l'intron dans le transcrit est situé entre les limites de la fenêtres
			if introns[j] <= windows[i][max_window_limit] and introns[j] >= windows[i][min_window_limit]:
				# On ajoute +1 au compteur de la fenêtre pour tous nos transcrits
				UTR_density[windows[i]]+=1
				# Et +1 au compteur de nos introns, celui-ci étant dans cette fenêtre, il ne peut pas être dans les autres, plus besoin de la tester
				j+=1
			# Si l'intron n'est pas situé entres les limites de la fenêtres
			else:
				# On met +1 au compteur des fenêtres, si le premier intron ne match pas, aucun autre de la liste ne correspondra
				i+=1
		# Sinon, c'est que le compteur est allé au dela du nombre d'éléments de la liste des introns, ce qui veut dire qu'il n'y a plus d'introns et donc que ça ne sert plus à rien de parcourir les fenêtres
		else:
			break
	return(UTR_density)

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
CDS_by_transcript = get_all_exon_for_transcript(liste_exon,seq_info)
# Récupération des introns
intron_by_transcript = get_all_intron_for_transcript(liste_intron,seq_info)
# Association des introns et des exons dans un transcrit
transcript_complete = association_between_intron_and_exon(intron_by_transcript,CDS_by_transcript)
# Pour chaque transcrit, on annote chaque intron du transcrit pour vérifier si il est NMD visible
NMD_dic = PTC_annotation(transcript_complete)

