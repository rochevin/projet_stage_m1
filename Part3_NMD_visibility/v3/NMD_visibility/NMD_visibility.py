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
		# Position génomique du transcrit :
		# Début = position du premier élément du transcrit
		self.start = self.exons[0][0]
		# Fin = position de fin du dernier élément du transcrit
		self.end = self.exons[-1][1]
		# Le chromosome du transcrit : chromosome de n'importe lequel de ses éléments :
		self.chr = self.exons[0][2].chr
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
		self.three_UTR_region = (self.cds_stop+1,len(self.seq)-1)
		# Fenetres de position pour les régions 5' et 3' UTR
		self.five_UTR_windows = [(c,c-29) for c in cumsum_for_5UTR(self.five_UTR_region[1],self.five_UTR_region[0],-30) if c-29 >= self.five_UTR_region[0]]
		self.three_UTR_windows = [(c,c+29) for c in cumsum_for_3UTR(self.three_UTR_region[0],self.three_UTR_region[1],30) if c+29 <= self.three_UTR_region[1]]
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
		if (self.CDS_seq.startswith("ATG")==False or self.CDS_seq.endswith(stop_codon)==False):
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

# Classes qui vont contenir nos annotations, afin qu'elles soient écrites dans des fichiers
class PTCInfo(IntronInfo):
	"Classe qui contiendra les annotations de chaque intron pour la NMD visibilité"
	def __init__(self,intron_id,trans_id,gene_id,coords,seq,CDS_status,dist_last_intron,dist_next_stop,dist_CDS_stop,PTC_status,intron_position,intron_phase,intron_start,ptc_in_intron):
		IntronInfo.__init__(self,intron_id,trans_id,gene_id,coords,seq)
		self.CDS_status = CDS_status
		self.dist_last_intron = dist_last_intron
		self.dist_next_stop = dist_next_stop
		self.dist_CDS_stop = dist_CDS_stop
		self.PTC_status = PTC_status
		self.position = intron_position
		self.phase = intron_phase
		self.trans_start = intron_start
		self.ptc_in_intron = ptc_in_intron

	def format_print(self):
		return self.id+"\t"+self.trans_id.split('(')[0]+"\t"+self.gene_id+"\t"+self.CDS_status+"\t"+str(self.dist_last_intron)+"\t"+str(self.dist_next_stop)+"\t"+str(self.ptc_in_intron)+"\t"+str(self.dist_CDS_stop)+"\t"+self.PTC_status+"\t"+str(self.position)+"\t"+str(self.phase)+"\n"


class TranscriptAnnotate(object):
	"Classe qui contiendra les annotations du transcrit"
	def __init__(self,trans_id,gene_id,status,trans_chr,trans_start,trans_end,intron_number,exon_number,intron_bf_stop):
		self.id = trans_id
		self.gene_id = gene_id
		self.chr = trans_chr
		self.start = trans_start
		self.end = trans_end
		self.status = status
		self.intron_number = intron_number
		self.exon_number = exon_number
		self.intron_bf_stop = intron_bf_stop

	def format_print(self):
		return self.id+"\t"+self.gene_id+"\t"+self.status+"\t"+str(self.intron_number)+"\t"+str(self.exon_number)+"\t"+str(self.intron_bf_stop)+"\n"


# Fonctions qui servent à l'objet :

# Fonction qui va calculer la position de l'intron dans le transcrit
def intron_start(exon_list,one_intron,strand):
		CDS_intron = []
		CDS_intron.extend(exon_list)
		CDS_intron.append(one_intron)
		CDS_intron.sort(key=lambda x: int(x[1]))
		CDS_upon_intron_pos = pos_in_transcript(CDS_intron,strand)
		intron_pos = [elmt for elmt in CDS_upon_intron_pos if elmt[2].get_type == "intron"][0]
		return(intron_pos)
# Fonction de calcul de position, utilise la fonction de cumule des valeurs
def pos_in_transcript(list_object,strand):
	elmt_list = list_object
	if strand == "-":
		elmt_list.reverse()
	# Va sortir une liste de positions d'élément dans un transcrit, on ajoute la taille du transcrit précédent en commençant par 0
	start_list = [c for c in cumsum(elmt_list)]
	return [(start_list[i],start_list[i]+len(elmt_list[i][2].seq)-1,elmt_list[i][2]) for i in range(0,len(elmt_list),1)]

# Fonction qui va servir à sortir les positions de chaque élément d'un transcrit en fonction du précédent, en commençant par 0
def cumsum(liste):
	s = 0
	yield s
	for elmt in liste[:-1]:
		s += len(elmt[2].seq)
		yield s

# Fonction basée sur le même principe que cumsum, sauf qu'il faut déterminer le nombres de fenêtres qu'il peut y avoir dans le transcrit
def cumsum_for_3UTR(s,limit,window):
	yield s
	while s <= limit-60:
		s+=window
		yield s
def cumsum_for_5UTR(s,limit,window):
	yield s
	while s >= limit+60:
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
		seq = seq.upper() #On met en majuscule au cas ou la séquence ne l'est pas
		seq_info[seq_id] = seq # On enregistre la séquence dans le dictionnaire
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
		# Chromosome
		chromosome = trans_object.chr
		# Début/fin du transcrit :
		transcript_start = trans_object.start
		transcript_end = trans_object.end
		# Status du CDS du transcrit :
		CDS_status = trans_object.type
		# Nombre d'intron :
		intron_number = len(trans_object.introns)
		# Nombre d'exons :
		exon_number = len(trans_object.exons)
		# Nombre d'introns avant le stop :
		intron_before_stop = len([elmt for elmt in trans_object.intron_pos if elmt[0]<= trans_object.cds_stop])
		# On enregistre les annotations dans un objet
		Transcript = TranscriptAnnotate(trans_id=transcript_id,gene_id=gene_id,trans_chr=chromosome,trans_start=transcript_start,trans_end=transcript_end,status=CDS_status,intron_number=intron_number,exon_number=exon_number,intron_bf_stop=intron_before_stop)
		CDS_dic[Transcript.id]=Transcript
	return(CDS_dic)


# Fonction principale, pour chaque transcrit, va annoter tous ses introns, et donner les informations suivantes :
# ANNOTATIONS INTRONS
# -Transcrit
# -Gene
# -Numéro de l'intron
# -Status du CDS
# -Phase de l'intron
# -Si PTC en cas de rétention
# -Distance en bp entre le stop du CDS et la position de l'intron en cas de rétention
# -Distance en bp du dernier intron en cas de rétention
# -Distance en bp du prochain stop en cas de rétention

def PTC_annotation(transcript_complete):
	PTC_dic = {}
	# On initialise les compteurs à 0
	total_intron = intron_NMD = no_NMD = five = three = 0
	for trans_id,trans_object in transcript_complete.items():
		#### ANNOTATIONS SUR LE TRANSCRIT
		# Identifiant du transcrit :
		transcript_id = trans_id
		print("Écriture des annotations de",transcript_id,"...")
		print("Introns traités :",total_intron)
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
			CDS_status = determination_CDS_status(intron_start=intron_start,start=cds_start,stop=cds_stop,CDS_type=trans_object.type,PTC_list=trans_object.PTC)
			if CDS_status =="5'UTR":
				five+=1
			elif CDS_status == "3'UTR":
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
			dist_intron_last_intron = dist_last_intron(intron,introns_for_transcript)

			# Calcul de la distance du stop le plus proche :
			if CDS_status != "5'UTR" and CDS_status != "3'UTR":
				dist_next_stop,ptc_in_intron = next_stop(seq=seq_for_transcript,intron_seq=intron_object.seq,intron_start=intron_start,start=cds_start)
			else:
				dist_next_stop = "NA"
				ptc_in_intron = "NA"

			# On enregistre toutes nos annotations dans un objet NMD, qu'on enregistre dans un dictionnaire
			PTC_info = PTCInfo(intron_object.id,intron_object.trans_id,intron_object.gene_id,intron_object.coords,intron_object.seq,CDS_status,dist_intron_last_intron,dist_next_stop,dist_intron_CDS_stop,PTC_status,intron_rank,intron_phase,intron_start,ptc_in_intron)
			PTC_dic[PTC_info.id] = PTC_info
	print ("Introns NMD visibles :",intron_NMD,"sur",total_intron)
	print("No NMD :",no_NMD)
	print ("5'UTR :",five)
	print("3'UTR :",three)
	return PTC_dic

# Si l'intron est contenu dans le CDS
def determination_CDS_status(intron_start,start,stop,CDS_type,PTC_list):
	if intron_start > start and intron_start <= stop:
		# Si le transcrit est annté comme OK, c'est que son CDS a réussi les tests de control qualité (voir fonction objet CDS_type)
		if CDS_type == "OK":
			CDS_status = "CDS_OK"
		# Si le transcrit possède un codon stop prématuré dans la séquence, on annote l'intron comme étant CDS_PTC
		elif CDS_type == "PTC":
			CDS_status = "CDS_PTC_5'" if intron_start < min(PTC_list) else "CDS_PTC_3'"
		# Si le status du transcrit est partiel, c'est à dire qu'il ne remplit pas une des conditions suivante :
		#	-Pas de codon start
		#	-Pas de codon stop
		#	-Pas multiple de 3
		# Alors on annote l'intron comme étant CDS_Partial
		elif CDS_type =="Partial":
			CDS_status = "CDS_Partial"
	# Si l'intron est dans le 5' UTR, on annote l'intron faisant parti de cette région et non détéctable par le système NMD
	elif intron_start <= start:
		CDS_status = "5'UTR"
	# Ou dans la région 3' UTR
	elif intron_start > stop:
		CDS_status = "3'UTR"
	return(CDS_status)

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
def dist_last_intron(one_intron,intron_list):
	if one_intron[2].strand == "-":
		if one_intron == intron_list[0]:
			return 0
		last_intron = intron_list[0]
		return int(last_intron[0])-int(one_intron[0])+len(one_intron[2].seq)-1

	else:
		if one_intron == intron_list[-1]:
			return 0
		last_intron = intron_list[-1]
		return int(last_intron[0])-int(one_intron[0])+len(one_intron[2].seq)-1


# Fonction qui va calculer le stop le plus proche de l'intron, si il n'y en a pas la fonction retourne -1
def next_stop(seq,intron_seq,intron_start,start):
	seq_intron_and_next_transcript = seq[start:intron_start]+intron_seq+seq[intron_start:]
	stop_position_list = stop_position_in_seq(seq_intron_and_next_transcript)
	intron_len = len(intron_seq)
	intron_position_cds = intron_start-start-2 # -2 car on veut aussi 2 bp en amont
	if stop_position_list != None:
		correct_stop_position_list = [position for position in stop_position_list if position>=intron_position_cds]
		dist_between_stop_and_intron_start = min(correct_stop_position_list)-intron_position_cds if len(correct_stop_position_list)>0 else "NA"
		if dist_between_stop_and_intron_start != "NA":
			ptc_in_intron = "Yes" if dist_between_stop_and_intron_start-intron_len <0 else "No"
		else :
			ptc_in_intron = "No"
		return dist_between_stop_and_intron_start,ptc_in_intron
	else:
		return "NA","No"


##########################################################
##################FONCTIONS BORNES UTR####################
##########################################################
def intron_density_in_UTR(transcript_complete):
	# On définit deux dictionnaires qui contiendrons les fenêtres des régions 5' et 3'
	five_UTR_density = {}
	transcript_for_windows_in_five_UTR = {}
	three_UTR_density = {}
	transcript_for_windows_in_three_UTR = {}
	# Pour chaque transcrit de notre jeu de données
	for trans_id,trans_object in transcript_complete.items():
		# UTR window :
		five_windows = trans_object.five_UTR_windows
		three_windows = trans_object.three_UTR_windows
		five_UTR_region = trans_object.five_UTR_region
		three_UTR_region = trans_object.three_UTR_region
		five_windows = [(elmt[0]-five_UTR_region[1],elmt[1]-five_UTR_region[1]) for elmt in five_windows]
		three_windows = [(elmt[0]-three_UTR_region[0],elmt[1]-three_UTR_region[0]) for elmt in three_windows]
		# On construit le dictionnaire en initialisant les valeurs à 0 :
		# On construit également le dictionnaire qui contiendra le nombre de transcrit possédant chaque fenêtres
		five_UTR_density,transcript_for_windows_in_five_UTR = dictionnary_construction(UTR_windows=five_windows,UTR_dictionnary=five_UTR_density,transcript_dictionnary=transcript_for_windows_in_five_UTR)
		three_UTR_density,transcript_for_windows_in_three_UTR = dictionnary_construction(UTR_windows=three_windows,UTR_dictionnary=three_UTR_density,transcript_dictionnary=transcript_for_windows_in_three_UTR)
		# Introns in UTR :
		five_introns = [elmt-five_UTR_region[1] for elmt in trans_object.intron_in_five_UTR]
		three_introns = [elmt-three_UTR_region[0] for elmt in trans_object.intron_in_three_UTR]
		# Assignation du nombre d'introns par fenêtres :
		# Pour le 5' UTR
		five_UTR_density = calcul_density(UTR_density=five_UTR_density,introns=five_introns,windows=five_windows,reverse=True)
		# Pour le 3' UTR
		three_UTR_density = calcul_density(UTR_density=three_UTR_density,introns=three_introns,windows=three_windows,reverse=False)
	return(five_UTR_density,transcript_for_windows_in_five_UTR,three_UTR_density,transcript_for_windows_in_three_UTR)

# Fonction qui étant donné deux dictionnaires pour le nombre de transcrits possédant la fenêtre et pour le nombres d'introns par fenêtres
# Va initialiser le dictionnaire pour les introns
# Et ajouter le nombre de transcrits pour chaque fenêtres
def dictionnary_construction(UTR_windows,UTR_dictionnary,transcript_dictionnary):
	# Pour chaque fenêtres du transcrit
	for elmt in UTR_windows:
		# Si cette fenêtre est deja annotée, on ajoute +1 (+1 transcrit possédant cette fenêtre)
		if elmt in transcript_dictionnary:
			transcript_dictionnary[elmt]+=1
		# Sinon on initialise le dictionnaire à 1
		else:
			transcript_dictionnary[elmt]=1
		# Si la fenêtre ne possède aucun intron d'annoté, on l'initialise à 0
		if not elmt in UTR_dictionnary:
			UTR_dictionnary[elmt]=0
	return(UTR_dictionnary,transcript_dictionnary)

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
	i = j = 0
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

# Fonction d'écriture du fichier de sortie pour les introns
def write_file_for_intron(PTC_dic,file_name_for_intron):

	file_out_intron=open(file_name_for_intron,"w")
	for key,value in PTC_dic.items():
		# Préparation du format BED
		chromosome = value.chr
		intron_start = value.start
		intron_end = value.end
		intron_strand = value.strand
		# Écriture des lignes BED obligatoires :
		Bed_format = chromosome+"\t"+intron_start+"\t"+intron_end+"\t"+intron_strand
		# On enregistre les annotations pour les transcrits et les introns
		annotations_intron=value.format_print()
		line_out_intron = Bed_format+"\t"+annotations_intron
		file_out_intron.write(line_out_intron)

	file_out_intron.close()

# Fonction d'écriture du fichier de sortie pour les transcrits
def write_file_for_transcript(CDS_dic,file_name_for_transcript):

	file_out_transcript=open(file_name_for_transcript,"w")
	for value in CDS_dic.values():
		transcript_object = value
		# Préparation du format BED
		chromosome = transcript_object.chr
		transcript_start = transcript_object.start
		transcript_end = transcript_object.end
		# Écriture des lignes BED obligatoires :
		Bed_format = chromosome+"\t"+transcript_start+"\t"+transcript_end
		# On enregistre les annotations pour les transcrits et les introns
		annotations_transcript=transcript_object.format_print()
		line_out_transcript = Bed_format+"\t"+annotations_transcript
		file_out_transcript.write(line_out_transcript)

	file_out_transcript.close()

# Fonction d'écriture du fichier des fenêtres
def write_file_for_windows_density(five_UTR_density,transcript_for_windows_in_five_UTR,three_UTR_density,transcript_for_windows_in_three_UTR,file_name):
	# On prépare le nom pour les deux fichiers
	file_for_five_prime = "/".join(file_name.split('/')[:-1])+"/5_"+file_name.split('/')[-1]
	file_for_three_prime = "/".join(file_name.split('/')[:-1])+"/3_"+file_name.split('/')[-1]
	# On créer les deux fichiers, un pour chaque région
	file_out_five = open(file_for_five_prime,"w")
	file_out_three = open(file_for_three_prime,"w")
	# On prépare et écris l'en tête
	head = "Begin\tEnd\tTranscripts for windows\tIntrons in windows\n"
	file_out_five.write(head)
	file_out_three.write(head)
	# On lance les deux fonctions d'écritures pour les deux fichiers
	for key,value in sorted(transcript_for_windows_in_five_UTR.items(), reverse=True):
		value_for_key_in_five_UTR_density = five_UTR_density[key]
		line = str(key[0])+"\t"+str(key[1])+"\t"+str(value)+"\t"+str(value_for_key_in_five_UTR_density)+"\n"
		file_out_five.write(line)
	for key,value in sorted(transcript_for_windows_in_three_UTR.items()):
		value_for_key_in_three_UTR_density = three_UTR_density[key]
		line = str(key[0])+"\t"+str(key[1])+"\t"+str(value)+"\t"+str(value_for_key_in_three_UTR_density)+"\n"
		file_out_three.write(line)
	# On ferme les deux fichiers
	file_out_five.close()
	file_out_three.close()

# Interface avec l'utilisateur :
opts, args = getopt.getopt(sys.argv[1:],'',['liste_exon=','liste_intron=','fasta=','output_transcript=','output_intron=','output_windows=',])
for elmts in opts:
	if elmts[0] == '--liste_exon':
		liste_exon = elmts[1] # Fichier contenant les coordonnées pour chaque exon : liste_exon.tab

	elif elmts[0] == '--liste_intron':
		liste_intron = elmts[1] # Fichier contenant les coordonnées pour chaque intron : list_intron.tab

	elif elmts[0] == '--fasta':
		fasta_file = elmts[1] # Fichier contenant les séquences au format fasta des exons et des introns : seq.fa

	elif elmts[0] == '--output_transcript': # Variable qui va contenir le nom de fichier de sortie
		output_file_transcript = elmts[1]

	elif elmts[0] == '--output_intron': # Variable qui va contenir le nom de fichier de sortie
		output_file_intron = elmts[1]

	elif elmts[0] == '--output_windows': # Variable qui va contenir le nom de fichier de sortie
		output_file_windows = elmts[1]

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
PTC_dic = PTC_annotation(transcript_complete)
# On annote également les transcrits
CDS_dic = CDS_annotation(transcript_complete)
# On lance la fonction qui va permettre de calculer la densité des fenêtres UTR
five_UTR_density,transcript_for_windows_in_five_UTR,three_UTR_density,transcript_for_windows_in_three_UTR = intron_density_in_UTR(transcript_complete)
# On lance les fonctions d'écritures :
# -Pour les introns :
write_file_for_intron(PTC_dic,output_file_intron)
# -Pour les transcrits :
write_file_for_transcript(CDS_dic,output_file_transcript)
# Pour les fenêtres des UTRs
write_file_for_windows_density(five_UTR_density,transcript_for_windows_in_five_UTR,three_UTR_density,transcript_for_windows_in_three_UTR,output_file_windows)



####Rajouté ensuite pour calculer la densité des introns dans le CDS
# RAJOUT CDS

# for key,value in transcript_complete.items():
# 	value.CDS_region = (value.cds_start+1,value.cds_stop)
# 	value.CDS_windows = [(c,c+29) for c in cumsum_for_3UTR(value.CDS_region[0],value.CDS_region[1],30) if c+29 <= value.CDS_region[1]]
# 	value.intron_in_CDS = [elmt[0] for elmt in value.intron_pos if (elmt[0] >= value.CDS_region[0] and elmt[0] <= value.CDS_region[1])]


def intron_density_in_CDS(transcript_complete):
	# On définit deux dictionnaires qui contiendrons les fenêtres des régions 5' et 3'
	CDS_density = {}
	transcript_in_CDS = {}
	# Pour chaque transcrit de notre jeu de données
	for trans_id,trans_object in transcript_complete.items():
		CDS_region = trans_object.CDS_region
		CDS_windows = [(elmt[0]-CDS_region[0],elmt[1]-CDS_region[0]) for elmt in trans_object.CDS_windows]
		CDS_introns = [elmt-CDS_region[0] for elmt in trans_object.intron_in_CDS]
		# On construit également le dictionnaire qui contiendra le nombre de transcrit possédant chaque fenêtres
		CDS_density,transcript_for_windows_in_five_UTR = dictionnary_construction(UTR_windows=CDS_windows,UTR_dictionnary=CDS_density,transcript_dictionnary=transcript_in_CDS)
		CDS_density = calcul_density(UTR_density=CDS_density,introns=CDS_introns,windows=CDS_windows,reverse=False)
	return(CDS_density,transcript_in_CDS)

def write_file_for_CDS_density(CDS_density,transcript_in_CDS,file_name):
	# On prépare le nom pour les deux fichiers
	file_CDS = "/".join(file_name.split('/')[:-1])+"/CDS_"+file_name.split('/')[-1]
	# On créer les deux fichiers, un pour chaque région
	file_out = open(file_CDS,"w")
	# On prépare et écris l'en tête
	head = "Begin\tEnd\tTranscripts for windows\tIntrons in windows\n"
	file_out.write(head)
	# On lance les deux fonctions d'écritures pour les deux fichiers
	for key,value in sorted(transcript_in_CDS.items()):
		value_for_key_in_CDS_density = CDS_density[key]
		line = str(key[0])+"\t"+str(key[1])+"\t"+str(value)+"\t"+str(value_for_key_in_CDS_density)+"\n"
		file_out.write(line)
	# On ferme les deux fichiers
	file_out.close()
