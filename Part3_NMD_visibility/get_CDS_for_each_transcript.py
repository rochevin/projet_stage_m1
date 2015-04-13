import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq
class ExonInfo(object):
	"On définit un classe qui contiendra les annotations de chaque exon"
	def __init__(self,exon_id,trans_id,gene_id,coords,seq):
		self.id = exon_id
		self.coords = coords
		self.gene_id = gene_id
		self.trans_id = trans_id
		self.seq = seq

		regex = re.compile("(chr[^\+-]+)([\+-])([0-9]+):([0-9]+)")
		result = regex.findall(self.coords)
		self.chr = result[0][0]
		self.strand = result[0][1]
		self.start = result[0][2]
		self.end = result[0][3]
	def get_seq(self):
		if self.strand == "-":
			return str(self.seq.reverse_complement())
		else:
			return str(self.seq)


class IntronInfo(object):
	"On définit un classe qui contiendra les annotations de chaque intron"
	def __init__(self,id_ens,id_braunch,id_trans,id_gene,coords,seq):
		self.id_ens = id_ens
		self.id_braunch = id_braunch
		self.id_trans = id_trans
		self.id_gene = id_gene
		self.coords = coords # On prend les coordonnées du fichier IntronAnnotate car celle de notre fichier intron_association contiennent parfois des coordonnées issues de chromosomes non assemblés
		self.seq = seq 

		regex = re.compile("^(chr[^\+-]+)([\+-])[0-9]+:([0-9]+)_([0-9]+):[0-9]+")
		result = regex.findall(self.coords)
		self.chr = result[0][0]
		self.strand = result[0][1]
		self.start = result[0][2]
		self.end = result[0][3]
	def get_seq(self):
		if self.strand == "-":
			return str(self.seq.reverse_complement())
		else:
			return str(self.seq)

	def formating_coord_for_fasta(self):
		return self.chr+":"+str(int(self.start)-1)+"-"+self.end

	def formatring_coord_for_sort(self):
		return self.chr+self.strand+self.start+":"+self.end

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


def get_all_cds_for_transcript(file_name,seq_info):

	CDS_by_transcript = {}
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête "Intron	Coordinates	type	txRegion	PTCstatus"

	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t")
		exon_id = content[0]
		exon_seq = seq_info[exon_id]
		exon = ExonInfo(exon_id,content[1],content[2],content[3].replace('\n', ''),exon_seq)
		if exon.trans_id in CDS_by_transcript:
			CDS_by_transcript[exon.trans_id].append((exon.coords,exon))
		else:
			CDS_by_transcript[exon.trans_id]=[]
			CDS_by_transcript[exon.trans_id].append((exon.coords,exon))
	file_in.close()
	return(CDS_by_transcript)

def extract_coord(file_name):
	coord_by_intron = {}
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête 

	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		coord_by_intron[content[0]]=content[1]
	file_in.close()
	return(coord_by_intron)


def get_all_intron_for_transcript(file_name,coord_by_intron,seq_info_for_intron):
	intron_by_transcript = {}
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête

	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") # On split le contenu dans une liste
		ens_id = content[0]
		braunch_id = content[1]
		if braunch_id in coord_by_intron:
			real_coords = coord_by_intron[braunch_id]
			regex = re.compile("^(chr[^\+-]+)[\+-][0-9]+:([0-9]+)_([0-9]+):[0-9]+")
			result = regex.findall(real_coords)
			coords_for_fasta = result[0][0]+":"+str(int(result[0][1])-1)+"-"+result[0][2]
			seq = seq_info_for_intron[coords_for_fasta]
			intron = IntronInfo(ens_id,braunch_id,content[2],content[3],real_coords,seq)
			if intron.id_trans in intron_by_transcript:
				intron_by_transcript[intron.id_trans].append((intron.formatring_coord_for_sort(),intron))
			else:
				intron_by_transcript[intron.id_trans]=[]
				intron_by_transcript[intron.id_trans].append((intron.formatring_coord_for_sort(),intron))
	file_in.close()
	return(intron_by_transcript)

def association_between_intron_and_exon(intron_by_transcript,CDS_by_transcript):
	transcript_complete = {}
	for trans_id, list_exon in CDS_by_transcript.items():
		if trans_id in intron_by_transcript:
			introns = intron_by_transcript[trans_id]
			full_list = []
			full_list.extend(list_exon)
			full_list.extend(introns)
			full_list.sort(key=lambda x: x[0].split(':')[1])
			transcript_complete[trans_id]=full_list
	return(transcript_complete)



opts, args = getopt.getopt(sys.argv[1:],'',['annotation_braunch=','exon_ens_list=','exon_fasta=','intron_association=','intron_fasta=',])
for elmts in opts:
	if elmts[0] == '--annotation_braunch':
		annotation_braunch = elmts[1] 

	elif elmts[0] == '--exon_ens_list':
		exon_ens_list = elmts[1]

	elif elmts[0] == '--exon_fasta':
		exon_fasta = elmts[1]

	elif elmts[0] == '--intron_association':
		intron_association = elmts[1]

	elif elmts[0] == '--intron_fasta':
		intron_fasta = elmts[1]

# On récupère les séquences fasta des exons codants (CDS) et des introns de chaque transcrit :
seq_info_for_cds = get_seq_fasta(exon_fasta) # Les clés seront les id des exons
seq_info_for_intron = get_seq_fasta(intron_fasta) # Les clés seront les coordonnées des introns
# On récupère l'annotation de chaque exon codant qu'on tri par transcrit :
CDS_by_transcript = get_all_cds_for_transcript(exon_ens_list,seq_info_for_cds)
# On récupère les coordonnées de chacun de nos introns de notre jeu de données
coord_by_intron = extract_coord(annotation_braunch)
# On récupère l'annotationn de chaque intron qu'on tri par transcrit
intron_by_transcript = get_all_intron_for_transcript(intron_association,coord_by_intron,seq_info_for_intron)
# On créer le transcrit complet avec les CDS et les introns :
transcript_complete = association_between_intron_and_exon(intron_by_transcript,CDS_by_transcript)
