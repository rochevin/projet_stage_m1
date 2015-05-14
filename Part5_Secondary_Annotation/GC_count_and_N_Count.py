import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq
#Définition des classes


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

def get_seq_of_interest(file_name):
    dico_info = {}
    handle = open(file_name, "rU") #On ouvre le fichier en mode lecture
    for record in SeqIO.parse(handle, "fasta") : #On parse le fichier avec seqIO
        seq_id = record.id #On enregistre le nom de la sequence, ici les coordonnées
        seq = record.seq #On enregistre la sequence fasta dans une variable
        seq = seq.upper() #On met en majuscule au cas ou la séquence ne l'est pas
        if seq_id in dico_info:
            dico_info[seq_id].append((seq,record.description.split(" ")[1]))
        else:
            dico_info[seq_id] = []
            dico_info[seq_id].append((seq,record.description.split(" ")[1]))
    handle.close()
    return(dico_info)

def GC_rate_for_seq_region(dico_info):
	GC_info = {}
	for key,value in dico_info.items():
		seq_region="".join([str(elmt[0]) for elmt in value])
		A= seq_region.count('A')
		C = seq_region.count('C')
		G = seq_region.count('G')
		T = seq_region.count('T')
		N = seq_region.count('N')
		# On détermine une valeur minimale de nucleotides pour calculer le taux de GC
		#Calcul du taux de GC
		GCrate = ((G+C)/(A+T+G+C))*100
		GCrate = round(GCrate,2)
		GC_info[key]=[GCrate,N]
	return GC_info

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

def paste_file(file_name,intron_by_transcript,GC_info):
	file_out=open(file_name,'w')
	header = "Chromosome\tstart\tend\tIntron_id\tTranscript_id\tGene_id\tN_count\tGC_flanking\tN_flanking\n"
	file_out.write(header)

	for key,value in intron_by_transcript.items():
		for intron in value:
			trans_id = key.split('(')[0]
			N_in_intron = intron[1].seq.count('N')
			N_and_GC_for_region = GC_info[trans_id]
			GC_for_region_seq = N_and_GC_for_region[0]
			N_for_region_seq = N_and_GC_for_region[1]
			line = intron[1].chr+"\t"+intron[1].start+"\t"+intron[1].end+"\t"+intron[1].id+"\t"+trans_id+"\t"+intron[1].gene_id+"\t"+str(N_in_intron)+"\t"+str(GC_for_region_seq)+"\t"+str(N_for_region_seq)+"\n"
			file_out.write(line)
	file_out.close()


