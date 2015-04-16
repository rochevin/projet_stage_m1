# Programme qui va récupérer la liste des exons codants et des introns pour chaque transcrit canonique ensembl
# A partir de ces annotations, va organiser chaque transcrit canonique de façon à ce qu'il contienne tous ses exons/introns
# A partir de cette liste, va construire le CDS de chaque transcrit, puis introduire un à un chaque intron 
# Pour chaque intron, on déterminera si celui-ci est NMD visible, et si il décale la phase, on indiquera donc la phase de l'intron


import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq
from collections import Counter # On importe Counter qui va construire un dictionnaire pour compter le nombre d'occurence dans nos codons
#Définition des classes

class TranscriptInfo(object):
	"Classe qui contiendra les annotations de chaque exon"
	def __init__(self,transcript_id,gene_id,pos_start,pos_stop,list_object_exon,list_object_intron):
		self.id = transcript_id
		self.cds_start = int(pos_start)-1
		self.cds_stop = int(pos_stop)
		self.exons = [(exon.start,exon.stop,exon) for exon in list_object_exon]
		self.introns = [(intron.start,intron.stop,intron) for intron in list_object_intron]
	
	def strand(self):
		return self.exons[0][2].strand

	def transcript_seq(self):
		return "".join([exon[2].seq for exon in self.exons])

	def CDS_seq(self):
		return transcript_seq()[self.cds_start:self.cds_stop]

	def CDS_control(self):
		complete_CDS = self.CDS_seq()
		PTC_pos = None
		# On définit un tuple contenant tous les codons stop possible
		stop_codon = ("TGA","TAG","TAA")
		# On vérifie qu'il ne contient pas un codon stop prématuré
		codons = [complete_CDS[i:i+3] for i in range(0,len(complete_CDS),3)]
		index_codons = [ i for i in range(0,len(codons),1) if codons[i] in stop_codon]
		stop_control = len(index_codons)
		if stop_control > 1 : 
			status = "PTC"
			PTC_pos = index_codons[:-1]
		# On vérifie qu'il commence par un start et fini par un stop :
		if not (complete_CDS.endswith(stop_codon) or complete_CDS.startswith("ATG")) : 
			status = "Partial"
		# On vérifie qu'il est bien multiple de 3 :
		if not len(complete_CDS)%3 == 0 : 
			status = "Partial"
		else:
			status = "OK"
		
		return(status,PTC_pos)

	def transcript_content(self):
		full_list = []
		full_list.extend(self.exons)
		full_list.extend(self.introns)
		# Fonction de tri, on a des coordonnées celon la conventien : chr[+/-]debut:fin
		# On va donc prendre en compte uniquement la fin, donc on split la chaine pour obetenir fin, qu'on converti en entier, et qu'on tri dans l'ordre
		full_list.sort(key=lambda x: int(x[0].split(':')[1]))
		return full_list

	def exon_position(self,one_exon):
		start_list = [c for c in cumsum(self.exons)]
		end_list = [start[i]+len(self.exons[i][2].seq)-1 for i in range(0,len(self.exons),1)]
		return [(start_list[i],end_list[i],self.exons[i][2]) for i in range(0,len(self.exons),1)]

	def intron_start(self,one_intron):
		full_transcript = self.transcript_content()
		return exon_position(full_transcript.index(one_intron))[1]


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


def cumsum(liste):
	s = 0
	yield s
	for elmt in liste:
		s += len(elmt[2].seq)
		yield s
