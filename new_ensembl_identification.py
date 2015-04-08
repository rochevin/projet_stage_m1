from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq
import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python

class ExonInfo(object):
	"Classe contenant les données de chaque exons issus du transcrit Ensembl (canonique dans la pluspart des cas)"
	def __init__(self,refseq,gene,transcript,exon,number,chromosome,strand,coords,canonical):
		self.ref_id = refseq
		self.gene_id = gene
		self.ens_id = transcript
		self.id = exon 
		self.chr = chromosome
		self.strand = strand
		self.coords = coords
		regex = re.compile('^([0-9]+)_([0-9]+)')
		#On enregistre les résultats dans une liste (deux dimensions)
		result = regex.findall(self.coords)
		self.start = result[0][0]
		self.end = result[0][1]

def extract_exon_annotation(file_name):
	Dataset_Ensembl = {}
	file_in=open(file_name,"r")
	header = file_in.readline() # On enregistre l'en-tête
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		line.replace('\n', '')
		content=line.split("\t")
		gene_id = content[1]
		exon = ExonInfo(content[0],content[1],content[2],content[3],content[4],content[5],content[6],content[7],content[8])
		if gene_id in Dataset_Ensembl:
			Dataset_Ensembl[gene_id].append(exon)
		else:
			Dataset_Ensembl[gene_id]=[]
			Dataset_Ensembl[gene_id].append(exon)