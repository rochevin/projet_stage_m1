import vcf # Parseur PyVCF
import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import time
class IntronInfo(object):
	"Classe qui contiendra les annotations de chaque intron"
	def __init__(self,intron_id,trans_id,gene_id,chromosome,strand,start,end,donor_site,acceptor_site):
		self.id = intron_id
		self.chr = chromosome
		self.strand = strand
		self.start = start
		self.end = end
		self.gene_id = gene_id
		self.trans_id = trans_id
		self.donor_site = donor_site
		self.acceptor_site = acceptor_site

		self.start_interval = (int(self.start)-20,int(self.start)+19) if self.strand == "+" else (int(self.end)-19,int(self.end)+20)
		self.end_interval = (int(self.end)-19,int(self.end)+20) if self.strand == "+" else (int(self.start)-20,int(self.start)+19)

class IntronInfoPoly(object):
	"Classe contenant l'information sur le polymorphisme"
	def __init__(self,intron_id,chromosome,strand,start,end,donor_poly,acceptor_poly):
		self.intron_id
		self.chromosome
		self.strand
		self.start
		self.end
		self.donor_poly = donor_poly
		self.acceptor_poly = acceptor_poly

class PolyPos(object):
	"Classe contenant l'annotation sur le polymorphisme d'une position"
	def __init__(self,position_type,pos,ref,alt,der,anc,af,daf):
		self.type = position_type
		self.pos = pos
		self.ref = ref
		self.alt = alt
		self.der = der
		self.anc = anc
		self.af = af
		self.daf = daf



def get_pos_and_coords_by_intron(file_name):
	chromosome_list={} # On définit un dictionnaire qui contiendra les positions des introns par chromosome
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête

	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") # On split la ligne dans une liste
		Intron_Object = IntronInfo(intron_id=content[3], trans_id=content[5], gene_id=content[6], chromosome=content[0], strand=content[4], start=content[1], end=content[2], donor_site=content[7], acceptor_site=content[8])
		if content[0] in chromosome_list:
			chromosome_list[content[0]].append(Intron_Object)
		else:
			chromosome_list[content[0]]=[Intron_Object]

	for key,value in chromosome_list.items():
			value.sort(key=lambda x: int(x.start))

	return(chromosome_list)





chromosome_list = get_pos_and_coords_by_intron('/Users/Vincent/Documents/STAGE_LYON/Projet_stage_RV/results/result_human/result_primary_v2.bed')


tps1 = time.clock()
file_out = open('result.txt','w')

for key,value in chromosome_list.items():
	file_name = key+".vcf.gz"
	try:
		vcf_reader = vcf.Reader(filename=file_name)
		print('Fichier',file_name,"ouvert")
	except:
		print('Fichier',file_name,"introuvable")
		continue
	chromosome = key[3:]

	for intron in value:
		intervals = (intron.start_interval,intron.end_interval)
		for interval in intervals:
			polymorphism = vcf_reader.fetch(chromosome, interval[0], interval[1])
			if polymorphism:
				for record in polymorphism:
					if record.is_snp and len(record.ALT)==1:
						string = str(intron.id)+"\t"+chromosome+"\t"+str(record.POS)+"\t"+str(record.REF)+"\t"+str(record.ALT[0])+"\t"+str(record.INFO['AF'][0])+"\t"+str(record.INFO['AA'])[0]+"\n"
						file_out.write(string)

file_out.close()
tps2 = time.clock()
print(tps2 - tps1)

