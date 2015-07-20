from Bio.Seq import Seq
from Bio import motifs

class IntronInfo(object):
	"Classe qui contiendra les annotations de chaque intron"
	def __init__(self,intron_id,trans_id,gene_id,donor_site,acceptor_site,chromosome,start,end,GC_rate):
		self.id = intron_id
		self.gene_id = gene_id
		self.trans_id = trans_id
		self.donor = donor_site
		self.acceptor = acceptor_site
		self.chr = chromosome
		self.start = start
		self.end = end
		self.GC_rate = GC_rate
	def get_donor(self):
		if self.donor !="NA":
			return(self.donor[18:26])
		else:
			return("NA")
	def get_acceptor(self):
		if self.acceptor !="NA":
			return(self.acceptor[7:21])
		else:
			return("NA")

def get_seq_for_each_intron(file_name):
	introns = {}
	donors = []
	acceptors = []
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête

	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t")
		intron_id = content[3]
		intron_chr = content[0]
		intron_start = content[1]
		intron_end = content[2]
		intron_gene_id = content[5]
		intron_trans_id = content[4]
		donor_seq = Seq(content[6])
		acceptor_seq = Seq(content[7])
		intron_GC_rate = content[11]
		introns[intron_id]=IntronInfo(intron_id,intron_trans_id,intron_gene_id,donor_seq,acceptor_seq,intron_chr,intron_start,intron_end,intron_GC_rate)

	file_in.close()

	return(introns)

def get_matrix_frequency(seqs):
	m = motifs.create(seqs)
	PFM = m.pwm
	PPM = {}
	for key,value in PFM.items():
		percentages = [elmt*100 for elmt in value]
		PPM[key]=percentages
	return(PPM)

def donor_score(matrix,seq):
	maxt=sum(max(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(8))
	mint=sum(min(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(8))
	t=sum(matrix[seq[i]][i] for i in range(len(seq)))
	score=100*((t-mint)/(maxt-mint))
	return(round(score,2))

def acceptor_score(matrix,seq):
	t1 = sum(sorted([matrix[seq[i]][i] for i in range(10)],reverse=True)[0:8])
	t2 = sum(matrix[seq[i]][i] for i in range(10,14))
	l1 = sum(sorted([min(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(10)])[0:8])
	l2 = sum(sorted([min(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(10,14)]))
	h1 = sum(sorted([max(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(10)],reverse=True)[0:8])
	h2 = sum(sorted([max(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(10,14)]))
	score = 100*((t1-l1)/(h1-l1)+(t2-l2)/(h2-l2))/2
	return(round(score,2))


def type_of_data(object_dictionnary,type_data="All"):
	if type_data == "All":
		donor_seq = [value.get_donor() for value in object_dictionnary.values() if (value.GC_rate!="NA" and value.get_donor().count('N')==0)]
		acceptor_seq = [value.get_acceptor() for value in object_dictionnary.values() if (value.GC_rate!="NA" and value.get_acceptor().count('N')==0)]
	if type_data == "low_gc":
		donor_seq = [value.get_donor() for value in object_dictionnary.values() if (value.GC_rate!="NA" and float(value.GC_rate)<=35 and value.get_donor().count('N')==0)]
		acceptor_seq = [value.get_acceptor() for value in object_dictionnary.values() if (value.GC_rate!="NA" and float(value.GC_rate)<=35 and value.get_acceptor().count('N')==0)]
	if type_data == "high_gc":
		donor_seq = [value.get_donor() for value in object_dictionnary.values() if (value.GC_rate!="NA" and float(value.GC_rate)>=60 and value.get_donor().count('N')==0)]
		acceptor_seq = [value.get_acceptor() for value in object_dictionnary.values() if (value.GC_rate!="NA" and float(value.GC_rate)>=60 and value.get_acceptor().count('N')==0)]

	return(donor_seq,acceptor_seq)

def score_for_each_intron(file_name,introns):
	file_out=open(file_name,"w")

	all_donors,all_acceptors = type_of_data(introns,type_data="All")
	low_gc_donors,low_gc_acceptors = type_of_data(introns,type_data="low_gc")
	high_gc_donors,high_gc_acceptors = type_of_data(introns,type_data="high_gc")


	all_donor_matrix = get_matrix_frequency(all_donors)
	all_acceptor_matrix = get_matrix_frequency(all_acceptors)
	low_gc_donor_matrix = get_matrix_frequency(low_gc_donors)
	low_gc_acceptor_matrix = get_matrix_frequency(low_gc_acceptors)
	high_gc_donor_matrix = get_matrix_frequency(high_gc_donors)
	high_gc_acceptor_matrix = get_matrix_frequency(high_gc_acceptors)


	head = "chromosome\tstart\tend\tIntron_id\tAll_Donor_score\tAll_Acceptor_score\tlow_gc_Donor_score\tlow_gc_Acceptor_score\thigh_gc_Donor_score\thigh_gc_Acceptor_score\n"
	file_out.write(head)
	for key,value in introns.items():
		intron=value
		if (intron.get_donor() != "NA" and intron.get_donor().count('N')==0):
			all_score_for_donor_seq = donor_score(all_donor_matrix,intron.get_donor())
			low_gc_score_for_donor_seq = donor_score(low_gc_donor_matrix,intron.get_donor())
			high_gc_score_for_donor_seq = donor_score(high_gc_donor_matrix,intron.get_donor())
		else:
			score_for_donor_seq = "NA"
		if (intron.get_acceptor() != "NA" and intron.get_acceptor().count('N')==0):
			all_score_for_acceptor_seq = acceptor_score(all_acceptor_matrix,intron.get_acceptor())
			low_gc_score_for_acceptor_seq = acceptor_score(low_gc_acceptor_matrix,intron.get_acceptor())
			high_gc_score_for_acceptor_seq = acceptor_score(high_gc_acceptor_matrix,intron.get_acceptor())
		else:
			score_for_acceptor_seq = "NA"
		line = intron.chr+"\t"+intron.start+"\t"+intron.end+"\t"+intron.id+"\t"+str(all_score_for_donor_seq)+"\t"+str(all_score_for_acceptor_seq)+"\t"+str(low_gc_score_for_donor_seq)+"\t"+str(low_gc_score_for_acceptor_seq)+"\t"+str(high_gc_score_for_donor_seq)+"\t"+str(high_gc_score_for_acceptor_seq)+"\n"
		file_out.write(line)

	file_out.close()
