from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq

class IntronInfo(object):
	"Classe qui contiendra les annotations de chaque intron"
	def __init__(self,intron_id,trans_id,gene_id,chromosome,strand,start,end,GC_rate,start_site,end_site,donor_site,acceptor_site):
		self.id = intron_id
		self.chr = chromosome
		self.strand = strand
		self.start = start
		self.end = end
		self.gene_id = gene_id
		self.trans_id = trans_id
		self.donor_site = donor_site
		self.acceptor_site = acceptor_site
		self.GC_rate = GC_rate
		#Les sequences start/end sont déjà inversée dans le fichier fasta, donc on part toujours de start, par contre, pour les positions, si brin -, le début de l'intron correspond à la position end dans le génome
		self.start_site = start_site
		self.end_site = end_site

		self.start_interval = (int(self.start)-20,int(self.start)+19) if self.strand == "+" else (int(self.end)-19,int(self.end)+20)
		self.end_interval = (int(self.end)-19,int(self.end)+20) if self.strand == "+" else (int(self.start)-20,int(self.start)+19)

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

class Polymorphism(object):
	"Classe qui contiendra un SNV à une position données"
	def __init__(self,chromosome,position,REF,ALT,Frq_ALT,ANC):
		self.chromosome = chromosome
		self.pos = position
		self.REF = REF
		self.ALT = ALT
		self.Frq_ALT = Frq_ALT
		self.ANC = ANC


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

def get_pos_and_coords_by_intron(file_name,dico_info):
	intron_list=[] # On définit une liste qui contiendra les introns
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête

	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") # On split la ligne dans une liste

		seq_of_interest = dico_info[content[3]]
		start_site = seq_of_interest[0][0]
		end_site = seq_of_interest[1][0]
		Intron_Object = IntronInfo(intron_id=content[3], trans_id=content[5], gene_id=content[6], chromosome=content[0], strand=content[4], start=content[1], end=content[2], donor_site=content[7], acceptor_site=content[8],GC_rate=content[12],start_site=start_site,end_site=end_site)
		intron_list.append(Intron_Object)


	return(intron_list)

def get_polymorphism(file_name):
	poly_info = {}
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture

	for line in file_in:
		content = line.split("\t")
		poly_object = Polymorphism(chromosome=content[1],position=content[2],REF=content[3],ALT=content[4],Frq_ALT=content[5],ANC=content[6].replace("\n",""))
		intron_id = content[0]
		poly_info[intron_id,poly_object.pos]=poly_object

	file_in.close
	return(poly_info)

###Fonctions de calcul du score

###Détermination de la séquence alternative avec le polymorphisme
def splice_site_preparation(alt,seq,alt_pos,strand):
	seq_alt = seq[:alt_pos]+alt+seq[alt_pos+1:]
	return seq_alt if strand =="+" else str(Seq(seq_alt).reverse_complement())

###Construction des séquences pour la matrice
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

###Construction de la matrice
def get_matrix_frequency(seqs):
	m = motifs.create(seqs)
	PFM = m.pwm
	PPM = {}
	for key,value in PFM.items():
		percentages = [elmt*100 for elmt in value]
		PPM[key]=percentages
	return(PPM)

###Calcul du score donneur
def donor_score(matrix,seq):
	sub_seq = seq[18:26]
	maxt=sum(max(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(8))
	mint=sum(min(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(8))
	t=sum(matrix[seq[i]][i] for i in range(len(sub_seq)))
	score=100*((t-mint)/(maxt-mint))
	return(round(score,2))

###Calcul du score accepteur
def acceptor_score(matrix,seq):
	sub_seq = seq[7:21]
	t1 = sum(sorted([matrix[sub_seq[i]][i] for i in range(10)],reverse=True)[0:8])
	t2 = sum(matrix[sub_seq[i]][i] for i in range(10,14))
	l1 = sum(sorted([min(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(10)])[0:8])
	l2 = sum(sorted([min(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(10,14)]))
	h1 = sum(sorted([max(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(10)],reverse=True)[0:8])
	h2 = sum(sorted([max(matrix['A'][i],matrix['C'][i],matrix['G'][i],matrix['T'][i]) for i in range(10,14)]))
	score = 100*((t1-l1)/(h1-l1)+(t2-l2)/(h2-l2))/2
	return(round(score,2))


def write_annotation_polymorphism(file_name,intron_list,poly_info):
	file_out = open(file_name,"w")

	nuc_type = {"A":"W","T":"W","C":"S","G":"S"} #On construit un dico qui contient l'information sur le type de nucléotide, on pourra alors déterminer le type de mutation (ex : si A vers T => WW, si A vers G => WS)

	###Détermination de la matrice totale
	all_donors,all_acceptors = type_of_data(introns,type_data="All")
	all_donor_matrix = get_matrix_frequency(all_donors)
	all_acceptor_matrix = get_matrix_frequency(all_acceptors)


	for intron in intron_list:
		intron_id = intron.id
		splice_donor_interval = intron.start_interval
		splice_acceptor_interval = intron.end_interval

		###Calcul du score des sites de l'intron sans polymorphisme
		score_for_donor_seq = donor_score(all_donor_matrix,intron.get_donor())
		score_for_acceptor_seq = acceptor_score(all_acceptor_matrix,intron.get_acceptor())


		position_in_splice_site = 0 #Compteur pour recupérer la sequence de référence quand pas annoté d'un SNP
		for i in range(splice_donor_interval[0],splice_donor_interval[1]+1):
			if (intron_id,str(i)) in poly_info:
				polymorphism = poly_info[intron_id,str(i)]
				seq_alt = splice_site_preparation(alt=polymorphism.ALT,seq=intron.start_site,alt_pos=position_in_splice_site,strand=intron.strand)
				alt_score_for_donor_seq = donor_score(all_donor_matrix,seq_alt)
				delta_score = score_for_donor_seq-alt_score_for_donor_seq
				if polymorphism.ANC not in (".","-","N"):
					ANC = polymorphism.ANC.upper()
					DER = polymorphism.ALT if ANC == polymorphism.REF else polymorphism.REF
					DAF = polymorphism.Frq_ALT if ANC == polymorphism.REF else 1-float(polymorphism.Frq_ALT)
					type_mut = nuc_type[ANC]+nuc_type[DER]
					ANC_Qual = "2" if ANC.isupper() else "1"
				else:
					ANC = "NA"
					DER = "NA"
					DAF = "NA"
					type_mut = "NA"
					ANC_Qual = "0"
				line = intron.chr+"\t"+str(i)+"\t"+"Donor_site"+"\t"+str(position_in_splice_site+1)+"\t"+intron_id+"\t"+polymorphism.REF+"\t"+polymorphism.ALT+"\t"+ANC+"\t"+DER+"\t"+str(polymorphism.Frq_ALT)+"\t"+str(DAF)+"\t"+type_mut+"\t"+ANC_Qual+"\t"+delta_score+"\n"
			else:
				line = intron.chr+"\t"+str(i)+"\t"+"Donor_site"+"\t"+str(position_in_splice_site+1)+"\t"+intron_id+"\t"+intron.start_site[position_in_splice_site]+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\n"
			file_out.write(line)
			position_in_splice_site +=1
		position_in_splice_site = 0 # On remet le compteur à 0 pour le parcours de la séquence du site accepteur
		for i in range(splice_acceptor_interval[0],splice_acceptor_interval[1]+1):
			if (intron_id,str(i)) in poly_info:
				polymorphism = poly_info[intron_id,str(i)]
				seq_alt = splice_site_preparation(alt=polymorphism.ALT,seq=intron.end_site,alt_pos=position_in_splice_site,strand=intron.strand)
				alt_score_for_acceptor_seq = acceptor_score(all_donor_matrix,seq_alt)
				delta_score = score_for_acceptor_seq-alt_score_for_acceptor_seq
				if polymorphism.ANC not in (".","-","N"):
					ANC = polymorphism.ANC.upper()
					DER = polymorphism.ALT if ANC == polymorphism.REF else polymorphism.REF
					DAF = polymorphism.Frq_ALT if ANC == polymorphism.REF else 1-float(polymorphism.Frq_ALT)
					type_mut = nuc_type[ANC]+nuc_type[DER]
					ANC_Qual = "2" if ANC.isupper() else "1"
				else:
					ANC = "NA"
					DER = "NA"
					DAF = "NA"
					type_mut = "NA"
					ANC_Qual = "0"
				line = intron.chr+"\t"+str(i)+"\t"+"Acceptor_site"+"\t"+str(position_in_splice_site+1)+"\t"+intron_id+"\t"+polymorphism.REF+"\t"+polymorphism.ALT+"\t"+ANC+"\t"+DER+"\t"+str(polymorphism.Frq_ALT)+"\t"+str(DAF)+"\t"+type_mut+"\t"+ANC_Qual+"\n"
			else:
				line = intron.chr+"\t"+str(i)+"\t"+"Acceptor_site"+"\t"+str(position_in_splice_site+1)+"\t"+intron_id+"\t"+intron.start_site[position_in_splice_site]+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\n"
			file_out.write(line)
			position_in_splice_site +=1



	file_out.close()
