from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq
import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python

#Création de la classe contenant les annotations d'ensembl sur les exons des transcrits canoniques
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

		def formating_coord_for_intron(self): # Converti les coordonnées de l'intron pour retrouver les mêmes coordonnées lors de la construction d'intron d'Ensembl
			return self.chr+":"+str(self.start)+"-"+str(self.end)

# Classe contenant les annotations sur l'intron d'un transcrit Braunch 
class IntronInfoAnnotate(object):
	"Classe contenant les données des introns annotées générées via script précédent"
	def __init__(self,Intron_id,Intron_coord,Intron_type,Intron_txRegion,Intron_PTCstatus,seq_annotation_begin,seq_annotation_end,intron_len,intron_position,total_intron,GCrate):
		self.id = Intron_id # Contient l'identifiant de l'intron
		self.coords = Intron_coord # Contient les coordonnées de l'intron
		#Expression régulière qui récupère le nom du chromosome, 
		#les coordonnées de la dernière base de l'exon au site d'épissage 5'
		#les coordonnées de la première base de l'exon au site d'épissage 3'
		regex = re.compile('^(chr[\w]+)([-\+])[0-9]+:([0-9]+)_([0-9]+):[0-9]+')
		#On enregistre les résultats dans une liste (deux dimensions)
		result = regex.findall(self.coords)
		self.chr, self.signe, self.debut, self.fin = [result[0][0],result[0][1],int(result[0][2]),int(result[0][3])-1]
		#L'annotation étant en référentiel 1, on laisse les coordonnées de la dernière base de l'exon flanquant pour déterminer le début de l'intron


		self.type = Intron_type # Contient le type de l'intron
		self.txregion = Intron_txRegion # Contient txRegion
		self.PTC = Intron_PTCstatus # Contient le status PTC
		self.begin_annotation = seq_annotation_begin # Contient la séquence d'intéret de début
		self.end_annotation = seq_annotation_end # Contient la séquence d'intéret de fin
		self.len = intron_len # Contient la taille totale de l'intron
		self.pos = intron_position # Contient la position de l'intron par rapport aux autres dans le transcrit
		self.total = total_intron # Nombre total d'intron dans le transcrit
		self.GCrate = GCrate # Taux de GC de l'intron sans les séquences d'intéret

	def formating_coord(self): # Converti les coordonnées de l'intron en format num_chr:start-stop
		return self.chr+":"+str(self.debut)+"-"+str(self.fin)

	def formating_coord_for_exon(self): # Converti les coordonnées de l'intron pour retrouver les mêmes coordonnées lors de la construction d'intron d'Ensembl
		return self.chr+":"+str(self.debut)+"-"+str(self.fin+1)

	def flank_exon_left(self): # Permet d'obtenir les coordonnées de l'exon flanquant gauche
		regex = re.compile('^(chr[\w]+)[-\+]([0-9]+):([0-9]+)_[0-9]+:[0-9]+')
		result = regex.findall(self.coords)
		coord_exon_left = result[0][0]+":"+result[0][1]+"-"+result[0][2]
		return coord_exon_left
	def flank_exon_right(self): # Permet d'obtenir les coordonnées de l'exon flanquant droite
		regex = re.compile('^(chr[\w]+)[-\+][0-9]+:[0-9]+_([0-9]+):([0-9]+)')
		result = regex.findall(self.coords)
		coord_exon_right = result[0][0]+":"+result[0][1]+"-"+result[0][2]
		return coord_exon_right
	def BED_coord(self): # Converti les coordonnées au format BED
		return self.chr+"\t"+str(self.debut)+"\t"+str(self.fin)+"\t"+self.id

#Fonction qui extrait les annotations des exons d'ensembl pour les enregistrer dans un objet exon
def extract_exon_annotation(file_name):
	Dataset_Ensembl = {}
	file_in=open(file_name,"r")
	header = file_in.readline() # On enregistre l'en-tête
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		line.replace('\n', '')
		content=line.split("\t")
		exon = ExonInfo(content[0],content[1],content[2],content[3],content[4],content[5],content[6],content[7],content[8])
		exon_id = exon.ref_id+":"+exon.ens_id
		if exon_id in Dataset_Ensembl:
			Dataset_Ensembl[exon_id].append(exon)
		else:
			Dataset_Ensembl[exon_id]=[]
			Dataset_Ensembl[exon_id].append(exon)

	file_in.close()
	return(Dataset_Ensembl)

def extract_braunch_annotation(file_name):
	Dataset_Braunch = {} # dictionnaire qui va contenir les objets issues de IntronInfoAnnotate
	Braunch_intron_coord = {}
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête "Intron	Coordinates	type	txRegion	PTCstatus"

	regex = re.compile('^[^:]+:.+:(N[RM]_.+):[0-9]{,3}') # On définit une expression régulière qui va capturer le nom du gène et l'id refseq pour définir un transcrit Braunch
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		intron_name = content[0] #On récupère le nom de l'intron
		#On vérifie que l'intron est annoté d'un gène, sinon cela ne nous sert à rien de l'enregistrer

		for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
			content=line.split("\t") #On split le contenu dans une liste
			#On enregistre les annotations dans un objet
			intron_name = IntronInfoAnnotate(content[0],content[1],content[2],content[3],content[4],content[5],content[6],content[7],content[8],content[9],content[10].replace('\n', ''))
			#On enregistre l'objet dans le dictionnaire avec l'identifiant de transcrit Braunch comme clé
			result = regex.findall(intron_name.id)
			if result[0] in Dataset_Braunch:
				Dataset_Braunch[result[0]].append(intron_name)
			else:
				Dataset_Braunch[result[0]]=[]
				Dataset_Braunch[result[0]].append(intron_name)
			#Enregistrement des coordonnées des introns par transcrit :
			if result[0] in Braunch_intron_coord:
				Braunch_intron_coord[result[0]].append(intron_name.formating_coord_for_exon())
			else:
				Braunch_intron_coord[result[0]]=[]
				Braunch_intron_coord[result[0]].append(intron_name.formating_coord_for_exon())
	file_in.close() 

	return(Dataset_Braunch,Braunch_intron_coord)


### Fonction qui va construire les introns des transcrits ensembl à partir des coordonnées des exons flanquants
### Afin d'associer un transcrit Braunch à un transcrit Ensembl
def transform_exon_to_intron(Dataset_Ensembl):
	Ensembl_intron_coord = {} # Dictionnaire qui va transformer les coordonnées exoniques en introns afin de comparer les introns
	retrieve_transcript_with_coord = {} # Dictionnaire de coordonnées qui va enregistrer une coordonnée et à quel(s) transcrit(s) elle appartient
	for transcript,liste_exon in Dataset_Ensembl.items():
		for exon in liste_exon:
			# Si notre élément n'est pas le dernier, on peut continuer
			if liste_exon.index(exon)+1 != len(liste_exon):
				# On récupère l'exon suivant celui qu'on parcours
				exon_left = exon
				exon_right = liste_exon[liste_exon.index(exon)+1]
				
				position_exon_left = exon_left.end
				position_exon_right = exon_right.start
				chromosome = exon_left.chr

				intron_coord = chromosome+":"+str(position_exon_left)+"-"+str(position_exon_right)

				if transcript in Ensembl_intron_coord:
					Ensembl_intron_coord[transcript].append(intron_coord)
				else:
					Ensembl_intron_coord[transcript]=[]
					Ensembl_intron_coord[transcript].append(intron_coord)

	return(Ensembl_intron_coord)


def verify_canonical_transcript(Ensembl_intron_coord,Dataset_Ensembl,Braunch_intron_coord,Dataset_Braunch):
	match=nomatch=0
	for transcript,Ensembl_liste_introns in Ensembl_intron_coord.items():
		regex = re.compile('^(.+):.+')
		result = regex.findall(transcript)
		if result[0] in Braunch_intron_coord:
			Braunch_liste_introns = Braunch_intron_coord[result[0]]
			Ensembl_liste_introns.sort(key=lambda x: x.split('-')[1]) 
			Braunch_liste_introns.sort(key=lambda x: x.split('-')[1])
			if Braunch_liste_introns == Ensembl_liste_introns:
				print('COUCOU')
				match+=1
			else:
				print('PAS COUCOU')
				print(transcript)
				print('Ensembl :',Ensembl_liste_introns)
				print('Braunch :',Braunch_liste_introns)
				nomatch+=1
	print('match =',match)
	print('no match =',nomatch)


file_name_dataset_ensembl = sys.argv[1]
file_name_dataset_braunch = sys.argv[2]

Dataset_Ensembl = extract_exon_annotation(file_name_dataset_ensembl)
Dataset_Braunch,Braunch_intron_coord = extract_braunch_annotation(file_name_dataset_braunch)

Ensembl_intron_coord = transform_exon_to_intron(Dataset_Ensembl)

verify_canonical_transcript(Ensembl_intron_coord,Dataset_Ensembl,Braunch_intron_coord,Dataset_Braunch)









	



