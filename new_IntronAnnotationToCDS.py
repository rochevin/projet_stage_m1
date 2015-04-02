####################################################################
#Auteur : ROCHER Vincent										   #
#																   #
#But : Programme principal										   #
#	   -Récupère les transcrits de référence					   #
#	   -Recherche les CDS de chaque genes associés aux introns	   #
#	   -Utilise les coordonnées produites pour créer des séquences #
#		au format fasta pour les CDS et introns					   #
#	   -Parcours chaque CDS et vérifie la présence ou non de nos   #
#		introns -> si oui vérifie si intron NMD visible			   #
####################################################################

from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq
import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
import tempfile # On importe tempfile pour créer des fichiers temporaires
import gffutils # On importe gffutils pour parser le fichier GTF

# Classe contenant les annotations sur les transcrits de référence
class CanonicalTranscript(object):
	"Classe contenant les annotations de chaque transcrit Braunch référencé par un transcrit Ensembl"
	def __init__(self,Braunch_id,Ensembl_id,Gene_id,Chr,Strand,Type,Intron_number,Exon_number):
		self.B_id = Braunch_id
		self.E_id = Ensembl_id
		self.gene_id = Gene_id
		self.chr = Chr
		self.strand = Strand
		self.type = Type
		self.intron_number = Intron_number
		self.exon_number = Exon_number

	def formating_id(self):
		return self.chr+self.strand+self.B_id+":"+self.E_id

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
		return self.chr+"\t"+str(self.debut-1)+"\t"+str(self.fin)+"\t"+self.id


# Classe contenant les annotations sur l'exon d'un transcrit Ensembl
class ExonInfo(object):
	"Classe contenant les données de chaque CDS (exon sans 5'-UTR/3' UTR)"
	def __init__(self,cds_chr,cds_start,cds_stop,cds_strand,cds_transcript,cds_exon_number,gene_id,feature_type,source):
		# Identifiant de l'exon construit à partir de l'id du transcrit, de sa position dans le transcrit, ses coordonnées et son type
		self.id = cds_transcript+":"+cds_exon_number+":"+str(cds_start)+"-"+str(cds_stop)+":"+feature_type 
		self.chr = cds_chr # Chromosome qui contient l'exon
		self.start = cds_start # Position de début de séquence dans le chromosome (ref 1)
		self.stop = cds_stop # Position de fin de séquence dans le chromosome (ref 1)
		self.strand = cds_strand # Brin auquel appartient l'exon (+/-)
		self.transcript = cds_transcript # Transcrit auquel fait parti l'exon
		self.exon_number = cds_exon_number # Position de l'exon dans le transcrit par rapport aux autres exons
		self.gene_id = gene_id # Identifiant du gène auquel appartient le transcrit de l'exon
		self.feature_type = feature_type # Type de l'exon (CDS, exon, stop_codon)
		self.source = source # Type du transcrit -> non_sense, protein_coding

	def formating_coord(self): # Fonction pour afficher les coordonnées de l'exon d'une façon spécifique num_chr:start-stop
		return "chr"+self.chr+":"+str(self.start)+"-"+str(self.stop)

	def BED_coord(self): # Converti les coordonnées au format BED
		return self.chr+"\t"+str(self.start-1)+"\t"+str(self.stop)+"\t"+self.id

### On récupère les transcrits de référence
def get_canonical_transcript(file_name):		
	canonical_transcript = {} # Dictionnaire qui va contenir les transcrits sous forme d'objet
	canonical_list = {} # Va contenir l'association entre l'identifiant Braunch et Ensembl
	file_in = open(file_name,"r")
	header = file_in.readline()

	for line in file_in:
		result = line.split("\t") #On split le contenu dans une liste
		canonical = CanonicalTranscript(result[0],result[1],result[2],result[3],result[4],result[5],result[6],result[7])
		canonical_transcript[canonical.B_id] = canonical # On enregistre l'objet dans le dictionnaire
		canonical_list[canonical.B_id] = canonical.E_id # On associe l'identifiant Ensembl à celui de Braunch pour chacun de nos transcrits

	file_in.close()

	return(canonical_transcript,canonical_list)

### On définit une fonction qui va extraire les annotations et enregistrer les coordonnées au format BED
def extract_coord(file_name,canonical_transcript):
	intron_content = {} # dictionnaire qui va contenir les objets issues de IntronInfoAnnotate

	file_in = open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête "Intron	Coordinates	type	txRegion	PTCstatus"

	regex = re.compile('^[^:]+:(.+):[0-9]{,3}') # On définit une expression régulière qui va capturer le nom du gène et l'id refseq pour définir un transcrit Braunch
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		intron_name = content[0] #On récupère le nom de l'intron
		transcript_name = regex.findall(intron_name)[0]
		#On vérifie que l'intron est présent dans un de nos transcrits, sinon on passe
		if transcript_name in canonical_transcript:
			#On enregistre les annotations dans un objet
			intron_name = IntronInfoAnnotate(content[0],content[1],content[2],content[3],content[4],content[5],content[6],content[7],content[8],content[9],content[10].replace('\n', ''))
			#On enregistre l'objet dans le dictionnaire avec l'identifiant de transcrit Braunch comme clé
			if transcript_name in intron_content:
				intron_content[transcript_name].append(intron_name)
			else:
				intron_content[transcript_name]=[]
				intron_content[transcript_name].append(intron_name)

	file_in.close() 

	return(intron_content)

### Fonction qui va parser le fichier GTf à partir de la base de données générée et récupérer les annotations dans un objet
def parsing_GTF(dbname,canonical_list):
	inversed_canonical_list = {v: k for k, v in canonical_list.items()} # On inverse le dictionnaire pour parser les id Ensembl
	#Connexion à la base de données du GTF
	db = gffutils.FeatureDB(dbname, keep_order=True)

	exon_content = {} # On définit un dico en dehors de la boucle qui contiendra tous les transcrits en clé et tous les exons du transcrit en valeur
	count = 0 # On définit un compteur pour chaque fois où un id ENS provenant de notre jeu de donnée match dans le GTF

	# Pour chaque identifiant de gène Ensembl associé à un de nos introns
	for id_ENS in inversed_canonical_list:
		try: # On tente de retrouver l'identifiant dans la base de donnée du fichier GTF
			g = db[id_ENS]
			count += 1 # On ajoute 1 au compteur 
		except: # Si l'identifiant ne match pas
			message = "Impossible de trouver "+id_ENS+" dans le fichier GTF" # On print un message d'alerte
			print(message)
			continue # Et on passe à l'itération suivante
		#Pour tous les éléments issues du transcrit, on récupère les exons et les codons stop uniquement
		for i in db.children(g, featuretype=['CDS','stop_codon']):
			# On récupère les annotations de l'élément (exon ou codon stop) qu'on enregistre dans un objet
			id_exon = i.attributes['exon_number'][0]+":"+i.attributes['transcript_id'][0]+":"+i.featuretype
			id_exon = ExonInfo(cds_chr=i.seqid,cds_start=i.start,cds_stop=i.stop,cds_strand=i.strand,cds_transcript=i.attributes['transcript_id'][0],cds_exon_number=i.attributes['exon_number'][0],gene_id=id_ENS,feature_type=i.featuretype,source=i.source)
			# Puis on enregistre l'élément dans le dictionnaire via son identifiant de transcrit
			# Ex : Transcrit_id : exon1,exon2,exon3,stop_codon
			if id_ENS in exon_content:
				exon_content[id_ENS].append(id_exon)
			else:
				exon_content[id_ENS]=[]
				exon_content[id_ENS].append(id_exon)


	print(count,'/',len(inversed_canonical_list),' matchs',sep='')
	# On retourne le dictionnaire contenant chaque transcrit et ses éléments (exon ou stop)
	return(exon_content)
def BED_writing(exon_content,intron_content):

	file_out = tempfile.NamedTemporaryFile(mode="w+",delete=False,suffix='.bedPos')
###Parsing des éléments en argument ####
create_db = False # Par defaut, gffutils va créer une base de donnée issue du fichier GTF pour le parser

opts, args = getopt.getopt(sys.argv[1:],'',['ref_trans=','annotation=','intron_ens=','genref=','gtf=','gtfdb=','output_fa='])
for elmts in opts:
	if elmts[0] == '--ref_trans':
		ref_trans_file = elmts[1] # Nom du fichier de sortie

	elif elmts[0] == '--annotation':
		annotation_file = elmts[1] # Nom du premier fichier d'annotation

	elif elmts[0] == '--intron_ens':
		intron_ens_file = elmts[1] # Nom du fichier liant l'identifiant de chacun de nos introns à un id de gene Ensembl

	elif elmts[0] == '--genref':
		gen_ref_file = elmts[1] # Nom du fichier contenant le genome de reference (format binaire)

	elif elmts[0] == '--gtf':
		gtf_file = elmts[1] # Nom du fichier au format GTF contenant toutes les annotations de la release (ici release 65)

	elif elmts[0] == '--gtfdb': # On demande si la base de donnée issue du fichier gtf à parser existe
		create_db = True
		gtf_db_file = elmts[1]
	elif elmts[0] == '--output_fa':
		output_fasta = elmts[1]



canonical_transcript,canonical_list = get_canonical_transcript(ref_trans_file)

intron_content = extract_coord(annotation_file,canonical_transcript)

exon_content = parsing_GTF(gtf_db_file,canonical_list)



