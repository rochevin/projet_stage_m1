####################################################################
#Auteur : ROCHER Vincent										   #
#																   #
#But : Programme principal										   #
#	   -Associe les données de Braunschweig + données annotées	   #
#		à l'id de chaque gene Ensembl							   #
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

####Définition des classes######

class CDSInfo(object):
	"Classe contenant les données de chaque CDS (exon sans 5'-UTR/3' UTR)"
	def __init__(self,cds_chr,cds_start,cds_stop,cds_strand,cds_transcript,cds_exon_number,gene_id,feature_type):

		self.id = cds_transcript+":"+cds_exon_number+":"+str(cds_start)+"-"+str(cds_stop)+":"+feature_type
		self.chr = cds_chr
		self.start = str(int(cds_start)-1) # On enlève 1 pour passer en référentiel 0
		self.stop = cds_stop
		self.strand = cds_strand
		self.transcript = cds_transcript
		self.exon_number = cds_exon_number
		self.gene_id = gene_id
		self.feature_type = feature_type

	def formating_coord(self):
		return "chr"+self.chr+":"+str(self.start)+"-"+str(self.stop)

	def BED_coord(self):
		return "chr"+self.chr+"\t"+str(self.start)+"\t"+str(self.stop)+"\t"+self.id


class CDSInfoComplete(CDSInfo):
	"Classe contenant les mêmes informations que CDSInfo avec la séquence en plus"
	def __init__(self,cds_chr,cds_start,cds_stop,cds_strand,cds_transcript,cds_exon_number,gene_id,feature_type,seq):
		CDSInfo.__init__(self,cds_chr,cds_start,cds_stop,cds_strand,cds_transcript,cds_exon_number,gene_id,feature_type)
		self.seq = seq


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

	def formating_coord(self):
		return self.chr+":"+str(self.debut)+"-"+str(self.fin)

	def formating_coord_for_exon(self):
		return self.chr+":"+str(self.debut)+"-"+str(self.fin+1)

	def flank_exon_left(self):
		regex = re.compile('^(chr[\w]+)[-\+]([0-9]+):([0-9]+)_[0-9]+:[0-9]+')
		result = regex.findall(self.coords)
		coord_exon_left = result[0][0]+":"+result[0][1]+"-"+result[0][2]
		return coord_exon_left
	def flank_exon_right(self):
		regex = re.compile('^(chr[\w]+)[-\+][0-9]+:[0-9]+_([0-9]+):([0-9]+)')
		result = regex.findall(self.coords)
		coord_exon_right = result[0][0]+":"+result[0][1]+"-"+result[0][2]
		return coord_exon_right
	def BED_coord(self):
		return self.chr+"\t"+str(self.debut)+"\t"+str(self.fin)+"\t"+self.id


class IntronInfoComplete(IntronInfoAnnotate):
		"Classe contenant les informations des introns de la classe parente + seq + ID Ensembl"

		def __init__(self,Intron_id,Intron_coord,Intron_type,Intron_txRegion,Intron_PTCstatus,seq_annotation_begin,seq_annotation_end,intron_len,intron_position,total_intron,GCrate,seq,ensembl):
			IntronInfoAnnotate.__init__(self,Intron_id,Intron_coord,Intron_type,Intron_txRegion,Intron_PTCstatus,seq_annotation_begin,seq_annotation_end,intron_len,intron_position,total_intron,GCrate)
			self.seq = seq
			self.ensembl = ensembl

###Fonction qui va extraire les identifiant de genes Ensembl associés à notre jeu d'intron et les associer à leur identifiants ainsi qu'être listés pour parser le GTF
def extract_ensembl_id(file_name):
	IDlist = {} # Dictionnaire contenant en clé l'id de l'intron et en valeur l'id Ensembl
	list_Ensembl_ids = {} # Dictionnaire contenant tous les ids Ensembl en clé pour éviter la redondance -> parsing GTF
	file_in=open(file_name,"r")
	header = file_in.readline() # On enregistre l'en-tête
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t")
		IDlist[content[0]] = content[1].replace('\n', '')
		list_Ensembl_ids[content[1].replace('\n', '')] =""
	file_in.close()
	return(IDlist,list_Ensembl_ids)

###On définit une fonction qui va extraire les annotations et enregistrer les coordonnées au format BED
def extract_coord(file_name):
	retrieve_name = {} # dictionnaire qui va contenir les objets issues de IntronInfoAnnotate
	coord_to_intron = {} # dictionnaire qui contient en clé les coordonnées de nos introns et en valeur l'identifiant de chaque intron

	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête "Intron	Coordinates	type	txRegion	PTCstatus"

	file_out = tempfile.NamedTemporaryFile(mode="w+",delete=False,suffix='.bedPos')
	file_out_name = file_out.name
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		intron_name = content[0] #On récupère le nom de l'intron
		#On enregistre les annotations dans un objet
		intron_name = IntronInfoAnnotate(content[0],content[1],content[2],content[3],content[4],content[5],content[6],content[7],content[8],content[9],content[10].replace('\n', ''))
		#On enregistre l'objet dans le dictionnaire avec son identifiant comme clé
		retrieve_name[intron_name.id]=intron_name
		#On enregistre les coordonnées de l'intron dans le dictionnaire qui va identifier un intron à ses coordonnées via son identifiant
		if intron_name.formating_coord_for_exon() in coord_to_intron:
			coord_to_intron[intron_name.formating_coord_for_exon()].append(intron_name.id)
		else:
			coord_to_intron[intron_name.formating_coord_for_exon()]=[]
			coord_to_intron[intron_name.formating_coord_for_exon()].append(intron_name.id)
		# Écriture du fichier BED
		file_out.write(intron_name.BED_coord())
		file_out.write("\n")
	file_in.close() 
	file_out.close()
	return(retrieve_name,file_out_name,coord_to_intron)

###Fonction qui va parser le fichier GTf à partir de la base de données générée et récupérer les annotations dans un objet
def parsing_GTF(temporary_file_name,db_name):
	#Connexion à la base de données du GTF
	db = gffutils.FeatureDB(db_name, keep_order=True)

	file_out=open(temporary_file_name,"a") # On ouvre le fichier temporaire en mode ajout
	liste_transcripts = {} # On définit un dico en dehors de la boucle qui contiendra tous les transcrits en clé et tous les CDS du transcrit en valeur
	count = 0 # On définit un compteur pour chaque fois où un id ENS provenant de notre jeu de donnée match dans le GTF

	gene_with_stop = 0 # Compteur du nombre de gene avec un codon stop annoté
	for id_ENS in list_Ensembl_ids:
		try:
			g = db[id_ENS]
			count += 1
		except:
			if id_ENS != "NA":
				message = "Impossible de trouver "+id_ENS+" dans le fichier GTF"
				print(message)
			continue
		compteur_stop = 0
		
		for i in db.children(g, featuretype=['exon','stop_codon']):
			if i.source == 'protein_coding':
				if i.featuretype =="stop_codon":
					compteur_stop =1
				id_CDS = i.attributes['exon_number'][0]+":"+i.attributes['transcript_id'][0]+":"+i.featuretype
				id_CDS = CDSInfo(cds_chr=i.seqid,cds_start=i.start,cds_stop=i.stop,cds_strand=i.strand,cds_transcript=i.attributes['transcript_id'][0],cds_exon_number=i.attributes['exon_number'][0],gene_id=id_ENS,feature_type=i.featuretype)
				if i.attributes['transcript_id'][0] in liste_transcripts:
					liste_transcripts[i.attributes['transcript_id'][0]].append(id_CDS)
				else:
					liste_transcripts[i.attributes['transcript_id'][0]]=[]
					liste_transcripts[i.attributes['transcript_id'][0]].append(id_CDS)
				file_out.write(id_CDS.BED_coord())
				file_out.write("\n")

			
		if compteur_stop ==1:
				gene_with_stop+=1

	print(gene_with_stop,'genes possèdent un codon stop sur',len(list_Ensembl_ids)-1)
	file_out.close()
	print(count,'/',len(list_Ensembl_ids)-1,' matchs',sep='')
	return(liste_transcripts)


###On définit une fonction pour récupérer la séquence fasta de chacun de nos introns et de nos CDS
def extract_fasta_info(filename,retrieve_name,liste_transcripts,IDlist): 
	
	seq_info = {} # On définit un dictionnaire qui contient la séquence en valeur et ses coordonnées en clé
	dataset_intron = {} # On définit un dictionnaire qui contiendra les données clé : id_intron, valeur = objet_intron
	dataset_CDS = {} # Dictionnaire qui contiendra les CDS(transcrits) avec leurs exons (sans 5' et 3' UTR)

	CDS_content = {} # On définit un dictionnaire qui contiendra les données clé : id_CDS, valeur = objet_CDS


	handle = open(filename, "rU") #On ouvre le fichier en mode lecture
	#Parsing du fichier fasta pour récupérer toutes les séquences
	for record in SeqIO.parse(handle, "fasta") : #On parse le fichier avec seqIO
		seq_id = record.id #On enregistre le nom de la sequence, ici les coordonnées
		seq = record.seq #On enregistre la sequence fasta dans une variable
		seq= seq.upper() #On met en majuscule au cas ou la séquence ne l'est pas
		seq_info[seq_id]=seq # On enregistre la séquence dans le dictionnaire avec ses coordonnées comme clé
		
	handle.close()
	#Parsing du fichier pour récupérer les séquences des introns
	for intron_id in retrieve_name.values():
		coords = intron_id.formating_coord() # On récupère les coordonées de façon formatés pour retrouver la sequence
		seq = seq_info[coords] # On va dans le dictionnaire qui contient les séquences et on récupère la séquence correspondant aux coordonnées de l'intron
		Ensembl_id = IDlist[intron_id.id] # On récupère l'id ensembl du gene de l'intron à partir de l'id de l'intron


		if intron_id.signe == "-":
			seq = seq.reverse_complement()
		intron_id = IntronInfoComplete(intron_id.id,intron_id.coords,intron_id.type,intron_id.txregion,intron_id.PTC,intron_id.begin_annotation,intron_id.end_annotation,intron_id.len,intron_id.pos,intron_id.total,intron_id.GCrate,str(seq),Ensembl_id) 

		dataset_intron[intron_id.id]=intron_id


	#Parsing du fichier pour récupérer les séquences des CDS
	for Trans_id in liste_transcripts.values():
		for CDS_id in Trans_id:
			coords = CDS_id.formating_coord() # On récupère les coordonées de façon formatés pour retrouver la sequence
			seq = seq_info[coords] # On va dans le dictionnaire qui contient les séquences et on récupère la séquence correspondant aux coordonnées de l'intron


			if CDS_id.strand == "-":
				seq = seq.reverse_complement()
			CDS_id = CDSInfoComplete(CDS_id.chr,CDS_id.start,CDS_id.stop,CDS_id.strand,CDS_id.transcript,CDS_id.exon_number,CDS_id.gene_id,CDS_id.feature_type,str(seq))
			if CDS_id.transcript in dataset_CDS:
				dataset_CDS[CDS_id.transcript].append(CDS_id.id)
			else:
				dataset_CDS[CDS_id.transcript] = [] # On définit la valeur comme étant une liste si elle n'était pas définie avant
				dataset_CDS[CDS_id.transcript].append(CDS_id.id)
			CDS_content[CDS_id.id] = CDS_id

	return(dataset_intron,dataset_CDS,CDS_content)

###Parsing des éléments en argument ####
create_db = False # Par defaut, gffutils va créer une base de donnée issue du fichier GTF pour le parser

opts, args = getopt.getopt(sys.argv[1:],'o:',['annotation=','intron_ens=','genref=','gtf=','gtfdb=','output_fa='])
for elmts in opts:
	if elmts[0] == '-o':
		output_file = elmts[1] # Nom du fichier de sortie

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


print('Parsing de',intron_ens_file,'...')
IDlist,list_Ensembl_ids = extract_ensembl_id(intron_ens_file)
print('OK, IDs Ens associés aux IDs Introns')
print('Parsing de',annotation_file,'...')
retrieve_name,file_out_name,coord_to_intron = extract_coord(annotation_file)
print('Ok, annotations sur les Introns enregistrés')
print('Parsing de',gtf_file,'via',gtf_db_file,'...')
###Création de la base de données du GTF si n'existe pas####
if create_db == False:
	print('Création de la database issu du fichier GTF')
	db = gffutils.create_db(result[0], dbfn=gtf_db_file, force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
	print('Database OK')
liste_transcripts = parsing_GTF(file_out_name,gtf_db_file)
print('Ok, annotations sur les CDS enregistrés')
print('Lancement de TwobitToFa -> Récupération des séquences fasta via les coordonnées ...')
if os.path.isfile(output_fasta):
	print("Fichier fasta déjà existant, passe à l'étape suivante.")
else:
	command = "twobitToFa -bedPos -bed="+file_out_name+" "+gen_ref_file+" "+output_fasta # Commande qui lance TwobitToFa pour récupérer la séquence de chaque coordonnée intronique
	os.system(command) # Execute la commande précédente
	print('Ok, séquences enregistrées')
os.unlink(file_out_name) # Suprimme le fichier BED temporaire
print('Extraction des séquences ...')
dataset_intron,dataset_CDS,CDS_content = extract_fasta_info(filename=output_fasta,retrieve_name=retrieve_name,liste_transcripts=liste_transcripts,IDlist=IDlist)
print('Ok, séquences enregistrés dans les annotations')	


#Parcours de tous les CDS #####
intron_by_transcripts = {} # Dictionnaire qui contiendra la liste des identifiant des introns dans chaque transcrit
number_CDS_no_stop = 0 # On définit un compteur pour chaque fois ou un transcrit n'a pas de codon stop -> assemblage incomplet 
for Transcript_id,CDS in dataset_CDS.items():
	stop_codon_control = 0 # Control pour savoir si le transcrit possède un codon stop

	CDS.sort(key=lambda x: x.split(':')[2]) # Fonction de tri qui met les exons du CDS dans l'ordre en fonction des coordonnées
	

	for elmt in CDS:
		if CDS.index(elmt)+1 != len(CDS): # Si notre élément n'est pas le dernier, on peut continuer
			position_exon_left = CDS.index(elmt)
			position_exon_right = int(position_exon_left)+1
			exon_left = CDS_content[elmt]
			exon_right = CDS_content[CDS[position_exon_right]]

			if (exon_left.feature_type == 'stop_codon' or exon_right.feature_type == 'stop_codon'):
				continue
			print(exon_left.id,"and",exon_right.id)
			intron_coord = "chr"+exon_left.chr+":"+str(int(exon_left.stop))+"-"+str(int(exon_right.start)+2)
			if intron_coord in coord_to_intron:
				intron_list = coord_to_intron[intron_coord]
				for intron in intron_list:
					if Transcript_id in intron_by_transcripts:
						intron_by_transcripts[Transcript_id].append(intron)
					else:
						intron_by_transcripts[Transcript_id]=[]
						intron_by_transcripts[Transcript_id].append(intron)
		else: # Sinon on arrete le parcours
			break
	if stop_codon_control == 0: # Si le control vaut 0, c'est qu'on a pas de codon stop dans la séquence
		number_CDS_no_stop +=1 # On ajoute donc +1 au compteur de transcrits sans codon stop



print(number_CDS_no_stop,'/',len(dataset_CDS),' transcrits n\'ont pas de codon stop (',round((number_CDS_no_stop/len(dataset_CDS))*100,2),'%)',sep='')
print(len(intron_by_transcripts),'/',len(dataset_CDS),' transcrits possèdent des introns de nos jeux de données (',round((len(intron_by_transcripts)/len(dataset_CDS))*100,2),'%)',sep='')


