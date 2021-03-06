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

class ExonInfo(object):
	"Classe contenant les données de chaque CDS (exon sans 5'-UTR/3' UTR)"
	def __init__(self,cds_chr,cds_start,cds_stop,cds_strand,cds_transcript,cds_exon_number,gene_id,feature_type,transcript_name):

		self.id = cds_transcript+":"+cds_exon_number+":"+str(cds_start)+"-"+str(cds_stop)+":"+feature_type
		self.chr = cds_chr
		self.start = str(int(cds_start)-1) # On enlève 1 pour passer en référentiel 0
		self.stop = cds_stop
		self.strand = cds_strand
		self.transcript = cds_transcript
		self.exon_number = cds_exon_number
		self.gene_id = gene_id
		self.feature_type = feature_type
		self.transcript_name = transcript_name

	def formating_coord(self):
		return "chr"+self.chr+":"+str(self.start)+"-"+str(self.stop)

	def BED_coord(self):
		return "chr"+self.chr+"\t"+str(self.start)+"\t"+str(self.stop)+"\t"+self.id


class ExonInfoComplete(ExonInfo):
	"Classe contenant les mêmes informations que ExonInfo avec la séquence en plus"
	def __init__(self,cds_chr,cds_start,cds_stop,cds_strand,cds_transcript,cds_exon_number,gene_id,feature_type,transcript_name,seq):
		ExonInfo.__init__(self,cds_chr,cds_start,cds_stop,cds_strand,cds_transcript,cds_exon_number,gene_id,feature_type,transcript_name)
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
def extract_coord(file_name,IDlist):
	intron_content = {} # dictionnaire qui va contenir les objets issues de IntronInfoAnnotate
	coord_to_intron = {} # dictionnaire qui contient en clé les coordonnées de nos introns et en valeur l'identifiant de chaque intron

	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête "Intron	Coordinates	type	txRegion	PTCstatus"

	file_out = tempfile.NamedTemporaryFile(mode="w+",delete=False,suffix='.bedPos')
	file_out_name = file_out.name
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		intron_name = content[0] #On récupère le nom de l'intron
		#On enregistre les annotations dans un objet
		if intron_name in IDlist: # On vérifie que l'intron est annoté d'un gène, sinon cela ne nous intérèsse pas
			intron_name = IntronInfoAnnotate(content[0],content[1],content[2],content[3],content[4],content[5],content[6],content[7],content[8],content[9],content[10].replace('\n', ''))
			#On enregistre l'objet dans le dictionnaire avec son identifiant comme clé
			intron_content[intron_name.id]=intron_name
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
	return(intron_content,file_out_name,coord_to_intron)


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
				id_CDS = ExonInfo(cds_chr=i.seqid,cds_start=i.start,cds_stop=i.stop,cds_strand=i.strand,cds_transcript=i.attributes['transcript_id'][0],cds_exon_number=i.attributes['exon_number'][0],gene_id=id_ENS,feature_type=i.featuretype,transcript_name=i.attributes['transcript_name'][0])
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
def extract_fasta_info(filename,intron_content,liste_transcripts,IDlist): 
	
	seq_info = {} # On définit un dictionnaire qui contient la séquence en valeur et ses coordonnées en clé
	dataset_intron = {} # On définit un dictionnaire qui contiendra les données clé : id_intron, valeur = objet_intron
	dataset_exon = {} # On définit un dictionnaire qui contiendra les données clé : id_exon, valeur = objet_exon
	exon_by_transcripts = {} # Dictionnaire qui contiendra les CDS(transcrits) avec leurs exons (sans 5' et 3' UTR)

	
	handle = open(filename, "rU") #On ouvre le fichier en mode lecture
	#Parsing du fichier fasta pour récupérer toutes les séquences
	for record in SeqIO.parse(handle, "fasta") : #On parse le fichier avec seqIO
		seq_id = record.id #On enregistre le nom de la sequence, ici les coordonnées
		seq = record.seq #On enregistre la sequence fasta dans une variable
		seq= seq.upper() #On met en majuscule au cas ou la séquence ne l'est pas
		seq_info[seq_id]=seq # On enregistre la séquence dans le dictionnaire avec ses coordonnées comme clé
		
	handle.close()
	#Parsing du fichier pour récupérer les séquences des introns
	for intron_id in intron_content.values():
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
			CDS_id = ExonInfoComplete(CDS_id.chr,CDS_id.start,CDS_id.stop,CDS_id.strand,CDS_id.transcript,CDS_id.exon_number,CDS_id.gene_id,CDS_id.feature_type,CDS_id.transcript_name,str(seq))
			if CDS_id.transcript in exon_by_transcripts:
				exon_by_transcripts[CDS_id.transcript].append([CDS_id.id,CDS_id.formating_coord()])
			else:
				exon_by_transcripts[CDS_id.transcript] = [] # On définit la valeur comme étant une liste si elle n'était pas définie avant
				exon_by_transcripts[CDS_id.transcript].append([CDS_id.id,CDS_id.formating_coord()])
			dataset_exon[CDS_id.id] = CDS_id

	return(dataset_intron,exon_by_transcripts,dataset_exon)

#Fonction va parcourir notre liste d'identifiants de transcrits pour chaque gène et associer les introns de notre jeu de données
def junction_introns_to_transcripts(exon_by_transcripts,dataset_exon,coord_to_intron):
	intron_by_transcripts = {} # Dictionnaire qui contiendra la liste des identifiant des introns dans chaque transcrit
	number_CDS_no_stop = 0 # On définit un compteur pour chaque fois ou un transcrit n'a pas de codon stop -> assemblage incomplet 
	for Transcript_id,exons in exon_by_transcripts.items():
		stop_codon_control = 0 # Control pour savoir si le transcrit possède un codon stop

		exons.sort(key=lambda x: x[0].split(':')[2]) # Fonction de tri qui met les exons du CDS dans l'ordre en fonction des coordonnées
		exons_without_stop = []
		for elmt in exons:
			#Vérification si transcrit possède un codon stop
			exon_left = dataset_exon[elmt[0]]
			if exon_left.feature_type == 'stop_codon':
				stop_codon_control += 1 #Si oui on met +1 au control
			else:
				exons_without_stop.append(elmt)

		for elmt in exons_without_stop:
			if exons.index(elmt)+1 != len(exons): # Si notre élément n'est pas le dernier, on peut continuer
				position_exon_left = exons.index(elmt)
				position_exon_right = int(position_exon_left)+1
				exon_left = dataset_exon[elmt[0]]
				exon_right = dataset_exon[exons[position_exon_right][0]]
				
				intron_coord = "chr"+exon_left.chr+":"+str(int(exon_left.stop))+"-"+str(int(exon_right.start)+2)
				#Vu qu'on peut avoir plusieurs introns pour les mêmes coordonnées, on doit parcourir cette liste et choisir le bon intron
				if intron_coord in coord_to_intron:
					intron_list = coord_to_intron[intron_coord]
					for intron in intron_list:
						name = intron.split(":")[1]
						if re.search(name , exon_left.transcript_name):
							if Transcript_id in intron_by_transcripts:
								intron_by_transcripts[Transcript_id].append(intron_coord)
							else:
								intron_by_transcripts[Transcript_id]=[]
								intron_by_transcripts[Transcript_id].append(intron_coord)

		if stop_codon_control == 0: # Si le control vaut 0, c'est qu'on a pas de codon stop dans la séquence
			number_CDS_no_stop +=1 # On ajoute donc +1 au compteur de transcrits sans codon stop
	
	return(intron_by_transcripts,number_CDS_no_stop)

##Fonction qui va trier les introns de sorte à les lister par transcrits, puis qui va les comparer à la liste des transcrits ensembm, et vérifier
##Si la liste d'intron dans les deux listes sont égales, si c'est le cas, on considerera le transcrit ensembl comme le canonique pour notre jeu de données
def retrieve_ensembl_transcript(dataset_intron,intron_by_transcripts):
	group_by_intron={} # Dictionnaire qui contiendra la liste des introns portant le même identifiant
	transcript_no_match = {} # Dictionnaire qui contiendra la liste des transcrits n'étant pas considéré comme référence
	#On parcours les identifiants pour déterminer le nombre d'introns présent dans notre jeu de données dans chaque gene
	regex = re.compile('^[^:]+:(.+):[0-9]{,3}')
	for key,value in dataset_intron.items():
		result = regex.findall(key)
		intron_and_coord = [key,value.formating_coord_for_exon()]
		if result[0] in group_by_intron:
			group_by_intron[result[0]].append(intron_and_coord)
		else:
			group_by_intron[result[0]]=[]
			group_by_intron[result[0]].append(intron_and_coord)
		group_by_intron[result[0]].sort(key=lambda x: x[1].split(':')[1])

	Canonical_Transcript = {} # Contient les identifiants des transcrits qui contiennent tous les introns d'un gène de notre jeu de données
	for key,value in intron_by_transcripts.items():

		introns_in_transcript = value # Liste des introns dans chaque transcrit
		result = regex.findall(introns_in_transcript[0])

		if result[0] in group_by_intron:
			introns_by_group = []
			[introns_by_group.append(elmt[0]) for elmt in group_by_intron[result[0]] ] # On prend uniquement l'identifiant pour comparer
		else:
			print('Le transcrit',key,'ne contient pas d\'introns de notre jeu de données')
			continue
		if introns_by_group == introns_in_transcript:
			Canonical_Transcript[key] = result[0]
		else:
			transcript_no_match[key]=value
	return(Canonical_Transcript,transcript_no_match)

def get_CDS_by_canonical_transcript(dataset_intron,dataset_exon,exon_by_transcripts,Canonical_Transcript):
	for transcript,introns in Canonical_Transcript.items():
		exon_list = exon_by_transcripts[transcript]
		intron_list = introns


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
intron_content,file_out_name,coord_to_intron = extract_coord(annotation_file,IDlist)
print('Ok, annotations sur les Introns enregistrés')
print('Parsing de',gtf_file,'via',gtf_db_file,'...')
###Création de la base de données du GTF si n'existe pas####
if create_db == False:
	print('Création de la database issu du fichier GTF')
	db = gffutils.create_db(result[0], dbfn=gtf_db_file, force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
	print('Database OK')
liste_transcripts = parsing_GTF(file_out_name,gtf_db_file)
print('Ok, annotations sur les transcrits enregistrés')
print('Lancement de TwobitToFa -> Récupération des séquences fasta via les coordonnées ...')
if os.path.isfile(output_fasta):
	print("Fichier fasta déjà existant, passe à l'étape suivante.")
else:
	command = "twobitToFa -bedPos -bed="+file_out_name+" "+gen_ref_file+" "+output_fasta # Commande qui lance TwobitToFa pour récupérer la séquence de chaque coordonnée intronique
	os.system(command) # Execute la commande précédente
	print('Ok, séquences enregistrées')
os.unlink(file_out_name) # Suprimme le fichier BED temporaire
print('Extraction des séquences ...')
dataset_intron,exon_by_transcripts,dataset_exon = extract_fasta_info(filename=output_fasta,intron_content=intron_content,liste_transcripts=liste_transcripts,IDlist=IDlist)
print('Ok, séquences enregistrés dans les annotations')	
intron_by_transcripts,number_CDS_no_stop = junction_introns_to_transcripts(exon_by_transcripts,dataset_exon,coord_to_intron)
print(number_CDS_no_stop,'/',len(exon_by_transcripts),' transcrits n\'ont pas de codon stop (',round((number_CDS_no_stop/len(exon_by_transcripts))*100,2),'%)',sep='')
print(len(intron_by_transcripts),'/',len(exon_by_transcripts),' transcrits possèdent des introns de nos jeux de données (',round((len(intron_by_transcripts)/len(exon_by_transcripts))*100,2),'%)',sep='')
Canonical_Transcript,transcript_no_match = retrieve_ensembl_transcript(dataset_intron,intron_by_transcripts)
print(len(Canonical_Transcript),'/',len(exon_by_transcripts),' transcrits sont identifiés comme étant les canoniques (',round((len(Canonical_Transcript)/len(exon_by_transcripts))*100,2),'%)',sep='')


