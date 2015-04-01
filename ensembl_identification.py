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
	def __init__(self,cds_chr,cds_start,cds_stop,cds_strand,cds_transcript,cds_exon_number,gene_id,feature_type):

		self.id = cds_transcript+":"+cds_exon_number+":"+str(cds_start)+"-"+str(cds_stop)+":"+feature_type
		self.chr = cds_chr
		self.start = cds_start
		self.stop = cds_stop
		self.strand = cds_strand
		self.transcript = cds_transcript
		self.exon_number = cds_exon_number
		self.gene_id = gene_id
		self.feature_type = feature_type

	def formating_coord(self):
		return "chr"+self.chr+":"+str(self.start)+"-"+str(self.stop)


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

	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline() #On enregistre la première ligne étant l'en tête "Intron	Coordinates	type	txRegion	PTCstatus"

	regex = re.compile('^[^:]+:(.+):[0-9]{,3}') # On définit une expression régulière qui va capturer le nom du gène et l'id refseq pour définir un transcrit Braunch
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		intron_name = content[0] #On récupère le nom de l'intron
		#On vérifie que l'intron est annoté d'un gène, sinon cela ne nous sert à rien de l'enregistrer
		if IDlist[intron_name] != "NA":
			#On enregistre les annotations dans un objet
			intron_name = IntronInfoAnnotate(content[0],content[1],content[2],content[3],content[4],content[5],content[6],content[7],content[8],content[9],content[10].replace('\n', ''))
			#On enregistre l'objet dans le dictionnaire avec l'identifiant de transcrit Braunch comme clé
			result = regex.findall(intron_name.id)
			if result[0] in intron_content:
				intron_content[result[0]].append(intron_name)
			else:
				intron_content[result[0]]=[]
				intron_content[result[0]].append(intron_name)

	file_in.close() 

	return(intron_content)


###Fonction qui va parser le fichier GTf à partir de la base de données générée et récupérer les annotations dans un objet
def parsing_GTF(dbname,list_Ensembl_ids):
	#Connexion à la base de données du GTF
	db = gffutils.FeatureDB(dbname, keep_order=True)

	exon_content = {} # On définit un dico en dehors de la boucle qui contiendra tous les transcrits en clé et tous les exons du transcrit en valeur
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
				id_CDS = ExonInfo(cds_chr=i.seqid,cds_start=i.start,cds_stop=i.stop,cds_strand=i.strand,cds_transcript=i.attributes['transcript_id'][0],cds_exon_number=i.attributes['exon_number'][0],gene_id=id_ENS,feature_type=i.featuretype)
				if i.attributes['transcript_id'][0] in exon_content:
					exon_content[i.attributes['transcript_id'][0]].append(id_CDS)
				else:
					exon_content[i.attributes['transcript_id'][0]]=[]
					exon_content[i.attributes['transcript_id'][0]].append(id_CDS)

			
		if compteur_stop ==1:
				gene_with_stop+=1

	print(gene_with_stop,'genes possèdent un codon stop sur',len(list_Ensembl_ids)-1)

	print(count,'/',len(list_Ensembl_ids)-1,' matchs',sep='')
	return(exon_content)

#Fonction qui va suprimmer les stop pour ne plus les prendres en compte
def exon_construction(exon_content):
	number_CDS_no_stop = 0
	exon_coord_in_ensembl_transcript = {} # Dictionnaire de transcrit qui contient chaque coordonnées de ses introns
	for transcript,exons in exon_content.items():
		stop_codon_control = 0 # Control pour savoir si le transcrit possède un codon stop

		
		exons_without_stop = []
		for elmt in exons:
			#Vérification si transcrit possède un codon stop
			if elmt.feature_type == 'stop_codon':
				stop_codon_control += 1 #Si oui on met +1 au control
			else:
				exons_without_stop.append(elmt)

		#Va associer pour chaque transcrits une liste de coordonnées et trier la liste en fonction des coordonnées
		for elmt in exons_without_stop:
			exon = elmt.formating_coord()
			if transcript in exon_coord_in_ensembl_transcript:
				exon_coord_in_ensembl_transcript[transcript].append(exon)
			else:
				exon_coord_in_ensembl_transcript[transcript] = []
				exon_coord_in_ensembl_transcript[transcript].append(exon)

			exon_coord_in_ensembl_transcript[transcript].sort(key=lambda x: x.split('-')[1]) # Fonction de tri qui met les exons du CDS dans l'ordre en fonction des coordonnées

		if stop_codon_control == 0: # Si le control vaut 0, c'est qu'on a pas de codon stop dans la séquence
			number_CDS_no_stop +=1 # On ajoute donc +1 au compteur de transcrits sans codon stop
	print(number_CDS_no_stop,'/',len(exon_content),' transcrits n\'ont pas de codon stop (',round((number_CDS_no_stop/len(exon_content))*100,2),'%)',sep='')
	return(exon_coord_in_ensembl_transcript)

def transform_exon_to_intron(exon_coord_in_ensembl_transcript,exon_content):
	intron_coord_in_ensembl_transcript = {} # Dictionnaire qui va identifier des introns à un identifiant de transcrit
	retrieve_transcript_with_coord = {} # Dictionnaire de coordonnées qui va enregistrer une coordonnée et à quel(s) transcrit(s) elle appartient
	for transcript,liste_exon in exon_coord_in_ensembl_transcript.items():
		liste_exon.sort(key=lambda x: x.split('-')[1])
		for elmt in liste_exon:
			if liste_exon.index(elmt)+1 != len(liste_exon): # Si notre élément n'est pas le dernier, on peut continuer
				position_exon_left = liste_exon.index(elmt)
				position_exon_right = int(position_exon_left)+1
				#exon de gauche
				regex = re.compile('^(chr[\w]+):[0-9]+-([0-9]+)')
				result = regex.findall(liste_exon[position_exon_left])
				exon_chr = result[0][0]
				stop_left = result[0][1]
				#exon de droite
				regex = re.compile('^(chr[\w]+):([0-9]+)-[0-9]+')
				result = regex.findall(liste_exon[position_exon_right])
				start_right = result[0][1]

				intron_coord = exon_chr+":"+str(stop_left)+"-"+str(start_right)
				if transcript in intron_coord_in_ensembl_transcript:
					intron_coord_in_ensembl_transcript[transcript].append(intron_coord)
				else:
					intron_coord_in_ensembl_transcript[transcript]=[]
					intron_coord_in_ensembl_transcript[transcript].append(intron_coord)

				if intron_coord in retrieve_transcript_with_coord:
					retrieve_transcript_with_coord[intron_coord].append(transcript)
				else:
					retrieve_transcript_with_coord[intron_coord] = []
					retrieve_transcript_with_coord[intron_coord].append(transcript)
	return(intron_coord_in_ensembl_transcript,retrieve_transcript_with_coord)

def retrieve_canonical_transcript(intron_coord_in_ensembl_transcript,intron_content,retrieve_transcript_with_coord):
	canonical_transcripts = {}

	for B_transcript,intron_list in intron_content.items():
		liste = [] # Liste qui va servir à augmenter la rapidité du parcours en cherchant les identifiants ensembl associés à nos coordonnées
		for elmt in intron_list:
			coords = elmt.formating_coord_for_exon()
			if coords in retrieve_transcript_with_coord:
				resultat = retrieve_transcript_with_coord[coords]
				liste.extend(resultat)
		#On suprimme les doublons
		liste = set(list(liste))
		B_intron_list = []
		#On transforme les objets intron en coordonnées
		[B_intron_list.append(elmt.formating_coord_for_exon()) for elmt in intron_list]
		B_intron_list.sort(key=lambda x: x.split('-')[1])
		success = 0
		for elmt in liste:
			E_intron_list = intron_coord_in_ensembl_transcript[elmt]
			if B_intron_list == E_intron_list:
				canonical_transcripts[elmt]=B_transcript
	return(canonical_transcripts)

#Fonction qui détermine les doublons de notre jeu de données
def determine_doublons(canonical_transcripts):
	dictionnary_without_redondance = {} # On associe un transcrit braunch a une liste de transcrits Ensembl qui sont potentiellements plusieurs

	for key,value in canonical_transcripts.items():
		if value in dictionnary_without_redondance:
			dictionnary_without_redondance[value].append(key)
		else:
			dictionnary_without_redondance[value]=[]
			dictionnary_without_redondance[value].append(key)
	doublon = {}
	for elmt,value in dictionnary_without_redondance.items():
		if len(value)>1:
			doublon[elmt]=value
	return(doublon,dictionnary_without_redondance)


def transform_intron_to_exon(intron_content):
	exon_coord_in_braunch_transcript = {} #Dictionnaire qui contient les coordonnées de nos exons pour chaque transcrit Braunch
	for B_transcript, B_intron_list in intron_content.items():
		exon_coord = [] # Liste qui va contenir les coordonnées de tous les exons flanquants nos introns
		for intron in B_intron_list:
			exon_coord.append(intron.flank_exon_left())
			exon_coord.append(intron.flank_exon_right())
		exon_coord = list(set(exon_coord))
		exon_coord.sort(key=lambda x: x.split('-')[1])
		exon_coord_in_braunch_transcript[B_transcript]=exon_coord

	return(exon_coord_in_braunch_transcript)

# Fonction qui va identifier le transcrit le plus proche parmis les doublons
def scoring_doublon(doublon,exon_coord_in_braunch_transcript,exon_coord_in_ensembl_transcript):
	best_score = {} # Va contenir pour chaque transcrit braunch quel transcrit Ensembl est le plus proche
	for B_transcript,liste_of_E_transcripts in doublon.items(): # On obtient un identifiant Braunch et une liste d'identifiants Ensembl
		B_intron_list = exon_coord_in_braunch_transcript[B_transcript] # On récupère la liste des introns de notre transcrit Braunch

		B_intron_list.sort(key=lambda x: x.split('-')[1]) # On trie les coordonnées pour bien prendre le premier et dernier exon
		#On détermine les coordonnées du premier exon de l'intron comme base 0
		regex_first = re.compile('^chr[\w]+:([0-9]+)-[0-9]+')
		result = regex_first.findall(B_intron_list[0])
		base_first_exon = int(result[0])
		#On détermine les coordonnées du dernier exon de l'intron comme base 0
		regex_last = re.compile('^chr[\w]+:[0-9]+-([0-9]+)')
		result = regex_last.findall(B_intron_list[-1])
		base_last_exon = int(result[0])

		best_score_for_transcript = {} # Va contenir le score en clé et l'identifiant du transcrit correspondant au score
		for E_transcript in liste_of_E_transcripts:
			E_intron_list = exon_coord_in_ensembl_transcript[E_transcript]
			#On détermine les coordonnées du premier exon de du transcrit ensembl
			result_first = regex_first.findall(E_intron_list[0])
			coord_first_exon = int(result_first[0])
			#On détermine les coordonnées du dernier exon de du transcrit ensembl
			result_last = regex_last.findall(E_intron_list[-1])
			coord_last_exon = int(result_last[0])
			#Calcul du score :
			score = abs(base_first_exon-coord_first_exon)+abs(base_last_exon-coord_last_exon)
			best_score_for_transcript[score]=E_transcript
		best_E_transcript = best_score_for_transcript[min(best_score_for_transcript)]
		best_score[B_transcript] = best_E_transcript
	return(best_score)

def join_best_score_with_canonical(best_score,dictionnary_without_redondance):
	new_canonical_transcript = {} # Dictionnaire final qui va contenir les transcrits Braunch et a quel Transcrit Ensembl il appartient sans doublon
	new_canonical_transcript = dictionnary_without_redondance.copy() # Ajoute les données du dictionnaire avec les doublons
	new_canonical_transcript.update(best_score) # Puis remplace les doublons par l'identifiant calculé par notre score
	return(new_canonical_transcript)

def write_file_for_canonical_transcript(new_canonical_transcript,output_file,exon_content,intron_content):
	file_out = open(output_file,"w")
	header = "Braunschweig id\tEnsembl id\t Gene id\tIntron number\tExon number\n"
	file_out.write(header)
	for B_id,E_id in new_canonical_transcript.items():
		if type(E_id) == list:
			E_id = E_id[0]
		exon_number = len(exon_content[E_id])
		intron_number = len(intron_content[B_id])
		gene_id = exon_content[E_id][0].gene_id
		message=B_id+"\t"+E_id+"\t"+gene_id+"\t"+str(intron_number)+"\t"+str(exon_number)+"\n"
		file_out.write(message)
	file_out.close()

def write_file_for_no_match_transcript(new_canonical_transcript,intron_content,exon_coord_in_braunch_transcript,no_match_file_name):
	file_out = open(no_match_file_name,"w")
	header = "Transcript id\tIntron number\n"
	file_out.write(header)
	for key,value in intron_content.items():
		if not key in new_canonical_transcript:
			message = key+"\t"+str(len(value))+"\n"
			file_out.write(message)
	file_out.close()

###Parsing des éléments en argument ####
create_db = True # Par defaut, gffutils va créer une base de donnée issue du fichier GTF pour le parser

opts, args = getopt.getopt(sys.argv[1:],'',['output_file=','annotation=','intron_ens=','gtf=','gtfdb=','no_match='])
for elmts in opts:
	if elmts[0] == '--output_file':
		output_file = elmts[1] # Nom du fichier de sortie

	elif elmts[0] == '--annotation':
		annotation_file = elmts[1] # Nom du premier fichier d'annotation

	elif elmts[0] == '--intron_ens':
		intron_ens_file = elmts[1] # Nom du fichier liant l'identifiant de chacun de nos introns à un id de gene Ensembl

	elif elmts[0] == '--gtf':
		gtf_file = elmts[1] # Nom du fichier au format GTF contenant toutes les annotations de la release (ici release 65)

	elif elmts[0] == '--gtfdb': # On demande si la base de donnée issue du fichier gtf à parser existe
		create_db = False
		gtf_db_file = elmts[1]

	elif elmts[0] == '--no_match': # On demande un nom de fichier pour la liste des transcrits qui ne concordent pas
		no_match_file_name = elmts[1]

print('Parsing de',intron_ens_file,'...')
IDlist,list_Ensembl_ids = extract_ensembl_id(intron_ens_file)
print('OK, IDs Ens associés aux IDs Introns')
print('Parsing de',annotation_file,'...')
intron_content= extract_coord(annotation_file,IDlist)
print('Ok, annotations sur les transcrits Braunschweig enregistrés')
print('Parsing de',gtf_file,'via',gtf_db_file,'...')
###Création de la base de données du GTF si n'existe pas####
if create_db == True:
	print('Création de la database issu du fichier GTF')
	db = gffutils.create_db(result[0], dbfn=gtf_db_file, force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
	print('Database OK')
exon_content = parsing_GTF(gtf_db_file,list_Ensembl_ids)
print('Ok, annotations sur les transcrits Ensembl enregistrés')
exon_coord_in_ensembl_transcript = exon_construction(exon_content)
print('Récupération des introns des transcrits Ensembl à partir des exons ...')
intron_coord_in_ensembl_transcript,retrieve_transcript_with_coord = transform_exon_to_intron(exon_coord_in_ensembl_transcript,exon_content)
print('OK, construction des introns pour chaque transcrit Ensembl terminé')
print('Détermination de l\'identifiant Ensembl correspondant au transcrit de Braunschweig ...')
canonical_transcripts = retrieve_canonical_transcript(intron_coord_in_ensembl_transcript,intron_content,retrieve_transcript_with_coord)
print('Détermination des doublons pour éviter la redondance ...')
doublon,dictionnary_without_redondance = determine_doublons(canonical_transcripts)
print(len(dictionnary_without_redondance),'/',len(intron_content),' transcrits de Braunschweig associés aux transcrits d\'Ensembl',sep='')
print('Calcul du score pour déterminer le transcrit Ensembl le plus proche du transcrit de Braunschweig dans les doublons ...')
exon_coord_in_braunch_transcript = transform_intron_to_exon(intron_content)
best_score = scoring_doublon(doublon,exon_coord_in_braunch_transcript,exon_coord_in_ensembl_transcript)
new_canonical_transcript = join_best_score_with_canonical(best_score,dictionnary_without_redondance)
print('Doublons éliminés, liste des transcrits Braunschweig associés à un transcrit Ensembl :',len(new_canonical_transcript))
print('Écriture du fichier',output_file,'=> Transcrits de référence')
write_file_for_canonical_transcript(new_canonical_transcript,output_file,exon_content,intron_content)
print('Écriture du fichier',no_match_file_name,'=> Transcrits sans identifiants Ensembl')
write_file_for_no_match_transcript(new_canonical_transcript,intron_content,exon_coord_in_braunch_transcript,no_match_file_name)
