import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python

class IntronBraunch(object):
	"Classe contenant l'id de l'intron, ses coordonnées et les id Gene et Transcript d'ensembl"
	def __init__(self,intron_id,gene_id,trans_id,coords):
		self.id = intron_id
		self.gene_id = gene_id
		self.trans_id = trans_id
		self.coords = coords


class IntronEns(object):
	"Classe contenant l'id de l'intron, ses coordonnées et les id Gene et Transcript d'ensembl"
	def __init__(self,intron_id,gene_id,trans_id,coords):
		self.id = intron_id
		self.gene_id = gene_id
		self.trans_id = trans_id
		self.coords = coords

def extract_Ens_intron_coord(file_name):
	Ens_coord = {}
	Ens_data = {}
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline()

	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t") #On split le contenu dans une liste
		intron_coord = content[3].replace('\n','')
		if intron_coord in Ens_coord:
			Ens_coord[intron_coord].append(content[0])
		else:
			Ens_coord[intron_coord]=[]
			Ens_coord[intron_coord].append(content[0])
		intron_ens = IntronEns(content[0],content[1],content[2],intron_coord)
		if content[2] in Ens_data:
			Ens_data[content[2]].append((intron_ens.id,intron_ens))
		else:
			Ens_data[content[2]]=[]
			Ens_data[content[2]].append((intron_ens.id,intron_ens))
	file_in.close()
	return(Ens_coord,Ens_data)

def extract_Braunch_intron(file_name):
	file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
	name = file_in.readline()
	Braunch_intron = {}
	Braunch_coord = {}
	regex1 = re.compile('^[^:]+:(.+):[0-9]{,3}') # On définit une expression régulière qui va capturer le nom du gène et l'id refseq pour définir un transcrit Braunch
	for line in file_in:
		content=line.split("\t") #On split le contenu dans une liste
		result = regex1.findall(content[0])
		braunch_id = result[0]
		regex2 = re.compile('^(chr[\w]+[-\+])[0-9]+:([0-9]+_[0-9]+):[0-9]+')
		result = regex2.findall(content[1])
		new_coord = result[0][0]+result[0][1]
		if braunch_id in Braunch_intron:
			Braunch_intron[braunch_id].append((content[0],new_coord))
		else:
			Braunch_intron[braunch_id]=[]
			Braunch_intron[braunch_id].append((content[0],new_coord))
		if new_coord in Braunch_coord:
			Braunch_coord[new_coord].append(content[0])
		else:
			Braunch_coord[new_coord]=[]
			Braunch_coord[new_coord].append(content[0])
	file_in.close
	return(Braunch_intron,Braunch_coord)
def association_between_Braunch_Ensembl_coords(Braunch_intron,Ens_data,Ens_coord):
	Dataset = {}
	for key,value in Braunch_intron.items():
		for elmt in value:
			if elmt[1] in Ens_coord:
				Ens_ids = Ens_coord[elmt[1]]
				Dataset[elmt]=Ens_ids
	return(Dataset)

def association_between_Ensembl_Braunch_coords(Braunch_coord,Ens_data,Ens_coord):
	Dataset = {}
	for key,value in Ens_data.items():
		for elmt in value:
			if elmt[1].coords in Braunch_coord:
				Braunch_ids = Braunch_coord[elmt[1].coords] 
				Dataset[elmt[1].id]=Braunch_ids
	return(Dataset)

def complete_annotation(Ens_data,Dataset_all_canonical):
	Dataset_complete = {}
	more_one_id = {}
	match = nomatch = no_complete = 0
	for key,value in Ens_data.items():
		verif = 0
		for intron_ens in value:
			if intron_ens[0] in Dataset_all_canonical:
				verif+=1
				braunch_ids = Dataset_all_canonical[intron_ens[0]]
				if len(braunch_ids) > 1:
					more_one_id[key]=[]
				if key in Dataset_complete:
					Dataset_complete[key].append((intron_ens[0],braunch_ids))
				else:
					Dataset_complete[key]=[]
					Dataset_complete[key].append((intron_ens[0],braunch_ids))
			else:
				if key in Dataset_complete:
					Dataset_complete[key].append((intron_ens[0],"None"))
				else:
					Dataset_complete[key]=[]
					Dataset_complete[key].append((intron_ens[0],"None"))
		if verif == 0:
			nomatch+=1
		elif verif >= 1:
			if verif == len(value):
				match+=1
			else:
				no_complete+=1
	print("Transcrits sans introns :",nomatch)
	print("Transcrits avec un ou plusieurs introns manquants :",no_complete)
	print("Transcrits avec tous ses introns annotés :",match)
	return(Dataset_complete,more_one_id)


def no_doublon(Dataset_complete,more_one_id):
	Dataset_without_doublon = more_one_id
	for key in more_one_id:
		intron_list = Dataset_complete[key]
		id_braunch_by_transcript = {}
		for intron_ids in intron_list:
			liste_id_braunch = intron_ids[1]
			if type(liste_id_braunch) != str:
				for elmt in liste_id_braunch:
					regex = re.compile('^[^:]+:(.+):[0-9]{,3}')
					result = regex.findall(elmt)
					if result[0] in id_braunch_by_transcript:
						id_braunch_by_transcript[result[0]]+=1
					else:
						id_braunch_by_transcript[result[0]]=1
		# Selection de l'id Braunch représentatif : 
		for intron_ids in intron_list:
			id_ens = intron_ids[0]
			liste_id_braunch = intron_ids[1]
			if type(liste_id_braunch) != str:
				# Par défaut, on considère que le meilleur id est le premier
				regex = re.compile('^[^:]+:(.+):[0-9]{,3}')
				result = regex.findall(liste_id_braunch[0])
				best_id = result[0]
				for one_id in liste_id_braunch:
					result = regex.findall(one_id)
					if id_braunch_by_transcript[result[0]] > id_braunch_by_transcript[best_id]:
						best_id = result[0]
				complete_id = [elmt for elmt in liste_id_braunch if re.search(best_id, elmt)]
				Dataset_without_doublon[key].append((id_ens,complete_id))
			else:
				Dataset_without_doublon[key].append((id_ens,"None"))
	return(Dataset_without_doublon)


Ens_coord,Ens_data = extract_Ens_intron_coord(sys.argv[1])
Braunch_intron,Braunch_coord = extract_Braunch_intron(sys.argv[2])
Dataset_all_Braunch = association_between_Braunch_Ensembl_coords(Braunch_intron,Ens_data,Ens_coord)
Dataset_all_canonical = association_between_Ensembl_Braunch_coords(Braunch_coord,Ens_data,Ens_coord)
Dataset_complete,more_one_id = complete_annotation(Ens_data,Dataset_all_canonical)
Dataset_without_doublon = no_doublon(Dataset_complete,more_one_id)
