import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python


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
					Dataset_complete[key].append((intron_ens[0],"NA"))
				else:
					Dataset_complete[key]=[]
					Dataset_complete[key].append((intron_ens[0],"NA"))
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
				Dataset_without_doublon[key].append((id_ens,"NA"))
	return(Dataset_without_doublon)

def fusion_dataset(Dataset_complete,Dataset_without_doublon):
	for trans_id,intron_list in Dataset_without_doublon.items():
		Dataset_complete[trans_id]=intron_list
	return(Dataset_complete)

#Fonction qui va créer un fichier contenant la liste des transcrits ensembl annotés, indiquer si celui-ci est complètement ou partiellement annoté
def transcript_information_file(file_name,Dataset_complete):
	file_in=open(file_name,"w")
	header = "Ensembl Transcript id\tEnsembl Gene id\tIntron Annotation\tAnnotation Status\n"
	file_in.write(header)
	for trans_id,intron_list in Dataset_complete.items():
		None_comptor = 0
		patchwork_id = {}
		for intron_annotation in intron_list:
			ens_id = intron_annotation[0]
			gene_id = ens_id.split(":")[0]
			if type(intron_annotation[1]) == list:
				association = intron_annotation[1][0]
			else:
				association = intron_annotation[1]
			if association != "NA" : 
				regex = re.compile('^[^:]+:(.+):[0-9]{,3}')
				result = regex.findall(association)
				patchwork_id[result[0]]="" 
			else : None_comptor+=1 
		annotation = "Complete"
		if None_comptor == len(intron_annotation) : annotation = "NA" 
		elif None_comptor >= 1 : annotation = "Incomplete" 
		if annotation == "Complete":
			status = "Full"
			if len(patchwork_id) > 1 : 
				status = "Patchwork"
		else:
			status = "NA"
		line = trans_id+"\t"+gene_id+"\t"+annotation+"\t"+status+"\n"
		file_in.write(line)
	file_in.close()


def intron_annotation_file(file_name,Dataset_complete,Ens_data):
	file_in=open(file_name,"w")
	header = "Ensembl intron id\tRefSeq intron id\tTranscript id\tGene id\tCoordinates\n"
	file_in.write(header)
	for trans_id,intron_list in Dataset_complete.items():
		for intron_annotation in intron_list:
			intron_id = intron_annotation[0]
			intron_list_object = Ens_data[trans_id]
			intron_object = [elmt[1] for elmt in intron_list_object if elmt[0] == intron_id]
			if type(intron_annotation[1]) == list:
				association = intron_annotation[1][0]
			else:
				association = intron_annotation[1]
			line = intron_id+"\t"+association+"\t"+trans_id+"\t"+intron_object[0].gene_id+"\t"+intron_object[0].coords+"\n"
			file_in.write(line)
	file_in.close()



opts, args = getopt.getopt(sys.argv[1:],'',['output_file_trans=','output_file_intron=','intron_braunch=','intron_ens='])
for elmts in opts:
	if elmts[0] == '--output_file_trans':
		output_file_trans = elmts[1]
	elif elmts[0] == '--output_file_intron':
		output_file_intron = elmts[1]
	elif elmts[0] == '--intron_braunch':
		annotation_braunch = elmts[1]
	elif elmts[0] == '--intron_ens':
		annotation_ens = elmts[1]


Ens_coord,Ens_data = extract_Ens_intron_coord(annotation_ens)
Braunch_intron,Braunch_coord = extract_Braunch_intron(annotation_braunch)
Dataset_all_Braunch = association_between_Braunch_Ensembl_coords(Braunch_intron,Ens_data,Ens_coord)
Dataset_all_canonical = association_between_Ensembl_Braunch_coords(Braunch_coord,Ens_data,Ens_coord)
Dataset_complete,more_one_id = complete_annotation(Ens_data,Dataset_all_canonical)
Dataset_without_doublon = no_doublon(Dataset_complete,more_one_id)
Dataset_complete = fusion_dataset(Dataset_complete,Dataset_without_doublon)
#Écriture des fichiers de sortie : 
transcript_information_file(output_file_trans,Dataset_complete)
intron_annotation_file(output_file_intron,Dataset_complete,Ens_data)
