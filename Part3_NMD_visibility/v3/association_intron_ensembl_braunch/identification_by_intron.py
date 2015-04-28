import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
from collections import Counter

# Définition des classes
class IntronEns(object):
    "Classe contenant l'id de l'intron, ses coordonnées et les id Gene et Transcript d'ensembl"
    def __init__(self,intron_id,gene_id,trans_id,coords):
        self.id = intron_id
        self.gene_id = gene_id
        self.trans_id = trans_id
        self.coords = coords

# Fonctions principales :
# Fonction qui va récupérer les coordonnées de chaque introns de chaque transcrit Ensembl, et l'enregistrer dans un dictionnaire
def get_ens_intron_coord(file_name):
    Ens_data = {}
    file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
    name = file_in.readline() #On enregistre la première ligne étant l'en tête
    for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
        content=line.split("\t")
        intron_id = content[0]
        trans_id = content[1]
        gene_id = content[2]
        coords = content[3].replace('\n', '')
        intron = IntronEns(intron_id,gene_id,trans_id,coords)
        if intron.trans_id in Ens_data:
            Ens_data[intron.trans_id].append((intron.id,intron))
        else:
            Ens_data[intron.trans_id]=[]
            Ens_data[intron.trans_id].append((intron.id,intron))
    return(Ens_data)

# Fonction qui va récupérer les coordonnées et les identifiants des introns de Braunchweig
def get_braunch_intron_coord(file_name):
    Braunch_coord = {}
    file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
    name = file_in.readline()
    for line in file_in:
        content=line.split("\t") #On split le contenu dans une liste
        regex = re.compile('^(chr[\w]+[-\+])[0-9]+:([0-9]+)_([0-9]+):[0-9]+')
        result = regex.findall(content[1])
        new_coord = result[0][0]+str(int(result[0][1])+1)+":"+str(int(result[0][2])-1)
        if new_coord in Braunch_coord:
            Braunch_coord[new_coord].append(content[0])
        else:
            Braunch_coord[new_coord]=[]
            Braunch_coord[new_coord].append(content[0])
    file_in.close
    return(Braunch_coord)

# Fonction d'association entre un identifiant Ens et Braunch pour chaque transcrit de notre jeu de données
def association_between_Ensembl_and_Braunch(Ens_data,Braunch_coord):
    association = {}
    for trans_id,intron_list in Ens_data.items():
        for intron in intron_list:
            intron_coords = intron[1].coords
            if intron_coords in Braunch_coord:
                braunch_id = Braunch_coord[intron_coords]
                association[intron[0]]=braunch_id
    return(association)

# Fonction qui va vérifier pour chaque transcrit si l'association pour chacun de ses introns est faites, et trier les transcrit en fonction de ce paramètre
def check_association_for_transcript(association,Ens_data):
    dataset = {}
    more_one_id = {}
    for trans_id,intron_list in Ens_data.items():
        dataset[trans_id]=[]
        for intron in intron_list:
            if intron[0] in association:
                braunch_ids = association[intron[0]]
                if len(braunch_ids)>1:
                    more_one_id[trans_id]=[]
                dataset[trans_id].append((intron[0],braunch_ids))
            else:
                dataset[trans_id].append((intron[0],"NA"))
    return(dataset,more_one_id)





# Fonction qui pour tous les transcrits possédant un intron avec plusieurs id braunch aossocié va déterminer lequel choisir
def no_doublon(dataset,more_one_id):
    Dataset_without_doublon = more_one_id
    for key in more_one_id:
        intron_list = dataset[key]
        # Calcul du nombre de fois ou chaque identifiant est retrouvé :
        regex = re.compile('^[^:]+:(.+):[0-9]{,3}')
        # Double liste comprehension on parcours chaque identifiant ensembl de la liste, pour chacun de ces identifiants, on parcours les identifiants braunch qui match pour celui-ci, et on extrait son identifiant de transcrit avec une expression régulière, on obtient une liste d'identifiant de transcrit pour chaque identifiant d'intron ensembl, puis on compte le nombre de fois ou chaque id de transcrit est retrouvé, ce qui nous donne un dictionnaire avec les identifiants présents
        id_braunch_by_transcript = Counter([regex.findall(sub_elmt)[0] for elmt in intron_list for sub_elmt in elmt[1] if type(elmt[1]) != str])
        # Selection de l'id Braunch représentatif :
        # Parcours de chacun des introns du transcrit
        for intron_ids in intron_list:
            id_ens = intron_ids[0]
            liste_id_braunch = intron_ids[1]
            if type(liste_id_braunch) != str:
                # On tri les identifiants pour qu'ils soient dans le même ordre d'apparition
                liste_id_braunch.sort(key=lambda x: x.split(':')[1])
                # Par défaut, on considère que le meilleur id est le premier
                best_id = regex.findall(liste_id_braunch[0])[0]
                # Puis on parcours chacun des identifiants braunch pour un identifiant ensembl
                for one_id in liste_id_braunch:
                    one_id_transcript = regex.findall(one_id)[0]
                    if id_braunch_by_transcript[one_id_transcript] > id_braunch_by_transcript[best_id]:
                        best_id = one_id_transcript
                complete_id = [elmt for elmt in liste_id_braunch if re.search(best_id, elmt)]
                Dataset_without_doublon[key].append((id_ens,complete_id))
            else:
                Dataset_without_doublon[key].append((id_ens,"NA"))
    return(Dataset_without_doublon)

def fusion_dataset(Dataset_complete,Dataset_without_doublon):
    for trans_id,intron_list in Dataset_without_doublon.items():
        Dataset_complete[trans_id]=intron_list
    return(Dataset_complete)



# opts, args = getopt.getopt(sys.argv[1:],'',['output_file_trans=','output_file_intron=','intron_braunch=','intron_ens='])
# for elmts in opts:
#     if elmts[0] == '--output_file_trans':
#         output_file_trans = elmts[1]
#     elif elmts[0] == '--output_file_intron':
#         output_file_intron = elmts[1]
#     elif elmts[0] == '--intron_braunch':
#         annotation_braunch = elmts[1]
#     elif elmts[0] == '--intron_ens':
#         annotation_ens = elmts[1]
Ens_data = get_ens_intron_coord('/Users/Vincent/Documents/STAGE_LYON/Projet_stage_RV/results/file_for_NMD/list_and_fasta_v2/list_intron.tab')
Braunch_coord = get_braunch_intron_coord('/Users/Vincent/Documents/STAGE_LYON/Projet_stage_RV/data/DataBraunschweig/TableS6_IntronsHuman.tab')
association = association_between_Ensembl_and_Braunch(Ens_data,Braunch_coord)
dataset,more_one_id = check_association_for_transcript(association,Ens_data)
Dataset_without_doublon = no_doublon(dataset,more_one_id)
Dataset_complete = fusion_dataset(dataset,Dataset_without_doublon)
