# Troisième version du script, réalise l'annotation sur Ensembl, et non plus sur les transcrit de Braunch, on a deja les séquences d'intérêt, pas besoin de les construires
####################################################################
#Auteur : ROCHER Vincent                                           #
#                                                                  #
#But : Programme principal(remplace new_main.py)                   #
#      Ajoute nos propres annotations au jeu de données            #
#      -Séquence d'intérêt en début et fin de séquence des introns #
#      -Taux de GC                                                 #
#      -Taille de l'intron                                         #
#      -Position de l'intron (par rapport aux autres)              #
#      -Nombre d'intron par transcrit                              #
####################################################################

# Importation des modules
import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
from Bio.Seq import Seq

# Définition des classes :

# Classe qui va contenir l'annotation Ensembl de chaque intron
class IntronAnnotation(object):
    "Classe qui contiendra les annotations de chaque intron"
    def __init__(self,intron_id,transcript_id,gene_id,coords,seq,minimal_len,seq_left,seq_right):
        # Expression régulière qui va récupérer l'id, le debut et la fin du CDS dans la séquence du transcrit
        regex = re.compile('([\w]+)\(([0-9]+):([0-9]+)\)')
        result = regex.findall(transcript_id)
        self.trans_id = result[0][0] # ID du transcrit
        self.id = intron_id
        self.coords = coords
        self.gene_id = gene_id
        self.get_type = "intron"

        regex = re.compile("(chr[^\+-]+)([\+-])([0-9]+):([0-9]+)")
        result = regex.findall(self.coords)
        self.chr = result[0][0]
        self.strand = result[0][1]
        self.start = result[0][2]
        self.end = result[0][3]

        self.minimal_len = minimal_len
        self.seq = self.get_seq(seq)
        self.len_seq = len(seq)

        self.seq_left = seq_left
        self.seq_right = seq_right

    def get_seq(self,sequence):
        if len(sequence)<self.minimal_len:
            return "NA"
        else:
            return sequence
    # Calcul du taux de GC à partir de la séquence fasta
    def GCrate(self,total_len,nb_GC):
        if self.seq != "NA":
            # Détermination des bp à prendre en compte pour faire le calcul de GC : les bords de la séquence minimale de l'intron moins les nucleotides d'intéret
            intron_seq = self.seq[total_len:-total_len]
            A = intron_seq.count('A')
            C = intron_seq.count('C')
            G = intron_seq.count('G')
            T = intron_seq.count('T')
            # On détermine une valeur minimale de nucleotides pour calculer le taux de GC
            if ((A+T+C+G)>=nb_GC):
                #Calcul du taux de GC
                GCrate = ((G+C)/(A+T+G+C))*100
                GCrate = round(GCrate,2)
            else:
                GCrate = "NA"
            return GCrate
        else:
            return "NA"
    def seq_interest_construction(self,seq):
        if len(self.seq)<self.minimal_len:
            return "NA"
        if self.strand == "-":
            return seq.reverse_complement()
        else:
            return seq


class TranscriptInfo(object):
    "Classe qui contiendra les annotations de chaque transcrits"
    def __init__(self,transcript_id,list_object):
        self.id = transcript_id # ID du transcrit
        self.introns = [(intron[1].start,intron[1].end,intron[1]) for intron in list_object if intron[1].get_type == "intron"] # Liste des introns sous forme d'objet
        self.strand = self.introns[0][2].strand
# Fonctions pour déroulement du programme
#Fonction qui va récupérer les séquences fasta pour chaque exon/intron et les enregistrer dans un dictionnaire
def get_seq_fasta(file_name):
    seq_info = {}
    handle = open(file_name, "rU") #On ouvre le fichier en mode lecture

    for record in SeqIO.parse(handle, "fasta") : #On parse le fichier avec seqIO
        seq_id = record.id #On enregistre le nom de la sequence, ici les coordonnées
        seq = record.seq #On enregistre la sequence fasta dans une variable
        seq = seq.upper() #On met en majuscule au cas ou la séquence ne l'est pas
        seq_info[seq_id] = seq # On enregistre la séquence dans le dictionnaire
    handle.close()
    return(seq_info)

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
# Fonction qui va enregistrer chaque intron dans un objet, puis l'ajouter dans un dictionnaire pour chaque transcrit
def get_all_intron_for_transcript(file_name,seq_info,dico_info,minimal_len_seq= 50):

    intron_by_transcript = {}
    file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
    name = file_in.readline() #On enregistre la première ligne étant l'en tête

    for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
        content=line.split("\t")
        intron_id = content[0]
        trans_id = content[1]
        gene_id = content[2]
        coords = content[3].replace('\n', '')
        intron_seq = seq_info[intron_id]
        seq_of_interest = dico_info[intron_id]
        for elmt in seq_of_interest:
            if elmt[1] == "first":
                seq_left=elmt[0]
            else:
                seq_right=elmt[0]
        intron = IntronAnnotation(intron_id,trans_id,gene_id,coords,intron_seq,minimal_len_seq,seq_left,seq_right)
        if intron.trans_id in intron_by_transcript:
            intron_by_transcript[intron.trans_id].append((intron.coords,intron))
        else:
            intron_by_transcript[intron.trans_id]=[]
            intron_by_transcript[intron.trans_id].append((intron.coords,intron))
    file_in.close()
    return(intron_by_transcript)

def get_dictionnary_of_transcript(intron_by_transcript):
    transcript_complete = {}
    for trans_id, list_intron in intron_by_transcript.items():
        real_trans_id = trans_id.split("(")[0]
        full_list = []
        full_list.extend(list_intron)
        # Fonction de tri, on a des coordonnées celon la conventien : chr[+/-]debut:fin
        # On va donc prendre en compte uniquement la fin, donc on split la chaine pour obetenir fin, qu'on converti en entier, et qu'on tri dans l'ordre
        full_list.sort(key=lambda x: int(x[0].split(':')[1]))
        # On créer l'objet transcrit, qui contiendra les introns et les exons du transcrit
        transcript_object = TranscriptInfo(trans_id,full_list)
        transcript_complete[real_trans_id]=transcript_object
    return(transcript_complete)

def print_data(file_name,transcript_complete, len_GC = 40):
    file_out = open(file_name,"w")
    header = "Intron id\tTranscript id\tGene id\tCoordinates\tSeq left\tSeq right\tIntron length\tIntron position\tTotal intron\tGC rate\n"
    file_out.write(header)
    for trans_id, trans_object in transcript_complete.items():
        intron_list = trans_object.introns
        if trans_object.strand == "-":
            intron_list.reverse()
        numbers_of_introns = len(intron_list)
        for intron_and_coords in intron_list:
            intron = intron_and_coords[2]
            # On récupère la taille de l'intron :
            intron_length = intron.len_seq
            # On récupère la position de l'intron par rapport aux autres :
            position = intron_list.index(intron_and_coords)+1
            # On récupère les séquences d'intérêts :
            seq_left = intron.seq_interest_construction(intron.seq_left)
            seq_right = intron.seq_interest_construction(intron.seq_right)
            # On calcule le taux de GC :
            total_len = len(seq_left) # Détermine les bords de la séquence minimale de l'intron moins les nucleotides d'intéret
            intron_GC = intron.GCrate(total_len,len_GC)
            line = intron.id+"\t"+trans_id+"\t"+intron.gene_id+"\t"+intron.coords+"\t"+str(seq_left)+"\t"+str(seq_right)+"\t"+str(intron.len_seq)+"\t"+str(position)+"\t"+str(numbers_of_introns)+"\t"+str(intron_GC)+"\n"
            file_out.write(line)
    file_out.close()



# Interface avec l'utilisateur :
opts, args = getopt.getopt(sys.argv[1:],'',['fasta_interest=','liste_intron=','fasta=','len_GC=','len_seq=','output_file='])
for elmts in opts:
    if elmts[0] == '--liste_intron':
        liste_intron = elmts[1] # Fichier contenant les coordonnées pour chaque intron : list_intron.tab

    elif elmts[0] == '--fasta':
        fasta_file = elmts[1] # Fichier contenant les séquences au format fasta des exons et des introns : seq.fa

    elif elmts[0] == '--fasta_interest':
        fasta_file_interest = elmts[1] # Fichier contenant les séquences au format fasta des exons et des introns : seq.fa

    elif elmts[0] == '--len_GC':
        nb_GC = int(elmts[1]) # détermine la longueur minimale de l'intron -(nb_nt_exon+nb_nt_intron) pour calculer le taux de GC

    elif elmts[0] == '--len_seq':
        len_seq = int(elmts[1]) # détermine la taille minimale de la séquence pour annoter

    elif elmts[0] == '--output_file':
        output_file_name = elmts[1] # Nom du fichier qui contiendra les anotations sur les introns

# Lancement des fonctions :
# Récupération des séquences fasta
seq_info = get_seq_fasta(fasta_file)
# Récuperation des séquences d'intérêts
dico_info = get_seq_of_interest(fasta_file_interest)
# Récupération des introns
intron_by_transcript = get_all_intron_for_transcript(liste_intron,seq_info,dico_info,len_seq)
# Association des introns et des exons dans un transcrit
transcript_complete = get_dictionnary_of_transcript(intron_by_transcript)
# Ecriture des annotations
print_data(transcript_complete=transcript_complete, len_GC = nb_GC,file_name=output_file_name)
