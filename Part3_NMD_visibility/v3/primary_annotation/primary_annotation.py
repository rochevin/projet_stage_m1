# Seconde version du script, réalise l'annotation sur Ensembl, et non plus sur les transcrit de Braunch
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
    def __init__(self,intron_id,transcript_id,gene_id,coords,seq,minimal_len):
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

        self.seq = self.get_seq(seq,minimal_len)
        self.len_seq = len(seq)

    def get_seq(self,sequence,minimal_len):
        if len(sequence)<minimal_len:
            return "NA"
        if self.strand == "-":
            return sequence.reverse_complement()
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
    # Détermine nos séquences d'intéret, par exemple (-20/+30) et (-30/+20)
    def seq_begin(self,len_intron_seq):
        if self.seq != "NA":
            return self.seq[0:len_intron_seq]
        else:
            return "NA"
    def seq_end(self,len_intron_seq):
        if self.seq != "NA":
            taille = self.len_seq
            return self.seq[:len_intron_seq]
        else:
            return "NA"


class ExonInfo(object):
    "Classe qui contiendra les annotations de chaque exon"
    def __init__(self,exon_id,transcript_id,gene_id,coords,seq):
        # Expression régulière qui va récupérer l'id, le debut et la fin du CDS dans la séquence du transcrit
        regex = re.compile('([\w]+)\(([0-9]+):([0-9]+)\)')
        result = regex.findall(transcript_id)
        self.trans_id = result[0][0] # ID du transcrit
        self.id = exon_id
        self.coords = coords
        self.gene_id = gene_id
        self.get_type = "exon"

        regex = re.compile("(chr[^\+-]+)([\+-])([0-9]+):([0-9]+)")
        result = regex.findall(self.coords)
        self.chr = result[0][0]
        self.strand = result[0][1]
        self.start = result[0][2]
        self.end = result[0][3]

        self.seq = self.get_seq(seq)
        self.len_seq = len(self.seq)

    def get_seq(self,sequence):
        if self.strand == "-":
            return sequence.reverse_complement()
        else:
            return sequence

    def flanking_seq_left(self,len_flanking_seq):
        return self.seq[-int(len_flanking_seq):]

    def flanking_seq_right(self,len_flanking_seq):
        return self.seq[:int(len_flanking_seq)]

class TranscriptInfo(object):
    "Classe qui contiendra les annotations de chaque transcrits"
    def __init__(self,transcript_id,list_object):
        self.id = transcript_id # ID du transcrit
        self.transcript_content = [(elmt[1].start,elmt[1].end, elmt[1]) for elmt in list_object] # Liste d'objets contenus dans le transcrits (introns + exons)
        self.exons = [(exon[1].start,exon[1].end,exon[1]) for exon in list_object if exon[1].get_type == "exon"] # Liste des exons sous forme d'objet
        self.introns = [(intron[1].start,intron[1].end,intron[1]) for intron in list_object if intron[1].get_type == "intron"] # Liste des introns sous forme d'objet
        self.strand = self.transcript_content[0][2].strand
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

# Fonction qui va enregistrer chaque intron dans un objet, puis l'ajouter dans un dictionnaire pour chaque transcrit
def get_all_intron_for_transcript(file_name,seq_info,minimal_len_seq= 50):

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
        intron = IntronAnnotation(intron_id,trans_id,gene_id,coords,intron_seq,minimal_len_seq)
        if intron.trans_id in intron_by_transcript:
            intron_by_transcript[intron.trans_id].append((intron.coords,intron))
        else:
            intron_by_transcript[intron.trans_id]=[]
            intron_by_transcript[intron.trans_id].append((intron.coords,intron))
    file_in.close()
    return(intron_by_transcript)

# Fonction qui va enregistrer chaque exon dans un objet, puis l'ajouter dans un dictionnaire pour chaque transcrit
def get_all_exon_for_transcript(file_name,seq_info):

    exon_by_transcript = {}
    file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
    name = file_in.readline() #On enregistre la première ligne étant l'en tête

    for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
        content=line.split("\t") # On split la ligne dans une liste
        # On récupère chaque élément de la liste dans une variable
        exon_id = content[0]
        trans_id = content[1]
        gene_id = content[2]
        coords = content[3].replace('\n', '')
        # On récupère la séquence au format BioSeq via le dictionnaire
        exon_seq = seq_info[exon_id]
        # On enregistre le tout dans l'objet
        exon = ExonInfo(exon_id,trans_id,gene_id,coords,exon_seq)
        # Puis on ajoute l'objet dans un dictionnaire, annoté par transcrit
        if exon.trans_id in exon_by_transcript:
            exon_by_transcript[exon.trans_id].append((exon.coords,exon))
        else:
            exon_by_transcript[exon.trans_id]=[]
            exon_by_transcript[exon.trans_id].append((exon.coords,exon))
    file_in.close()
    return(exon_by_transcript)

def association_between_intron_and_exon(intron_by_transcript,exon_by_transcript):
    transcript_complete = {}
    for trans_id, list_exon in exon_by_transcript.items():
        real_trans_id = trans_id.split("(")[0]
        if trans_id in intron_by_transcript:
            introns = intron_by_transcript[trans_id]
            full_list = []
            full_list.extend(list_exon)
            full_list.extend(introns)
            # Fonction de tri, on a des coordonnées celon la conventien : chr[+/-]debut:fin
            # On va donc prendre en compte uniquement la fin, donc on split la chaine pour obetenir fin, qu'on converti en entier, et qu'on tri dans l'ordre
            full_list.sort(key=lambda x: int(x[0].split(':')[1]))
            # On créer l'objet transcrit, qui contiendra les introns et les exons du transcrit
            transcript_object = TranscriptInfo(trans_id,full_list)
            transcript_complete[real_trans_id]=transcript_object
    return(transcript_complete)

def print_data(file_name,transcript_complete, exon_len = 20, intron_len = 30, len_GC = 40):
    file_out = open(file_name,"w")
    header = "Intron id\tTranscript id\tGene id\tCoordinates\tSeq left\tSeq right\tIntron length\tIntron position\tTotal intron\tGC rate\n"
    file_out.write(header)
    for trans_id, trans_object in transcript_complete.items():
        intron_list = trans_object.introns
        transcript_content = trans_object.transcript_content
        if trans_object.strand == "-":
            transcript_content.reverse()
            intron_list.reverse()
        numbers_of_introns = len(intron_list)
        for intron_and_coords in intron_list:
            intron = intron_and_coords[2]
            # On récupère la séquence de l'intron :
            intron_seq = intron.seq
            # On récupère la taille de l'intron :
            intron_length = intron.len_seq
            # On récupère la position de l'intron par rapport aux autres :
            position = intron_list.index(intron_and_coords)+1
            # On récupère les séquences d'intérêts :
            if intron_seq != "NA":
                # On récupère les séquences bordantes des exons :
                seq_left_exon,seq_right_exon = get_flanking_exon_for_intron(transcript_content,intron_and_coords,exon_len)
                # On récupère le début et la fin de la séquence de l'intron :
                intron_seq_left = intron.seq_begin(intron_len)
                intron_seq_right = intron.seq_end(intron_len)
                # On construit les deux séquences d'intérêt :
                seq_left = sequence_of_interest(intron_seq_left,seq_left_exon,"left")
                seq_right = sequence_of_interest(intron_seq_right,seq_right_exon,"right")
            else:
                seq_left = seq_right = "NA"
            # On calcule le taux de GC :
            total_len = int(intron_len) # Détermine les bords de la séquence minimale de l'intron moins les nucleotides d'intéret
            intron_GC = intron.GCrate(total_len,len_GC)
            line = intron.id+"\t"+trans_id+"\t"+intron.gene_id+"\t"+intron.coords+"\t"+str(seq_left)+"\t"+str(seq_right)+"\t"+str(intron.len_seq)+"\t"+str(position)+"\t"+str(numbers_of_introns)+"\t"+str(intron_GC)+"\n"
            file_out.write(line)
    file_out.close()

def get_flanking_exon_for_intron(full_list,one_intron,number_of_nucleotides = 20):

    intron_position_in_list = full_list.index(one_intron)
    # On récupère l'exon flanquant gauche :
    exon_left = full_list[intron_position_in_list-1]
    # Si la séquence de l'exon est inférieure au nombre de nucleotides qu'on veut récupérer :
    if len(exon_left[2].seq) < number_of_nucleotides:
        # On calcule le nombre de base qu'il nous manque
        if exon_left == full_list[0]:
            seq_exon_left = exon_left[2].flanking_seq_left(number_of_nucleotides)
        else:
            missing_length = number_of_nucleotides-exon_left[2].len_seq
            exon_position_in_list = full_list.index(exon_left)
            intron_before_left_exon = full_list[exon_position_in_list-1]
            # On récupère de l'intron les nucleotides qu'il nous manque
            missing_nucleotides = intron_before_left_exon[2].seq[-missing_length:]
            seq_exon_left = missing_nucleotides+exon_left[2].flanking_seq_left(number_of_nucleotides)
    else:
        seq_exon_left = exon_left[2].flanking_seq_left(number_of_nucleotides)
    # On récupère l'exon flanquant droit :
    exon_right = full_list[intron_position_in_list+1]
    # Si la séquence de l'exon est inférieure au nombre de nucleotides qu'on veut récupérer :
    if len(exon_right[2].seq) < number_of_nucleotides:
        # On calcule le nombre de base qu'il nous manque
        if exon_right == full_list[-1]:
            seq_exon_right = exon_right[2].flanking_seq_right(number_of_nucleotides)
        else:
            missing_length = number_of_nucleotides-exon_right[2].len_seq
            exon_position_in_list = full_list.index(exon_right)
            intron_after_right_exon = full_list[exon_position_in_list+1]
            # On récupère de l'intron les nucleotides qu'il nous manque
            missing_nucleotides = intron_after_right_exon[2].seq[-missing_length:]
            seq_exon_right = missing_nucleotides+exon_right[2].flanking_seq_right(number_of_nucleotides)
    else:
        seq_exon_right = exon_right[2].flanking_seq_right(number_of_nucleotides)
    return (seq_exon_left,seq_exon_right)

def sequence_of_interest(intron_seq,exon_seq,position):
    return exon_seq+intron_seq if position =="left" else intron_seq+exon_seq



# Interface avec l'utilisateur :
opts, args = getopt.getopt(sys.argv[1:],'',['liste_exon=','liste_intron=','fasta=','exon_len=','intron_len=','len_GC=','len_seq=','output_file='])
for elmts in opts:
    if elmts[0] == '--liste_exon':
        liste_exon = elmts[1] # Fichier contenant les coordonnées pour chaque exon : liste_exon.tab

    elif elmts[0] == '--liste_intron':
        liste_intron = elmts[1] # Fichier contenant les coordonnées pour chaque intron : list_intron.tab

    elif elmts[0] == '--fasta':
        fasta_file = elmts[1] # Fichier contenant les séquences au format fasta des exons et des introns : seq.fa

    elif elmts[0] == '--exon_len':
        nb_nt_exon = int(elmts[1]) # détermine le nombre de nucleotides à prendre en compte dans les exons flanquants

    elif elmts[0] == '--intron_len':
        nb_nt_intron = int(elmts[1]) # détermine le nombre de nucleotides à prendre en compte aux extrémités de l'intron

    elif elmts[0] == '--len_GC':
        nb_GC = int(elmts[1]) # détermine la longueur minimale de l'intron -(nb_nt_exon+nb_nt_intron) pour calculer le taux de GC

    elif elmts[0] == '--len_seq':
        len_seq = int(elmts[1]) # détermine la taille minimale de la séquence pour annoter

    elif elmts[0] == '--output_file':
        output_file_name = elmts[1] # Nom du fichier qui contiendra les anotations sur les introns

# Lancement des fonctions :
# Récupération des séquences fasta
seq_info = get_seq_fasta(fasta_file)
# Récupération des exons
exon_by_transcript = get_all_exon_for_transcript(liste_exon,seq_info)
# Récupération des introns
intron_by_transcript = get_all_intron_for_transcript(liste_intron,seq_info)
# Association des introns et des exons dans un transcrit
transcript_complete = association_between_intron_and_exon(intron_by_transcript,exon_by_transcript)
# Ecriture des annotations
print_data(transcript_complete=transcript_complete, exon_len = nb_nt_exon, intron_len = nb_nt_intron, len_GC = nb_GC,file_name=output_file_name)
