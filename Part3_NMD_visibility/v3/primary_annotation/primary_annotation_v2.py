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
import tempfile # On importe tempfile pour créer des fichiers temporaires
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
    def seq_begin(self,total_len_of_interest):
        if self.seq != "NA":
            return self.seq[:total_len_of_interest]
        else:
            return "NA"

    def seq_end(self,total_len_of_interest):
        if self.seq != "NA":
            return self.seq[-total_len_of_interest:]
        else:
            return "NA"

    def seq_without_exon(self,exon_len):
        return len(self.seq[exon_len:-exon_len])

# Fonctions pour déroulement du programme
# Fonction qui va générer un fichier BED en temporaire pour récupérer les séquences au format fasta avec TwoBitTofa
def BED_generation(file_name):
    file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
    name = file_in.readline() #On enregistre la première ligne étant l'en tête "Intron  Coordinates type    txRegion    PTCstatus"

    file_out = tempfile.NamedTemporaryFile(mode="w+",delete=False,suffix='.bedPos')
    file_out_name = file_out.name
    #file_out=open(file_out_name,"w") # On créer le fichier BED qui va contenir les données
    for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
        content=line.split("\t")
        intron_id = content[0]
        coords = content[3].replace('\n', '')
        regex = re.compile("(chr[^\+-]+)[\+-]([0-9]+):([0-9]+)")
        result = regex.findall(coords)
        intron_chr = result[0][0]
        # Les positions doivent être données en base 0 et pour la fin en non inclusif, c'est à dire la première position après la séquence, on doit donc enlever 1 au start pour passer en base 0 mais on ne touche pas à la position de fin
        intron_start = int(result[0][1])-1
        intron_end = result[0][2]
        bed_line = intron_chr+"\t"+str(intron_start)+"\t"+intron_end+"\t"+intron_id+"\n"
        # Écriture du fichier BED
        file_out.write(bed_line)
    file_in.close()
    file_out.close()
    return file_out_name

def get_fasta_seq(file_name):
    seq_info = {} # On définit un dictionnaire qui contient la séquence en valeur et ses coordonnées en clé

    handle = open(file_name, "rU") #On ouvre le fichier en mode lecture

    for record in SeqIO.parse(handle, "fasta") : #On parse le fichier avec seqIO
        seq_id = record.id #On enregistre le nom de la sequence, ici les coordonnées
        seq = record.seq #On enregistre la sequence fasta dans une variable
        seq= seq.upper() #On met en majuscule au cas ou la séquence ne l'est pas
        seq_info[seq_id]=seq # On enregistre la séquence dans le dictionnaire avec ses coordonnées comme clé

    handle.close()
    return seq_info

# Fonction qui va enregistrer chaque intron dans un objet, puis l'ajouter dans un dictionnaire pour chaque transcrit
def get_all_intron_for_transcript(file_name,seq_info,minimal_len_seq = 50):

    intron_by_transcript = {}
    intron_coords = {}
    file_in=open(file_name,"r") #On ouvre le fichier qui contient les données en mode lecture
    name = file_in.readline() #On enregistre la première ligne étant l'en tête

    for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
        content=line.split("\t")
        intron_id = content[0]
        trans_id = content[1]
        gene_id = content[2]
        coords = content[3].replace('\n', '')
        coords_for_fasta_file = formating_coord(full_coordinates=coords)
        intron_seq = seq_info[coords_for_fasta_file]
        intron = IntronAnnotation(intron_id,trans_id,gene_id,coords,intron_seq,minimal_len_seq)
        if intron.trans_id in intron_by_transcript:
            intron_by_transcript[intron.trans_id].append((intron.start,intron))
        else:
            intron_by_transcript[intron.trans_id]=[]
            intron_by_transcript[intron.trans_id].append((intron.start,intron))
    file_in.close()
    return(intron_by_transcript)

# Fonction qui va écrire les annotations dont nous avons besoin dans un fichier :
def print_data(intron_by_transcript,len_of_interest_for_exon = 20,len_of_interest_for_intron = 30, minimal_len_for_GC = 40):
    total_len_of_interest = len_of_interest_for_exon+len_of_interest_for_intron
    for trans_id, intron_list in intron_by_transcript.items():
        # Première annotation : on calcule le nombre d'introns pour le transcrit
        numbers_of_introns = len(intron_list)
        # On tri la liste pour obtenir la position exacte de l'intron par rapport aux autres dans le transcrit :
        intron_list.sort(key=lambda x: int(x[0]))
        # Et on reverse la liste dans le cas ou la séquence se situe sur le brin -
        if intron_list[0][1].strand == "-":
            intron_list.reverse()
        # On fait le parcours de la liste des introns pour le transcrit :
        for intron_with_coords in intron_list:
            intron = intron_with_coords[1]
            # Deuxième annotation : on calcule la taille de la séquence de l'intron, en enlevant les nucléotides appartenant aux exons
            intron_len = intron.seq_without_exon(len_of_interest_for_exon)
            # Troisième annotation : on détermine la position de l'intron par rapport aux autres dans le transcrit
            intron_pos = intron_list.index(intron)+1
            # Quatrième annotation : on récupère les séquences d'intérêt de l'intron à gauche et à droite :
            intron_seq_left = intron.seq_begin(total_len_of_interest)
            intron_seq_right = intron.seq_end(total_len_of_interest)
            # Cinquième annotation : on récupère le GC rate
            intron_GC = intron.GCrate(total_len_of_interest,minimal_len_for_GC)
            line = intron.id+"\t"+intron.trans_id+"\t"+intron.gene_id+"\t"+intron.coords+"\t"+intron_seq_left+"\t"+intron_seq_right+"\t"+intron_len+"\t"+intron_pos+"\t"+numbers_of_introns+"\t"+intron_GC+"\n"
            print(line)
# Fonctions secondaires
# Formate les coordonnées pour retrouver la séquence fasta (le nom de la séquence est sous ce format)
def formating_coord(full_coordinates):
    regex = re.compile("(chr[^\+-]+)[\+-]([0-9]+):([0-9]+)")
    result = regex.findall(full_coordinates)
    intron_chr = result[0][0]
    # Les positions doivent être données en base 0 et pour la fin en non inclusif, c'est à dire la première position après la séquence, on doit donc enlever 1 au start pour passer en base 0 mais on ne touche pas à la position de fin
    intron_start = int(result[0][1])-1
    intron_end = result[0][2]
    return intron_chr+":"+str(intron_start)+"-"+str(intron_end)

# Interaction utilisateur
opts, args = getopt.getopt(sys.argv[1:],'',['liste_intron=','fasta=','exon_len=','intron_len=','len_GC=',"len_seq=","gen_ref="])
for elmts in opts:
    if elmts[0] == '--liste_intron':
        liste_intron = elmts[1] # Fichier contenant les coordonnées pour chaque intron : list_intron.tab

    elif elmts[0] == '--fasta':
        fasta_file = elmts[1] # Fichier qui contiendra les séquences fasta des introns et des exons flanquants

    elif elmts[0] == '--exon_len':
        nb_nt_exon = int(elmts[1]) # détermine le nombre de nucleotides à prendre en compte dans les exons flanquants

    elif elmts[0] == '--intron_len':
        nb_nt_intron = int(elmts[1]) # détermine le nombre de nucleotides à prendre en compte aux extrémités de l'intron

    elif elmts[0] == '--len_GC':
        nb_GC = int(elmts[1]) # détermine la longueur minimale de l'intron -(nb_nt_exon+nb_nt_intron) pour calculer le taux de GC

    elif elmts[0] == '--len_seq':
        len_seq = int(elmts[1]) # détermine la taille minimale de la séquence pour annoter

    elif elmts[0] == '--gen_ref':
        gen_ref = elmts[1] # On récupère le chemin d'acces du génome de référence en format binaire hg19.2bit


file_out_name = BED_generation(liste_intron)
command = "twobitToFa -bedPos -bed="+file_out_name+" "+gen_ref+" "+fasta_file
####Execution de la commande en ligne de commande
os.system(command)
####On suprimme le fichier temporaire .bedPos généré
os.unlink(file_out_name)
seq_info = get_fasta_seq(fasta_file)
intron_by_transcript = get_all_intron_for_transcript(liste_intron,seq_info,len_seq)
