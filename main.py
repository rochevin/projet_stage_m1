####################################################################
#Auteur : ROCHER Vincent										   #
#																   #
#But : Programme principal										   #
#	   Execute tous les autres									   #
####################################################################


from Bio import SeqIO # On importe seqIO pour parser le fichier fasta
import sys # Sert à récupérer les noms de fichiers en arguments
import re # On importe re pour expression régulière
import os # On importe os pour éxecuter des commandes terminal dans python
# On importe chaque fonction issue des autres scripts
from extract_from_introntab import extract_coord
from extract_fasta_info import extract_fasta_info
from paste_data_file import paste_data

#On lance le premier script : récupère les coordonnées et enregistre au format BED
retrieve_name = extract_coord(sys.argv[1],sys.argv[1]+".bedPos")

#On lance twobitToFa pour récupérer les séquences fasta à partir des coordonnées
os.system("twobitToFa -bedPos -bed=testintronsHuman.tab.bedPos /Users/vincentrocher/Documents/STAGE_LYON/Projet_stage_RV/data/hg19.2bit output.fa")
#On lance le deuxième script pour récupérer les séquences fasta et capturer les informations : GC rate (-20/+30) (-30/+20)
seq_info=extract_fasta_info("output.fa",retrieve_name)


# On lance le troisième script qui va coller les nouvelles données obtenues à celle de base
paste_data(sys.argv[1],"result_"+sys.argv[1],seq_info) 



