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