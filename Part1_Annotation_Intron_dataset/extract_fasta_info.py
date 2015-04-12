####################################################################
#Auteur : ROCHER Vincent										   #
#																   #
#But : Extrait les séquences du fichier fasta					   #
#	   Enregistre (-20/+30) de l'intron de chaque côté			   #
####################################################################
# On importe seqIO pour parser le fichier fasta et Seq pour le brin complémentaire
from Bio import SeqIO 
from Bio.Seq import Seq


def extract_fasta_info(filename,retrieve_name): #On définit une fonction pour récupérer les informations : GC rate (-20/+30) (-30/+20)
	seq_info={} # On définit un dictionnaire qui contiendra les données
	handle = open(filename, "rU") #On ouvre le fichier en mode lecture
	for record in SeqIO.parse(handle, "fasta") : #On parse le fichier avec seqIO

		seq_id = record.id #On enregistre le nom de la sequence, ici les coordonnées
		seq = record.seq #On enregistre la sequence fasta dans une variable
		seq= seq.upper() #On met en majuscule au cas ou la séquence ne l'est pas
		value = retrieve_name[seq_id] # On récupère la liste d'identifiants intronique pour chaque position
		for intron_id in value: # On parcours cette liste
			if len(seq)>=100: #On vérifie que la séquence est bien supérieure à 100 nucléotides
				taille = len(intron_id)
				sens = intron_id[taille-1] #On détermine le signe
				if sens == "-":
					#print("brin +",seq,sep=" : ")
					seq = seq.reverse_complement()
					#print("brin -",seq,sep=" : ")
				area_begin = seq[0:50] # séquence d'intéret : (-20/+30)
				taille = len(seq)
				area_end = seq[taille - 50 :taille] # séquence d'intéret : (-30/+20)
				intron_seq = seq[50:-50] # séquence intronique pour GC rate sans les régions d'intéret
				intron_length = len(seq[20:-20])
				list_info = [seq_id, str(area_begin) , str(area_end) , intron_length ,GCrate(intron_seq)]
				#On enregistre la liste dans le dictionnaire
				seq_info[intron_id[0:-1]] = list_info # On enregistre l'identifiant comme clé sans le signe
			else:
				#Si la taille est inférieure à 100 nucléotides, on met NA
				list_info = [seq_id, "NA" , "NA" , "NA" , "NA"]
				seq_info[intron_id[0:-1]] = list_info
	#On ferme le fichier et on renvoit le dictionnaire
	handle.close()
	return(seq_info)

def GCrate(intron_seq):
	#On calcule le taux de ACTG
	A = intron_seq.count("A") 
	C = intron_seq.count("C")
	G = intron_seq.count("G")
	T = intron_seq.count("T")
	#On vérifie si on a suffisament de ATGC pour calculer un taux de GC
	if ((A+T+C+G)>=40):
		#Calcul du taux de GC
		GCrate = ((G+C)/(A+T+G+C))*100
		GCrate = round(GCrate,2)
	else:
		GCrate = "NA"
	return GCrate


#Calcul du taux des GC et des bout d'intéret :
# area_begin = seq[0:50] # séquence d'intéret : (-20/+30)
# 					area_end = seq[-51:-1] # séquence d'intéret : (-30/+20)
# 					intron_seq = seq[50:-50] # séquence intronique pour GC rate
# 					intron_length = len(seq[20:-20])
# 					#On calcule le taux de ACTG
					# A = intron_seq.count("A") 
					# C = intron_seq.count("C")
					# G = intron_seq.count("G")
					# T = intron_seq.count("T")
					#On vérifie si on a suffisament de ATGC pour calculer un taux de GC
					# if ((A+T+C+G)>=40):
					# 	#Calcul du taux de GC
					# 	GCrate = ((G+C)/(A+T+G+C))*100
					# 	GCrate = round(GCrate,2)
# 						# On enregistre sous le format liste
# 						list_info = [str(area_begin) , str(area_end) , intron_length , GCrate]
# 					#Sinon on met NA
# 					else:
# 						list_info = [str(area_begin) , str(area_end) , intron_length , "NA"]