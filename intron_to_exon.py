#Fonction qui va associer une coordonnée d'exon à une liste de transcrits ensembl -> sert a up la rapidité du parcours
def exon_coord_by_transcript(exon_content):
	coord_to_transcript = {} # Contient une coordonnée et un identifiant Ensembl se rapportant à ces coordonnées d'exons
	for transcript,exon_list in exon_content.items():
		for elmt in exon_list:
			coords = elmt.formating_coord()
			if coords in coord_to_transcript:
				coord_to_transcript[coords].append(transcript)
			else:
				coord_to_transcript[coords]=[]
				coord_to_transcript[coords].append(transcript)
	return(coord_to_transcript)


def retrieve_transcripts_with_exon_coord(exon_coord_in_braunch_transcript,exon_coord_in_ensembl_transcript,coord_to_transcript):
	canonical_with_exon = {}
	for B_transcript,B_intron_list in exon_coord_in_braunch_transcript.items():
		liste = [] # Liste qui va servir à augmenter la rapidité du parcours en cherchant les identifiants ensembl associés à nos coordonnées
		for elmt in B_intron_list:
			if elmt in coord_to_transcript:
				resultat = coord_to_transcript[elmt]
				liste.extend(resultat)
		#On suprimme les doublons
		liste = list(set(liste))
		for E_transcript in liste:
			E_intron_list = exon_coord_in_ensembl_transcript[E_transcript]
			if B_intron_list == E_intron_list:
				canonical_with_exon[E_transcript]=B_transcript
			else:
				if len(E_intron_list) == len(B_intron_list):
					sE = set(E_intron_list)
					sB = set(B_intron_list)
					resultE = [x for x in E_intron_list if x not in sB]
					positionE = [(E_intron_list.index(x)+1) for x in E_intron_list if x not in sB]
					resultB = [x for x in B_intron_list if x not in sE]
					positionB = [(B_intron_list.index(x)+1) for x in B_intron_list if x not in sE]

					listeE=[]
					for elmt in resultE:
						position = resultE.index(elmt)
						jointure = (elmt,positionE[position])
						listeE.append(jointure)
					listeB=[]
					for elmt in resultB:
						position = resultB.index(elmt)
						jointure = (elmt,positionB[position])
						listeB.append(jointure)
					if (len(resultE) != len(E_intron_list) and len(resultB) != len(B_intron_list)):
						print(E_transcript," : ",len(resultE),'/',len(B_intron_list),sep='')
						print(listeE)
						print(B_transcript," : ",len(resultB),'/',len(B_intron_list),sep='')
						print(listeB)
					print("##################################################")
	return(canonical_with_exon)


