#Parcours de tous les CDS #####
intron_by_transcripts = {} # Dictionnaire qui contiendra la liste des identifiant des introns dans chaque transcrit
number_CDS_no_stop = 0 # On définit un compteur pour chaque fois ou un transcrit n'a pas de codon stop -> assemblage incomplet 
for Transcript_id,CDS in dataset_CDS.items():
	stop_codon_control = 0 # Control pour savoir si le transcrit possède un codon stop

	CDS.sort(key=lambda x: x.split(':')[2]) # Fonction de tri qui met les exons du CDS dans l'ordre en fonction des coordonnées
	#Assemblage de la séquence du transcrit
	#On prend le premier exon du transcrit et on vérifie son signe :
	if CDS_content[CDS[0]].strand == "-": # Si c'est un signe "-", on veut la liste dans l'ordre décroissant
		CDS.reverse()
	sequence_CDS = "" # On définit une chaine de caractères vide qui contiendra la séquence du CDS
	first_exon = CDS_content[CDS[0]]
	début_sequence = int(first_exon.start) # On définit une variable qui sera la position 0 du CDS -> soit les coordonnées de début du premier exon
	for elmt in CDS:
		exon = CDS_content[elmt]
		#Vérification si transcrit possède un codon stop
		if exon.feature_type == 'stop_codon':
			stop_codon_control += 1 #Si oui on met +1 au control
		sequence_CDS += exon.seq # On ajoute la séquence de l'exon à notre CDS total
	if stop_codon_control == 0: # Si le control vaut 0, c'est qu'on a pas de codon stop dans la séquence
		number_CDS_no_stop +=1 # On ajoute donc +1 au compteur de transcrits sans codon stop


print('#############')
print(Transcript_id,":",CDS_content[CDS[0]].start,"-",CDS_content[CDS[-1]].stop,sep='')
print(sequence_CDS)