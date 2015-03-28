import gffutils
import sys

class CDSInfo(object):
	"Classe contenant les données de chaque CDS (exon sans 5'-UTR/3' UTR)"
	def __init__(self,cds_chr,cds_start,cds_stop,cds_strand,cds_transcript,cds_exon_number,gene_id):

		self.id = cds_transcript+":"+cds_exon_number
		self.chr = cds_chr
		self.start = cds_start
		self.stop = cds_stop
		self.strand = cds_strand
		self.transcript = cds_transcript
		self.exon_number = cds_exon_number

	def formating_coord(self):
		return self.chr+":"+str(self.start)+"-"+str(self.stop)

	def BED_coord(self):
		return "chr"+self.chr+"\t"+str(self.start)+"\t"+str(self.stop)+"\t"+self.id


def extract_ensembl_id(file_name):
	IDlist = {} # Dictionnaire contenant en clé l'id de l'intron et en valeur l'id Ensembl
	list_Ensembl_ids = {} # Dictionnaire contenant tous les ids Ensembl en clé pour éviter la redondance
	file_in=open(file_name,"r")
	header = file_in.readline() # On enregistre l'en-tête
	for line in file_in: #On parcours chaque ligne du fichier à partir de la deuxième
		content=line.split("\t")
		IDlist[content[0]] = content[1].replace('\n', '')
		list_Ensembl_ids[content[1].replace('\n', '')] ="None"
	file_in.close()
	return(IDlist,list_Ensembl_ids)



# db = gffutils.create_db(sys.argv[1], dbfn=sys.argv[1]+'.db', force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
# print('database ok')
IDlist,list_Ensembl_ids = extract_ensembl_id('data/DataBraunschweig/Expression/Human_EnsemblToIntron.tab')

db = gffutils.FeatureDB(sys.argv[1]+'.db', keep_order=True)

liste_transcripts = {} # On définit un dico en dehors de la boucle qui contiendra tous les transcrits du gènes en clé et tous les CDS du transcrit en valeur
count = 0 # On définit un compteur pour chaque fois où un id ENS provenant de notre jeu de donnée match dans le GTF
for id_ENS in list_Ensembl_ids:
	try:
		g = db[id_ENS]
		count += 1
	except:
		if id_ENS != "NA":
			message = "Impossible de trouver "+id_ENS+" dans le fichier GTF"
			print(message)
		pass
	
	for i in db.children(g, featuretype='CDS'):
		if i.source == 'protein_coding':
			id_CDS = i.attributes['exon_number'][0]+":"+i.attributes['transcript_id'][0]
			id_CDS = CDSInfo(cds_chr=i.seqid,cds_start=i.start,cds_stop=i.stop,cds_strand=i.strand,cds_transcript=i.attributes['transcript_id'][0],cds_exon_number=i.attributes['exon_number'][0],gene_id=id_ENS)
			if i.attributes['transcript_id'][0] in liste_transcripts:
				liste_transcripts[i.attributes['transcript_id'][0]].append(id_CDS)
			else:
				liste_transcripts[i.attributes['transcript_id'][0]]=[]
				liste_transcripts[i.attributes['transcript_id'][0]].append(id_CDS)



for key,value in liste_transcripts.items():
	print('##########################')
	print(key)
	[print(elmt.id,elmt.formating_coord(),sep=" : ") for elmt in value]
print(count,"/",len(list_Ensembl_ids)-1,"matchs")
		
		





