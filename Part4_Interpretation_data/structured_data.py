####################################################################
#Auteur : ROCHER Vincent                                           #
#                                                                  #
#But : Utiliser une structure de données pratique en forme d'objet #
#      Les annotations seront regroupés en objet par fichier       #
#      Puis regroupés en une seule structure globale               #
####################################################################

# Définition des classes
class PrimaryItem(object):
    "Classe contenant les annotations issues du script primary_annotation.py"
    def __init__(self,intron_id,trans_id,gene_id,coordinates,seq_left,seq_right,intron_len,intron_pos,total_intron,GCrate):
    self.id = intron_id
    self.trans_id = trans_id
    self.gene_id = gene_id
    self.coords = coordinates
    self.seq_left = seq_left
    self.seq_right = seq_right
    self.len = intron_len
    self.pos = intron_pos
    self.total_intron = total_intron
    self.GCrate = GCrate


class PIRInfo(object):
    "Classe contenant la PIR pour chaque intron dans chaque tissu"
    def __init__(self,intron_id,PIR_dictionnary):
        self.id = intron_id
        self.PIR_by_tissue = PIR_dictionnary

class NMDInfo(object):
    "Classe contenant les annotations relatives au système NMD"
    def __init__(self,intron_id,CDS_status,dist_last_intron,dist_next_stop,dist_CDS_stop,PTC_status,position,phase):
        self.id = intron_id
        self.CDS_status = CDS_status
        self.dist_last_intron = dist_last_intron
        self.dist_next_stop = dist_next_stop
        self.dist_CDS_stop = dist_CDS_stop
        self.PTC_status = PTC_status
        self.position = position
        self.phase = phase


# Définition des fonctions
# Fonction qui récupère les annotations tel que taux de GC, seq_left, seq_right
def get_primary_annotation(file_name):
    primary_annotation = {}
    file_in = open(file_name,"r")
    head = file_in.readline()

    for line in file_in:
        content = line.split('\t')
        primary_intron = PrimaryItem(intron_id=content[0],trans_id=content[1],gene_id=content[2],coordinates=content[3],seq_left=content[4],seq_right=content[5],intron_len=content[6],intron_pos=content[7],total_intron=content[8],GCrate=content[9].replace('\n',''))
        primary_annotation[primary_intron.id]=primary_intron
    file_in.close()
    return(primary_annotation)

# Fonction qui va récupérer l'association entre les id refseq/braunch et Ensembl
def get_association_between_intron_id(file_name):
    association_ids = {}
    file_in = open(file_name,"r")
    head = file_in.readline()

    for line in file_in:
        content = line.split('\t')
        intron_id_ens = content[0]
        intron_id_braunch = contenr[1]
        association_ids[intron_id_ens]=intron_id_braunch
    file_in.close()
    return(association_ids)

# Fonction qui recupère le Pourcentage d'Introns Retenus pour chaque introns
def get_PIR_by_braunch_id(file_name):
    Pir_annotation = {}
    file_in = open(file_name,"r")
    head = file_in.readline()
    tissues = head.split('\t')[1:]
    for line in file_in:
        content = line.split('\t')
        intron_id = content[0]
        PIR = content[1:]
        PIR_by_tissue = {tissues[i]:PIR[i] for i in range(0,len(PIR),1)}
        intron_pir = PIRInfo(intron_id=intron_id,PIR_dictionnary=PIR_by_tissue)
        Pir_annotation[intron_pir.id]=intron_pir
    file_in.close()
    return(Pir_annotation)

# Fonction qui va récupérer les annotations sur les introns relatif à la NMD visibilité
def get_annotation_for_NMD(file_name):
    NMD_annotation = {}

    file_in = open(file_name,"r")
    head = file_in.readline()

    for line in file_in:
        content = line.split('\t')
        intron_NMD = NMDInfo(intron_id=content[3],CDS_status=content[4],dist_last_intron=content[5],dist_next_stop=content[6],dist_CDS_stop=content[7],PTC_status=content[8],position=content[9],phase=content[10])
        NMD_annotation[intron_NMD.id]=intron_NMD
    file_in.close()
    return(NMD_annotation)

def get_windows_for_
