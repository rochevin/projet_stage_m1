####################################################################
#Auteur : ROCHER Vincent                                           #
#                                                                  #
#But : Convertir le taux d'expression et la PIR                    #
#      -Avec le bon identifiant (Ensembl)                          #
#      -Au format BED : chr debut   fin id                         #
####################################################################

# Importation des modules
import sys, getopt # Sert à récupérer les noms de fichiers en arguments
import re

# Définition des classes :

class IntronInfo(object):
    "Classe contenant les différentes annotations de l'intron, soit son PIR en fonction du tissu et son taux d'expression"
    def __init__(self,ens_id,braunch_id,gene_id,trans_id,coords,tissues,PIR,expression):
        self.id = ens_id
        self.braunch_id = braunch_id
        self.gene_id = gene_id
        self.trans_id = trans_id
        self.coords = coords
        self.tissues = tissues
        self.PIR = PIR
        self.expression = expression

        regex = re.compile("(chr[^\+-]+)([\+-])([0-9]+):([0-9]+)")
        result = regex.findall(self.coords)

        self.chr = result[0][0]
        self.strand = result[0][1]
        self.start = result[0][2]
        self.end = result[0][3]


# Fonctions :
# Fonction qui recupère le Pourcentage d'Introns Retenus pour chaque introns
def get_PIR_by_braunch_id(file_name):
    Pir_annotation = {}
    file_in = open(file_name,"r")
    head = file_in.readline()
    tissues = [elmt.replace('\n','') for elmt in head.split('\t')[1:]]
    for line in file_in:
        content = line.split('\t')
        intron_id = content[0]
        PIR = content[1:]
        PIR_by_tissue = {tissues[i]:PIR[i].replace('\n','') for i in range(0,len(PIR),1)}
        Pir_annotation[intron_id]=PIR_by_tissue
    # On créer un id NA qui ne contient que des NA
    PIR_by_tissue = {tissues[i]:"NA" for i in range(0,len(tissues),1)}
    Pir_annotation["NA"]=PIR_by_tissue
    file_in.close()
    return(Pir_annotation,tissues)

# Fonction qui recupère le taux d'expression pour chaque introns
def get_expression_by_braunch_id(file_name):
    exp_annotation = {}
    file_in = open(file_name,"r")
    head = file_in.readline()
    tissues = [elmt.replace('\n','') for elmt in head.split('\t')[3:]]
    for line in file_in:
        content = line.split('\t')
        intron_id = content[0]
        exp = content[3:]
        exp_by_tissue = {tissues[i]:exp[i].replace('\n','') for i in range(0,len(exp),1)}
        exp_annotation[intron_id]=exp_by_tissue
    # On créer un id NA qui ne contient que des NA
    exp_by_tissue = {tissues[i]:"NA" for i in range(0,len(tissues),1)}
    exp_annotation["NA"]=exp_by_tissue
    file_in.close()
    return(exp_annotation)

# Fonction qui va récupérer l'association entre les id refseq/braunch et Ensembl
def get_association_between_intron_id(file_name,dico_expression,dico_pir,tissues):
    association_ids = {}
    file_in = open(file_name,"r")
    head = file_in.readline()

    for line in file_in:
        # On récupère les annotations avec un split
        content = line.split('\t')
        intron_id_ens = content[0]
        intron_id_braunch = content[1]
        transcript_id = content[2]
        gene_id = content[3]
        coordinates = content[4]
        # On récupère la PIR et le taux d'expression relatif à l'id braunch
        if intron_id_braunch in dico_expression:
            list_expression = dico_expression[intron_id_braunch]
        else:
            list_expression = dico_expression["NA"]

        if intron_id_braunch in dico_pir:
            list_PIR = dico_pir[intron_id_braunch]
        else:
            list_PIR = dico_pir["NA"]
        association_ids[intron_id_ens] = IntronInfo(ens_id=intron_id_ens,braunch_id=intron_id_braunch,gene_id=gene_id,trans_id=transcript_id,coords=coordinates,tissues=tissues,PIR=list_PIR,expression=list_expression)
    file_in.close()
    return(association_ids)

# Fonction d'écriture
def write_file_in_BED_format(pathway,association_ids):
    file_name_PIR = pathway+"PIR_by_intron.bed"
    file_name_expression = pathway+"expression_by_intron.bed"
    file_out_PIR = open(file_name_PIR,"w")
    file_out_expression = open(file_name_expression,"w")

    for key,value in sorted(association_ids.items()):
        chromosome = value.chr
        debut = value.start
        fin = value.end
        intron_id = value.id
        PIR = value.PIR
        expression = value.expression
        line_PIR = chromosome+"\t"+debut+"\t"+fin+"\t"+intron_id+"\t"+"\t".join([PIR[tissu] for tissu in sorted(tissues)])+"\n"
        line_expression = chromosome+"\t"+debut+"\t"+fin+"\t"+intron_id+"\t"+"\t".join([expression[tissu] for tissu in sorted(tissues)])+"\n"
        file_out_PIR.write(line_PIR)
        file_out_expression.write(line_expression)
    file_out_PIR.close()
    file_out_expression.close()

# Interaction utilisateur
opts, args = getopt.getopt(sys.argv[1:],'',['PIR_file=','exp_file=','output_directory=','association_file='])
for elmts in opts:
    if elmts[0] == '--output_directory':
        output_directory = elmts[1] # Nom du fichier de sortie

    elif elmts[0] == '--association_file':
        association_file = elmts[1] # Nom du fichier contenant l'association entre les ids

    elif elmts[0] == '--exp_file':
        exp_file = elmts[1] # Nom du fichier contenant le taux d'expression

    elif elmts[0] == '--PIR_file':
        PIR_file = elmts[1] # Nom du fichier contenant le taux de PIR

# Lancement des fonctions :
Pir_annotation,tissues = get_PIR_by_braunch_id(PIR_file)
exp_annotation = get_expression_by_braunch_id(exp_file)

association_ids = get_association_between_intron_id(file_name=association_file, dico_expression=exp_annotation, dico_pir=Pir_annotation, tissues=tissues)

write_file_in_BED_format(output_directory,association_ids)
