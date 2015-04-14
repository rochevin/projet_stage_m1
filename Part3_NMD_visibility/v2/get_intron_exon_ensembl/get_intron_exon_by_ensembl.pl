# Programme d'utilisation de l'API Ensembl
# Pour chaque transcrit canonique d'Ensembl, va récupérer tous les introns et les exons(codants) ainsi que leur séquences
# Va créer trois fichiers : -fichier fasta contenant chaque séquence d'intron et d'exon
# 							-fichier contenant la liste des exons codants pour chaque transcrit
#							-fichier contenant la liste des introns pour chaque transcrit
#Utilise Gene_ens_75 qui contient toute la liste des gènes codants d'ensembl



use Bio::EnsEMBL::Registry;
use Bio::SeqIO;
use Data::Dumper;

# Définition du socket MYSQL pour que DBI utilise la bonne base de données
$ENV{MYSQL_UNIX_PORT} = "/Applications/MAMP/tmp/mysql/mysql.sock";


# Ouverture de la liste des gènes codant Ensembl : Gene_ens_75
unless ( open(file_list, $ARGV[0]) ) {
    print STDERR "Impossible de trouver $ARGV[0] ...\n\n";
    exit;
}

my %gene_list;

foreach my$line (<file_list>) {
	if ( $line =~ /^\s*$/ ){
    	next;
    }
    else {
    	my @line_content = split(" ", $line);
    	$gene_list{$line_content[0]}++;
    }
}
close file_list;
# Création du fichier fasta qui va contenir la liste des séquences des introns et des exons codants
my $fasta_file_out = ">".$ARGV[1];
my $out = Bio::SeqIO->new(-file => $fasta_file_out, -format => "Fasta");

# Création du fichier de sortie pour les exons codants
unless ( open(file_out_exon, ">".$ARGV[2]) ) {
    print STDERR "Impossible de trouver $ARGV[0] ...\n\n";
    exit;
}
print file_out_exon "ID exon\tID transcript\tID gene\tCoordinates\n";
# Création du fichier de sortie pour les introns
unless ( open(file_out_intron, ">".$ARGV[3]) ) {
    print STDERR "Impossible de trouver $ARGV[0] ...\n\n";
    exit;
}
print file_out_intron "ID intron\tID transcript\tID gene\tCoordinates\n";
# On utilise l'API d'ensembl
my $r = "Bio::EnsEMBL::Registry";
# On se connecte à la base de donnée locale ensembl 75
$r->load_registry_from_db(-host => "localhost", -user => "root", -pass => "root", -port => "8889", -verbose => "0");
# On utilise gene_adaptor pour récupérer les informations sur nos gènes
my $gene_adaptor=$r->get_adaptor("Human" ,"core", "Gene");
# Puis on parcours notre liste de gènes :
# Pour chaque gène =
foreach my$id (keys %gene_list) {
	# On récupère l'objet gène correspondant à notre identifiant
	my $gene = $gene_adaptor->fetch_by_stable_id($id);
	# On récupère le transcrit canonique de chaque gène
	my $canonical_transcript=$gene->canonical_transcript();
	# On enregistre son identifiant
	my $canonical_id = $canonical_transcript->stable_id();
	# On détermine son chromosome
	my $seq_region = $canonical_transcript->seq_region_name();
	# Ainsi que son brin : +/-
	my $strand;
	if ($canonical_transcript->strand() < 1) {$strand = "-";}
	else {$strand = "+";}
	

	###########################################################
	############Phase de récupération des introns##############
	###########################################################
	# On récupère tous nos introns dans une liste
	my @introns = @{ $canonical_transcript->get_all_Introns() };
	# On ne veut récupérer que les transcrits qui possèdent au moins un intron
	next unless (scalar(@introns) >= 1);
	# On parcours chacun de nos introns
	while (my $intron = shift @introns) {
		# On construit l'identifiant à partir des exons flanquants
		my $intron_id = $intron->prev_Exon->display_id()."_".$intron->next_Exon->display_id();
		# On récupère ses coordonnées 
		my $coords_for_one_intron = $intron->start().":".$intron->end();
		# On récupère la séquence et on l'enregistre au format BioSeq
		my $seq = $intron->seq();
		my $fasta_seq = Bio::Seq->new(-display_id => $intron_id, -seq => $seq);
		# Puis on l'écris dans le fichier
		$out->write_seq($fasta_seq);
		# On écris l'identifiant dans notre fichier d'annotation d'intron
		print file_out_intron $intron_id."\t".$canonical_id."\t".$id."\t"."chr".$seq_region.$strand.$coords_for_one_intron."\n";
	}


	# On récupère tous nos exons CDS dans une liste 
	my @exons = @{ $canonical_transcript->get_all_translateable_Exons() };
	# On parcours chacun de nos exons
	while (my $exon = shift @exons) {
		# On récupère ses coordonnées
		my $coords_for_one_exon = $exon->start().":".$exon->end();
		# On récupère sa séquence au format BioSeq
		my $seq = $exon->seq();
		# Puis on l'écris dans le fichier
		$out->write_seq($seq);
		# On écris l'identifiant dans notre fichier d'annotation d'exon
		print file_out_exon $exon->stable_id()."\t".$canonical_id."\t".$id."\t"."chr".$seq_region.$strand.$coords_for_one_exon."\n";
	}

}

$out->close();
close file_out_exon;
close file_out_intron;