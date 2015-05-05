use Bio::EnsEMBL::Registry;
use Bio::SeqIO;


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



# On utilise l'API d'ensembl
my $r = "Bio::EnsEMBL::Registry";
# On se connecte à la base de donnée locale ensembl 75
$r->load_registry_from_db(-host => "localhost", -user => "root", -pass => "root", -port => "8889", -verbose => "0");



my $gene_adaptor=$r->get_adaptor("Human" ,"core", "Gene");
my $slice_adaptor=$r->get_adaptor("Human" ,"core", "Slice");
foreach my$id (keys %gene_list) {
	my $gene = $gene_adaptor->fetch_by_stable_id($id);

	my $canonical_transcript=$gene->canonical_transcript();

	my $seq_region = $canonical_transcript->seq_region_name();

	my @introns = @{ $canonical_transcript->get_all_Introns() };
		# On ne veut récupérer que les transcrits qui possèdent au moins un intron
		next unless (scalar(@introns) >= 1);
		# On parcours chacun de nos introns
	while (my $intron = shift @introns) {
		my $begin_seq = "";
		my $end_seq = "";
		# On détermine l'identifiant de l'intron selon notre méthode
		my $intron_id = $intron->prev_Exon->display_id()."_".$intron->next_Exon->display_id();
		# On enregistre ses coordonnées dans le chromosome
		my $debut_intron = $intron->start();
		my $fin_intron = $intron->end();
		print "Identifiant de l'intron :".$intron_id."\nCoordonnées : ".$debut_intron."-".$fin_intron."\n";
		print "Région de l'intron : ".$seq_region."\n";
		print "Coordonnées des séquences d'intérêt : \n";
		# On détermine nos séquences d'intérêts et on récupère la séquence :
		# Pour la séquence de début
		my $debut = $debut_intron-20;
		my $fin = $debut_intron+29;
		print "Séquence début : ".$debut."-".$fin."\n";
		my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $seq_region, $debut, $fin );
		if (defined $slice){
			$begin_seq = $slice->seq();
		}else {
			($seq_exon_left,$seq_exon_right) = get_seq_with_exon($intron);
			$begin_seq = $seq_exon_left.substr($intron->seq(),0,30);
		}
		# On définit un objet Bio::Seq pour la séquence
		my $begin_fasta_seq = Bio::Seq->new(-display_id => $intron_id, -seq => $begin_seq,-description => "first");
		# Puis on l'écris dans le fichier
		$out->write_seq($begin_fasta_seq);
		# Pour la séquence de fin 
		my $debut = $fin_intron-29;
		my $fin = $fin_intron+20;
		print "Séquence fin : ".$debut."-".$fin."\n";
		my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $seq_region, $debut, $fin );
		if (defined $slice){
			$end_seq = $slice->seq();
		}else {
			($seq_exon_left,$seq_exon_right) = get_seq_with_exon($intron);
			$end_seq = substr($intron->seq(),length($intron->seq())-30).$seq_exon_right;
		}
		# On définit un objet Bio::Seq pour la séquence
		my $end_fasta_seq = Bio::Seq->new(-display_id => $intron_id, -seq => $end_seq,-description => "last");
		# Puis on l'écris dans le fichier
		$out->write_seq($end_fasta_seq);
	}
}
$out->close();

sub get_seq_with_exon {
	(my$intron) = @_;
	$exon_left_seq = $intron->prev_Exon->seq->seq;
	$exon_right_seq = $intron->next_Exon->seq->seq;
	$sub_seq_left = substr($exon_left_seq,length($exon_left_seq)-20);
	$sub_seq_right = substr($exon_right_seq,0,20);
	print "Slice marche pas, donc on test les séquences des exons : \n";
	print $sub_seq_left."\n";
	print $sub_seq_right."\n";
	return($sub_seq_left,$sub_seq_right);
}