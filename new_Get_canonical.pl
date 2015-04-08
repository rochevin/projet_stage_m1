use Bio::EnsEMBL::Registry;
use Bio::SeqIO;
use Data::Dumper;
# On possède une liste de transcrits Ensembl associés à un id gene et un id de transcrit Refseq, tous les transcrits RefSeq présents sont annotés dans notre jeu de données
# Contenant les introns de la publication de Braunch, on cherche parmis ces IDs à identifier le transcrit canonique

$ENV{MYSQL_UNIX_PORT} = "/Applications/MAMP/tmp/mysql/mysql.sock";

# Ouverture du fichier d'annotation Refseq Gene Ensembl
unless ( open(file_list, $ARGV[0]) ) {
    print STDERR "Impossible de trouver $ARGV[0] ...\n\n";
    exit;
}

my %transcript_list;

foreach my$line (<file_list>) {
	if ( $line =~ /^\s*$/ ){
    	next;
    }
    else {
    	my @line_content = split(" ", $line);
    	$transcript_list{$line_content[1]}{$line_content[2]}=$line_content[0];
    }
}
close file_list;

# Ouverture du fichier d'annotation des introns Braunch
unless ( open(annotation_braunch, $ARGV[1]) ) {
    print STDERR "Impossible de trouver $ARGV[1] ...\n\n";
    exit;
}

my %annotation_braunch;

foreach my$line (<annotation_braunch>) {
	if ( $line =~ /^\s*$/ ){
    	next;
    }
    else {
    	my @line_content = split("\t", $line);
    	my $intron_braunch = $line_content[0];
    	my ($transcript_braunch) = ($intron_braunch =~ /^[^:]+:.+:(.+):[0-9]+/);
    	my $intron_coord = $line_content[1];
    	$annotation_braunch{$transcript_braunch}{$intron_braunch}=$intron_coord;
    }
}
close annotation_braunch;


my $r = "Bio::EnsEMBL::Registry";

$r->load_registry_from_db(-host => "localhost", -user => "root", -pass => "root", -port => "8889", -verbose => "0");


my $gene_adaptor=$r->get_adaptor("Human" ,"core", "Gene");

foreach my$id (keys %transcript_list) {




	my $gene = $gene_adaptor->fetch_by_stable_id($id);

	my $canonical_transcript=$gene->canonical_transcript();
	$canonical_id = $canonical_transcript->stable_id();
	my @list_transcripts = @{ $gene->get_all_Transcripts };
	
	foreach my $one_transcript (@list_transcripts){
		my $transcript_id = $one_transcript->stable_id();
		next unless (exists $transcript_list{$id}{$transcript_id});
		# On récupère nos introns Braun correspondant au transcrit refseq/Ensembl associé
		my %intron_liste_by_refseq_id = %{$annotation_braunch{$transcript_list{$id}{$transcript_id}}};

	}

}



sub get_all_exon_for_transcript {
	my($transcript) = @_;
	my @exons = @{ $transcript->get_all_Exons };	
	my $exon_number = 1;
	my @content;
	while (my$exon = shift @exons) {
		my $stable_id  = $exon->stable_id();
    	my $seq_region = $exon->seq_region_name();
    	my $start      = $exon->start();
    	my $end        = $exon->end();
    	my $strand     = $exon->strand();
    	my $content = sprintf ("%s\t%s\t%s\t%s\t%d\tchr%s\t%+d\t%s_%s\t%s\n",$refseq,$gene_id,$canonical_transcript->stable_id(),$stable_id,$exon_number,$seq_region,$strand,$start,$end,$canonical);
		push (@content, $content);
		$exon_number++;
	}
	return @content;
}

