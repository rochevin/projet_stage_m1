use Bio::EnsEMBL::Registry;
use Bio::SeqIO;
use Data::Dumper;

$ENV{MYSQL_UNIX_PORT} = "/Applications/MAMP/tmp/mysql/mysql.sock";


# Ouverture du fichier des transcrits canonique ensembl : transcript_list1050.tab
unless ( open(file_list, $ARGV[0]) ) {
    print STDERR "Impossible de trouver $ARGV[0] ...\n\n";
    exit;
}

my %transcript_list;

foreach my$line (<file_list>) {
	if ( $line =~ /^\s*$/ ){
    	next;
    }
    elsif ( $line =~ /^Ensembl/ ) {
    	next;
    }
    else {
    	my @line_content = split("\t", $line);
    	$transcript_list{$line_content[0]}=$line_content[1];
    }
}
close file_list;

my $fasta_file_out = ">".$ARGV[1];
my $out = Bio::SeqIO->new(-file => $fasta_file_out, -format => "Fasta");
unless ( open(file_out, ">".$ARGV[2]) ) {
    print STDERR "Impossible de trouver $ARGV[0] ...\n\n";
    exit;
}

my $r = "Bio::EnsEMBL::Registry";

$r->load_registry_from_db(-host => "localhost", -user => "root", -pass => "root", -port => "8889", -verbose => "0");

my $transcript_adaptor = $r->get_adaptor( 'Human', 'Core','Transcript' );

foreach my$transcript_id (keys %transcript_list) {

	#On récupère tout d'abbord nos identifiants de chaque exon, pour pouvoir associer les exons CDS qui sortiront à nos données
	# my %coord_for_each_exon = reverse %{$exon_by_transcript{$transcript}} if exists($exon_by_transcript{$transcript});
	

	my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);

	my $seq_region = $transcript->seq_region_name();
	my $strand;
	if ($transcript->strand() < 1) {$strand = "-";}
	else {$strand = "+";}

	my @exons = @{ $transcript->get_all_translateable_Exons() };

	while (my $exon = shift @exons) {
		my $coords_for_one_exon = $exon->start().":".$exon->end();
		my $seq = $exon->seq();
		$out->write_seq($seq);
		print file_out $exon->stable_id()."\t".$transcript_id."\t".$transcript_list{$transcript_id}."\t"."chr".$seq_region.$strand.$coords_for_one_exon."\n";
	}

}

$out->close();
close file_out;