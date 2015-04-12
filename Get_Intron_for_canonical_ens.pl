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
unless ( open(exon_file_output, ">".$ARGV[1]) ) {
    print STDERR "Impossible de trouver $file_accession_output ...\n\n";
    exit;
}
print exon_file_output "Intron_id\tGene_id\tTrans_id\tCoordinates\n";
my $r = "Bio::EnsEMBL::Registry";

$r->load_registry_from_db(-host => "localhost", -user => "root", -pass => "root", -port => "8889", -verbose => "0");


my $gene_adaptor=$r->get_adaptor("Human" ,"core", "Gene");

foreach my$id (keys %gene_list) {

	my $gene = $gene_adaptor->fetch_by_stable_id($id);

	my $canonical_transcript=$gene->canonical_transcript();
	my $canonical_id = $canonical_transcript->stable_id();
	my $seq_region = $canonical_transcript->seq_region_name();
	my $strand;
	if ($canonical_transcript->strand() < 1) {$strand = "-";}
	else {$strand = "+";}


	my @exons = @{ $canonical_transcript->get_all_Exons };
	my @exons_coords;
	while (my$exon = shift @exons) {
		my $coords_for_one_exon = $exon->start().":".$exon->end();
		push (@exons_coords,$coords_for_one_exon);
	}
	@exons_coords = sort @exons_coords;
	my @Ensembl_intron_coord;
	for(my$pos = 0; $pos<scalar(@exons_coords);$pos++) {
		if ($pos+1 != scalar(@exons_coords)){
			my($intron_pos_left) = $exons_coords[$pos] =~ /^[0-9]+:([0-9]+)/;
			my($intron_pos_right) = $exons_coords[$pos+1] =~ /^([0-9]+):[0-9]+/;
			my $intron_coord = $intron_pos_left."_".$intron_pos_right;
			push(@Ensembl_intron_coord, $intron_coord);
		} 
	}
	@Ensembl_intron_coord = reverse(@Ensembl_intron_coord) if ($canonical_transcript->strand() <1);
	my @intron_line;
	my $intron_num = 1;
	foreach my $intron (@Ensembl_intron_coord) {
		my $intron_id= $id.":".$canonical_id.":".$intron_num;
		my $line = $intron_id."\t".$id."\t".$canonical_id."\t"."chr".$seq_region.$strand.$intron."\n";
		push(@intron_line,$line);
		$intron_num++;
	}
	@intron_line = reverse @intron_line if ($canonical_transcript->strand() <1);
	print exon_file_output $_ for @intron_line;
	
}