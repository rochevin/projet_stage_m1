use Bio::EnsEMBL::Registry;
use Bio::SeqIO;


$ENV{MYSQL_UNIX_PORT} = "/Applications/MAMP/tmp/mysql/mysql.sock";

unless ( open(file_list, $ARGV[0]) ) {
    print STDERR "Impossible de trouver $file_list ...\n\n";
    exit;
}

my %transcript_list;

foreach my$line (<file_list>) {
	if ( $line =~ /^\s*$/ ){
    	next;
    }
    else {
    	my @line_content = split(" ", $line);
    	$transcript_list{$line_content[0]}{"Gene_ID"}=$line_content[1];
    	$transcript_list{$line_content[0]}{"Trans_ID"}=$line_content[2];
    }
}
close file_list;

# Fichier qui va servir à annoter les exons de chaques transcrits canoniques, pour pouvoir comparer avec nos transcrits, si ils ont autant d'introns
unless ( open(exon_file_output, ">".$ARGV[1]) ) {
    print STDERR "Impossible de trouver $file_accession_output ...\n\n";
    exit;
}
print exon_file_output "RefSeq\tEnsembl Transcript ID\tEnsembl Gene ID\tEnsembl Exon ID\tExon number\tChr\tStrand\tCoordinates\n";

my $fasta_file_out = ">".$ARGV[2];
my $out = Bio::SeqIO->new(-file => $fasta_file_out, -format => "Fasta");

unless ( open(cds_file_output, ">".$ARGV[3]) ) {
    print STDERR "Impossible de trouver $file_accession_output ...\n\n";
    exit;
}


my $r = "Bio::EnsEMBL::Registry";

$r->load_registry_from_db(-host => "localhost", -user => "root", -pass => "root", -port => "8889", -verbose => "0");


my $gene_adaptor=$r->get_adaptor("Human" ,"core", "Gene");

foreach my$id (keys %transcript_list) {
	my $gene = $gene_adaptor->fetch_by_stable_id($transcript_list{$id}{"Gene_ID"});

	my $canonical_transcript=$gene->canonical_transcript();
	my $transcript = $canonical_transcript->stable_id();

	print $gene->stable_id()." => ".$canonical_transcript->stable_id()."\n";
	my @exons = @{ $canonical_transcript->get_all_Exons };
	
	my $exon_number = 1;
	while (my$exon = shift @exons) {
		my $stable_id  = $exon->stable_id();
    	my $seq_region = $exon->seq_region_name();
    	my $start      = $exon->start();
    	my $end        = $exon->end();
    	my $strand     = $exon->strand();
    	my $content = sprintf ("%s\t%s\t%s\t%s\t%d\tchr%s\t%+d\t%s_%s\n",$id,$transcript,$transcript_list{$id}{"Gene_ID"},$stable_id,$exon_number,$seq_region,$strand,$start,$end);
    	print exon_file_output $content;
		$exon_number++;
	}

	my @translateable_exons = @{ $canonical_transcript->get_all_translateable_Exons };
	while (my$exon = shift @translateable_exons) {

		my $stable_id  = $exon->stable_id();
    	my $seq_region = $exon->seq_region_name();
    	my $start      = $exon->start();
    	my $end        = $exon->end();
    	my $strand     = $exon->strand();
    	my $content = sprintf ("%s\t%s\t%s\t%s\t%d\tchr%s\t%+d\t%s_%s\n",$id,$transcript,$transcript_list{$id}{"Gene_ID"},$stable_id,$exon_number,$seq_region,$strand,$start,$end);
    	print cds_file_output $content;


		$exonseqobj = $exon->seq();
		$exonseqobj->description($id.":".$transcript_list{$id}{"Gene_ID"});
		$out->write_seq($exonseqobj);
	}
	my @introns = @{ $canonical_transcript->get_all_Introns };
	while (my$intron = shift @introns) {
		$tag = $intron->prev_Exon->stable_id()."_".$intron->next_Exon->stable_id();
		$intronseqobj = Bio::Seq->new( -display_id => $tag,
                             -seq => $intron->seq(),
                             -description => $id.":".$transcript_list{$id}{"Gene_ID"});
		$out->write_seq($intronseqobj);
	}

}
$out->close();
close exon_file_output;
close cds_file_output;