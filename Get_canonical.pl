use Bio::EnsEMBL::Registry;
use Bio::SeqIO;
use Data::Dumper;
# On possède une liste de transcrits Ensembl associés à un id gene et un id de transcrit Refseq, tous les transcrits RefSeq présents sont annotés dans notre jeu de données
# Contenant les introns de la publication de Braunch, on cherche parmis ces IDs à identifier le transcrit canonique

$ENV{MYSQL_UNIX_PORT} = "/Applications/MAMP/tmp/mysql/mysql.sock";

# Ouverture du fichier d'annotation Refseq Gene Ensembl
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
    	$transcript_list{$line_content[1]}{$line_content[2]}=$line_content[0];
    }
}
close file_list;
# Fichier qui va servir à annoter les exons de chaques transcrits canoniques, pour pouvoir comparer avec nos transcrits, si ils ont autant d'introns
unless ( open(exon_file_output, ">".$ARGV[1]) ) {
    print STDERR "Impossible de trouver $file_accession_output ...\n\n";
    exit;
}
print exon_file_output "RefSeq\tEnsembl Transcript ID\tEnsembl Gene ID\tEnsembl Exon ID\tExon number\tChr\tStrand\tCoordinates\tCanonical\n";

unless ( open(exon_no_match, ">".$ARGV[2]) ) {
    print STDERR "Impossible de trouver $file_accession_output ...\n\n";
    exit;
}

# my $fasta_file_out = ">".$ARGV[2];
# my $out = Bio::SeqIO->new(-file => $fasta_file_out, -format => "Fasta");

# unless ( open(cds_file_output, ">".$ARGV[3]) ) {
#     print STDERR "Impossible de trouver $file_accession_output ...\n\n";
#     exit;
# }


my $r = "Bio::EnsEMBL::Registry";

$r->load_registry_from_db(-host => "localhost", -user => "root", -pass => "root", -port => "8889", -verbose => "0");


my $gene_adaptor=$r->get_adaptor("Human" ,"core", "Gene");

foreach my$id (keys %transcript_list) {
	my $gene = $gene_adaptor->fetch_by_stable_id($id);

	my $canonical_transcript=$gene->canonical_transcript();
	my $transcript = $canonical_transcript->stable_id();
	if (exists $transcript_list{$id}{$transcript}) {
		# print $gene->stable_id()." => ".$canonical_transcript->stable_id()."\n";
		my @content = get_all_exon_for_transcript($canonical_transcript,$transcript_list{$id}{$transcript},$id,"YES");
		print exon_file_output $_ for @content;
	}
	else {
		print sprintf ("%s n'est pas dans la liste des transcrits de %s\n", $transcript,$id);
		$nb_trans_by_gene = scalar( keys %{$transcript_list{$id}});
		if ($nb_trans_by_gene <= 1) {
			my @list_transcripts = @{ $gene->get_all_Transcripts };
			foreach my$one_transcript (@list_transcripts) {
				my $transcript_id = $one_transcript->stable_id();
				next unless (exists $transcript_list{$id}{$transcript_id});
				print $id." : On a trouvé un remplaçant au canonique\n";
				my @content = get_all_exon_for_transcript($one_transcript,$transcript_list{$id}{$transcript_id},$id,"NO");
				print exon_file_output $_ for @content;
			}
		}
		else {
			print exon_no_match $id."\t".$nb_trans_by_gene."\n";
		}
	}

}

close exon_file_output;
close exon_no_match;
# close cds_file_output;
# $out->close();


sub get_all_exon_for_transcript {
	my($canonical_transcript,$refseq,$gene_id,$canonical) = @_;
	my @exons = @{ $canonical_transcript->get_all_Exons };	
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

	# my $exon_number = 1;
	# my @translateable_exons = @{ $canonical_transcript->get_all_translateable_Exons };
	# while (my$exon = shift @translateable_exons) {
	# 	my $stable_id  = $exon->stable_id();
 #    	my $seq_region = $exon->seq_region_name();
 #    	my $start      = $exon->start();
 #    	my $end        = $exon->end();
 #    	my $strand     = $exon->strand();
 #    	my $content = sprintf ("%s\t%s\t%s\t%s\t%d\tchr%s\t%+d\t%s_%s\n",$transcript_list{$id}{$transcript},$id,$transcript,$stable_id,$exon_number,$seq_region,$strand,$start,$end);
 #    	print cds_file_output $content;
	# 	$exonseqobj = $exon->seq();
	# 	$exonseqobj->description($id.":".$transcript_list{$id}{"Gene_ID"});
	# 	$out->write_seq($exonseqobj);
	# }


	# my @introns = @{ $canonical_transcript->get_all_Introns };
	# while (my$intron = shift @introns) {
	# 	$tag = $intron->prev_Exon->stable_id()."_".$intron->next_Exon->stable_id();
	# 	$intronseqobj = Bio::Seq->new( -display_id => $tag,
 #                             -seq => $intron->seq(),
 #                             -description => $id.":".$transcript_list{$id}{"Gene_ID"});
	# 	$out->write_seq($intronseqobj);
	# }





