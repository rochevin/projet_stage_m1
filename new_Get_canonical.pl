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

my %gene_match;
my $gene_with_more_transcript = 0;
my $total_gene_with_more_transcript = 0;
print "Gènes dans notre jeu de données : ".scalar(keys %transcript_list)."\n";
my $passe_au_travers = scalar(keys %transcript_list);
my $nomatch = 0;
foreach my$id (keys %transcript_list) {

	my $gene = $gene_adaptor->fetch_by_stable_id($id);

	my $canonical_transcript=$gene->canonical_transcript();
	my $canonical_id = $canonical_transcript->stable_id();
	my @list_transcripts = @{ $gene->get_all_Transcripts };
	my $chromosome = $gene->seq_region_name();

	my @transcript_match; # Liste de transcrit pour un gène qui a les mêmes introns que notre jeu de données
	foreach my $one_transcript (@list_transcripts){
		my $transcript_id = $one_transcript->stable_id();
		next unless (exists $transcript_list{$id}{$transcript_id});
		# On récupère nos introns Braun correspondant au transcrit refseq/Ensembl associé
		# Et on formate les coordonnées des introns
		my %intron_liste_by_refseq_id = %{$annotation_braunch{$transcript_list{$id}{$transcript_id}}};
		my @Braunch_intron_coord;
		foreach $elmt (keys %intron_liste_by_refseq_id) {
			my($coord_intron) = $intron_liste_by_refseq_id{$elmt} =~ /chr[\w]+[-\+][0-9]+:([0-9]+_[0-9]+):[0-9]+/; 
			push (@Braunch_intron_coord, $coord_intron);
		}

		my @exons = @{ $one_transcript->get_all_Exons };
		my @exons_coords;
		while (my $exon = shift @exons) {
			my $coords_for_one_exon = $exon->start().":".$exon->end();
			push (@exons_coords,$coords_for_one_exon);
		}

		my @exons_coords = sort @exons_coords;
		my @Ensembl_intron_coord = intron_by_exon(@exons_coords);
		my @Ensembl_intron_coord = sort @Ensembl_intron_coord;
		my @Braunch_intron_coord = sort @Braunch_intron_coord;
		#On ajoute tous les transcrits de notre gene qui ont les mêmes introns que notre jeu de données dans la liste transcript match
		if (@Braunch_intron_coord eq @Ensembl_intron_coord) {
			push(@transcript_match,$transcript_id);
		}
	}
	#Si on a qu'un seul transcrit par gène, on l'ajoute directement
	if (scalar(@transcript_match) == 1) {
		push(@{$gene_match{$id}},$transcript_match[0]);
		$passe_au_travers += -1;
	}
	#Sinon on cherche le canonique parmis ceux là
	elsif (scalar(@transcript_match) > 1) {
		$gene_with_more_transcript++;
		$total_gene_with_more_transcript++;
		foreach my $each_transcript (@transcript_match) {
			next unless($each_transcript eq $canonical_id);
			$gene_with_more_transcript+= -1;
			$passe_au_travers += -1;
			push(@{$gene_match{$id}},$each_transcript);
		}
	}
	if (scalar(@transcript_match) == 0) {$nomatch++;}
}
print "Gènes avec transcrits identifiés : ".scalar(keys %gene_match)."\n";
print "Gènes avec plus d'un transcrit identifié : ".$total_gene_with_more_transcript."\n";
print "Gènes avec plus d'un transcrit identifié qui n'a pas le canonique dans la liste : ".$gene_with_more_transcript."\n";
print "Gènes identifiés qui n'ont aucun transcrits equivalents à notre jeu de données : ".$nomatch."\n";
my $combien_canonical = 0;
foreach my$elmt (keys %gene_match){
	my $transcript_id = @{$gene_match{$elmt}}[0];
	my $gene = $gene_adaptor->fetch_by_stable_id($elmt);
	my $canonical_transcript=$gene->canonical_transcript();
	my $canonical_id = $canonical_transcript->stable_id();
	if ($canonical_id eq $transcript_id) {$combien_canonical++;}
}
print "Gènes avec un transcrit canonique identifié : ".$combien_canonical."\n";
print "Gènes qui n'ont pas de transcrits identifiés : ".$passe_au_travers."\n";

unless ( open(exon_file_output, ">".$ARGV[2]) ) {
    print STDERR "Impossible de trouver $file_accession_output ...\n\n";
    exit;
}
print exon_file_output "RefSeq\tEnsembl Transcript ID\tEnsembl Gene ID\tEnsembl Exon ID\tExon number\tChr\tStrand\tCoordinates\tCanonical\n";
my $transcript_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( 'Human', 'Core',
    'Transcript' );
foreach my$elmt (keys %gene_match){
	my $transcript_id = @{$gene_match{$elmt}}[0];
	my $tr = $transcript_adaptor->fetch_by_stable_id($transcript_id);
	my $refseq = $transcript_list{$elmt}{$transcript_id};
	my @content = get_all_exon_for_transcript($tr,$refseq,$elmt);
	print exon_file_output $_ for @content;
}

sub get_all_exon_for_transcript {
	my($transcript,$refseq,$gene_id) = @_;
	my @exons = @{ $transcript->get_all_Exons };	
	my $exon_number = 1;
	my @content;
	while (my$exon = shift @exons) {
		my $stable_id  = $exon->stable_id();
    	my $seq_region = $exon->seq_region_name();
    	my $start      = $exon->start();
    	my $end        = $exon->end();
    	my $strand     = $exon->strand();
    	my $content = sprintf ("%s\t%s\t%s\t%s\t%d\tchr%s\t%+d\t%s_%s\t%s\n",$refseq,$gene_id,$transcript->stable_id(),$stable_id,$exon_number,$seq_region,$strand,$start,$end,"YES");
		push (@content, $content);
		$exon_number++;
	}
	return @content;
}


sub intron_by_exon {
	my(@exon_list)=@_;
	my @intron_list;
	for(my$pos = 0; $pos<scalar(@exon_list);$pos++) {
		if ($pos+1 != scalar(@exon_list)){
			my($intron_pos_left) = $exon_list[$pos] =~ /^[0-9]+:([0-9]+)/;
			my($intron_pos_right) = $exon_list[$pos+1] =~ /^([0-9]+):[0-9]+/;
			my $intron_coord = $intron_pos_left."_".$intron_pos_right;
			push(@intron_list, $intron_coord)
		} 
	}
	return @intron_list;
}


sub no_doublon {
	my(@liste) = @_;
	my %unique;
	my @unique;
	foreach my $elmt (@liste) {
		unless ($unique{$elmt}){
			$unique{$elmt}++;
			push(@unique,$elmt);
		}
	}
	return @unique;
}



