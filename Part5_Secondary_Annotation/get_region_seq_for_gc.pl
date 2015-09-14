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

# On utilise l'API d'ensembl
my $r = "Bio::EnsEMBL::Registry";
# On se connecte à la base de donnée locale ensembl 75
$r->load_registry_from_db(-host => "localhost", -user => "root", -pass => "root", -port => "8889", -verbose => "0");
# On utilise gene_adaptor pour récupérer les informations sur nos gènes
my $gene_adaptor=$r->get_adaptor("Mouse" ,"core", "Gene");
my $slice_adaptor=$r->get_adaptor("Mouse" ,"core", "Slice");

foreach my$id (keys %gene_list) {
    my $gene = $gene_adaptor->fetch_by_stable_id($id);

    my $gene_strand = $gene->strand();

    my $canonical_transcript=$gene->canonical_transcript();

    my $transcript_id = $canonical_transcript->stable_id();

    my $seq_region = $canonical_transcript->seq_region_name();
    #On définit les tailles des régions
    my $debut = $canonical_transcript->start()-20000;
    my $fin = $canonical_transcript->end()+20000;
    #Et on récupère les informations sur la région
    my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $seq_region, $debut, $fin );
    if (defined $slice){
        $full_seq = $slice->seq();
    }else {
        # $full_seq = get_seq_with_exon($intron);
        my $slice = $slice_adaptor->fetch_by_region( 'supercontig', $seq_region, $debut, $fin );
        $full_seq = $slice->seq();
    }
    # Si le brin est -, on fait le reverse complement de la séquence
    if ($gene_strand<0) {
        $full_seq = Bio::Seq->new(-seq => $full_seq);
        $full_seq->revcom;
        $full_seq= $full_seq->seq();
    }
    # On détermine nos séquences d'intérêt
    my $begin_seq = substr($full_seq,0,20000);
    my $end_seq = substr($full_seq,-20000);
    ($begin_seq,$end_seq) = ($end_seq,$begin_seq) if ($gene_strand<0);
    # On définit un objet Bio::Seq pour les séquences
    my $begin_fasta_seq = Bio::Seq->new(-display_id => $transcript_id, -seq => $begin_seq,-description => "first");
    my $end_fasta_seq = Bio::Seq->new(-display_id => $transcript_id, -seq => $end_seq,-description => "last");
    # Puis on écris les séquences dans les fichier
    $out->write_seq($begin_fasta_seq);
    $out->write_seq($end_fasta_seq);
}
$out->close();