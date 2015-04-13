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
    else {
    	my @line_content = split("\t", $line);
    	$transcript_list{$line_content[0]}=$line_content[1];
    }
}
close file_list;


my $r = "Bio::EnsEMBL::Registry";

$r->load_registry_from_db(-host => "localhost", -user => "root", -pass => "root", -port => "8889", -verbose => "0");

$transcript_adaptor = r->get_adaptor( 'Human', 'Core','Transcript' );

$transcript = $transcript_adaptor->fetch_by_stable_id('ENST00000201961');