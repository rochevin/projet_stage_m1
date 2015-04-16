use Bio::EnsEMBL::Registry;
use Bio::SeqIO;
use Data::Dumper;

$ENV{MYSQL_UNIX_PORT} = "/Applications/MAMP/tmp/mysql/mysql.sock";

my $r = "Bio::EnsEMBL::Registry";

$r->load_registry_from_db(-host => "localhost", -user => "root", -pass => "root", -port => "8889", -verbose => "0");

my $transcript_adaptor = $r->get_adaptor( 'Human', 'Core','Transcript' );

my $transcript = $transcript_adaptor->fetch_by_stable_id($ARGV[0]);

print "$ARGV[0] :".$transcript->seq_region_start."-".$transcript->seq_region_end."\n";

print "CDS de $ARGV[0] :\n".$transcript->translateable_seq(), "\n";

print "Exons du CDS : \n";
my @exons = @{ $transcript->get_all_translateable_Exons() };
while (my $exon = shift @exons) {
	print $exon->stable_id()."\n";
	print $exon->seq->seq."\n";
}


print "Tous les exons : \n";
my @exons = @{ $transcript->get_all_Exons() };
while (my $exon = shift @exons) {
	print $exon->stable_id()."\n";
	print $exon->start()."-".$exon->end()."\n";
	print $exon->seq->seq."\n";
}


print "sequence de $ARGV[0] :\n".$transcript->seq->seq, "\n";

print $transcript->cdna_coding_start."\n";
print $transcript->cdna_coding_end."\n";



# @features = @{$feature->get_all_alt_locations()};
# foreach $f (@features) {
# 	print $f->slice->seq_region_name,' ',$f->start, $f->end,"\n";
# }