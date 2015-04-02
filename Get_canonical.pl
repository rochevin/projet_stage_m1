use Bio::EnsEMBL::Registry;


my $r = "Bio::EnsEMBL::Registry";

$r->load_registry_from_db(-host => "localhost", -user => "root", -pass => "root", -port => "8889", -verbose => "0");

my $gene_adaptor=$r->get_adaptor("Human" ,"core", "Gene");

my $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000000003');

my $canonical_transcript=$gene->canonical_transcript();

print $canonical_transcript->stable_id();
# 	my $r ="Bio::EnsEMBL::Registry";

 
# 	# Pour chaque id de la liste gene ...
# 	foreach my $main_id (@gene_id) {
# 		# On définit une nouvelle connexion à la base de donnée ensembl
# 		$r->load_registry_from_db(-host => "ensembldb.ensembl.org", -user => "anonymous", -verbose => "0");
# 		# On se connecte à la bonne base de donnée
# 		my $orga = $dataset{$main_id.'_Taxonomy_Name'};
# 		my $gene_adaptor=$r->get_adaptor($orga ,"core", "Gene");

# 		# On récupère l'information grâce au gene symbol
# 		my $gene = $gene_adaptor->fetch_by_display_label($dataset{$main_id.'_Gene_Symbol'});

# 		$dataset{$main_id.'_ENS_geneID'}=$gene->stable_id(); #On enregistre l'id du gene 
# 		# On enregistre ses coordonnées
# 		$dataset{$main_id.'_ENS_browserLOC'}=$gene->seq_region_name().":".$gene->seq_region_start()."-".$gene->seq_region_end();
# 		# On obtient la liste de tous les transcrits de ce gene
# 		my@trans=@{$gene->get_all_Transcripts()};
# 		my@liste_prot=();
# 		# On parcours ces transcrits...
# 		foreach $trans (@trans) {
# 			if (defined $trans->translation()) { #Si il existe une proteine pour ce transcrit
# 				$prot=$trans->translation();
# 				# Alors on enregistre dans dataset un tableau associatif avec comme clé l'id du transcrit et sa valeur l'id de la protéine
# 				$dataset{$main_id.'_ENS_transID'}{$trans->stable_id()}=$prot->stable_id;
# 			}
# 			else {
# 				# Sinon on enregistre juste l'id du transcrit comme clé et "none" comme valeur
# 				$dataset{$main_id.'_ENS_transID'}{$trans->stable_id()}="none";
# 			}
# 		}
# 		# On récupère le transcrit canonique définit par ensembl -> nous servira à récupérer l'id uniprot "principal"
# 		my $canonical_transcript=$gene->canonical_transcript();
# 		$canonical_prot=$canonical_transcript->translation();
# 		# De la même façon on enregistre la clé id transcrit canonique et sa valeur la protéine
# 		$dataset{$main_id.'_ENS_canonical_transID'}{$canonical_transcript->stable_id()}=$canonical_prot->stable_id();


# 		# À utiliser uniquement si on veut récupérer les termes GO via la base de donnée EnsEMBL

# 		# my @dblinks = @{$canonical_transcript->get_all_DBLinks("GO")};
# 		# my %verif_GO_term;
# 		# while (my $dblinks = shift @dblinks) {
# 			# my $GO_id = $dblinks->display_id;
# 			# unless (exists $verif_GO_term{$GO_id}) {
# 	      		# $verif_GO_term{$GO_id}=1;
# 	      		# push @{$dataset{$main_id.'_ENS_GOTERM'}},$GO_id;
# 			# }
# 		# }
# 	print "EnsEMBL => ".$main_id." OK !\n";
# 	}
# 	# On renvoit le tableau associatif avec les nouvelles données dedans
# 	return %dataset;
# }
# return 1;