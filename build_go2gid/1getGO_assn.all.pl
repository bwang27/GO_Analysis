#get all GO childs and parents for any GO term  

#!/usr/bin/perl -w

require "../../perl.fun.inc.pl";
use GO::Parser;
$goVersion="190105";
$go_defi_f="../../../data/geneOntology/go-basic.obo";
my $parser = new GO::Parser({handler=>'obj'}); # create parser object
$parser->parse("$go_defi_f"); # parse file -> objects

my $output = "GO_assn/$goVersion/go2go.all.out";

##############################################################
create_dir_ifNotExist($output);
my %cnt;

open (OUT, "> $output");	
print OUT "#Category\tGO_ID\tGO_term\tNum_Parent\tNum_Child\tALL_Children\n";

my $graph = $parser->handler->graph;  # get L<GO::Model::Graph> object
my $all_terms = $graph->get_all_nodes();

foreach my $term (sort {$a->acc cmp $b->acc} @$all_terms) {

	next if ($term->acc !~ /^GO:/);

	$cnt{'all'}++;

	my $nsp = $term->namespace;
	if ($term->is_obsolete) {
		$cnt{'obs'}++;
		next;
	}

	my $parent_terms = $graph->get_recursive_parent_terms($term->acc);
	my $child_terms = $graph->get_recursive_child_terms($term->acc);

	my @a = &proc($parent_terms);   #remove duplicate terms
	my @b = &proc($child_terms);	#remove duplicate terms
	my $child_list = join ("|",@b);
	if ($child_list eq "") { $child_list = "-"; }

	my $line = sprintf "%s\t%s\t%s\t%s\t%s", $term->acc, $term->name, scalar @a, scalar @b, $child_list;

	printf OUT "$nsp\t$line\n";
	$cnt{"OUTPUT: $output"}++;
}
close OUT;

sub proc {
	my $ref = $_[0];

	my @arrays;
	foreach my $term (@$ref) {
		push (@arrays, $term->acc);
	}

	my %h   = map { $_, 1 }		@arrays; #remove duplicates
	my @arrays2 = sort keys %h;
	return @arrays2;
}

foreach my $temp (sort keys %cnt) {
	print "$temp\t$cnt{$temp}\n";
}