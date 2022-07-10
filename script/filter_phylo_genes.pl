#!/bin/perl
## Filter gene fasta based on averge missing rate across samples and generate a gene list for downstream analysis
## run as: cat {gene_names} | perl filter_phylo_genes.pl > {list}
## Kaichi Huang 2022 May

use List::Util qw(sum);
use List::MoreUtils qw(none);

my @gene_list;

while(<STDIN>){
	chomp;
	my $gene_name = $_;
	my @missing;
	open GENE, "cache/$gene_name.fasta";
	while(<GENE>){
		chomp;
		if(/^>/){
			next;
		} else {
			my $seq = $_;
			my $count = ($seq=~tr/N/N/);
			my $missing = $count / length($seq);
			push @missing, $missing;
		}
	}
	close GENE;
	if ((sum(@missing)/($#missing+1) < 0.3) and (none {$_>0.5} @missing)) {
		push @gene_list, "$gene_name";
	}
}

# Output the good ones
foreach my $gene (@gene_list) {
	print "$gene\n";
}
