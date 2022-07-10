#!/bin/perl
## Generate fasta from GVCF file
## run as: cat {g.vcf} | perl gvcf2fasta_nogaps.pl > {fa}
## modified from G. L. Owens
## Kaichi Huang 2022 May

use warnings;
use strict;

my $min_depth = 2; # minimum depth to call a genotype
my $sample;
my $header;
my $seq;
my $position = 1;

while(<STDIN>){
	chomp;
	if ($_ =~ /^##/) {next;} # skip header rows
	if ($_ =~ /^#/) {
		# extract sample name from the GVCF file
		my @a = split(/\t/,$_);
		$sample = $a[9];
		$header .= ">$sample";
		next;
	}
	my @a = split(/\t/,$_);
	my $ref = $a[3];
	my $alt = $a[4];
	my @alts = split(/,/,$alt);
	my $genotype;
	my @field_names = split(/:/,$a[8]);
	# Skip weird rows
	unless (exists($field_names[2])) {
		$genotype = "N";
		goto PRINTOUT;
	}
	if ($field_names[2] ne "DP") {
		$genotype = "N";
		goto PRINTOUT;
	}
	my $info = $a[9];
	my @infos = split(/:/,$info);
	my $depth = $infos[3];
	my $call = $infos[0];
	my @calls;
	if ($call =~ /\//) {
		@calls = split(/\//,$call);
	} elsif ($call =~ /\|/) {
		@calls = split(/\|/,$call);
	}
	my $rand_call = int(rand(2));
	my $called_value = $calls[$rand_call];
	if ($depth < $min_depth) {
		$genotype = "N";
		goto PRINTOUT;
	}
	if ($called_value == 0) {
		$genotype = $ref;
	if (length($genotype) > 1) {
		$genotype = "N";
	}
	} else {
		$genotype = $alts[($called_value - 1)];
		if (length($genotype) > 1){
			$genotype = "N";
		}
	}
	if ($genotype eq '*') {
		$genotype = "N";
	}
	PRINTOUT:
	$seq .= "$genotype";
	$position++;
}

# Output the sequence in FASTA format
print "$header\n$seq";
