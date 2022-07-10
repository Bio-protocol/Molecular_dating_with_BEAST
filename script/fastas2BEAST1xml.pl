#!/bin/perl
## Take a list of FASTA files and convert to BEASTv1 XML
## Hard coded for: unlinked substitution model - HKY, base frequencies estimated, gamma model of rate heterogeneity with 4 categories, unlinked molecular clock model - Strict Clock (uniform rates across branches), linked tree prior - Constant Size, as well as most default priors
## Kaichi Huang 2022 May

use warnings;
use strict;
use List::Util qw/sum max min/;
use Getopt::Long;
use Pod::Usage;

my $usage = "\nThis script takes a list of input files and converts them to a BEASTv1 XML. Usage:
	perl fastas2BEAST1xml.pl [options]
		--fastaList	[File]	A file list of all FASTAs of the partitions, one FASTA each line.
		--treeFile	[File]	Optional staring tree in NEWICK format. If specified, topology will be fixed.
		--out	[str]	Output prefix. Default: out.
		--chainLength	[int]	MCMC chain length. Default: 1000000.
		--screenEcho	[int]	Echo state to screen every. Default: 1000.
		--fileLog	[int]	Log MCMC parameters and trees to output every. Default: 1000.
		--help|-h	Display this help info\n";

# pre-defined
my $fastaList = '';
my $treeFile = '';
my $out = "out";
my $chainLength = 1000000;
my $screenEcho = 1000;
my $fileLog = 100;
my $help = undef;

if ( !GetOptions( 'fastaList=s'   => \$fastaList,
			   'treeFile=s'    => \$treeFile,
			   'out=s'         => \$out,
			   'chainLength=i' => \$chainLength,
			   'screenEcho=i'  => \$screenEcho,
			   'fileLog=i'     => \$fileLog,
			   'help|h!'       => \$help ) )
{
	pod2usage( { -verbose => 1,
			   -exitval => 1,
			   -message => 'Failed to parse command line.' } );
}
if ($help) {
	pod2usage( { -verbose => 0,
			   -exitval => 0,
			   -message => "$usage\n" } );
}
if (! -s $fastaList) {
	pod2usage( { -verbose => 0,
			   -exitval => 2,
			   -message => "Mandatory files not found:\n --fastaList\n"."$usage\n" } );
}

# Read sequences partition by partition
print STDERR "Reading FASTA files...\n";
my @taxa;
my %seqs;
my %seqs_len;
my $part_num = 1;
open LIST, $fastaList;
while (<LIST>) {
	chomp;
	my $file = $_;
	open FASTA, $file;
	my $sample;
	my $seq;
	$seqs{$part_num} = {};
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		if ($_ =~ m/>/) {
			$_ =~ s/>//g;
			$sample = $_;
			if ($part_num == 1) {
				push @taxa, $sample;
			}
		} else {
			$seq = $_;
			$seqs{$part_num}{$sample} = $seq;
		}
	}
	close FASTA;
	$seqs_len{$part_num} = length($seq);
	$part_num = $part_num + 1;
}
close LIST;

# Read (fixed) starting tree
my $fixTree;
if ($treeFile) {
	print STDERR "Reading NEWICK tree $treeFile...\n";
	open TRE, $treeFile;
	while (<TRE>) {
		chomp;
		$fixTree = $_;
	}
	close TRE;
}

my $n_part = keys %seqs;
my $total_len = sum(values(%seqs_len));

print STDERR "Start to write the XML file...\n";
open OUT, ">$out.xml";
# Head
print OUT
"<beast version=\"1.10.4\">
\n\n";
# The list of taxa to be analysed
print OUT
"	<taxa id=\"taxa\">
";
foreach my $sample (@taxa) {
print OUT "<taxon id=\"$sample\"/>\n";
}
print OUT
"	</taxa>
\n";
# The sequence alignments
foreach my $j (sort {$a<=>$b} keys %seqs) {
print OUT
"	<alignment id=\"alignment$j\" dataType=\"nucleotide\">
";
foreach my $sample (sort keys %{$seqs{$j}}) {
print OUT
"		<sequence>
			<taxon idref=\"$sample\"/>
			$seqs{$j}{$sample}
		</sequence>
";
}
print OUT
"	</alignment>
\n";
}
# The unique patterns from 1 to end
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"	<patterns id=\"part$a.patterns\" from=\"1\" strip=\"false\">
		<alignment idref=\"alignment$a\"/>
	</patterns>
\n";
}
# Initial tree model (Coalescent: Constant Size)
print OUT
"	<constantSize id=\"constant\" units=\"years\">
		<populationSize>
			<parameter id=\"constant.popSize\" value=\"0.04\" lower=\"0.0\"/>
		</populationSize>
	</constantSize>
\n";

# Construct a starting tree
if ($treeFile) {
print OUT
"<newick id=\"startingTree\">
	$fixTree
</newick>
\n"; # (*-w-*)
} else {
print OUT
"	<coalescentSimulator id=\"startingTree\">
		<taxa idref=\"taxa\"/>
		<constantSize idref=\"constant\"/>
	</coalescentSimulator>
\n";
}
# Generate a tree model (Link Trees)
print OUT
"	<treeModel id=\"treeModel\">
		<coalescentTree idref=\"startingTree\"/>
		<rootHeight>
			<parameter id=\"treeModel.rootHeight\"/>
		</rootHeight>
		<nodeHeights internalNodes=\"true\">
			<parameter id=\"treeModel.internalNodeHeights\"/>
		</nodeHeights>
		<nodeHeights internalNodes=\"true\" rootNode=\"true\">
			<parameter id=\"treeModel.allInternalNodeHeights\"/>
		</nodeHeights>
	</treeModel>
\n"; # /*'-'*/
print OUT
"	<treeLengthStatistic id=\"treeLength\">
		<treeModel idref=\"treeModel\"/>
	</treeLengthStatistic>
\n";
# Tree prior (Coalescent: Constant Size)
print OUT
"	<coalescentLikelihood id=\"coalescent\">
		<model>
			<constantSize idref=\"constant\"/>
		</model>
		<populationTree>
			<treeModel idref=\"treeModel\"/>
		</populationTree>
	</coalescentLikelihood>
\n";

# Clock model (Link Clock Models; Strict clock)
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"	<strictClockBranchRates id=\"part$a.branchRates\">
		<rate>
			<parameter id=\"part$a.clock.rate\" value=\"1.0\"/>
		</rate>
	</strictClockBranchRates>
	<rateStatistic id=\"part$a.meanRate\" name=\"part$a.meanRate\" mode=\"mean\" internal=\"true\" external=\"true\">
		<treeModel idref=\"treeModel\"/>
		<strictClockBranchRates idref=\"part$a.branchRates\"/>
	</rateStatistic>
\n";
}

# Nucleotide substitution model (Unlink Subst. Models) frequency estimated)
for (my $a=1; $a<=$n_part; $a=$a+1) {
my $weight = $total_len/$seqs_len{$a};
my $nu = 1/$n_part;
print OUT
"	<HKYModel id=\"part$a.hky\">
		<frequencies>
			<frequencyModel dataType=\"nucleotide\">
				<frequencies>
					<parameter id=\"part$a.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>
				</frequencies>
			</frequencyModel>
		</frequencies>
		<kappa>
			<parameter id=\"part$a.kappa\" value=\"2.0\" lower=\"0.0\"/>
		</kappa>
	</HKYModel>

	<siteModel id=\"part$a.siteModel\">
		<substitutionModel>
			<HKYModel idref=\"part$a.hky\"/>
		</substitutionModel>
		<relativeRate weight=\"$weight\">
			<parameter id=\"part$a.nu\" value=\"$nu\" lower=\"0.0\" upper=\"1.0\"/>
		</relativeRate>
		<gammaShape gammaCategories=\"4\">
			<parameter id=\"part$a.alpha\" value=\"0.5\" lower=\"0.0\"/>
		</gammaShape>
	</siteModel>
	
	<statistic id=\"part$a.mu\" name=\"mu\">
		<siteModel idref=\"part$a.siteModel\"/>
	</statistic>
\n";
}

# Likelihood for tree given sequence data
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"	<treeDataLikelihood id=\"part$a.treeLikelihood\" useAmbiguities=\"false\">
		<partition>
			<patterns idref=\"part$a.patterns\"/>
			<siteModel idref=\"part$a.siteModel\"/>
		</partition>
		<treeModel idref=\"treeModel\"/>
		<strictClockBranchRates idref=\"part$a.branchRates\"/>
	</treeDataLikelihood>\n";
}

# Define operators
print OUT
"	<operators id=\"operators\" optimizationSchedule=\"log\">
";
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"	<scaleOperator scaleFactor=\"0.75\" weight=\"1\">
		<parameter idref=\"part$a.kappa\"/>
	</scaleOperator>
	<deltaExchange delta=\"0.01\" weight=\"1\">
		<parameter idref=\"part$a.frequencies\"/>
	</deltaExchange>
	<scaleOperator scaleFactor=\"0.75\" weight=\"1\">
		<parameter idref=\"part$a.alpha\"/>
	</scaleOperator>
\n";
}
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"	<scaleOperator scaleFactor=\"0.75\" scaleAll=\"true\" ignoreBounds=\"true\" weight=\"3\">
		<parameter idref=\"treeModel.allInternalNodeHeights\"/>
	</scaleOperator>\n";
}
# Copyright: Kaichi Huang
if ( ! $treeFile ) {
print OUT
"		<subtreeSlide size=\"1.0\" gaussian=\"true\" weight=\"30\">
			<treeModel idref=\"treeModel\"/>
		</subtreeSlide>
		<narrowExchange weight=\"30\">
			<treeModel idref=\"treeModel\"/>
		</narrowExchange>
		<wideExchange weight=\"3\">
			<treeModel idref=\"treeModel\"/>
		</wideExchange>
		<wilsonBalding weight=\"3\">
			<treeModel idref=\"treeModel\"/>
		</wilsonBalding>
";
}
print OUT
"		<scaleOperator scaleFactor=\"0.75\" weight=\"3\">
			<parameter idref=\"treeModel.rootHeight\"/>
		</scaleOperator>
		<uniformOperator weight=\"30\">
			<parameter idref=\"treeModel.internalNodeHeights\"/>
		</uniformOperator>
		<scaleOperator scaleFactor=\"0.75\" weight=\"3\">
			<parameter idref=\"constant.popSize\"/>
		</scaleOperator>
	</operators>
\n";

# Define MCMC
print OUT
"	<mcmc id=\"mcmc\" chainLength=\"$chainLength\" autoOptimize=\"true\" operatorAnalysis=\"$out.ops\">
";
print OUT
"		<joint id=\"joint\">
			<prior id=\"prior\">
";
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"				<logNormalPrior mu=\"1.0\" sigma=\"1.25\" offset=\"0.0\">
					<parameter idref=\"part$a.kappa\"/>
				</logNormalPrior>
				<dirichletPrior alpha=\"1.0\" sumsTo=\"1.0\">
					<parameter idref=\"part$a.frequencies\"/>
				</dirichletPrior>
				<exponentialPrior mean=\"0.5\" offset=\"0.0\">
					<parameter idref=\"part$a.alpha\"/>
				</exponentialPrior>
";
}
print OUT
"				<gammaPrior shape=\"10.0\" scale=\"0.004\" offset=\"0.0\">
					<parameter idref=\"constant.popSize\"/>
				</gammaPrior>
				<coalescentLikelihood idref=\"coalescent\"/>
";
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"				<strictClockBranchRates idref=\"part$a.branchRates\"/>
";
}
print OUT
"			</prior>
			<likelihood id=\"likelihood\">
";
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"				<treeDataLikelihood idref=\"part$a.treeLikelihood\"/>
";
}
print OUT
"			</likelihood>
		</joint>
		<operators idref=\"operators\"/>
\n";
# Write log to screen
print OUT
"		<log id=\"screenLog\" logEvery=\"$screenEcho\">
			<column label=\"Joint\" dp=\"4\" width=\"12\">
				<joint idref=\"joint\"/>
			</column>
			<column label=\"Prior\" dp=\"4\" width=\"12\">
				<prior idref=\"prior\"/>
			</column>
			<column label=\"Likelihood\" dp=\"4\" width=\"12\">
				<likelihood idref=\"likelihood\"/>
			</column>
			<column label=\"rootHeight\" sf=\"6\" width=\"12\">
				<parameter idref=\"treeModel.rootHeight\"/>
			</column>
		</log>
\n";
# Write log to file
print OUT
"		<log id=\"fileLog\" logEvery=\"$fileLog\" fileName=\"$out.log\" overwrite=\"false\">
			<joint idref=\"joint\"/>
			<prior idref=\"prior\"/>
			<likelihood idref=\"likelihood\"/>
			<parameter idref=\"treeModel.rootHeight\"/>
			<treeLengthStatistic idref=\"treeLength\"/>
			<parameter idref=\"constant.popSize\"/>
";
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"			<parameter idref=\"part$a.kappa\"/>
			<parameter idref=\"part$a.frequencies\"/>
			<parameter idref=\"part$a.alpha\"/>
";
}
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"			<parameter idref=\"part$a.clock.rate\"/>
";
}
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"			<rateStatistic idref=\"part$a.meanRate\"/>
";
}
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"			<treeDataLikelihood idref=\"part$a.treeLikelihood\"/>
";
}
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"			<strictClockBranchRates idref=\"part$a.branchRates\"/>
";
}
print OUT
"			<coalescentLikelihood idref=\"coalescent\"/>
		</log>
\n";
# Write tree log to file
print OUT
"		<logTree id=\"treeFileLog\" logEvery=\"$fileLog\" nexusFormat=\"true\" fileName=\"$out.trees\" sortTranslationTable=\"true\">
			<treeModel idref=\"treeModel\"/>
";
for (my $a=1; $a<=$n_part; $a=$a+1) {
print OUT
"			<trait name=\"rate\" tag=\"gene$a.rate\">
				<strictClockBranchRates idref=\"part$a.branchRates\"/>
			</trait>
";
}
print OUT
"			<joint idref=\"joint\"/>
		</logTree>
";
# <zzz>
print OUT
"	</mcmc>
\n";
print OUT
"	<report>
		<property name=\"timer\">
			<mcmc idref=\"mcmc\"/>
		</property>
	</report>
\n\n";
# Tail
print OUT
"</beast>
\n";
close OUT;
print STDERR "Done!\n";
