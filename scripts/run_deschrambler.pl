#!/usr/bin/env perl
use warnings;
use strict;

use FindBin qw($Bin);
use File::Basename;
use Getopt::Long;
use Cwd 'abs_path';
use Switch;

## Required parameters
my $help;				## Help page
my $threads=1;			## Number of threads
my $params;				## Parameter file
my $resolution=10000;	## Resolution
my $outDir="./";		## Output directory
my $path_conf = abs_path("$Bin/../path.conf");

my $tree;				## path of newick tree
my $reference;			## name of reference
my @outgroup=();		## name of outgroup
my $reference_fa;		## path of reference genome sequence
my @outgroup_fa=();		## path of outgroup genome sequence

my $deschrambler;		## path of deschrambler
my $kent;				## path of kent directory

my $options = GetOptions (
			"threads|t=i" => \$threads,
			"params|p=s" => \$params,
			"outDir|o=s" => \$outDir,
			"help|h" => \$help,
		);

sub HELP{
	my $src = basename($0);
	print "\nUsage: $src [options] <filename>\n\n";
	print "\t-t <integer> Number of threads (defulat: 1)\n";
	print "\t-p <filename> Parameter file\n";
	print "\t-o <filename> Results of output directory (default: ./)\n";
	print "\t-h|help print this page\n\n";
	exit;
}

### Help Prints ###
if (!-e $path_conf){    
	print STDERR "[ERROR] Not exists 'path.conf' file\n";
	HELP();
}
if (!defined $params){
	print STDERR "[ERROR] Please check path of parameter file\n";
	HELP();
}

$outDir = abs_path ("$outDir");


open PATH, "$path_conf";
while (<PATH>){
	chomp;
	next if /^#/;
	next if /""/;
	my ($program, $path) = split (/=/, $_);
	switch ($program){
		case("DESCHRAMBLER"){	$deschrambler = $path;	}
		case("kent"){	$kent = $path;	}
	}
}
close PATH;

open PARAM, "$params";
while (<PARAM>){
	chomp;
	if ($_ =~ /^#/ || $_ eq ""){	next;	}
	my @split = split (/\s+/, $_);
	if ($split[0] eq "Reference"){
		if (!-e $split[2]){  die "[ERROR] Please check path of the reference genome sequence\n"; }
		$reference = $split[1];
		$reference_fa = $split[2];
	}
	elsif ($split[0] eq "Outgroup"){
		if (!-e $split[2]){  die "[ERROR] Please check path of the outgroup genome sequence\n";  }
		push (@outgroup, $split[1]);
		push (@outgroup_fa, $split[2]);
	}
	elsif ($split[0] eq "Resolution"){	$resolution = $split[1];	}
	elsif ($split[0] eq "TREE"){	$tree = $split[1];	}
}
close PARAM;


## Running the whole genome alignment ##
open LOG, ">>$outDir/chainNet/log";

$reference_fa = "$outDir/chainNet/$reference\.fa";
for (my $i=0; $i<=$#outgroup; $i++){	$outgroup_fa[$i] = "$outDir/chainNet/$outgroup[$i]\.fa";	}

print STDERR "### 3. Step of meta-assembly ###\n";
print STDERR "    (1) Whole genome alignment of RACA assemblies with $reference\n";
print STDERR "       (1-1) Whole genome alignment of RACA_spades to $reference\n";
print LOG "### Whole genome alignment of RACA_spades to $reference ###\n";
WGA ($reference_fa, "$outDir/RACA/RACA_spades.fa");

print STDERR "       (1-2) Whole genome alignment of RACA_masurca to $reference\n";
print LOG "\n### Whole genome alignment of RACA_masurca to $reference ###\n";
WGA ($reference_fa, "$outDir/RACA/RACA_masurca.fa");

print STDERR "       (1-3) Whole genome alignment of RACA_soap2 to $reference\n";
print LOG "\n### Whole genome alignment of RACA_soap2 to $reference ###\n";
WGA ($reference_fa, "$outDir/RACA/RACA_soap2.fa");
close LOG;

## Making DESCHRMALBER inputs ##
`mkdir -p $outDir/meta-assembly`;
print STDERR "    (2) Making the inputs for meta-assembly\n";
DESCHINPUT ("spades", "masurca", "soap2");

## Running DESCHRAMBLER ##
chdir ("$outDir/meta-assembly");
print STDERR "    (3) Running the inputs for meta-assembly\n";
exeDESCH ("$outDir/meta-assembly/params", "$outDir/RACA/spades/rec_chrs.denovo_spades.unused.fa");


#=======================================================================================================#
#												Functions												#
#=======================================================================================================#

### Whole genome alignment ###
sub WGA {
	my ($ref_fa, $tar_fa) = @_;
	`$Bin/whole_genome_alignment.pl -p $threads -m $resolution -r $ref_fa -t $tar_fa -o $outDir/chainNet 2>> $outDir/chainNet/log`;
}


### Making the DESCHRAMBLER inputs ###
sub DESCHINPUT {
	my ($ass_1, $ass_2, $ass_3) = @_;
	`mkdir -p $outDir/meta-assembly/params`;

	## config.file
	open CONFIG, ">$outDir/meta-assembly/params/config.file";
	print CONFIG ">netdir\n$outDir\/chainNet\n\n";
	print CONFIG ">chaindir\n$outDir\/chainNet\n\n";
	print CONFIG ">species\n$reference 0 0\nRACA_$ass_1 1 0\nRACA_$ass_2 1 0\nRACA_$ass_3 1 0\n";
	for (my $i=0; $i<=$#outgroup; $i++){	print CONFIG "$outgroup[$i] 2 0\n";	}
	print CONFIG "\n>resolution\n<resolutionwillbechanged>\n";
	close CONFIG;

	## Makefile
	`cp $Bin/DESCHRAMBLER_Makefile.SFs $outDir/meta-assembly/params/Makefile.SFs`;

	## tree.txt
	open TREE, "$tree";
	open TREEF, ">$outDir/meta-assembly/params/tree.txt";
	my $format = <TREE>;
	if ($format =~ (s/TARGET@\:(\d+\.*\d*)/((RACA_$ass_1:0.000001,RACA_$ass_2:0.000001):0.000001,RACA_$ass_3:0.000002)@\:$1/)){	print TREEF "$format";	}
	close TREEF;
	close TREE;

	## deschrambler_param.txt
	open PARAM, ">$outDir/meta-assembly/params/deschrambler_param.txt";
	print PARAM "REFSPC=$reference\n\n";
	print PARAM "OUTPUTDIR=$outDir/meta-assembly\n\n";
	print PARAM "RESOLUTION=$resolution\n\n";
	print PARAM "TREEFILE=$outDir/meta-assembly/params/tree.txt\n\n";
	print PARAM "MINADJSCR=0.0001\n\n";
	print PARAM "CONFIGSFSFILE=$outDir/meta-assembly/params/config.file\n";
	print PARAM "MAKESFSFILE=$outDir/meta-assembly/params/Makefile.SFs";
	close PARAM;
}


### DESCHRAMBLER ###
sub exeDESCH {
	my ($paramDir, $unused_scf) = @_;
	`$deschrambler $paramDir/deschrambler_param.txt 2> $paramDir/../log`;
	`$Bin/DESCHRAMBLER_create_fasta.block.pl $paramDir/../APCF_RACA_spades.merged.map $outDir/RACA/RACA_spades.fa $paramDir/../APCF_RACA_spades.fa`;
	`cp $paramDir/../APCF_RACA_spades.fa $outDir/meta_assembled_scaf.fa`;
	`cat $unused_scf >> $outDir/meta_assembled_scaf.fa`;
}

`rm -f $outDir/chainNet/$reference\.fa`;
for (my $i=0; $i<=$#outgroup; $i++){	`rm -f $outDir/chainNet/$outgroup[$i]\.fa`;	}
