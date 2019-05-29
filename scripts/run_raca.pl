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
my $resolution=1;	## Resolution
my $outDir="./";		## Output directory
my $path_conf = abs_path("$Bin/../path.conf");

my $tree;				## path of newick tree
my $reference;			## name of reference
my $numchr;			## no. of chromosomes
my @outgroup=();		## name of outgroup
my $reference_fa;		## path of reference genome sequence
my @outgroup_fa=();		## path of outgroup genome sequence

my $bwa;				## path of bwa
my $raca;				## path of raca
my $kent;				## path of kent directory
my $samtools;				## path of samtools

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
if (!-e $path_conf){    print STDERR "[ERROR] Not exists 'path.conf' file\n";   }
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
		case("bwa"){	$bwa = $path;	}
		case("RACA"){	$raca = $path;	}
		case("kent"){	$kent = $path;	}
		case("samtools"){	$samtools = $path;	}
	}
}
close PATH;


my @reads = ();
open PARAM, "$params";
while (<PARAM>){
	chomp;
	if ($_ =~ /^#/ || $_ eq ""){	next;	}
	my @split = split (/\s+/, $_);
	if ($split[0] eq "Reference"){
		if (!-e $split[2]){	die "[ERROR] Please check path of the reference genome sequence\n";	}
		$reference = $split[1];
		$reference_fa = $split[2];
	}
	elsif ($split[0] eq "Outgroup"){
		if (!-e $split[2]){	die "[ERROR] Please check path of the outgroup genome sequence\n";	}
		push (@outgroup, $split[1]);
		push (@outgroup_fa, $split[2]);
	}
	elsif ($split[0] eq "p1" || $split[0] eq "p2" || $split[0] eq "m1" || $split[0] eq "m2"){
		push (@reads, $split[1]);
	}
	elsif ($split[0] eq "Resolution"){	$resolution = $split[1];	}
	elsif ($split[0] eq "TREE"){	$tree = $split[1];	}
}
close PARAM;


### Running the whole genome alignment ###
`mkdir -p $outDir/chainNet`;

`$Bin/fasta_id_convert.pl $reference_fa $outDir/chainNet/$reference\.fa`;
$reference_fa = "$outDir/chainNet/$reference\.fa";

$numchr = `grep ">" $reference_fa | wc -l`;
chomp($numchr);

for (my $i=0; $i<=$#outgroup; $i++){
	`$Bin/fasta_id_convert.pl $outgroup_fa[$i] $outDir/chainNet/$outgroup[$i]\.fa`;
	$outgroup_fa[$i] = "$outDir/chainNet/$outgroup[$i]\.fa";
}

open LOG, ">$outDir/chainNet/log";
print STDERR "### 2. Step of RACA ###\n";
print STDERR "    (1) Whole genome alignment of denovo assemblies with $reference\n";
print STDERR "       (1-1) Whole genome alignment of denovo_spades to $reference\n";
print LOG "### Whole genome alignment of denovo_spades to $reference ###\n";
WGA ($reference_fa, "$outDir/denovo_assemblies/denovo_spades.fa");

print STDERR "       (1-2) Whole genome alignment of denovo_masurca to $reference\n";
print LOG "\n### Whole genome alignment of denovo_masurca to $reference ###\n";
WGA ($reference_fa, "$outDir/denovo_assemblies/denovo_masurca.fa");

print STDERR "       (1-3) Whole genome alignment of denovo_soap2 to $reference\n";
print LOG "\n### Whole genome alignment of denovo_soap2 to $reference ###\n";
WGA ($reference_fa, "$outDir/denovo_assemblies/denovo_soap2.fa");

for (my $i=0; $i<=$#outgroup; $i++){
	my $j=$i+4;
	print STDERR "       (1-$j) Whole genome alignment of $outgroup[$i] to $reference\n";
	print LOG "\n### Whole genome alignment of $outgroup[$i] to $reference ###\n";
	WGA ($reference_fa, $outgroup_fa[$i]); 
}
close LOG;

### Read Mapping ###
`mkdir -p $outDir/mapping`;
print STDERR "    (2) Mapping the reads to denovo assemblies\n";

print STDERR "       (2-1) Mapping the reads to denovo_spades\n";
`mkdir -p $outDir/mapping/denovo_spades`;
chdir ("$outDir/mapping/denovo_spades");
READMAP ("$outDir/denovo_assemblies/denovo_spades.fa", @reads);

print STDERR "       (2-2) Mapping the reads to denovo_masurca\n";
`mkdir -p $outDir/mapping/denovo_masurca`;
chdir ("$outDir/mapping/denovo_masurca");
READMAP ("$outDir/denovo_assemblies/denovo_masurca.fa", @reads);

print STDERR "       (2-3) Mapping the reads to denovo_soap2\n";
`mkdir -p $outDir/mapping/denovo_soap2`;
chdir ("$outDir/mapping/denovo_soap2");
READMAP ("$outDir/denovo_assemblies/denovo_soap2.fa", @reads);

### Making the RACA inputs ###
`mkdir -p $outDir/RACA`;
print STDERR "    (3) Making the inputs for RACA\n";

`mkdir -p $outDir/RACA/spades`;
print STDERR "       (3-1) Making the inputs for RACA of denovo_spades\n";
RACAINPUT ("spades", "$outDir/denovo_assemblies/denovo_spades.fa");

`mkdir -p $outDir/RACA/masurca`;
print STDERR "       (3-2) Making the inputs for RACA of denovo_masurca\n";
RACAINPUT ("masurca", "$outDir/denovo_assemblies/denovo_masurca.fa");

`mkdir -p $outDir/RACA/soap2`;
print STDERR "       (3-3) Making the inputs for RACA of denovo_soap2\n";
RACAINPUT ("soap2", "$outDir/denovo_assemblies/denovo_soap2.fa");

### Running RACA ###
chdir ("$outDir/RACA/spades");
print STDERR "    (4) Running the RACA\n";

print STDERR "       (4-1) Running the RACA of denovo_spades\n";
exeRACA ("$outDir/RACA/spades/params", "denovo_spades", "$outDir/denovo_assemblies/denovo_spades.fa");
`cp $outDir/RACA/spades/rec_chrs.denovo_spades.RACA.fa $outDir/RACA/RACA_spades.fa`;

chdir ("$outDir/RACA/masurca");
print STDERR "       (4-2) Running the RACA of denovo_masurca\n";
exeRACA ("$outDir/RACA/masurca/params", "denovo_masurca", "$outDir/denovo_assemblies/denovo_masurca.fa");	
`cp $outDir/RACA/masurca/rec_chrs.denovo_masurca.RACA.fa $outDir/RACA/RACA_masurca.fa`;

chdir ("$outDir/RACA/soap2");
print STDERR "       (4-3) Running the RACA of denovo_soap2\n";
exeRACA ("$outDir/RACA/soap2/params", "denovo_soap2", "$outDir/denovo_assemblies/denovo_soap2.fa");
`cp $outDir/RACA/soap2/rec_chrs.denovo_soap2.RACA.fa $outDir/RACA/RACA_soap2.fa`;

#=======================================================================================================#
#												Functions												#
#=======================================================================================================#

### Whole genome alignment ###
sub WGA {
	my ($ref_fa, $tar_fa) = @_;
	`$Bin/whole_genome_alignment.pl -p $threads -m $resolution -r $ref_fa -t $tar_fa -o $outDir/chainNet 2>> $outDir/chainNet/log`;
}


### Mapping ###
sub READMAP {
	my @input = @_;
	## Indexing ##
	`$bwa index $input[0] 2> log`;

	## Mapping ##
	for (my $i=1; $i<=$#input; $i+=2){
		my $j=2*$i-1;
		`$bwa mem $input[0] -t $threads $input[$i] $input[$i+1] 2> log | $samtools view -h -@ $threads -o denovo_mapping_$j.bam -`;
	}
}


### Making the RACA inputs ###
sub RACAINPUT {
	my ($ass_name, $target) = @_;
	`mkdir -p $outDir/RACA/$ass_name/params`;

	my $mapDir = "$outDir/mapping/denovo_$ass_name";

	## config.file 
	open CONFIG, ">$outDir/RACA/$ass_name\/params/config.file";
	print CONFIG ">netdir\n$outDir/chainNet\n\n";
	print CONFIG ">chaindir\n$outDir/chainNet\n\n";
	print CONFIG ">species\n$reference 0\ndenovo_$ass_name 1\n";
	for (my $i=0; $i<=$#outgroup; $i++){	print CONFIG "$outgroup[$i] 2\n";	}
	print CONFIG "\n>resolution\n$resolution\n\n";
	print CONFIG ">numchr\n$numchr\n";
	close CONFIG;

	## library_info.txt & library_mapping.txt
	my @mapping_files = ();
	opendir MAPDIR, "$mapDir";
	while (my $files = readdir (MAPDIR)){	@mapping_files = glob ("$mapDir/*.bam");	}
	closedir MAPDIR;
	
	`touch $outDir/RACA/$ass_name/params/library_info.txt`;
	
	open MAPLIB, ">$outDir\/RACA/$ass_name\/params/library_mapping.txt";
	for (my $i=0; $i<=$#mapping_files; $i++){
		my @split = split (/\//, $mapping_files[$i]);
		print MAPLIB "$split[-1]\.mapping\t$split[-1]\n";

		my $bname = basename($mapping_files[$i]);
		`$Bin/sambam2raca_files.pl $bname $mapping_files[$i] $mapDir/$bname.mapping $samtools >> $outDir/RACA/$ass_name/params/library_info.txt`;
	}
	close MAPLIB;

	`mkdir -p $outDir/readmap/denovo_$ass_name`;
	`mv $outDir/mapping/denovo_$ass_name/*.mapping $outDir/readmap/denovo_$ass_name`;

	## reliable_adjs.txt
	if (-e "$outDir/RACA/$ass_name/params/reliable_adjs.txt"){	`rm -f "$outDir/RACA/$ass_name/params/reliable_adjs.txt`;	}
	`touch $outDir/RACA/$ass_name/params/reliable_adjs.txt`;

	## tree.nwk
	open TREE, "$tree";
	open TREEF, ">$outDir\/RACA/$ass_name\/params/tree.nwk";
	my $format = <TREE>;
	if ($format =~ (s/TARGET/denovo_$ass_name/)){	print TREEF "$format";	}
	close TREEF;
	close TREE;

	## Makefile
	`cp $Bin/RACA_Makefile.SFs $outDir/RACA/$ass_name/params/Makefile.SFs`;

	## raca_param.txt
	open PARAM, ">$outDir/RACA/$ass_name\/params/raca_param.txt";
	print PARAM "INSERTLIBFILE=$outDir\/RACA/$ass_name\/params/library_info.txt\n";
	print PARAM "INSERTSIZETHR=1000\n\n";
	print PARAM "READMAPPINGDIR=$outDir\/readmap/denovo_$ass_name\n";
	print PARAM "READMAPPINGLIB=$outDir\/RACA/$ass_name\/params/library_mapping.txt\n\n";
	print PARAM "NCPUS=$threads\n\n";
	print PARAM "SCFSIZEFILE=$target\.sizes\n";
	print PARAM "SCFPREFIX=scaf\n";
	print PARAM "SCFSEQFILE=$target\n\n";
	print PARAM "REFSPC=$reference\n";
	print PARAM "TARSPC=denovo_$ass_name\n\n";
	print PARAM "WINDOWSIZE=1000\n\n";
	print PARAM "OUTPUTDIR=$outDir\/RACA\/$ass_name\n";
	print PARAM "RESOLUTION=$resolution\n\n";
	print PARAM "MIN_INTRACOV_PERC=5\n";
	print PARAM "IGNORE_ADJS_WO_READS=0\n\n";
	print PARAM "TREEFILE=$outDir\/RACA/$ass_name\/params/tree.nwk\n";
	print PARAM "BENADJFILE=$outDir\/RACA/$ass_name\/params/reliable_adjs.txt\n\n";
	print PARAM "CONFIGSFSFILE=$outDir\/RACA/$ass_name\/params/config.file\n";
	print PARAM "MAKESFSFILE=$outDir\/RACA/$ass_name\/params/Makefile.SFs\n";
	close PARAM;
}


### RACA ###
sub exeRACA {
	my ($paramDir, $ass_name, $target_scf_f) = @_;
	`$raca $paramDir/raca_param.txt 2> $paramDir/../log`;
	`$Bin/RACA_create_fasta.block.pl $paramDir/../rec_chrs.refined.txt $target_scf_f $paramDir/../rec_chrs.$ass_name.fa $paramDir/../rec_chrs.$ass_name.RACA.fa $paramDir/../rec_chrs.$ass_name.unused.fa`;
	`$kent/faSize -detailed $paramDir/../rec_chrs.$ass_name.fa > $paramDir/../rec_chrs.$ass_name.fa.sizes`;
	`$kent/faSize -detailed $paramDir/../rec_chrs.$ass_name.RACA.fa > $paramDir/../rec_chrs.$ass_name.RACA.fa.sizes`;
	`$kent/faSize -detailed $paramDir/../rec_chrs.$ass_name.unused.fa > $paramDir/../rec_chrs.$ass_name.unused.fa.sizes`;
}
