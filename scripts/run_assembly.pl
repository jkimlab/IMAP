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
my $outDir="./";		## Output directory
my $path_conf = abs_path("$Bin/../path.conf");

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
if (!-e $path_conf){	print STDERR "[ERROR] Not exists 'path.conf' file\n";	}
if (!defined $params){
	print STDERR "[ERROR] Please check path of parameter file\n";
	HELP();
}

$outDir = abs_path("$outDir");
`mkdir -p $outDir/denovo_assemblies`;


### Parameter ###
my $kmer;				## Number of kmers
my $minContigLength;	## Minimum length of contigs
my $readLength;			## Length of reads
my $n=0;				## Check if a file exists

my $masurca;			## path of masurca
my $spades;				## path of spades
my $soap2_127;			## path of soapdenovo2-127mer
my $soap2_63;			## path of soapdenovo2-63mer
my $gapcloser;			## path of gapcloser
my $kent;				## path of kent directory

my $lib=0;
my @insertSize = ();	## Average of Insert sizes
my @insertSizeSD = ();	## Standard deviation of Insert sizes

my %paired = ();		## $paired{$lib}{1}	/	$paired{$lib}{2}
my %mate = ();			## $mate{$lib}{1}	/	$mate{$lib}{2}


open PATH, "$path_conf";
while (<PATH>){
	chomp;
	next if /^#/;
	next if /""/;
	my ($program, $path) = split (/=/, $_);
	switch ($program){
		case("masurca"){	$masurca = $path;	}
		case("spades"){	$spades = $path;	}
		case("soapdenovo2_127"){	$soap2_127 = $path;	}
		case("soapdenovo2_63"){	$soap2_63 = $path;	}
		case("gapcloser"){	$gapcloser = $path;	}
		case("kent"){	$kent = $path;	}
	}
}
close PATH;

open PARA, "$params";
while (<PARA>){
	chomp;
	next if /^#/;
	next if /""/;
	my ($param, $val) = split (/\s+/, $_);
	switch ($param){
		case ("MinContigLength"){	$minContigLength = $val;	}
		case ("MaxReadLength"){	$readLength = $val;	}
		case ("Kmer"){	$kmer = $val;	}
		case ("[LIB]"){	$lib++;	}
		case ("p1"){
			if (-e $val){	push (@{$paired{$lib}{1}}, $val);	$n=0;	}
			else {	die "[ERROR] Please check the path of paired-end read(s)\n";	}
		}
		case ("p2"){
			if (-e $val){	push (@{$paired{$lib}{2}}, $val);	}
			else {	die "[ERROR] Please check the path of paired-end read(s)\n";	}
		}
		case ("m1"){
			if (-e $val){	push (@{$mate{$lib}{1}}, $val);	$n=0;	}
			else {	die "[ERROR] Please check the path of mate pair read(s)\n";	}
		}
		case ("m2"){
			if (-e $val){	push (@{$mate{$lib}{2}}, $val);	}
			else {	die "[ERROR] Please check the path of mate pair read(s)\n";	}
		}
		case ("insertSize"){	
			if ($val eq "" || $val < 0){	die "[ERROR] Please check insert size\n";	}
			else {	push (@insertSize, $val);	}
		}
		case ("insertSizeSD"){	
			if ($val eq "" || $val == 0){	$val = int($insertSize[-1]) * 0.15;	}
			elsif ($val < 0){	die "[ERROR] Please check stdev of insert size\n";	}
			push (@insertSizeSD, $val);
		}
	}
}
close PARA;


### Running the Assemblies ###
## SPAdes ##
print STDERR "### 1. Step of de novo assembly ###\n";
print STDERR "    (1) SPAdes assembly\n";
`mkdir -p $outDir/denovo_assemblies/spades`;
chdir ("$outDir/denovo_assemblies/spades");
print STDERR "       (1-1) Making the configuration file for SPAdes assembly\n";
makeSPACMD();
print STDERR "       (1-2) Starting the assembly\n";
if (-e "$outDir/denovo_assemblies/spades/log"){	`rm -f "$outDir/denovo_assemblies/spades/log"`;	}
exeSpades();
changeID("$outDir/denovo_assemblies/spades/scaffolds.fasta", "$outDir/denovo_assemblies/denovo_spades.fa");
`$kent/faSize -detailed $outDir/denovo_assemblies/denovo_spades.fa > $outDir/denovo_assemblies/denovo_spades.fa.sizes`;

## MaSuRCA ##
print STDERR "    (2) MaSuRCA assembly\n";
`mkdir -p $outDir/denovo_assemblies/masurca`;
chdir ("$outDir/denovo_assemblies/masurca");
print STDERR "       (2-1) Making the configuration file for MaSuRCA assembly\n";
makeMasCMD();
print STDERR "       (2-2) Starting the assembly\n";
if (-e "$outDir/denovo_assemblies/masurca/log"){	`rm -f "$outDir/denovo_assemblies/masurca/log"`;	}
exeMasurca();
changeID("$outDir/denovo_assemblies/masurca/CA/final.genome.scf.fasta", "$outDir/denovo_assemblies/denovo_masurca.fa");
`$kent/faSize -detailed $outDir/denovo_assemblies/denovo_masurca.fa > $outDir/denovo_assemblies/denovo_masurca.fa.sizes`;

## SOAPdenovo2 ##
print STDERR "    (3) SOAPdenovo2 assembly\n";
`mkdir -p $outDir/denovo_assemblies/soap2`;
chdir ("$outDir/denovo_assemblies/soap2");
print STDERR "       (3-1) Making the configuration file for SOAPdenovo2 assembly\n";
makeSOAPCMD();
print STDERR "       (3-2) Starting the assembly\n";
if (-e "$outDir/denovo_assemblies/soap2/log"){	`rm -f "$outDir/denovo_assemblies/soap2/log"`;	}
exeSoapdenovo2();
changeID("$outDir/denovo_assemblies/soap2/soap2.gc.scafSeq", "$outDir/denovo_assemblies/denovo_soap2.fa");
`$kent/faSize -detailed $outDir/denovo_assemblies/denovo_soap2.fa > $outDir/denovo_assemblies/denovo_soap2.fa.sizes`;


#=======================================================================================================#
#												Functions												#
#=======================================================================================================#

### cmd of SPAdes ###
sub makeSPACMD {
	open SPACMD, ">./spades.yaml";
	print SPACMD "[\n";
	for (my $i=0; $i<=$#insertSize; $i++){
		my $l = $i+1;
		if (exists $paired{$l}{1}){
			my @paired_1 = @{$paired{$l}{1}};
			my @paired_2 = @{$paired{$l}{2}};
			
			if ($#paired_1 != -1 || $#paired_2 != -1){	
				print SPACMD " {\n";
				print SPACMD "  orientation: \"fr\",\n";
				print SPACMD "  type: \"paired-end\",\n";
				print SPACMD "  right reads: [\n";
				print SPACMD "  $paired_1[0]";
				if ($#paired_1 != 0){
					for (my $i=1; $i<=$#paired_1; $i++){	print SPACMD ",\n  $paired_1[$i]";	}
				}
				print SPACMD "\n  ],\n";
				print SPACMD "  left reads: [\n";
				print SPACMD "  $paired_2[0]";
				if ($#paired_2 != 0){
					for (my $i=1; $i<=$#paired_2; $i++){	print SPACMD ",\n  $paired_2[$i]";	}
				}
				if ($l == $#paired_2){	print SPACMD "\n  ]\n }\n";	}
				else {	print SPACMD "\n  ]\n },\n";	}
			}
		}
		if (exists $mate{$l}{1}){
			my @mate_1 = @{$mate{$l}{1}};
			my @mate_2 = @{$mate{$l}{2}};
			if ($#mate_1 != -1 || $#mate_2 != -1){
				print SPACMD " {\n";
				print SPACMD "  orientation: \"rf\",\n";
				print SPACMD "  type: \"mate-pairs\",\n";
				print SPACMD "  right reads: [\n";
				print SPACMD "  $mate_1[0]";
				if ($#mate_1 != 0){
					for (my $i=1; $i<=$#mate_1; $i++){	print SPACMD ",\n  $mate_1[$i]";	}
				}
				print SPACMD "\n  ],\n";
				print SPACMD "  left reads: [\n";
				print SPACMD "  $mate_2[0]";
				if ($#mate_2 != 0){
					for (my $i=1; $i<=$#mate_2; $i++){	print SPACMD ",\n  $mate_2[$i]";	}
				}
				if ($l == $#mate_2){	print SPACMD "\n  ]\n }\n";	}
				else {	print SPACMD "\n  ]\n },\n";	}
			}
		}
	}
	print SPACMD "]";
	close SPACMD;
}

### cmd of MaSuRCA ###
sub makeMasCMD {
	open MasCMD, ">./masurca.cmd";
	print MasCMD "DATA\n";
	for (my $i=0; $i<=$#insertSize; $i++){
		my $l = $i+1;
		if (exists $paired{$l}{1}){
			print MasCMD "PE= p$l $insertSize[$i] $insertSizeSD[$i]";
			my @paired_1 = @{$paired{$l}{1}};
			my @paired_2 = @{$paired{$l}{2}};
			for (my $j=0; $j<=$#paired_1; $j++){	print MasCMD " $paired_1[$j] $paired_2[$j]";	}
		}
		if (exists $mate{$l}{1}){
			print MasCMD "JUMP= p$l $insertSize[$i] $insertSizeSD[$i]";
			my @mate_1 = @{$mate{$l}{1}};
			my @mate_2 = @{$mate{$l}{2}};
			for (my $j=0; $j<=$#mate_1; $j++){	print MasCMD " $mate_1[$j] $mate_2[$j]";	}
		}
		print MasCMD "\n";
	}
	print MasCMD "END\nPARAMETERS\n";
	print MasCMD "GRAPH_KMER_SIZE=$kmer\n";
	print MasCMD "NUM_THREADS=$threads\n";
	print MasCMD "JF_SIZE=130000000\n";
	print MasCMD "END\n";
	close MasCMD;
}

### cmd of SOAPdenovo2 ###
sub makeSOAPCMD {
	open SOAPCMD, ">./soap2.cmd";
	print SOAPCMD "max_rd_length=$readLength\n";
	for (my $i=0; $i<=$#insertSize; $i++){
		print SOAPCMD "[LIB]\n";
		print SOAPCMD "avg_ins=$insertSize[$i]\n";
		my $l = $i+1;
		if (exists $paired{$l}{1}){
			my @paired_1 = @{$paired{$l}{1}};
			my @paired_2 = @{$paired{$l}{2}};
			if ($#paired_1 != -1 || $#paired_2 != -1){
				print SOAPCMD "reverse_seq=0\n";
				print SOAPCMD "asm_flags=3\n";
				print SOAPCMD "rank=$l\n";
				for (my $i=0; $i<=$#paired_1; $i++){	print SOAPCMD "q1=$paired_1[$i]\nq2=$paired_2[$i]\n";	}
			}
		}
		if (exists $mate{$l}{1}){
			my @mate_1 = @{$mate{$l}{1}};
			my @mate_2 = @{$mate{$l}{2}};
				
			if ($#mate_1 != -1 || $#mate_2 != -1){
				print SOAPCMD "reverse_seq=1\n";
				print SOAPCMD "asm_flags=3\n";
				print SOAPCMD "rank=$l\n";
				for (my $i=0; $i<=$#mate_1; $i++){	print SOAPCMD "q1=$mate_1[$i]\nq2=$mate_2[$i]\n";	}
			}
		}
	}
	close SOAPCMD;
}


### SPAdes ###
sub exeSpades {
	`$spades --careful --dataset ./spades.yaml -k $kmer -t $threads -o . 2> log`;
}

### MaSuRCA ###
sub exeMasurca {
	`$masurca ./masurca.cmd 2> log`;
	`./assemble.sh 2> log`;
}

### SOAPdenovo2 ###
sub exeSoapdenovo2 {
	if ($kmer > 63){
		`$soap2_127 all -K $kmer -F -R -E -w -u -p $threads -L $readLength -s ./soap2.cmd -o ./soap2 2> log`;
		`$gapcloser -b soap2.cmd -a soap2.scafSeq -o ./soap2.gc.scafSeq -t $threads 2> log`;
	}
	else {
		`$soap2_63 all -K $kmer -F -R -E -w -u -p $threads -L $readLength -s ./soap2.cmd -o ./soap2 2> log`;
		`$gapcloser -b soap2.cmd -a soap2.scafSeq -o ./soap2.gc.scafSeq -t $threads 2> log`;
	}
}

sub changeID {
	my ($fasta_F, $out_fasta_F) = @_;
	my $id=1;
	open(F,"$fasta_F");
	open(W,">$out_fasta_F");
	while(<F>){
		chomp;
		if($_ =~ /^>/){
			print W ">scaf$id\n";
			$id++;
		} else {
			print W "$_\n";
		}
	}
	close(W);
	close(F);
}
