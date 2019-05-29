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
my $iterNum=0;			## Number of iteraction
my $outDir="./";		## Output directory
my $path_conf = abs_path("$Bin/../path.conf");

my $metaFa;
my @p_reads = ();
my @m_reads = ();

my $samtools;			## path of samtools
my $bwa;				## path of bwa
my $pilon;				## path of pilon
my $gatk;				## path of gatk
my $kent;				## path of kent

my $options = GetOptions (
			"threads|t=i" => \$threads,
			"iterNum|i=i" => \$iterNum,
			"metaFa|m=s" => \$metaFa,
			"param|p=s" => \$params,
			"outDir|o=s" => \$outDir,
			"help|h" => \$help,
		);

sub HELP{
	my $src = basename($0);
	print "\nUsage: $src [options]\n\n";
	print "\t-t <integer> Number of threads (default: 1)\n";
	print "\t-i <integer> Number of iteration (default: 1)\n";
	print "\t-m <filename> Meta-assembly file\n";
	print "\t-p <filename> Parameter file\n";
	print "\t-o <filename> Results of output directory (default: ./)\n";
	print "\t-h|help print this page\n\n";
	exit;
}

### Help Prints ###
if (!-e $path_conf){	print STDERR "[ERROR] Not exists 'path.conf' file\n";	}
if (!defined $params){	print STDERR "[ERROR] Please check path of parameter file\n";	HELP();	}

$outDir = abs_path("$outDir");

open PARAM, "$params";
while (<PARAM>){
	chomp;
	if ($_ =~ /^#/ || $_ eq "" || $_ =~ /^\[/){	next;	}
	my @split = split (/\s+/, $_);
	if ($split[0] eq "p1" || $split[0] eq "p2"){	push (@p_reads, "$split[0]\t$split[1]");	}
	elsif ($split[0] eq "m1" || $split[1] eq "m2"){	push (@m_reads, "$split[0]\t$split[1]");	}
}
close PARAM;


print STDERR "### 4. Step of error correction ###\n";
for (my $i=1; $i<=$iterNum; $i++){
	`mkdir -p $outDir/iter_$i`;
	print STDERR "    ($i) Iteration of error correction\n";
	open INPUT, ">$outDir/iter_$i/input.txt";
	if ($i == 1){	print INPUT "fasta\t$metaFa\n";	}
	else {	my $l=$i-1;	print INPUT "fasta\t$outDir/iter_$l/finalEdit/meta_assembled_scaf.pilon.n50.edit.norm.fa\n";	}
	if ($#p_reads > 0){	
		for (my $j=0; $j<=$#p_reads; $j++){	print INPUT "$p_reads[$j]\n";	}
	}
	if ($#m_reads > 0){
		for (my $j=0; $j<=$#m_reads; $j++){	print INPUT "$m_reads[$j]\n";	}
	}
	close INPUT;

	`bash $Bin/error_correction.sh -t $threads -p $path_conf -o $outDir/iter_$i -i $outDir/iter_$i/input.txt -s $Bin`;
	if($i > 1){
		`ln -s $outDir/iter_$i/finalEdit/meta_assembled_scaf.pilon.n50.edit.norm.pilon.n50.edit.norm.fa $outDir/iter_$i/finalEdit/meta_assembled_scaf.pilon.n50.edit.norm.fa`;
	}
}
