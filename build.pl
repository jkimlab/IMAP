#!/usr/bin/env perl

use strict;
use warnings;
use FindBin '$Bin';
use Getopt::Long qw(:config no_ignore_case);
use lib "$Bin/scripts/lib";
use Check::Modules;

my $install;
my $uninstall;
my $check;
my $example;

GetOptions (
	"example" => \$example,
	"install" => \$install,
	"uninstall" => \$uninstall,
	"check" => \$check,
);


if(!$example && !$check && !$install && !$uninstall){
	print STDERR "\$./build.pl --check\n";
	print STDERR "\$./build.pl --install\n";
	print STDERR "\$./build.pl --uninstall\n";
	print STDERR "\$./build.pl --example\n";
}

my $thirdparty_path = "$Bin/scripts/thirdparty";
my $source_dir = "$thirdparty_path/sources";

if($check){
	if(!check_modules()){}
	else{ print STDERR "All perl modules exists!!\n";}
} elsif($example){
	`git clone https://github.com/jkimlab/IMAP_EX.git`;
	`gunzip $Bin/IMAP_EX/S288C.fa.gz`;
	`gunzip $Bin/IMAP_EX/Saccharomyces_dairensis.fa.gz`;
	print STDERR "Download read dataset...";
	`$Bin/scripts/fastq-dump --split-files --gzip SRR1569791 -O $Bin/IMAP_EX/`;
	print STDERR "Done.\n";
} elsif($uninstall){
	`rm -f $Bin/scripts/snp_edit`;
	chdir("$thirdparty_path");
	`rm -rf *.log bwa-0.7.17 GapCloser kent lastz-distrib-1.04.00 DESCHRAMBLER RACA-0.9.1.1 SOAPdenovo2-bin-LINUX-generic-r240 SPAdes-3.7.1-Linux MaSuRCA-3.3.1 samtools-1.9`;
	chdir("$Bin");
	`rm -f $Bin/path.conf`;
} elsif($install){
`g++ $Bin/scripts/snp_edit.cpp -o $Bin/scripts/snp_edit`;
print STDERR ">> Unzip kent utilities... ";
#### Unzip Kent utilities
`tar xf $source_dir/kent.tar.gz -C $thirdparty_path`;
print STDERR "Done\n";

print STDERR ">> Unzip SPAdes... ";
#### Unzip SPAdes
`tar xf $source_dir/SPAdes-3.7.1-Linux.tar.gz -C $thirdparty_path`;
print STDERR "Done\n";

print STDERR ">> Unzip SOAPdenovo2 & GapCloser... ";
#### Unzip SOAPdenovo2 & GapCloser
`tar xf $source_dir/SOAPdenovo2-bin-LINUX-generic-r240.tar.gz -C $thirdparty_path`;
`mkdir -p $thirdparty_path/GapCloser`;
`tar xf $source_dir/GapCloser-master.tgz -C $thirdparty_path/GapCloser`;
print STDERR "Done\n";

print STDERR ">> Compile RACA... ";
#### Compile RACA
`tar xf $source_dir/RACA-0.9.1.1.tar.gz -C $thirdparty_path`;
`make -C $thirdparty_path/RACA-0.9.1.1 2> $thirdparty_path/raca.log`;
print STDERR "Done\n";

print STDERR ">> Compile DESCHRAMBLER... ";
#### Compile DESCHRAMBLER
`tar xf $source_dir/DESCHRAMBLER.tar.gz -C $thirdparty_path`;
`make -C $thirdparty_path/DESCHRAMBLER 2> $thirdparty_path/meta.log`;
print STDERR "Done\n";

print STDERR ">> Compile BWA... ";
#### Compile BWA
`tar xf $source_dir/bwa-0.7.17.tar.bz2 -C $thirdparty_path`;
`make -C $thirdparty_path/bwa-0.7.17 2> $thirdparty_path/bwa.log`;
if(-f "$thirdparty_path/bwa-0.7.17/bwa"){
	print STDERR "Done\n";
} else {
	print STDERR "ERROR\n";
	exit(1);
}

print STDERR ">> Compile lastz... ";
#### Compile lastz
`tar xf $source_dir/lastz-distrib-1.04.00.tar.gz -C $thirdparty_path`;
`make -C $thirdparty_path/lastz-distrib-1.04.00/src 2> $thirdparty_path/lastz.log`;
if(-f "$thirdparty_path/lastz-distrib-1.04.00/src/lastz"){
	print STDERR "Done\n";
} else {
	print STDERR "ERROR\n";
	exit(1);
}

print STDERR ">> Compile samtools... ";
#### Compile samtools
`tar xf $source_dir/samtools-1.9.tar.bz2 -C $thirdparty_path`;
`make -C $thirdparty_path/samtools-1.9 2> $thirdparty_path/samtools.log`;
if(-f "$thirdparty_path/samtools-1.9/samtools"){
	print STDERR "Done\n";
} else {
	print STDERR "ERROR\n";
	exit(1);
}

print STDERR ">> Compile MaSuRCA... ";
#### Compile MaSuRCA
`tar xf $source_dir/MaSuRCA-3.3.1.tar.gz -C $thirdparty_path`;
chdir("$thirdparty_path/MaSuRCA-3.3.1");
`./install.sh 2> $thirdparty_path/masurca.log`;
chdir("$Bin");
if(-f "$thirdparty_path/MaSuRCA-3.3.1/bin/masurca"){
	print STDERR "Done\n";
} else {
	print STDERR "ERROR\n";
	exit(1);
}

open(W,">$Bin/path.conf");
print W "bwa=$thirdparty_path/bwa-0.7.17/bwa\n";
print W "lastz=$thirdparty_path/lastz-distrib-1.04.00/src/lastz\n";
print W "masurca=$thirdparty_path/MaSuRCA-3.3.1/bin/masurca\n";
print W "spades=$thirdparty_path/SPAdes-3.7.1-Linux/bin/spades.py\n";
print W "soapdenovo2_127=$thirdparty_path/SOAPdenovo2-bin-LINUX-generic-r240/SOAPdenovo-127mer\n";
print W "soapdenovo2_63=$thirdparty_path/SOAPdenovo2-bin-LINUX-generic-r240/SOAPdenovo-63mer\n";
print W "gapcloser=$thirdparty_path/GapCloser/GapCloser\n";
print W "RACA=$thirdparty_path/RACA-0.9.1.1/Run_RACA.pl\n";
print W "DESCHRAMBLER=$thirdparty_path/DESCHRAMBLER/DESCHRAMBLER.pl\n";
print W "pilon=$thirdparty_path/pilon-1.22.jar\n";
print W "gatk=$thirdparty_path/GenomeAnalysisTK.jar\n";
print W "picard=$thirdparty_path/picard.jar\n";
print W "samtools=$thirdparty_path/samtools-1.9/samtools\n";
print W "kent=$thirdparty_path/kent\n";
close(W);
}
