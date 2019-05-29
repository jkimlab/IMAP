#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd 'abs_path';
use FindBin '$Bin';
use Bio::TreeIO;
use Switch;

my $param_f;
my $threads = 1;
my $out_dir = "./IMAP_RESULT";
my $help    = 0;

### Program path information
my $path_conf = "$Bin/path.conf";
my $assembly_cmd = "$Bin/scripts/run_assembly.pl";
my $raca_cmd = "$Bin/scripts/run_raca.pl";
my $meta_cmd = "$Bin/scripts/run_deschrambler.pl";
my $errCorr_cmd = "$Bin/scripts/iter_errorCorrection.pl";

### Parameters
my %params = ();
my $n      = 0;
my $lib    = 0;
my $outgroup_num = 0;
my $ref_name = "";
my %paired = ();
my %mate   = ();
my $iterNum = 1;

GetOptions(
	"param|p=s"   => \$param_f,
	"threads|t=i" => \$threads,
	"outdir|o=s"  => \$out_dir,
	"help|h"      => \$help,
);

if ( $help == 1 || !$param_f ) {
	print STDERR "\n=====================================================\n\n";
	print STDERR "IMAP (Interactive Meta-Assembly Pipeline\n";
	print STDERR "Version 1.0.0\n\n";
	print STDERR "Usage: IMAP -t [threads] -p [parameter file] -o [out directory]\n\n";
	print STDERR "Options:\n";
	print STDERR "\t--threads|-t\t<integer> Number of threads (defulat: 1)\n";
    print STDERR "\t--params|-p\t<filename> Parameter file\n";
    print STDERR "\t--outdir|-o\t<filename> Output directory (default: ./IMAP_RESULT)\n";
    print STDERR "\t--help|-h\tPrint usages\n";
	print STDERR "\n=====================================================\n\n";
	exit;
}

$out_dir = abs_path("$out_dir");

### Validating path configures
open( F, $path_conf )
  or exit print STDERR "Failed to read path configure file ($path_conf)\n";
while (<F>) {
	chomp;
	if ( $_ =~ /^#/ || $_ eq "" ) { next; }
	my @split_arr = split(/=/);
	if ( !$split_arr[0] ne "kent" ) {
		if ( !-e $split_arr[1] ) { exit print STDERR "Required: $split_arr[0] ($split_arr[1])\n";
		}
	} else {
		my @kent_required = (
			"axtChain",   "chainMergeSort",  "chainSwap", "faSplit",
			"lavToPsl",   "chainAntiRepeat", "chainNet",  "faSize",
			"faToTwoBit", "netSyntenic"
		);
		foreach my $k (@kent_required) {
			if ( !-e $split_arr[1] ) {
				exit print STDERR "Required: $split_arr[0] ($split_arr[1])\n";
			}
		}
	}
}
close(F);

### Parsing & validating parameter file
my $matepair_read_flag = 0;
open( PARA, "$param_f" );
while (<PARA>) {
	chomp;
	if ( $_ =~ /^#/ ) { next; }
	my @param_arr = split(/\s+/);
	switch ($param_arr[0]) {
		case ("MinContigLength") { $params{"MinContigLength"} = $param_arr[1]; }
		case ("MaxReadLength")   { $params{"MaxReadLength"}   = $param_arr[1]; }
		case ("Kmer")            { $params{"Kmer"}            = $param_arr[1]; }
		case ("[LIB]")           { $lib++; }
		case ("p1") {
			if ( -e $param_arr[1] ) { push(@{$params{"libp"}{$lib}{1}},abs_path($param_arr[1])); }
			else           { exit print STDERR "Not existing ($param_arr[1])\n"; }
		}
		case ("p2") {

			if ( -e $param_arr[1] ) { push(@{$params{"libp"}{$lib}{2}},abs_path($param_arr[1])); }
			else           { exit print STDERR "Not existing ($param_arr[1])\n"; }

		}
		case ("m1") {
			if ( -e $param_arr[1] ) { push(@{$params{"libm"}{$lib}{1}},abs_path($param_arr[1])); $matepair_read_flag = 1;}
			else           { exit print STDERR "Not existing ($param_arr[1])\n"; }
		}
		case ("m2") {

			if ( -e $param_arr[1] ) { push(@{$params{"libm"}{$lib}{2}},abs_path($param_arr[1])); }
			else           { exit print STDERR "Not existing ($param_arr[1])\n"; }

		}
		case ("insertSize")   { $params{"lib"}{$lib}{"insertSize"}   = $param_arr[1]; }
		case ("insertSizeSD") { $params{"lib"}{$lib}{"insertSizeSD"} = $param_arr[1]; }
		case ("Reference") { 
			if ( -e $param_arr[2] ) {
				$ref_name = $param_arr[1];
				$params{"Reference"}{"Name"} = $param_arr[1];
				$params{"Reference"}{"Seq"} = abs_path($param_arr[2]);
			} else { exit print STDERR "Not existing ($param_arr[2])\n"; }
		}
		case ("Outgroup") { 
			if ( -e $param_arr[2] ) {
				$outgroup_num++;
				$params{"Outgroup"}{$outgroup_num}{"Name"} = $param_arr[1];
				$params{"Outgroup"}{$outgroup_num}{"Seq"} = abs_path($param_arr[2]);
			} else { exit print STDERR "Not existing ($param_arr[2])\n"; }
		}
		case ("Resolution") { $params{"Resolution"} = $param_arr[1]; }
		case ("TREE") { 
			if ( -e $param_arr[1] ) {
				my $tree = `cat $param_arr[1]`;
				chomp($tree);
				if($tree =~ /\((TARGET)\@:\d+\.*\d*,(\w+):\d+\.*\d*\)/){
					if($1 ne "TARGET"){
						exit print STDERR "Wrong formatted tree! ($param_arr[1]): TARGET@\n";
					}

					if($2 ne $ref_name){
						print STDERR "$2\t$ref_name\n";
						exit print STDERR "Wrong formatted tree! ($param_arr[1]): Reference\n";
					}
				} else {
					exit print STDERR "Wrong formatted tree! ($param_arr[1])\n";
				}
				$params{"TREE"} = abs_path($param_arr[1]); 
			} else { exit print STDERR "Not existing ($param_arr[1])\n"; }
		}
		case ("IterationNumber"){	$iterNum = $param_arr[1];	}
	}
}
close(PARA);

### Building parameter file
`mkdir -p $out_dir`;
my $running_param_f = "$out_dir/params.txt";
open( W, ">$out_dir/params.txt" );
print W "Kmer\t$params{'Kmer'}\n";
print W "MaxReadLength\t$params{'MaxReadLength'}\n";
print W "MinContigLength\t$params{'MinContigLength'}\n";
foreach my $lib_num (sort{$a <=> $b} keys %{$params{'libm'}}){
	print W "\n[LIB]\n";
	print W "insertSize\t$params{'lib'}{$lib_num}{'insertSize'}\n";
	print W "insertSizeSD\t$params{'lib'}{$lib_num}{'insertSizeSD'}\n";
	my @m1 = @{$params{'libm'}{$lib_num}{1}};
	my @m2 = @{$params{'libm'}{$lib_num}{2}};
	for(my $i = 0; $i <= $#m1;$i++){
		print W "m1\t$m1[$i]\nm2\t$m2[$i]\n";	
	}
}
foreach my $lib_num (sort{$a <=> $b} keys %{$params{'libp'}}){
	print W "\n[LIB]\n";
	print W "insertSize\t$params{'lib'}{$lib_num}{'insertSize'}\n";
	print W "insertSizeSD\t$params{'lib'}{$lib_num}{'insertSizeSD'}\n";
	my @p1 = @{$params{'libp'}{$lib_num}{1}};
	my @p2 = @{$params{'libp'}{$lib_num}{2}};
	for(my $i = 0; $i <= $#p1;$i++){
		print W "p1\t$p1[$i]\np2\t$p2[$i]\n";
	}
}
print W "\nReference\t$params{'Reference'}{'Name'}\t$params{'Reference'}{'Seq'}\n";
foreach my $o_num (sort{$a <=> $b} keys %{$params{'Outgroup'}}){
	print W "Outgroup\t$params{'Outgroup'}{$o_num}{'Name'}\t$params{'Outgroup'}{$o_num}{'Seq'}\n";
}
print W "TREE\t$params{'TREE'}\n";
print W "Resolution\t$params{'Resolution'}\n\n";
print W "IterationNumber\t$iterNum\n";
close(W);



print STDERR "\n";
print STDERR "===============   Parameters used in IMAP   ===============\n";
print STDERR "[General parameters]\n";
print STDERR "Threads:\t$threads\n";
print STDERR "Reference:\t$params{\"Reference\"}{\"Name\"}\t$params{\"Reference\"}{\"Seq\"}\n";
print STDERR "Outgroup:\t$params{\"Outgroup\"}{1}{\"Name\"}\t$params{\"Outgroup\"}{1}{\"Seq\"}\n";

if ($outgroup_num > 1){
	for (my $i=2; $i<=$outgroup_num; $i++){
		print STDERR "\t\t$params{\"Outgroup\"}{$i}{\"Name\"}\t$params{\"Outgroup\"}{$i}{\"Seq\"}\n";
	}
}

for (my $i=1; $i<=$lib; $i++){
	print STDERR "\n[lib_$i]\n";

	my @paired_1 = ();	my @paired_2 = ();

	if (exists $params{"libp"}{$i}{1} && exists $params{"libp"}{$i}{2}){
		@paired_1 = @{$params{"libp"}{$i}{1}};
		@paired_2 = @{$params{"libp"}{$i}{2}};
	}

	if ($#paired_1 != -1 && $#paired_2 != -1){
		print STDERR "Target (paired-end):\t$paired_1[0]\n";
		print STDERR "\t\t\t$paired_2[0]\n";
		if ($#paired_1 > 0 && $#paired_2 > 0){
			for (my $j=1; $j<=$#paired_1; $j++){
				print STDERR "\t\t\t$paired_1[$j]\n";
				print STDERR "\t\t\t$paired_2[$j]\n";
			}
		}
	}

	if($matepair_read_flag == 1){
		my @mate_1 = ();	my @mate_2 = ();
		if (exists $params{"libm"}{$i}{1} && exists $params{"libm"}{$i}{2}){
			@mate_1 = @{$params{"libm"}{$i}{1}};
			@mate_2 = @{$params{"libm"}{$i}{2}};
		}
	
		if ($#mate_1 != -1 && $#mate_2 != -1){
			print STDERR "Target (mate-pair):\t$mate_1[0]\n";
			print STDERR "\t\t\t$mate_2[0]\n";
			if ($#mate_1 > 0 && $#mate_2 > 0){
				for (my $j=1; $j<=$#mate_1; $j++){
					print STDERR "\t\t\t$mate_1[$j]\n";
					print STDERR "\t\t\t$mate_2[$j]\n";
				}
			}
		}
	}
}
print STDERR "\n";
print STDERR "[De novo assembly]\n";
print STDERR "kmer:\t$params{\"Kmer\"}\n\n";
print STDERR "[RACA & meta-assembly]\n";
print STDERR "Synteny resolution:\t$params{\"Resolution\"}\n\n";
print STDERR "[Error correction]\n";
print STDERR "Iteration number:\t$iterNum\n";
print STDERR "===========================================================\n\n\n";

### Running denovo assembly
`$assembly_cmd -t $threads -p $running_param_f -o $out_dir`;

### Running RACA
`$raca_cmd -t $threads -p $running_param_f -o $out_dir`;

### Running meta-assembly
`$meta_cmd -t $threads -p $running_param_f -o $out_dir`;

### Correcting errors
`mkdir -p $out_dir/errorCorrection`;
`$errCorr_cmd -t $threads -i $iterNum -m $out_dir/meta_assembled_scaf.fa -p $running_param_f -o $out_dir/errorCorrection`;

`cp $out_dir/errorCorrection/iter_$iterNum/meta_assembled_scaf.fa $out_dir/final_assembly.fa`;

print STDERR "Done.\n";
