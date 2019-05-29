#!/usr/bin/env perl

use strict;
use warnings;
use Cwd "abs_path";
use File::Basename;

my $fasta_f = abs_path(shift);
my $out_fa = shift;

my $num = 0;
open(W,">$out_fa");
open(F,$fasta_f);
while(<F>){
	chomp;
	if($_ =~ /^>/){
		$num++;
		print W ">seq$num\n";
	} else {
		my $uc = uc($_);
		print W "$uc\n";
	}
}
close(F);
