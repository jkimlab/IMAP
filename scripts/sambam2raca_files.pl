#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);

my $fname = shift;
my $f = shift;			# input SAM or BAM file
my $map_f = shift;		# output RACA mapping file
my $samtools = shift;

my ($sum, $len2sum, $cnt) = (0, 0, 0);
open(O,">$map_f");

my $fh;
my $prevline = "";
if ($fname =~ /\.bam$/) {
	open($fh, "$samtools view -h $f |");
} else {
	open($fh, "$f");
}
while(<$fh>) {
	chomp;
	if (length($_) == 0 || $_ =~ /^#/ || $_ =~ /^@/) { next; }
	my @ar = ();
	if($prevline ne ""){
		@ar = split(/\s+/,$prevline);
	} else {
		@ar = split(/\s+/);
	}
	my ($rname, $pos, $rnext, $pnext) = ($ar[2],$ar[3],$ar[6],$ar[7]);
	my $flag = $ar[1];
	my $seq = $ar[9];
	my $len = length($seq);
	my $tlen = $ar[8];

	# skip the second line
	my $secline = "";
	if($prevline ne ""){
		$secline = $_;
	} else {
		$secline = <$fh>;
		if(!$secline){last;}
	}
	my @arsecline = split(/\s+/, $secline);
	if ($ar[0] ne $arsecline[0]) {
		$prevline = $secline;
#print STDERR "warning: no read pair $ar[0]\n";

		next;
	} else {
		$prevline = "";
	}

	if ($rname eq "*" || $rnext eq "*") { next; }

	if ($rnext eq "=") { 
		$rnext = $rname; 
	}

	if ($rname eq $rnext && $pos == $pnext) { next; }

	if ($flag & 0x0002) {
		$sum += abs($tlen);
		$len2sum += (abs($tlen) * abs($tlen));
		$cnt++;
	}

	my $d = "+";
	if ($flag & 0x0010) { $d = "-"; }
	my $dnext = "+";
	if ($flag & 0x0020) { $dnext = "-"; }	
	
	print O "1\t$len\t$d\t$rname\t$pos\n";
	print O "2\t$len\t$dnext\t$rnext\t$pnext\n";
}
close($fh);
close(O);

my $mean = $sum/$cnt;
my $stdev = sqrt($len2sum/$cnt - $mean*$mean);

print "$fname\t$mean\t$mean\t$stdev\t$stdev\n";
