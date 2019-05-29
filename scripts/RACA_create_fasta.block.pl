#!/usr/bin/env perl

use strict;
use warnings;

my $str_f = shift;	# rec_chrs.refined.txt
my $genome_f = shift;	# genome.scf.fasta
my $total_f = shift;	# new: pseudo_chr (RACA) + unused_scf
my $pseudo_f = shift;	# new: pseudo_chr (RACA)
my $unused_f = shift;	# new: unused_scf (RACA)


my $gaps = "N"x100;

# read genome sequences
my $totalbp = 0;
my $totalscfs = 0;
my %hs_seqs = ();
my $spc = "";
my $seq = "";
open(F,"$genome_f");
while(<F>) {
	chomp;
	if ($_ =~ /^>/) {
		$totalscfs++;
		if (length($seq) > 0) {
			$hs_seqs{$spc} = $seq;
			my $tmpseq = $seq;
			$tmpseq =~ s/[Nn]//g;
			$totalbp += length($tmpseq);	
		}

		$seq = "";
		my @ars = split(/\s+/);
		$spc = substr($ars[0],1);
	} else {
		$seq .= $_;
	}
}
close(F);

if (length($seq) > 0) {
	$hs_seqs{$spc} = $seq;
	my $tmpseq = $seq;
	$tmpseq =~ s/[Nn]//g;
	$totalbp += length($tmpseq);	
}

my %hs_used_start = ();
my %hs_used_end = ();

open(TOTALF,">$total_f");
open(PSEUDOF, ">$pseudo_f");
open(UNUSEDF, ">$unused_f");

# read chr structure file
my $usedracabp = 0;
my $chrcnt = 1;
my $chrseq = "";
my $pchrid = "";

my %hs_usedscfs = ();

open(F,"$str_f");
while(<F>) {
	chomp;
	my @ar = split(/\s+/);
	my ($chrid,$start,$end) = ($ar[0],$ar[1],$ar[2]);
	if ($pchrid ne "" && $chrid ne $pchrid) {
		print TOTALF ">RACA_$chrcnt\n";
		print PSEUDOF ">RACA_$chrcnt\n";
		for (my $p = 0; $p < length($chrseq); $p += 80) {
			my $subseq = substr($chrseq, $p, 80);
			print TOTALF "$subseq\n";
			print PSEUDOF "$subseq\n";
		}	
		$chrcnt++;
		$chrseq = "";
	}

	if ($ar[3] eq "GAPS") {
		$chrseq .= $gaps;
	} else {
my $csize = $end - $start;
my $ssize = $ar[5] - $ar[4];
#if ($csize != $ssize) { print STDERR "DIFF $csize vs $ssize\n"; }

		my ($scfid,$scfstart,$scfend,$scfdir) = ($ar[3],$ar[4],$ar[5],$ar[6]);
		$hs_usedscfs{$scfid} = 1;
		my $scffasta = $hs_seqs{$scfid};
		if (defined($hs_used_start{$scfid})) {
			if ($hs_used_start{$scfid} > $scfstart) {
				$hs_used_start{$scfid} = $scfstart;
			}	
			if ($hs_used_end{$scfid} < $scfend) {
				$hs_used_end{$scfid} = $scfend;
			}	
		} else {
			$hs_used_start{$scfid} = $scfstart;
			$hs_used_end{$scfid} = $scfend;
		}
		my $scfseq = substr($scffasta,$scfstart,$scfend-$scfstart);
		if ($scfdir eq "-") {
			my $stmp = reverse($scfseq);
			$stmp =~ tr/NACGTacgt/NTGCAtgca/;		
			$scfseq = $stmp;
		}

		my $tmpscfseq = $scfseq;
		$tmpscfseq =~ s/[Nn]//g;	
		$usedracabp += length($tmpscfseq);
		$chrseq .= $scfseq;	
	}	

	$pchrid = $chrid;
}
close(F);
	
print TOTALF ">RACA_$chrcnt\n";
print PSEUDOF ">RACA_$chrcnt\n";
for (my $p = 0; $p < length($chrseq); $p += 80) {
	my $subseq = substr($chrseq, $p, 80);
	print TOTALF "$subseq\n";
	print PSEUDOF "$subseq\n";
}	
$chrcnt++;
$chrseq = "";

my $cntadded = 0;
foreach my $scfid (sort keys %hs_seqs) {
	if (defined($hs_usedscfs{$scfid})) { next; }
	my $seq = $hs_seqs{$scfid};
	print TOTALF ">$scfid\n";
	print UNUSEDF ">$scfid\n";

	for (my $p = 0; $p < length($seq); $p += 80) {
		my $subseq = substr($seq, $p, 80);
		print TOTALF "$subseq\n";
		print UNUSEDF "$subseq\n";
	}	
	$cntadded++;
}

close(TOTALF);
close(PSEUDOF);
close(UNUSEDF);
