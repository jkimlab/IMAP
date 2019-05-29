#!/usr/bin/env perl

use strict;
use warnings;

my $map_f = shift;	# APCF_xxx.merged.map 
my $genome_f = shift;	# assembly file
my $out_f = shift;

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

#my %hs_used_start = ();
#my %hs_used_end = ();

open(ON,">$out_f");

# read chr structure file
#my $usedracabp = 0;
my $chrcnt = 1;
my $chrseq = "";
my $pchrid = "";

my %hs_usedscfs = ();

open(F,"$map_f");
my @lines = <F>;
close(F);
chomp(@lines);

for (my $i = 0; $i < scalar(@lines); $i+=4) {
	my $ref_line = $lines[$i+1];
	my $tar_line = $lines[$i+2];

	$ref_line =~ /\S+\.(\S+):/;
	my $chrid = $1;

	#my @ar = split(/\s+/);
	#my ($chrid,$start,$end) = ($ar[0],$ar[1],$ar[2]);
	if ($pchrid ne "" && $chrid ne $pchrid) {
		print ON ">sscf_$chrcnt\n";
		for (my $p = 0; $p < length($chrseq); $p += 80) {
			my $subseq = substr($chrseq, $p, 80);
			print ON "$subseq\n";
		}	
		$chrcnt++;
		$chrseq = "";
	}

	if ($chrseq ne "") { $chrseq .= $gaps; } 

	$tar_line =~ /\S+\.(\S+):(\S+)\-(\S+) (\S+)/;
	my ($scfid, $scfstart, $scfend, $scfdir) = ($1, $2, $3, $4);

		#my ($scfid,$scfstart,$scfend,$scfdir) = ($ar[3],$ar[4],$ar[5],$ar[6]);
	$hs_usedscfs{$scfid} = 1;
	my $scffasta = $hs_seqs{$scfid};
	#if (defined($hs_used_start{$scfid})) {
	#	if ($hs_used_start{$scfid} > $scfstart) {
	#		$hs_used_start{$scfid} = $scfstart;
	#	}	
	#	if ($hs_used_end{$scfid} < $scfend) {
	#		$hs_used_end{$scfid} = $scfend;
	#	}	
	#} else {
	#	$hs_used_start{$scfid} = $scfstart;
	#	$hs_used_end{$scfid} = $scfend;
	#}
	my $scfseq = substr($scffasta,$scfstart,$scfend-$scfstart);
	if ($scfdir eq "-") {
		my $stmp = reverse($scfseq);
		$stmp =~ tr/NACGTacgt/NTGCAtgca/;		
		$scfseq = $stmp;
	}

	#my $tmpscfseq = $scfseq;
	#$tmpscfseq =~ s/[Nn]//g;	
	#$usedracabp += length($tmpscfseq);
	my $ssize = length($scfseq);
#print STDERR "added $ssize\n";
	$chrseq .= $scfseq;	
	$pchrid = $chrid;
}
close(F);
	
print ON ">sscf_$chrcnt\n";
for (my $p = 0; $p < length($chrseq); $p += 80) {
	my $subseq = substr($chrseq, $p, 80);
	print ON "$subseq\n";
}	
$chrcnt++;
$chrseq = "";

#my $cntadded = 0;
#foreach my $scfid (sort keys %hs_seqs) {
#	if (defined($hs_usedscfs{$scfid})) { next; }
#	my $seq = $hs_seqs{$scfid};
#	print ON ">$scfid\n";
#
#	for (my $p = 0; $p < length($seq); $p += 80) {
#		my $subseq = substr($seq, $p, 80);
#		print ON "$subseq\n";
#	}	
#	$cntadded++;
#}

close(ON);

#my $usedscfs = scalar(keys %hs_usedscfs);
#my $usedscfbp = 0;
#foreach my $scfid (keys %hs_usedscfs) {
#	my $seq = $hs_seqs{$scfid};
#	$usedscfbp += length($seq);
#} 
#
#my $scf_coverage = $usedscfbp/$totalbp;
#my $raca_coverage = $usedracabp/$totalbp;
#print "#total_scfs\tgenome_length\tused_scf_num\tused_scf_length\tused_scf_cov\tused_raca_length\tused_raca_cov\n";
#print "$totalscfs\t$totalbp\t$usedscfs\t$usedscfbp\t$scf_coverage\t$usedracabp\t$raca_coverage\n";
#print "$cntadded unused scaffold sequences were added\n";
