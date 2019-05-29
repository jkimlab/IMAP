#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;
use Switch;
use FindBin '$Bin';
use Cwd 'abs_path';

my $seq_list = shift; # size file
my $in_bam = shift;
my $outD = shift;
my $thread = shift;

my $path_conf = abs_path("$Bin/../path.conf");
my $samtools_cmd = "";
open(PATH,$path_conf);
while(<PATH>){
	chomp;
	next if /^#/;
	next if /""/;
	my ($program, $path) = split (/=/, $_);
	switch ($program){
		case("samtools") {$samtools_cmd = $path;}
	}
}
close(F);



my @coms = ();
open F,"$seq_list";
while(<F>){
        chomp;
        my @t = split(/\s+/);   # split pattern length > 2
        my @s = split(/\|/,$t[0]);
        my $cmd = "$samtools_cmd view -bh $in_bam \"$t[0]\" > $outD/$s[0].bam && $samtools_cmd index $outD/$s[0].bam";
        push(@coms,$cmd);
}
close F;

my $FM_lo = new Parallel::ForkManager($thread);
foreach my $com(@coms){
        my $lo_qid = $FM_lo -> start($com) and next;
        print "execute==> $com\n";
        system("$com");
        $FM_lo -> finish($com);
}
$FM_lo -> wait_all_children;
