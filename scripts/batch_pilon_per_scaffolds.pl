#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';
use Parallel::ForkManager;
use Switch;
use FindBin '$Bin';

my $mappingD = shift;
my $targetFaD = shift;
my $outD = shift;
my $thread = shift;
my $path_conf = abs_path("$Bin/../path.conf");
$mappingD = abs_path($mappingD);

my $pilon_jar = "";
open(PATH, $path_conf);
while(<PATH>){
	chomp;
	next if /^#/;
	next if /""/;
	my ($program, $path) = split (/=/, $_);
	switch ($program) {
		case("pilon") {$pilon_jar = $path;}
	}
}
close(PATH);

`mkdir -p $outD`;

my @f = <$targetFaD/*.fa>;
my @coms = ();
for(my $i = 0; $i <= $#f; $i++){
        my $scf = basename($f[$i],".fa");
        my $cmd = "java -jar $pilon_jar --genome $targetFaD/$scf.fa --frags $mappingD/$scf.bam --output $outD/$scf";
    push(@coms,$cmd);
}

my $FM_lo = new Parallel::ForkManager($thread);
foreach my $com(@coms){
        my $lo_qid = $FM_lo -> start($com) and next;
        print "execute==> $com\n";
        system("$com");
        $FM_lo -> finish($com);
}
$FM_lo -> wait_all_children;
