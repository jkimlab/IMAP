#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use Cwd;
use Cwd 'abs_path';
use File::Basename;
use Parallel::ForkManager;
use Getopt::Long qw(:config no_ignore_case);
use FindBin '$Bin';
use Switch;

my $out_dir;
my $core = 1;
my $resolution;
my $ref_fa;
my $tar_fa;
my $path_conf = abs_path("$Bin/../path.conf");

GetOptions (
		"core|p=i" => \$core,
		"res|m=i" => \$resolution,
		"ref|r=s" => \$ref_fa,
		"tar|t=s" => \$tar_fa,
        "outdir|o=s" => \$out_dir,
);

### ForkManager ###
my $pm = new Parallel::ForkManager($core);

### PATH  ###
my $lastz;
my $kent_src;
$ref_fa = abs_path($ref_fa);
$tar_fa = abs_path($tar_fa);
$out_dir = abs_path($out_dir);
my $bref_name = basename($ref_fa,".fa");
my $btar_name = basename($tar_fa,".fa");
my $chainNet_dir = "$out_dir/$bref_name/$btar_name";
$out_dir = "$out_dir/building/$bref_name.$btar_name";
`mkdir -p $out_dir`;
### Alignment parameter ###
my $chainPar = "-minScore=5000 -linearGap=$Bin/lastz_params/near.linearGap";
my $lastzPar = "E=150 H=2000 K=4500 L=2200 M=254 O=600 Q=$Bin/lastz_params/near.q T=2 Y=15000";

open PATH, "$path_conf";
while (<PATH>){
	chomp;
	next if /^#/;
	next if /""/;
	my ($program, $path) = split (/=/, $_);
	switch ($program){
		case("lastz"){	$lastz = $path;	}
		case("kent"){	$kent_src = $path;	}
	}
}
close PATH;


### create 2bit files ###
`mkdir -p $out_dir/ref`;
`mkdir -p $out_dir/tar`;
print STDERR "1. Create 2bit files & size files...\n";
my $ref_2bit = "$out_dir/ref/$bref_name.2bit";
my $tar_2bit = "$out_dir/tar/$btar_name.2bit";
my $ref_size = "$out_dir/ref/$bref_name.sizes";
my $tar_size = "$out_dir/tar/$btar_name.sizes";

my @two_size_cmds = ("$kent_src/faToTwoBit $ref_fa $ref_2bit","$kent_src/faToTwoBit $tar_fa $tar_2bit","$kent_src/faSize $ref_fa -detailed > $ref_size","$kent_src/faSize $tar_fa -detailed > $tar_size");
foreach my $c (0..$#two_size_cmds)
{
	my $pid = $pm->start($two_size_cmds[$c]) and next;
	`$two_size_cmds[$c]`;
	$pm->finish($c);
}
$pm -> wait_all_children;

### split fasta sequences ###
print STDERR "2. Split fasta files...\n";
`$kent_src/faSplit byName $ref_fa $out_dir/ref/`;
`$kent_src/faSplit byName $tar_fa $out_dir/tar/`;

my $dcnt = 0;
print STDERR "3.1 Removing small (< $resolution) target files...\n";
open(F,"$ref_size");
while(<F>) {
	chomp;
	my ($name, $size) = split(/\s+/);
	if (-f "$out_dir/ref/$name.fa" && $size < $resolution) {
		`rm -f $out_dir/ref/$name.fa`; 
		$dcnt++;
	}
}
close(F);
print STDERR "\t$dcnt target scaffolds were removed.\n";

$dcnt = 0;
print STDERR "3.2 Removing small (< $resolution) and merging query files...\n";
open(F,"$tar_size");
while(<F>) {
	chomp;
	my ($name, $size) = split(/\s+/);
	if (-f "$out_dir/tar/$name.fa" && $size < $resolution) {
			`rm -f $out_dir/tar/$name.fa`; 
			$dcnt++;
	}
}
close(F);
print STDERR "\t$dcnt query scaffolds were removed.\n";

### Make alignment jobs ###
my @ref_fas = <$out_dir/ref/*.fa>;
my @tar_fas = <$out_dir/tar/*.fa>;
my $num_ref = scalar(@ref_fas);
my $num_tar = scalar(@tar_fas);
my $total_jobs = $num_ref * $num_tar;
my $lav_dir = "$out_dir/lav";
`mkdir -p $lav_dir`;
my $psl_dir = "$out_dir/psl";
`mkdir -p $psl_dir`;
my $chain_dir = "$out_dir/chain";
`mkdir -p $chain_dir`;

# create a job file
my $job_cnt = 1;
foreach my $ref_file (@ref_fas) {
	my $ref_name = basename($ref_file,".fa");
	foreach my $tar_file (@tar_fas) {
		my $tar_name = basename($tar_file,".fa");
		`mkdir -p $lav_dir/$ref_name`;
		`mkdir -p $psl_dir/$ref_name`;
		`mkdir -p $chain_dir/$ref_name`;
		my $lav_f = "$lav_dir/$ref_name/$ref_name.$tar_name.lav";
		my $psl_f = "$psl_dir/$ref_name/$ref_name.$tar_name.psl";
		my $chain_f = "$chain_dir/$ref_name/$ref_name.$tar_name.chain";
		`mkdir -p $out_dir/jobs`;
		my $job_f = "$out_dir/jobs/$job_cnt.job"; 
		open(O,">$job_f");
		print O "#!/usr/bin/env perl\n";
		print O "use strict;\n";
		print O "use warnings;\n";
		print O "`$lastz $ref_file $tar_file $lastzPar > $lav_f`;\n";
		print O "`$kent_src/lavToPsl $lav_f $psl_f`;\n";
		print O "`$kent_src/axtChain -psl -verbose=0 $chainPar $psl_f $ref_2bit $tar_2bit stdout | $kent_src/chainAntiRepeat $ref_2bit $tar_2bit stdin $chain_f`;\n"; 
		close(O);
		$job_cnt++;
	}
}

### Submit jobs & wait for all job completion ###
my @fork_jobs = <$out_dir/jobs/*.job>;
my @align_cmds;
foreach my $job_f (@fork_jobs){
	push(@align_cmds, "perl $job_f");
}

foreach my $c (0..$#align_cmds){
	my $pid = $pm -> start($align_cmds[$c]) and next;
	my $job_num = $c+1;
	sleep(1);
	print STDERR "4. Lastz ($job_num/$total_jobs)\n";
	`$align_cmds[$c]`;
	$pm -> finish($c);
}
print STDERR "5. Wait for all job completion...\n";
$pm -> wait_all_children;

print STDERR "6. Chain merge sort...\n";
`$kent_src/chainMergeSort $chain_dir/*/*.chain > $out_dir/all.chain`;

print STDERR "7. Chain net...\n";
`$kent_src/chainNet $out_dir/all.chain -minSpace=1 $ref_size $tar_size stdout /dev/null | $kent_src/netSyntenic stdin $out_dir/all.net`;

print STDERR "8. Building ChainNet...\n";
`mkdir -p $chainNet_dir`;
`$Bin/splitChainNet.pl --chain $out_dir/all.chain --net $out_dir/all.net -o $chainNet_dir`;
