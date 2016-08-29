#!/usr/bin/perl

use warnings; use strict; use Getopt::Std; use foot;
use vars qw($opt_i $opt_r $opt_x $opt_y $opt_n);
getopts("i:r:x:y:n:");

# Usage {
my $usage = "

Usage: $YW$0 [Options] -r$CY <Fastq.fq>$N -n$CY <Output Folder>$N -i$CY <GeneIndex.bed>$N 

-r: [FastQ] FastQ file from Pacbio ccs
-n: [String] Folder of output
-i: [BED4+] file of all gene coordinate, tab separated.
    chr1	500	10000	RPL13A

${GN}Padding [default]:$N
-x: [0] Add this number from the$GN start$N of the index (strand taken into account)
-y: [0] Add this number from the$GN end$N of the index
    e.g. seq is:$GN chr1 200 500 SEQ 0 +$N and$GN -x -100 -y 50$N becomes: chr1 100 550 SEQ 0 +$N
    e.g. seq is:$GN chr1 200 500 SEQ 0 -$N and$GN -x -100 -y 50$N becomes: chr1 150 600 SEQ 0 -$N
    Strand is taken into account but if there's no strand info, it'll be treated as if it's + strand
";

# }

my $runID = join("", time(localtime));
open (OUTLOG, ">", "$runID.log") or die "Failed to create $runID.log: $!\n";
sanity_check();

my $geneIndexes = "$opt_i\_$x\_$y\_bp.bed";

system("bedtools_bed_change.pl -m -x $x -y $y -i $geneIndexes -o $geneIndexes\_$x\_$y\_bp.bed") == 0 or die "Failed to get +/- 100bp of $geneIndexes!\n";

#runs bismark (output file will be pacbio.fastq_bismark_bt2.sam) only if it hasn't been ran previously)
sub sanity_check {

	# Check if -i and -f exist
	my $error = "";
	$error .= err() . " -i Gene index [BED4+] isn't defined!\n" if not defined($opt_i);
	$error .= err() . " -r Pacbio Read [FastQ] isn't defined!\n" if not defined($opt_r);
	$error .= err() . " -i Gene index [BED4+] doesn't exist!\n" if defined($opt_i) and not -e $opt_i;
	$error .= err() . " -r Pacbio Read [FastQ] doesn't exist!\n" if defined($opt_r) and not -e $opt_r;
	my $die = 0;
	# Check -n output directory {
	$error .= err() . " -n Directory isn't defined!\n" if not defined($opt_n);
	if (defined($opt_n) and -e $opt_n) {
		my @mydir = split("/", getFullpath("$opt_n") . "/");
		$mydir;
		for (my $i = 0; $i < @mydir; $i++) {
			$mydir .= $mydir[$i];
			($error, $log, $die) = sys("mkdir $mydir", $error, $log) if not -d $mydir;
			last if $die == 1;
		}
	}
	if ($die == 1) {
		print OUTLOG $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
		die $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
	}
	# }
	
	# Check -i gene indexes and directory
	my $geneIndexes = $opt_i;
	my $geneIndexesFa = "$bismark_folder\/geneIndexes.fa";
	## Check the gene index itself {
	my @line = `cat $opt_i`;
	my $check = 0;
	for (my $i = 0; $i < @line; $i++) {
		chomp($line[$i]);
		my @arr = split("\t", $line[$i]);
		next if $line[$i] =~ /^track/; next if $line[$i] =~ /^#/;
		if (@arr >= 4 and $arr[2] =~ /^\d+$/ or $arr[3] =~ /^\d+/) {
			$check = 1;
		}
	}
	$check == 0 ? $error .= err() . " $opt_i doesn't seem to be a bed file (MUST contain$YW 4$N column of$YW tab$N separated:$LGN chromosome, start, end, and gene name$N!\n\n" : $log .= dat() . " $opt_i is a BED4+ file :)\n";
	if ($check == 0) {
		print OUTLOG $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
		die $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
	}
	## }
	## Check if folder exist {
	my $index = 1; # boolean to check create bismark index or not
	my ($md1) = `md5sum $geneIndexes` =~ /^(\w+)[ ]?/;
	my $bismark_folder = "~/Bismark_indexes/footLoop/$md1";
	my @folders = split("/", $bismark_folder);
	for (my $i = 0; $i < @folders; $i++) {
		if (not -d $folders[$i]) {
			$index = 0;
			($error, $log, $die) .= sys("mkdir $folders[$i]", $error, $log);
			last if $die == 1;
		}
	}
	if ($die == 1) {
		print OUTLOG $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
		die $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
	}
	($error, $log, $die) = sys("/bin/cp $geneIndexes $bismark_folder/")
	# }
	## If index doesn't exist, create
	print STDERR "\n${YW}1. Getting fasta sequence from $geneIndexes into $CY$geneIndexesFa$N\n";
	$log .= dat() . " ${YW}1. Getting fasta sequence from $geneIndexes into $CY$geneIndexesFa$N\n";
	if ($index == 0) {
		$log .= dat() . "\t- Running ${YW}bedtools getfasta$N -fi $opt_g -bed $geneIndexes -fo $geneIndexesFa -name\n";
		($err, $log, $die) = sys("fastaFromBed -fi $opt_g -bed $geneIndexes -fo $geneIndexesFa -name", $err, $log);
		if ($die == 1) {
			print OUTLOG $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
			die $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
		}
		($err, $log, $die) = sys("perl -pi -e 'tr/acgt/ACGT/' $geneIndexesFa", $err, $log);
		if ($die == 1) {
			print OUTLOG $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
			die $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
		}

		$log .= dat() . "\t${YW}1b. Running$CY bismark_genome_preparation$YW\n";
	  ($err, $log, $die) = sys("bismark_genome_preparation --bowtie2 $bismark_folder && md5sum $geneIndexesFa > $bismark_folder/Bisulfite_Genome/md5sum.txt");
	  if ($die == 1) {
			print OUTLOG $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
			die $usage . "\nLOG:\n$log\n\nERRORS:$error\n";
		}
	  $log .= dat() . "\t${GN}SUCCESS$N: $CY$bismark_folder\/Bisulfite_Genome$N already exist and is used!\n"
}

__END__

#runs bismark (output file will be pacbio.fastq_bismark_bt2.sam) only if it hasn't been ran previously
unless(-e "./pacbio.fastq_bismark_bt2.sam")
{
   if( -e "./geneIndexes.bed")
   {
     system("bedtools getfasta -fi $opt_g -bed geneIndexes.bed -fo geneIndexes.fa -name");
   }
   else
   {
      die "$YW Error:$N geneIndexes.bed does not exist in this directory";
   }
   system("bismark_genome_preparation --bowtie2 ./ > /dev/null");

   system("bismark --bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 ./ $opt_r > /dev/null") == 0 or die "Failed to run bismark: $!\n";
}
