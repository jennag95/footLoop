#!/usr/bin/perl

use warnings; use strict; use Getopt::Std; use foot; require foot;
use vars qw($opt_n $opt_t $opt_l $opt_r $opt_q $opt_s $opt_g $opt_c $opt_S $opt_i $opt_e $opt_f $opt_H $opt_L $opt_F $opt_p $opt_x $opt_y);
getopts("n:t:l:r:q:s:g:ci:S:e:fHL:Fpx:y:");

my $usage = "\nUsage: $YW$0$N ${GN}[options]$N -n $CY<Output name>$N -t $CY<conversion threshold>$N -l $CY<minimum R-loop length>$N -q $CY<minimum mapping quality for each read to be considered>$N -i$CY <geneIndexes.bed> -g $CY<path to reference genome for indexing>$N -r $CY<read file>$N -S $CY<sam file>

${GN}Extra information (required):$N
1) -n: Output name
2) Either:
${CY}-r$N: must be path to reads file itself named pacbio.fastq
${CY}-S$N: supply a sam file
3) -i: index bed file (bed4 with gene name on each)
4) -g: must be path to .fa file containing reference genome used for indexing (e.g. hg19.fa)

${GN}Optional [default]:$N
-x: [0] Add this number from the$GN start$N of the index (strand taken into account)
-y: [0] Add this number from the$GN end$N of the index 
e.g. seq is:$GN chr1 200 500 SEQ 0 +$N and$GN -x -100 -y 50$N becomes: chr1 100 550 SEQ 0 +$N
e.g. seq is:$GN chr1 200 500 SEQ 0 -$N and$GN -x -100 -y 50$N becomes: chr1 150 600 SEQ 0 -$N
-p: Run bismark script in paralel
-H: Run HMM peak caller
-c: <default: don't include CpG> consider Cs in CpG context
-t: [0.65] percentage (0.0-1.0) of Cs that must be converted to be considered an R-loop
-l: [100] minimum length in base pairs to be considered an R-loop peak (also for HMM)
-L: [500] minimum ${CY}read$N length in bp to be considered valid read
-F: Force re-create SAM file
-q: [0] q = map quality, where probability of (wrong) = 10^(q/-10). 
However for something like R-loop footprint from pacbio, where most stuff will be weird (indel+converted), this score is kind of meaningless so use default (0)
Some examples:
-  0 = 100     %
- 10 =  10     %
- 20 =   1     %
- 30 =    0.1  %
- 40 =    0.01 %

${GN}Example:$N

If you have .fq file but no SAM file:
$CY$0$N -H -p -n ${CY}myoutput$N -t ${CY}0.65$N -l$CY 250$N -e$CY /home/mitochi/shared/footLoop_example/hg19_appris_name.exon$N -i$CY /home/mitochi/shared/footLoop_example/pacbio12index_original.bed$N -g$CY /home/mitochi/shared/footLoop_example/hg19.fa.fa$N ${YW}-r /home/mitochi/shared/footLoop_example/pacbio12ccs.fq$N 

If you have SAM file, use -S instead of -r (everything is the same as above except$YW yellow$N):
$CY$0$N -H -p -n ${CY}myoutput$N -t ${CY}0.65$N -l$CY 250$N -e$CY /home/mitochi/shared/footLoop_example/hg19_appris_name.exon$N -i$CY /home/mitochi/shared/footLoop_example/pacbio12index_original.bed$N -g$CY /home/mitochi/shared/footLoop_example/hg19.fa.fa$N ${YW}-S /home/mitochi/shared/footLoop_example/pacbio12.sam$N

";

die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -n <output directory> not defined\n\n${LRD}   ########## USAGE ##########\n$N$usage" if not defined($opt_n);
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -r <read.fq> and $opt_s <read.sam> both not defined\n\n${LRD}   ########## USAGE ##########\n$N$usage" if not defined($opt_r) and not defined($opt_S);
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -i <geneindex.bed> not defined\n\n${LRD}   ########## USAGE ##########\n$N$usage" if not defined($opt_i);
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -l <min peak length in bp> must be positive integer!\n\n${LRD}   ########## USAGE ##########\n$N$usage" if defined($opt_l) and $opt_l !~ /^\d+$/;
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -L <min read length in bp> must be positive integer!\n\n${LRD}   ########## USAGE ##########\n$N$usage" if defined($opt_L) and $opt_L !~ /^\d+$/;
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -g <ref_genome.fa [hg19.fa]> not defined\n\n${LRD}   ########## USAGE ##########\n$N$usage" if not defined($opt_g);
my $read0 = defined($opt_r) ? $opt_r : "FALSE";
my $sam0 = defined($opt_S) ? $opt_S : "FALSE";
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -r $read0 and -S $sam0 both DOES NOT EXIST (provide at least one!)\n\n${LRD}########## FATAL ERROR ##########\n$N$usage" if (defined($opt_r) and not -e ($opt_r)) or (defined($opt_S) and not -e ($opt_S));
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -i $opt_i DOES NOT EXIST\n\n${LRD}########## FATAL ERROR ##########\n$N$usage" if not -e ($opt_i);
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -g $opt_g DOES NOT EXIST\n\n${LRD}########## FATAL ERROR ##########\n$N$usage" if not -e ($opt_g);
die "\n$LRD########## FATAL ERROR ##########\n\n!!!$N -S $CY$opt_S$N exists but seems to be empty!\n$N$usage" if defined($opt_S) and (not -e $opt_S or (-s $opt_S < 10));

#### Input ####
$opt_i = foot::getFullpath($opt_i);
$opt_g = foot::getFullpath($opt_g);
$opt_r = foot::getFullpath($opt_r) if defined($opt_r); my ($readname) = foot::getFullpath("$opt_n") . "/" . foot::getFilename($opt_r,'full') . "\_bismark_bt2.sam" if defined($opt_r);
my $mysam   = defined($opt_S) ? foot::getFullpath($opt_S) : $readname;
#-e $readname ? $readname : foot::getFullpath($opt_r,'full') . "_bismark_bt2.sam";
$opt_t = defined($opt_t) ? $opt_t : 0.65;
my $minPeakLength = defined($opt_l) ? $opt_l : 100;
my $minReadLength = defined($opt_L) ? $opt_L : 500;
$opt_q = defined($opt_q) ? $opt_q : 0;
my $myread = defined($opt_r) ? $opt_r : "${LRD}FALSE$N (Sam File Given)";
my $usecpg = $opt_c ? "${LGN}TRUE$N" : "${LRD}FALSE$N (Don't Use C from CpG)";
my $exonFile  = defined($opt_e) ? foot::getFullpath($opt_e) : "${LRD}FALSE$N (exon file not given)";
my ($x) = defined($opt_x) ? $opt_x : 0;
my ($y) = defined($opt_y) ? $opt_y : 0;

# Make directory
my $mydir 		= foot::getFullpath("$opt_n") . "/";
my $outFolder  = $mydir;

if (-e $outFolder) {
	print "\n\n$LRD------------------------ WARNING --------------------------------$N\n";
	print "\nOutput Folder ($CY-n$N) $CY$outFolder$N exists, continuing will overwrite this completely. Continue?\n\n(press$CY ENTER$N to continue or ${LGN}ctrl+c to cancel$N)
	\n";
	<STDIN>;
	#	if (-e $mysam and -s $mysam >= 10) {print "Expected samfile result $CY$mysam$N Exists! Are you sure you're overwriting results and not using SAM input option ($CY-S $mysam$N) instead?$N\n";<STDIN>;}
}
mkdir "$opt_n" if not -d "$opt_n";
mkdir $mydir if not -d $mydir;
chdir $mydir;
my $logFile = "$mydir/logFile.txt";
open(my $outLog, '>', $logFile);
my $mysam2 = "$mysam"; $mysam2 .= " $GN(will be generated)$N" if not defined($opt_S);
my $exonFileExists = defined($opt_e) ? "${GN}$opt_e$N" : "${LRD}FALSE$N";
my @time = localtime(time); my $time = join("", @time);
my $date = getDate();
my ($md) = `md5sum $opt_i` =~ /^(\w+)[ ]+/;
print $outLog "\n${YW}0. Initializing... output directory = $CY$mydir$N\n";
print $outLog "
Date   : $date
Run ID : $time

${YW}Input Parameters$N
1.  -n ${CY}Out Dir$N   :$mydir
2.  -r ${CY}Read$N      :$myread
3.  -S ${CY}SAM$N       :$mysam2
4.  -i ${CY}Index$N     :$opt_i (if padding added, will be $opt_i\_$x\_$y\_bp.bed)
5.  -g ${CY}Genome$N    :$opt_g
6.  -s ${CY}Seq$N       :/data/Bismark_indexes/footLoop/$md/geneIndexes.fa
7.  -c ${CY}UseCpG?$N   :$usecpg
8.  -t ${CY}Threshold$N :$opt_t
9.  -l ${CY}MinPkLen$N  :$minPeakLength
10. -L ${CY}MinRdLen$N  :$minReadLength
11. -q ${CY}MinMapQ$N   :$opt_q
12. -e ${CY}Exon$N      :$exonFileExists

";
print STDERR "
Date   : $date
Run ID : $time

${YW}Input Parameters$N
1.  -n ${CY}Out Dir$N   :$mydir
2.  -r ${CY}Read$N      :$myread
3.  -S ${CY}SAM$N       :$mysam2
4.  -i ${CY}Index$N     :$opt_i (if padding added, will be $opt_i\_$x\_$y\_bp.bed)
5.  -g ${CY}Genome$N    :$opt_g
6.  -s ${CY}Seq$N       :/data/Bismark_indexes/footLoop/$md/geneIndexes.fa
7.  -c ${CY}UseCpG?$N   :$usecpg
8.  -t ${CY}Threshold$N :$opt_t
9.  -l ${CY}MinPkLen$N  :$minPeakLength
10. -L ${CY}MinRdLen$N  :$minReadLength
11. -q ${CY}MinMapQ$N   :$opt_q
12. -e ${CY}Exon$N      :$exonFileExists

";
#runs bismark (output file will be pacbio.fastq_bismark_bt2.sam) only if it hasn't been ran previously


my $geneIndexes   = $opt_i; die "\n$LRD!!!\t$N$opt_i doesn't exist!\n\n" if not defined($geneIndexes);
system("bedtools_bed_change.pl -m -x $x -y $y -i $geneIndexes -o $geneIndexes\_$x\_$y\_bp.bed") == 0 or die "Failed to get +/- 100bp of $geneIndexes!\n";
$geneIndexes = "$geneIndexes\_$x\_$y\_bp.bed";
my ($bismark_folder_exist, $bismark_folder) = checkBismarkIndex($geneIndexes, $outLog);
my $geneIndexesFa = "$bismark_folder\/geneIndexes.fa";

print STDERR "\n${YW}1. Getting fasta sequence from $geneIndexes into $CY$geneIndexesFa$N\n";
print $outLog "\n${YW}1. Getting fasta sequence from $geneIndexes into $CY$geneIndexesFa$N\n";
if ($bismark_folder_exist == 0) {
	print $outLog "\t- Running ${YW}bedtools getfasta$N -fi $opt_g -bed $geneIndexes -fo $geneIndexesFa -name\n";
	system("fastaFromBed -fi $opt_g -bed $geneIndexes -fo $geneIndexesFa -name") == 0 ? print $outLog "\t${GN}SUCCESS$N: Output: $CY$geneIndexesFa$N\n" : die "Failed to run bedtools: $!\n";
	system("perl -pi -e 'tr/acgt/ACGT/' $geneIndexesFa") == 0 or die "\n$LRD!!!$N\tFailed to convert atgc into ATGC\n";
}

print $outLog "\n${YW}1b. Running$CY bismark$YW (output file will be$CY pacbio.fastq_bismark_bt2.sam$YW) only if it hasn't been ran previously$N\n";
print $outLog "\t- Running$YW bismark_genome_preparation$N --bowtie2 $outFolder\n";
my $check = 0;
#if (-d "$outFolder/Bisulfite_Genome/" and -e "$outFolder/Bisulfite_Genome/md5sum.txt") {
if ($bismark_folder_exist == 1) {#-d "$bismark_folder/Bisulfite_Genome/" and -e "$bismark_folder/Bisulfite_Genome/md5sum.txt") {
	my ($md5sum) = `cat $bismark_folder/Bisulfite_Genome/md5sum.txt`; chomp($md5sum);
	$md5sum = "NA" if not defined($md5sum); 
	($md5sum) = $md5sum =~ /^(\w+)[ ]+/;
	my $md5sum2 = `md5sum $geneIndexesFa`; chomp($md5sum2);
	($md5sum2) = $md5sum2 =~ /^(\w+)[ ]+/;
	$check = 1 if $md5sum eq $md5sum2;
	print $outLog "\tOlder bisulfite genome found but$LRD different$N (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n" if $md5sum ne $md5sum2;
	print $outLog "\tOlder bisulfite genome found and they're the$LGN same$N! (md5sum old = $CY$md5sum$N, new = $CY$md5sum2)$N\n" if $md5sum eq $md5sum2;
}
if ($bismark_folder_exist == 0) {#not -d "$bismark_folder/Bisulfite_Genome" or $check == 0) {
	system("bismark_genome_preparation --bowtie2 $bismark_folder && md5sum $geneIndexesFa > $bismark_folder/Bisulfite_Genome/md5sum.txt") == 0 or die "Failed to run bismark genome preparation: $!\n";
} 
else { 
	print $outLog "\t${GN}SUCCESS$N: $CY$bismark_folder\/Bisulfite_Genome$N already exist and is used!\n"
}
if ($opt_F or (not -e "$mysam" or -s $mysam <= 10)) {
	print "\t- Running$YW bismark$N --bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 $bismark_folder $opt_r\n";
	print $outLog "\t- Running$YW bismark$N --bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 $bismark_folder $opt_r\n";
	if ($opt_p) {
		print $outLog "\t  Splitting $CY$opt_r$N by 1000 sequences!\n";
		my $splitresult = `SplitFastq.pl -i $opt_r -o $outFolder -n 1000`; print $outLog "$splitresult\n";
		print "\t  Running bismark in paralel!\n";
		my $result = system("run_script_in_paralel2.pl -v \"bismark --bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 $bismark_folder FILENAME > FILENAME.bismark.log 2>&1\" $outFolder .part 20");
		my ($readFilename) = foot::getFilename($opt_r, "full");
		my @partSam = <$outFolder/*.part_bismark*.sam>; my $totalPartSam = @partSam;
		print "\t  All part.sam has been made (total = $totalPartSam). Now making $CY$mysam$N and then removing the part sam\n";
		for (my $p = 0; $p < @partSam; $p++) {
			my $partSam = $partSam[$p];
			print "\t\tPutting $partSam into $mysam and removing it!\n";
			system("cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >  $mysam") == 0 or die "Failed to cat $partSam: $!\n" if $p == 0;
			system("cat $partSam| awk '\$2 == 0 || \$2 == 16 {print}' >> $mysam") == 0 or die "Failed to cat $partSam: $!\n" if $p != 0;
			print "\t- Removing $CY$partSam$N: trash.pl $partSam\n";
			system("trash.pl $partSam") == 0 or die "Failed to trash.pl $partSam: $!\n";
		}
	}
	else {
		my $result = system("bismark --bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 $bismark_folder $opt_r");
		if ($result != 0) {
			print $outLog "\t${LRD}Bisulfte_Genome seems to be corrupted so re-running:\n\t${YW}-bismark_genome_preparation$N --bowtie2 $bismark_folder\n";
			system("bismark_genome_preparation --bowtie2 $bismark_folder") == 0 or die "Failed to run bismark genome preparation: $!\n";
			system("bismark --bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 $bismark_folder $opt_r") == 0 or die "$LRD!!!$N\tFailed to run bismark: $!\n";
		}
		print $outLog "\t${GN}SUCCESS$N: Output $mysam\n";
		print STDERR "\t${GN}SUCCESS$N: Output $mysam\n";
	}
}
else {
	print $outLog "\t${GN}SUCCESS$N: Output already exist: $CY$mysam$N\n";
	print STDERR "\t${GN}SUCCESS$N: Output already exist: $CY$mysam$N\n";
}

#takes sequence of gene and splits into array of individual bases
my $seqFile = $geneIndexesFa;
my %seq;
print STDERR "\n${YW}2. Parsing in sequence for genes from sequence file $CY$seqFile$N\n";
print $outLog "\n${YW}2. Parsing in sequence for genes from sequence file $CY$seqFile$N\n";
open(SEQ, $seqFile) or die "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!";
while (my $line = <SEQ>) {
	chomp($line);
	if ($line =~ />/) {
		my ($gene2) = $line =~ /^>(.+)$/;
		my $gene = uc($gene2);
		$line = <SEQ>;
		die "\n$LRD!!!$N\tFATAL ERROR: Corrupted fasta file $seqFile! $gene doesn't have sequence!\n\n" if $line =~ /^>/;
		@{$seq{$gene}{seq}} = split("", $line);
		$seq{$gene}{len} = @{$seq{$gene}{seq}};
		$seq{$gene}{total} = 0;
		$seq{$gene}{badlength} = 0;
		$seq{$gene}{lowq} = 0;
		$seq{$gene}{used} = 0;
		$seq{$gene}{pos} = 0;
		$seq{$gene}{neg} = 0;
		$seq{$gene}{orig} = $gene2;
	}
}
print STDERR "\t${GN}SUCCESS$N: Sequence has been parsed from fasta file $CY$seqFile$N\n";
print $outLog "\t${GN}SUCCESS$N: Sequence has been parsed from fasta file $CY$seqFile$N\n";

if ($opt_e) {
	#0 = white
	#1 = intron line
	#2 = exon line
	print $outLog "\t${BR}2a. Parsing in exon intron data from $CY$exonFile$N:\n";
	foreach my $gene (sort keys %seq) {
		next if -e "$outFolder/exon/$gene.exon"; # DELETE THIS
		my ($genepos) = `grep -i -P \"\\t$gene\\t\" $geneIndexes`; chomp($genepos);
		my $length_seq = $seq{$gene}{len};
		print $outLog "\n$LRD!!!$N\tWARNING: No sequence for gene $CY$gene$N is found in -s $CY$seqFile$N!\n\n" if $length_seq == 0;
		foot::parseExon($exonFile, $genepos, $gene, $outFolder, $length_seq);
		print STDERR "\t${GN}SUCCESS$N: Sequence of exon $CY$gene$N has been parsed from fasta file $CY$seqFile$N (length = $CY" . $length_seq . "$N bp)\n";
		print $outLog "\t${GN}SUCCESS$N: Sequence of exon $CY$gene$N has been parsed from fasta file $CY$seqFile$N (length = $CY" . $length_seq . "$N bp)\n";
	}
}

close SEQ;


# File description, open files, etc
my $samFile = $mysam;
my $positive = "$mydir/Positive.3";
my $negative = "$mydir/Negative.3";
my $notusedfile = "$mydir/NOTUSED.3";
print STDERR "${YW}3. Parsing sam file $CY$samFile$YW and getting only high quality reads\n";
print $outLog "${YW}3. Parsing sam file $CY$samFile$YW and getting only high quality reads\n";
open(my $positiveReads, ">", $positive) or die "Could not open $positive: $!";
open(my $negativeReads, ">", $negative) or die "Could not open $negative: $!";
open(my $notused, ">", $notusedfile) or die "Cannot open $notusedfile: $!\n";
open(my $sam, $samFile) or die "$LRD!!!$N\tFATAL ERROR: Could not open $samFile: $!";



## Some stats
my $linecount = 0; my %count; ($count{total}, $count{used}, $count{diffgene}, $count{lowq}, $count{badlength}) = (0,0,0,0,0);
my %readz;

#loops through each read and writes high quality reads into two separate files <gene>Positive.txt and <gene>Negative.txt
while(my $line = <$sam>) {
	$linecount ++;
	chomp($line);
	
	my @fields = split("\t", $line); #line is tab separated, first column is read name, the rest is value
	my $gene = uc($fields[2]);
	next if not defined($seq{$gene});
	#discounts the first 4 lines of information at the top of the .sam file
	next if(@fields < 6);
	
	# (EXTRA) This below is to shorten read name. Unique read name is the number between last number or behind ccs. If can't be found just use the whole read name.
	my ($read) = $fields[0] =~ /^.+\/(\d+)\/\d+/i;
	($read) = $fields[0] =~ /^.+\/(\d+)\/ccs/i if not defined($read);
	($read) = $fields[0] if not defined($read);
	$read = "SEQ_$read";
	
	# (EXTRA) This below is to show the user an example of parsed line and tells user if it's parsed every 20k read.
	my ($readname) = $fields[0] =~ /(.{20})$/; $readname = "..." . $readname;
	$count{total} ++;
	$seq{$gene}{total} ++;
	print $outLog "\tExample at read $count{total}: name=$CY$readname$N\tstrand=$CY$fields[1]$N\tchr/gene=$CY$fields[2]$N\tpos=$CY$fields[3]$N\tmapQ=$CY$fields[4]$N\n\n" if $count{total} == 1;
	print $outLog "\tDone $GN$count{total}$N\n" if $count{total} % 20000 == 0;
	
	#discounts if mapping quality is not at least phred score specified
	if($fields[4] < $opt_q)
	{
		print $notused "\t$CY$readname$N quality ($CY$fields[4]$N) is less than $CY$opt_q$N\n";
		$count{lowq} ++;
		$seq{$gene}{lowq} ++;
		next;
	}
	
	#discounts any read that is not the gene of interest, the proper length, or at the proper location in genome(accounting for indexing)
	elsif(length($fields[9]) < $minReadLength) # buffer is obsolete
	{
	# (EXTRA) This below is just for statistics
	$count{badlength} ++;
	$seq{$gene}{badlength} ++;
	die "READNAME undef\n" if not defined($readname);
	die "fields9 undef\n" if not defined($fields[9]);
	die "length of gene $fields[2] undef\n" if not defined($seq{$fields[2]}{len});
	
	print $notused "\t$CY$readname$N length of seq ($CY" . length($fields[9]). "$N) is less than 500bp (length of original sequence is ($CY" . $seq{$fields[2]}{len} . "$N)!\n";
	next
	}
	
	#counts number of CT conversions (takes into account whether or not the user wants to include Cs in CpG context)
	my $CT = ($fields[13] =~ tr/xhu/xhu/);
	if($opt_c)
	{
		$CT = ($fields[13] =~ tr/zxhu/zxhu/);
	}
	
	#writes positive reads into <gene>Positive.txt and negative reads into <gene>Negative.txt
	if($fields[1] == 0 || $fields[1] == 16)
	{
		my $to_be_printed = "$CT\t$line\n";
		$fields[1] == 0 ? print $positiveReads $to_be_printed : print $negativeReads $to_be_printed;
		$readz{$read} ++;
		$count{used} ++;
		$seq{$gene}{used} ++;
		$seq{$gene}{pos} ++ if $fields[1] == 0;
		$seq{$gene}{neg} ++ if $fields[1] == 16;
	}
}

my $passedFilterP = `wc -l < $positive`; chomp($passedFilterP);
my $passedFilterN = `wc -l < $negative`; chomp($passedFilterN);

print STDERR "\t${GN}SUCCESS$N: Total=$count{total}, used=$count{used}, Low Map Quality=$count{lowq}, Too short=$count{badlength}\n";
print $outLog "\n\t${GN}SUCCESS$N: Total=$count{total}, used=$count{used}, Low Map Quality=$count{lowq}, Too short=$count{badlength}\n";
print $outLog "\tOutputs are two separate files:\n\t\t- $CY$positive$N ($passedFilterP reads used)\n\t\t- $CY$negative$N ($passedFilterN reads used)\n\n";

print $outLog "
Reads that passed filters:
Positive: $passedFilterP
Negative: $passedFilterN
Total   : $count{total};

Per Gene:
";

foreach my $gene (keys %seq)
{
my $gene2 = $seq{$gene}{orig};
print $outLog "
- $gene (original name = $gene2):
Positive    = $seq{$gene}{pos}
Negative    = $seq{$gene}{neg}
Used        = $seq{$gene}{used}
Total       = $seq{$gene}{total}
Too Short   = $seq{$gene}{badlength}
Low Quality = $seq{$gene}{lowq}
";
}

#sorts by number of CT conversions, takes top 1000, and removes the number of CT conversions in prepartion for methylation extractor
print STDERR "${YW}4. Sorting $positive and $negative by number of CT conversions and removes the number of CT conversions in prepartion for methylation extractor\n";
print $outLog "${YW}4. Sorting $positive and $negative by number of CT conversions and removes the number of CT conversions in prepartion for methylation extractor\n";
my $finalPositive = "$mydir/PositiveFinal.txt";
system("sort -k1,1rn $positive | cut -f 2- > $finalPositive");
my ($finalPositiveLine) = `wc -l $finalPositive` =~ /^(\d+) /;

my $finalNegative = "$mydir/NegativeFinal.txt";
system("sort -k1,1rn $negative | cut -f 2- > $finalNegative");
my ($finalNegativeLine) = `wc -l $finalNegative` =~ /^(\d+) /;
print $outLog "\t${GN}SUCCESS$N: Output:\n\t\t- $CY$finalPositive$N ($finalPositiveLine reads used)\n\t\t- $CY$finalNegative$N ($finalNegativeLine reads used)\n\n";

#runs bismark methylation extractor on top 1000 reads of each strand
my $CPGpos = $mydir . "CpG_context_" . "PositiveFinal.txt";
my $CPGneg = $mydir . "CpG_context_" . "NegativeFinal.txt";
my $CHGpos = $mydir . "CHG_context_" . "PositiveFinal.txt";
my $CHGneg = $mydir . "CHG_context_" . "NegativeFinal.txt";
my $CHHpos = $mydir . "CHH_context_" . "PositiveFinal.txt";
my $CHHneg = $mydir . "CHH_context_" . "NegativeFinal.txt";
my @bismarkOutput = ($CPGpos, $CPGneg, $CHGpos, $CHGneg, $CHHpos, $CHHneg);
my $bismarkOutput = $CHGpos;

print STDERR "${YW}5. Running bismark_methylation_extractor on $CY$finalPositive$YW and $CY$finalNegative$N\n";
print $outLog "${YW}5. Running bismark_methylation_extractor on $CY$finalPositive$YW and $CY$finalNegative$N\n";
if (not $opt_f and (not -e $bismarkOutput or -s $bismarkOutput <= 10))
{
	system("bismark_methylation_extractor -s --comprehensive $finalPositive");
	system("bismark_methylation_extractor -s --comprehensive $finalNegative");
}
print STDERR "\t${GN}SUCCESS$N: Output: 4-6 files of <CpG/CHG/CHH>$CY\_context_$finalPositive$N:\n";
print $outLog "\t${GN}SUCCESS$N: Output: 4-6 files of <CpG/CHG/CHH>$CY\_context_$finalPositive$N:\n";
for (my $i = 0; $i < @bismarkOutput; $i++) {
	next if $opt_f;
	if (not -e $bismarkOutput[$i]) {
		print $outLog "\t\t$bismarkOutput[$i]: has$LRD 0$N reads!\n"; system("touch $bismarkOutput[$i]") == 0 or die; next;
	}
	my ($linecount) = `wc -l $bismarkOutput[$i]` =~ /^(\d+) /;
	my ($readnumber) = `unique_column.pl $bismarkOutput[$i] 1 | wc -l` =~ /^(\d+)$/;
	print STDERR "\t\t- $bismarkOutput[$i]: has$LRD 0$N reads!\n" if $linecount == 0;
	print STDERR "\t\t- $bismarkOutput[$i]: $CY$linecount$N total line and $CY$readnumber$N reads\n" if $linecount > 0;
	print $outLog "\t\t- $bismarkOutput[$i]: has$LRD 0$N reads!\n" if $linecount == 0;
	print $outLog "\t\t- $bismarkOutput[$i]: $CY$linecount$N total line and $CY$readnumber$N reads\n" if $linecount > 0;
}

#pulls the CHH, CHG, and CpG sites together into methylationPos<gene>.txt and methylationNeg<gene>.txt (takes into account -c option)
print STDERR "${YW}6. Combine CHH and CHG (and CpG if -c) sites together into$CY methylationPos<gene>.txt$YW and$CY methylationNeg<gene>.txt$N\n";
print $outLog "\n${YW}6. Combine CHH and CHG (and CpG if -c) sites together into$CY methylationPos<gene>.txt$YW and$CY methylationNeg<gene>.txt$N\n";
my $methylationPos = $mydir . "methylationPos" . ".txt";
my $methylationNeg = $mydir . "methylationNeg" . ".txt";
$methylationPos = $mydir . "methylationPos" . "CG.txt" if($opt_c);
$methylationNeg = $mydir . "methylationNeg" . "CG.txt" if($opt_c);

if($opt_c)
{
	system("cat $CPGpos $CHGpos $CHHpos | sort -n > $methylationPos") if not -e $methylationPos or -s $methylationPos < 10;
	system("cat $CPGneg $CHGneg $CHHneg | sort -n > $methylationNeg") if not -e $methylationNeg or -s $methylationNeg < 10;
}
else {#if (not defined($opt_f) or not -e $methylationPos or not -e $methylationNeg)
	system("cat $CHGpos $CHHpos | sort -n > $methylationPos") if not -e $methylationPos or -s $methylationPos < 10;
	system("cat $CHGneg $CHHneg | sort -n > $methylationNeg") if not -e $methylationNeg or -s $methylationNeg < 10;
}
	#gets the position of each conversion (conversions=1; not converted=0)
my %read; my %info; my %infoGene;
for (my $i = 0; $i < 2; $i++) {
	my $analyzeFile = $i == 0 ? $methylationPos : $methylationNeg;
	
	print $outLog "\t6.$i. Parsing $CY$analyzeFile$N\n";
	print "\t6.$i. Parsing $CY$analyzeFile$N\n";
	open(my $analyze, $analyzeFile) or die "$LRD!!!$N\tFATAL ERROR: Could not open $analyzeFile: $!";
	
	# Initialize array @read which are zeroes with length of @seq
	my $linecount = 0;
	while(my $line = <$analyze>) {
		chomp($line);
		print "Done $linecount\n" if $linecount % 100000 == 0;
		$linecount ++;
		next if $line  =~ /^Bismark/;
		my @fields = split("\t", $line);
		my ($read) = $fields[0] =~ /^.+\/(\d+)\/\d+/i;
		($read) = $fields[0] =~ /^.+\/(\d+)\/ccs/i if not defined($read);
		($read) = $fields[0] if not defined($read);
		$read = "SEQ_$read";
		my ($readname) = length($fields[0]) > 20 ? $fields[0] =~ /(.{20})$/ : $fields[0]; $readname = "..." . $readname;
		my ($name, $junk, $gene, $pos, $conv) = @fields; $pos = $pos - 1;
		print $outLog "\tExample: name=$CY$readname$N, junk=$CY$fields[1]$N, gene=$CY$fields[2]$N, pos=$CY$fields[3]$N, conv=$CY$fields[4]$N\n" if not defined($read{$gene}{$i});
		$read{$gene}{$i}{$read}{$pos} = $conv =~ /^[xzhu]$/ ? 1 : $conv =~ /^[XZHU]$/ ? 0 : die "$LRD!!!$N\tFATAL ERROR: conversion isn't x/z/h/u (case ins) ($CY$conv$N)in $CY$analyzeFile$N line:\n$line\n\n";
		$infoGene{$gene}{$read}{min} = $pos if not defined($infoGene{$gene}{$read}{min}) or $infoGene{$gene}{$read}{min} > $pos;
		$infoGene{$gene}{$read}{max} = $pos if not defined($infoGene{$gene}{$read}{max}) or $infoGene{$gene}{$read}{max} < $pos;
	}
}
	
print STDERR "\t${GN}SUCCESS$N: Done parsing methylationPos and methylationNeg files!\n";
print $outLog "\n\t${GN}SUCCESS$N: Done parsing methylationPos and methylationNeg files!\n";
	
print STDERR "${YW}7. Converting each gene's sequence position and methylation data\n$N";
print $outLog "${YW}7. Converting each gene's sequence position and methylation data\n$N";

open (my $SEQPOSOUT, ">", "$mydir/SEQ_GENE_POS.tsv") or die;
open (my $SEQNEGOUT, ">", "$mydir/SEQ_GENE_NEG.tsv") or die;
foreach my $gene (sort keys %seq) {
#	next if defined("$mydir/SEQ_GENE_POS.tsv") and defined("$mydir/SEQ_GENE_NEG.tsv");
	print $outLog "\t- Processing gene $CY$gene$N\n";
	my $filePos = $opt_c ? "$mydir/$gene\_CG_POS" : "$mydir/$gene\_POS";
	my $fileNeg = $opt_c ? "$mydir/$gene\_CG_NEG" : "$mydir/$gene\_NEG";
	
	my @seq = @{$seq{$gene}{seq}}; #@seq is the gene real sequence in + direction (even if the gene is - strand direction)
	my %ref = %{convert_seq(\@seq)};
	
	open(my $FILEPOS, ">", "$filePos.tsv") or die "Could not open $filePos.tsv: $!";
	open(my $FILENEG, ">", "$fileNeg.tsv") or die "Could not open $fileNeg.tsv: $!";
	open(my $FILEPOSFA, ">", "$filePos.customfa") or die "Could not open $filePos.customfa: $!";
	open(my $FILENEGFA, ">", "$fileNeg.customfa") or die "Could not open $fileNeg.customfa: $!";
	print $FILEPOS "#READ   ";
	print $FILENEG "#READ   ";
	# print ref header
	foreach my $strand (sort keys %{$read{$gene}}) {
		for (my $i = 0; $i < @seq; $i++) {
			my $ref = $ref{$strand}[$i];
			print $FILEPOS "\t$ref" if $strand == 0;
			print $FILENEG "\t$ref" if $strand == 1;
			print $FILEPOS "\n" if $i == @seq - 1 and $strand == 0;
			print $FILENEG "\n" if $i == @seq - 1 and $strand == 1;
		}
	}
	foreach my $strand (sort keys %{$read{$gene}}) {
		foreach my $name (keys %{$read{$gene}{$strand}}) {
			
			my @value; my @fasta; my @hmmval; #valueJ=joined
			for (my $i = 0; $i < @seq; $i++) {
				my $currvalue = $read{$gene}{$strand}{$name}{$i};
				my $ref = $ref{$strand}[$i];
				
				# B = border, C = C used, P = CpG used, N = Non C or CpG if not opt_c
				# if at border, just put as no data (in case of C of CpG)
				if ($i < $infoGene{$gene}{$name}{min} or $i > $infoGene{$gene}{$name}{max} or $ref eq "B") {
					$value[$i] = 6; $hmmval[$i] = $value[$i];
				}
				elsif ($ref eq "N") {
					$value[$i] = 2; $hmmval[$i] = $value[$i];
				}
				elsif ($ref eq "C") {
					$value[$i] = defined($currvalue) ? $currvalue : 6; $hmmval[$i] = $value[$i];
				}
				elsif ($ref eq "P") {
					$value[$i] = defined($currvalue) ? $currvalue + 3 : 6; 
					$hmmval[$i] = defined($currvalue) ? $currvalue : 6;
				}
				else {
					die "\n\nSomething very wrong happened at read=$name, i=$i\n\n";
				}
			}
			print $FILEPOS   "$name\t" . join("\t", @value) . "\n" if $strand == 0;
			print $FILENEG   "$name\t" . join("\t", @value) . "\n" if $strand == 1;
			print $FILEPOSFA ">$name\n" . join("", @hmmval) . "\n" if $strand == 0;
			print $FILENEGFA ">$name\n" . join("", @hmmval) . "\n" if $strand == 1;
		}
	}
	close $FILEPOS;
	close $FILENEG;
	
	# To check: fastaFromBed -fi $geneIndexesFa -bed test.bed -fo test.fa && cat test.fa
	# ALL MUST BE C
	my ($filePosLine) = `wc -l $filePos.tsv` =~ /^(\d+) /;
	my ($fileNegLine) = `wc -l $fileNeg.tsv` =~ /^(\d+) /;
	print STDERR "\t${GN}SUCCESS$N: Done parsing $YW$gene$N Output:\n\t\t- $CY$filePos.tsv$N ($filePosLine reads)\n\t\t- $CY$fileNeg.tsv$N ($fileNegLine reads)\n";
	print $outLog "\t${GN}SUCCESS$N: Done parsing $YW$gene$N Output:\n\t\t- $CY$filePos.tsv$N ($filePosLine reads)\n\t\t- $CY$fileNeg.tsv$N ($fileNegLine reads)\n\n";
	#system("mv $filePos.tsv $fileNeg.tsv $filePos.customfa $fileNeg.customfa $filePos.trans $fileNeg.trans $mydir/");
	my ($finalmd5) = -e "$filePos.tsv" ? "$filePos.tsv:" .`md5sum $filePos.tsv` : "md5sum $filePos.tsv DOES NOT EXIST!\n";
	($finalmd5)    = -e "$fileNeg.tsv" ? "$fileNeg.tsv:" .`md5sum $fileNeg.tsv` : "md5sum $fileNeg.tsv DOES NOT EXIST!\n";
	($finalmd5)    = -e "$filePos.customfa" ? "$filePos.customfa:" .`md5sum $filePos.customfa` : "md5sum $filePos.customfa DOES NOT EXIST!\n";
	($finalmd5)    = -e "$fileNeg.customfa" ? "$fileNeg.customfa:" .`md5sum $fileNeg.customfa` : "md5sum $fileNeg.customfa DOES NOT EXIST!\n";
	print $outLog "$finalmd5\n";
}

if ($opt_e) {
	foreach my $gene (sort keys %seq) {
		print $outLog "$mydir/exon/$gene.exon:" . `md5sum $mydir/exon/$gene.exon`;
	}
}
my @seq;

if ($opt_H) {
	print "${YW}8a. Running HMM$N\n\n";
	my $countGene = 1; my $totalGene = (keys %seq);
	foreach my $gene (sort keys %seq) {
		print "- (${LPR}$countGene/$totalGene$N)\tRunning HMM for gene: $CY$gene$N\n";
		my $filePos = $opt_c ? "$mydir/$gene\_CG_POS" : "$mydir/$gene\_POS";
		my $fileNeg = $opt_c ? "$mydir/$gene\_CG_NEG" : "$mydir/$gene\_NEG";
		print "\tCommand:$CY Pacbio_Heatmap.pl -l $minPeakLength $filePos.tsv $time$N\n"; 
		my $log1 = `Pacbio_Heatmap.pl -l 40 $filePos.tsv $time`; print $outLog "\nHMM:\n$log1\n";
		$countGene ++;
		if (-e "$filePos.HMMPEAK") {
			system("fastaFromBed -fi $geneIndexesFa -bed $filePos.HMMPEAK.bed -fo $filePos.HMMPEAK.fa -name");
		}
		if (-e "$fileNeg.HMMPEAK") {
			system("fastaFromBed -fi $geneIndexesFa -bed $fileNeg.HMMPEAK.bed -fo $fileNeg.HMMPEAK.fa -name");
		}
	}
}
print $outLog "${YW}8. determines which regions of conversion are R-loops based on conversion threshold$N\n";
print STDERR "${YW}8. determines which regions of conversion are R-loops based on conversion threshold$N\n";

foreach my $gene (sort keys %seq) {
	my @seq = @{$seq{$gene}{seq}};
	my $finalPos = $mydir . "$gene\_Pos" . ($opt_t*100) . ".txt";
	my $finalNeg = $mydir . "$gene\_Neg" . ($opt_t*100) . ".txt";
	$finalPos = $mydir . "$gene\_Pos" . ($opt_t*100) . "CG.txt" if($opt_c);
	$finalNeg = $mydir . "$gene\_Neg" . ($opt_t*100) . "CG.txt" if($opt_c);
	my $finalPos_NOPEAK = $mydir . "$gene\_Pos" . ($opt_t*100) . "_NOPEAK.txt";
	my $finalNeg_NOPEAK = $mydir . "$gene\_Neg" . ($opt_t*100) . "_NOPEAK.txt";
	$finalPos_NOPEAK = $mydir . "$gene\_Pos" . ($opt_t*100) . "_NOPEAK_CG.txt" if($opt_c);
	$finalNeg_NOPEAK = $mydir . "$gene\_Neg" . ($opt_t*100) . "_NOPEAK_CG.txt" if($opt_c);
	open(FINALPOS, ">", $finalPos) or die "Could not open $finalPos: $!";
	open(FINALNEG, ">", $finalNeg) or die "Could not open $finalNeg: $!";
	open(FINALPOS_NOPEAK, ">", $finalPos) or die "Could not open $finalPos: $!";
	open(FINALNEG_NOPEAK, ">", $finalNeg) or die "Could not open $finalNeg: $!";
	for(my $strand=0; $strand<2; $strand++)
	{
		my $filePos = defined($opt_c) ? "$mydir/$gene\_CG_POS" : "$mydir/$gene\_POS";
		my $fileNeg = defined($opt_c) ? "$mydir/$gene\_CG_NEG" : "$mydir/$gene\_NEG";
		my $fileLast = $strand == 0 ? "$filePos.tsv" : "$fileNeg.tsv";
		print $outLog "\t- Doing gene $CY$gene$N ($fileLast)\n";
		print $outLog "$fileLast doesn't exist!\n" and next if not -e $fileLast;
		print $outLog "$fileLast has no line!\n" and next if `wc -l < $fileLast` == 0;
		open (my $fileLastIn, "<", $fileLast) or print $outLog "Cannot open $fileLast\n" and die "Cannot open $fileLast: $!\n";
		my ($fileLastCount) = `wc -l $fileLast` =~ /^(\d+) /;
		
		my %peak;
		my $start = 1;
		my $end = $minPeakLength;
		my $linecount = 0;
		while (my $lineFinal = <$fileLastIn>)
		{
			$linecount ++;
			chomp($lineFinal);
			print $outLog "\tDone $GN$linecount$N\n" if $linecount % 500 eq 0;
			print "\t$fileLast: Done $GN$linecount$N / $fileLastCount\n" if $linecount % 500 eq 0;
			next if $lineFinal =~ /^\#READ/;
			my ($name, @fields) = split("\t", $lineFinal);
			$peak{$name}{length} = 0;
			#0= not converted (grey) # same
			#1= converted (green)   # same
			#2= non-C (white)       # same
			#3= non-converted CpG (blue) # is non converted
			#4= converted CpG  (black) # is converted
			#5= converted CpG in R-loop (purple) # conv CPG
			#6= no data
			#9= converted C in R-loop (red) # not exist
			# conversion of mine into Jenna's: 3 (no data) becomes 2 (non C)
			
			# peak call based on threshold
			my @peak; my $ispeak = 0;
			for (my $i = 0; $i < @fields - $minPeakLength; $i++)
			{
				$peak{$name}{length} ++ if $fields[$i] =~ /^[1459]$/;
				my ($conC, $nonConC, $conPer) = (0,0,0);
				# calculate % conv from start to end, where end is minPeakLength (e.g. 100 then 1 to 100)
				
				# calculate conC (converted C) and nonConC (non converted)
				for(my $j = $i; $j < $i + $minPeakLength; $j++)
				{
					$conC    ++ if($fields[$j] =~ /^[1459]$/);
					$nonConC ++ if($fields[$j] =~ /^[03]$/);
				}
				# calculate % conv
				if ($conC + $nonConC < 5) 
				{ # if too few C, (Alex Meissner uses 5 as threshold) we can't call it as anything ... Meissner et al (Nature 2012)
					$conPer = 0;
				}
				elsif($conC != 0 || $nonConC != 0)
				{
					$conPer = (($conC)/($conC+$nonConC));
				}
				else
				{
					$conPer = 0;
				}
				
				# if % conv is more than thershold, then turn 1 (conv C) into 9 and 4 (conv C in CpG) into 5
				if($conPer >= $opt_t)
				{
					$ispeak = 1;
					$peak[$i] = 1; # from i to i + minPeakLength is a peak
					for(my $j = $i; $j < $i + $minPeakLength; $j++)
					{
						$fields[$j] = 9 if $fields[$j] == 1;
						$fields[$j] = 5 if $fields[$j] == 4;
					}
				}
				else {
					$peak[$i] = 0;
				}
			}
			$peak{$name}{field} = join("\t", @fields);
			
			# if it's a read that has peak, then we need to put the peak length etc so it'll be presorted for R
			if ($ispeak == 1) {
				my ($maxLength) = 0; my $length = 0; my $lastPos = -2;
				for (my $i = 0; $i < @peak; $i++) {
					if ($lastPos != -2 and $lastPos == $i - 1 and $peak[$i] == 1) {
						$length ++; $lastPos = $i;
						$maxLength = $length if $maxLength < $length;
					}
					else {
						$length = 0; $lastPos = -2;
					}
				}
				$peak{$name}{length} = $maxLength;
				$peak{$name}{peak} = 1;
			}
			else {$peak{$name}{peak} = 0;}
		}
		my $peakcount = 0;
		foreach my $name (sort {$peak{$b}{peak} <=> $peak{$a}{peak} || $peak{$b}{length} <=> $peak{$a}{length}} keys %peak)
		{
			last if $peak{$name}{peak} == 0 and $peakcount > 1000;
			print FINALPOS "$name\t$peak{$name}{field}\n" if $strand == 0 and $peak{$name}{peak} == 1;
			print FINALPOS_NOPEAK "$name\t$peak{$name}{field}\n" if $strand == 0 and $peak{$name}{peak} == 0;
			print FINALNEG "$name\t$peak{$name}{field}\n" if $strand == 1 and $peak{$name}{peak} == 1;
			print FINALNEG_NOPEAK "$name\t$peak{$name}{field}\n" if $strand == 1 and $peak{$name}{peak} == 0;
			$peakcount ++;
		} 
	}
	
	my ($finalPosLine) = `wc -l $finalPos` =~ /^(\d+) /;
	my ($finalNegLine) = `wc -l $finalNeg` =~ /^(\d+) /;
	
	open (my $outReport, ">", "$mydir/foot.report") or die "Cannot write to $mydir/foot.report: $!\n";
	print $outLog "\n\t${GN}SUCCESS$N: Gene $CY$gene$N Output:\n\t\t- $CY$finalPos$N ($finalPosLine reads)\n\t\t- $CY$finalNeg$N ($finalNegLine reads)\n\n";
	print $outReport "\n\t${GN}SUCCESS$N: Gene $CY$gene$N Output:\n\t\t- $CY$finalPos$N ($finalPosLine reads)\n\t\t- $CY$finalNeg$N ($finalNegLine reads)\n\n";
	print STDERR "\t${GN}SUCCESS$N: Gene $CY$gene$N Output:\n\t\t- $CY$finalPos$N ($finalPosLine reads)\n\t\t- $CY$finalNeg$N ($finalNegLine reads)\n\n";
	
	#makes heatmaps for positive and negative strand
	
	#0= not converted (white) # same
	#1= converted (green)   # same
	#2= non-C (lightyellow)       # same
	#3= non-converted CpG (black) # is non converted
	#4= converted CpG  (blue) # is converted
	#5= converted CpG in R-loop (purple) # conv CPG
	#6= no data (grey)
	#9= converted C in R-loop (red) # not exist
	my $finalPosPDF = $mydir . "$gene\_Pos" . ($opt_t*100) . ".pdf";
	my $finalNegPDF = $mydir . "$gene\_Neg" . ($opt_t*100) . ".pdf";
	$finalPosPDF = $mydir . "$gene\_Pos" . ($opt_t*100) . "CG.pdf" if($opt_c);
	$finalNegPDF = $mydir . "$gene\_Neg" . ($opt_t*100) . "CG.pdf" if($opt_c);
	my $finalPosPDF_NOPEAK = $mydir . "$gene\_Pos" . ($opt_t*100) . "_NOPEAK.pdf";
	my $finalNegPDF_NOPEAK = $mydir . "$gene\_Neg" . ($opt_t*100) . "_NOPEAK.pdf";
	$finalPosPDF_NOPEAK = $mydir . "$gene\_Pos" . ($opt_t*100) . "_NOPEAK_CG.pdf" if($opt_c);
	$finalNegPDF_NOPEAK = $mydir . "$gene\_Neg" . ($opt_t*100) . "_NOPEAK_CG.pdf" if($opt_c);
	my $Rscript = "$mydir/$gene\_MakeHeatmap.R";
	open(my $out, ">", $Rscript) or die "Can't print to $Rscript: $!\n";
	my $breaks = defined($opt_c) ? "c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 9.5)" : "c(-0.5, 0.5, 1.5, 2.5, 6.5, 9.5)";
	my $colors = defined($opt_c) ? "c(\"white\",\"fec44f\",\"ffffb3\",\"deebf7\",\"9ecae1\",\"c51b8a\",\"grey\",\"red\")" : "c(\"white\",\"fec44f\",\"ffffb3\",\"grey\",\"red\")";
	if($opt_c)
	{
		print $out "

		.libPaths()
		library(\"GMD\")

		# print peak
		df = read.table(\"$finalPos\", sep=\"\t\", row.names=1)
		if (length(df) > 0) {
			png(\"$finalPosPDF\",height=2000,width=2000,pointsize=50)
			heatmap.3(
				x=df,
				dendrogram=\"none\",
				Rowv=TRUE, Colv=FALSE,
				labRow=FALSE,labCol=FALSE,
				breaks=$breaks,
				color.FUN=function(x) $colors
			)
			dev.off()
		}

		df = read.table(\"$finalNeg\", sep=\"\t\", row.names=1)

		if (length(df) > 0) {
			png(\"$finalNegPDF\",height=2000,width=2000,pointsize=50)
			heatmap.3(
				x=df,
				dendrogram=\"none\",
				Rowv=TRUE, Colv=FALSE,
				labRow=FALSE,labCol=FALSE,
				breaks=$breaks,
				color.FUN=function(x) $colors
			)
			dev.off()
		}

		df = read.table(\"$finalPos_NOPEAK\", sep=\"\t\", row.names=1)
		if (length(df) > 0) {
			png(\"$finalPosPDF_NOPEAK\",height=2000,width=2000,pointsize=50)
			heatmap.3(
				x=df,
				dendrogram=\"none\",
				Rowv=TRUE, Colv=FALSE,
				labRow=FALSE,labCol=FALSE,
				breaks=$breaks,
				color.FUN=function(x) $colors
			)
			dev.off()
		}

		df = read.table(\"$finalNeg_NOPEAK\", sep=\"\t\", row.names=1)

		if (length(df) > 0) {
			png(\"$finalNegPDF\",height=2000,width=2000,pointsize=50)
			heatmap.3(
				x=df,
				dendrogram=\"none\",
				Rowv=TRUE, Colv=FALSE,
				labRow=FALSE,labCol=FALSE,
				breaks=$breaks,
				color.FUN=function(x) $colors
			)
			dev.off()
		}
		";
	}
	close $out;
	my $cmd = "R --no-save < $Rscript >> $logFile\n";
	print ($cmd);
	system($cmd);
	system("rm Rplots.pdf") if -e "Rplots.pdf";
}

print STDERR "\n${LPR}If there is not PDF made from step (8) then it's due to too low number of read/peak$N\n\n";
sub getDate {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time); $year += 1900;
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my $date = "$mday $months[$mon] $year $hour:$min:$sec";
	my $timenow = $hour * 3600 + $min * 60 + $sec;
	return($date);
}

sub convert_seq {
	my @seq = @{$_[0]};
	my %new;
	for (my $strand = 0; $strand < 2; $strand++) {
		my @new;
		for (my $i = 0; $i < @seq; $i++) {
			my $lastpos = $strand == 0 ? 0 : @seq-1;
			my $add     = $strand == 0 ? 1 : -1; 
			my $nuc1 = $strand == 0 ? "C" : "G";
			my $nuc2 = $strand == 0 ? "G" : "C";
			
			if ($i == 0 || $i == @seq - 1) {
				$new[$i] = "B";
			}
			elsif ($seq[$i] ne $nuc1) {
				$new[$i] = "N";
			}
			else {
				if ($seq[$i+$add] eq $nuc2) {
					$new[$i] = "P" if defined($opt_c);
					$new[$i] = "N" if not defined($opt_c);
				}
				else {
					$new[$i] = "C";
				}
			}
		}
		@{$new{$strand}} = @new;
	}
	return(\%new);
}

