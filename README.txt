command: ./footLoop.pl -n <gene name> -t <minimum conversion threshold> -l <minimum length> -r <path to reads pacbio.fastq file> -q <minimum quality> -s <path to gene seq> -g <path to reference genome for mapping>

example command: ./footLoop.pl -n CALM3 -t .75 -l 100 -r ./pacbio.fastq -q 0 -s ./CALM3seq.txt -g /home/mitochi/Bowtie2_indexes/hg19/hg19.fa.fa -c

required options:
-n : name of the gene you want to generate heatmaps for (must be in all CAPS)
-t : minimum threshold percentage (expressed as a decimal, 0.0-1.0) of cytosines that must be converted in the window for it to be considered a footprint
-l : minimum length of footprint (this will be the size of the sliding window)
-r : path to fastq file containing the CCS reads *must be named pacbio.fastq* (This will only be used once to initially map all the reads. Once mapping has occurred, the script will recognized that a .sam file already exists when generating subsequent heatmaps and bismark will not be ran again.)
-q : minimum phred score quality threshold (e.g. score of 40 = 99.99% base call accuracy. Since we started using CCS reads, we have just been setting this value to 0.)
-s: path to .txt file containing gene sequence that was PCRed
-g : path to reference genome used for indexing (I've just been using the reference genome Stella has at /home/mitochi/Bowtie2_indexes/hg19/hg19.fa.fa)

optional option:
-c : consider cytosines in CpG context


***Before running the script, a bed file must be present in the same directory as the script. This bed file should contain all loci. The starting position must be 50 bp upstream of 5' PCR primer and the ending position must be 50 bp downstream of 3' PCR primer. This bed file must be named geneIndexes.bed***


If a .sam file doesn't already exist, the user is prompted to enter the chromosome number, starting and ending index of the PCR product, and the gene name (all separated by tabs). This is then used to create a .bed file that bedtools uses to create a .fa file. This geneIndexes.fa file is used by bismark to map the reads in pacbio.fastq. The mapped reads are contained in the .sam file. 

Reads from the .sam file are filtered by phred score (-q user entered) and gene of interest. Reads passing this filter are written into <gene>Positive.txt and <gene>Negative.txt (separated by +/- strand). Bismark methylation extractor is ran on these two files separately. 

The CHH, CHG, and CpG (if the -c option is used) outputs from the bismark methylation extractor are pulled together into methylationPos<gene>.txt and methylationNeg<gene>.txt. These methylation files list the position of each conversion in each read.

The methylation files are converted into a TSV file, with rows containing reads. A zero indicates the position was not converted and a one indicates the position was converted. 

Positions are then marked as follows with the color in parenthesis being its representation on the heatmap:
0 = not converted (grey)
1 = converted (green)
2 = non-cytosine base (white)
3 = converted CpG (blue) *only if -c option is used*
4 = non-converted CpG (black) *only if -c option is used*
5 = coverted CpG in R-loop (purple) *only if -c option is used*
9 = converted cytosine in R-loop (red)

Footprints are determined based on a sliding window method. The window size is the -l option entered by the user. If the percentage of converted cytosines is above the threshold specified by the user (-t option), all the converted cytosines are labelled as a footprint. 

The final heatmaps are saved in the same folder as the script and labelled as <gene name> <pos/neg> <conversion threshold> <CG if -c option is used>.pdf
A log file is also saved in the same folder as the script and named <gene name>logFile.txt. It contains statistics for mapping and peak calling.
