#!/usr/bin/env nextflow

//docker run -ti -v /home/etienne/masque/masque.nf:/masque.nf -v /home/etienne/masque/databases:/databases -v /home/etienne/masque/test/data:/mydata shaman_nextflow




// Colors for strings
ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";
def tored = {  str -> ANSI_RED + str + ANSI_RESET }
def toblack = {  str -> ANSI_BLACK + str + ANSI_RESET }
def togreen = {  str -> ANSI_GREEN + str + ANSI_RESET }
def toyellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def toblue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def tocyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def topurple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def towhite = {  str -> ANSI_WHITE + str + ANSI_RESET }



if (params.test){
	// FOR TEST
	
	exit 1
}


// HELP DISPLAY
params.help=false
if (params.size()==1 || params.help) {
	println("")
	println("Usage:")
	println("16S/18S:   nextflow masque.nf     -i </path/to/input/> -o </path/to/result/>")
	println("23S/28S:   nextflow masque.nf  -l -i </path/to/input/> -o </path/to/result/>")
	println("ITS:       nextflow masque.nf  -f -i </path/to/input/> -o </path/to/result/>")
	println("Amplicon:  nextflow masque.nf     -a <amplicon file>   -o </path/to/result/>")
	println("")
	println("All parameters:")
	println("--i                       Provide </path/to/input/directory/>")
	println("--a                       Provide <amplicon file>")
	println("--o                       Provide </path/to/result/directory/>")
	println("--n                       Indicate <project-name>")
	println("                          (default: use the name of the input directory or meta)")
	println("--t                       Number of <thread>")
	println("                          (default: all cpu will be used)")
	println("--c                       Contaminant filtering [danio,human,mouse,mosquito,phi]")
	println("                          (default: human,phi)")
	println("--s                       Perform OTU clustering with swarm")
	println("--b                       Perform taxonomical annotation with blast")
	println("                          (default: vsearch)")
	println("--l                       Perform taxonomical annotation")
	println("                          against LSU databases: Silva/RDP")
	println("--f                       Perform taxonomical annotation")
	println("                          against ITS databases: Unite/Findley/Underhill/RDP")
	println("--minreadlength           Minimum read length take in accound in the study")
	println("                          (default: 35nt)")
	println("--minphred                Qvalue must lie between [0-40]")
	println("                          (default: minimum qvalue 20)")
	println("--minphredperc            Minimum allowed percentage of correctly called")
	println("                          nucleotides [0-100] (default: 80)")
	println("--nbMismatchMapping       Maximum number of mismatch when mapping end-to-end")
	println("                          against Human genome and Phi174 genome")
	println("                          (default: 1 mismatch is accepted)")
	println("--maxoverlap              Maximum overlap when paired reads are considered")
	println("                          (default: 200)")
	println("--minoverlap              Minimum overlap when paired reads are considered")
	println("                          (default: 50)")
	println("--minampliconlength       Minimum amplicon length (default: 64)")
	println("--minotusize              Indicate minimum OTU size (default: 4)")
	println("--prefixdrep              Perform prefix dereplication")
	println("                          (default: full length dereplication)")
	println("--chimeraslayerfiltering  Use ChimeraSlayer database for chimera filtering")
	println("                          (default: Perform a de novo chimera filtering)")
	println("--otudiffswarm            Number of difference accepted in an OTU with swarm")
	println("                          (default: 1)")
	println("--evalueTaxAnnot          Evalue threshold for taxonomical annotation with blast")
	println("                          (default: evalue=1E-5)")
	println("--maxTargetSeqs           Number of hit per OTU with blast (default: 1)")
	println("--identityThreshold       Identity threshold for taxonomical annotation with")
	println("                          vsearch (default: 0.75)")
	println("--conservedPosition       Percentage of conserved position in the multiple")
	println("                          alignment considered for phylogenetic tree")
	println("                          (default: 0.8)")
	println("--accurateTree            Accurate tree calculation with IQ-TREE instead of")
	println("                          FastTree (default: FastTree)")
	println("--help                    Print this help")
    exit 1
}



// ASSEMBLY PARAMETERS
params.i=""
params.a=""
params.o=""
// params.n : see below
params.t=1
params.c=["human", "phi"]
params.s=0
params.b=0
params.l=0
// params.f ??
params.minreadlength=35
params.minphred=20
params.minphredperc=80
params.nbMismatchMapping=1
params.maxoverlap=550
params.minoverlap=10
params.minampliconlength=64
params.minotusize=4
params.prefixdrep=0
params.chimeraslayerfiltering=0
params.otudiffswarm=1
params.evalueTaxAnnot="1E-5"
params.maxTargetSeqs=1
params.identityThreshold=0.75
params.conservedPosition=0.5
params.accurateTree=0
fungi=0
paired=0
// OLD parameters
// accurateTree=0
// amplicon=""
// blast_tax=0
// chimeraslayerfiltering=0
// conservedPosition=0.5
// contaminant=("human" "phi")
// evalueTaxAnnot="1E-5"
// fungi=0
// identityThreshold=0.75
// input_dir=""
// lsu=0
// maxTargetSeqs=1
// maxoverlap=550
// minampliconlength=64
// minotusize=4
// minoverlap=10
// minphred=20
// minphredperc=80
// minreadlength=35
// NbMismatchMapping=1
// NbProc=$(grep -c ^processor /proc/cpuinfo)
// otudiffswarm=1
// paired=0
// prefixdrep=0
// ProjectName=""
// swarm_clust=0



// ARGUMENTS CHECK
// Mandatory parameters
if (params.i=="" && params.a=="") {println tored('Please indicate an input directory or an amplicon file.'); exit 1}
if (params.i!="" && params.a!="") {println tored('Please indicate an input directory or an amplicon file, but not both'); exit 1}
if (params.o=="") {println tored('Please indicate the output directory.'); exit 1}
// Definition of params.n
if (params.i =~ /\/*([a-zA-Z-_\s]+)\/*$/) {params.n=(params.i =~ /\/*([a-zA-Z-_\s]+)\/*$/)[0][1]}
else if (params.a =~ /\/*([a-zA-Z-._]+\.[a-zA-Z~]+)$/) {params.n=(params.a =~ /\/*([a-zA-Z-._]+\.[a-zA-Z~]+)$/)[0][1]}
else {params.n='masque_run'}



// DATABASES
// Databases directory
db = '/databases'
// ChimeraSlayer reference database
// http://drive5.com/uchime/uchime_download.html
gold="$db/gold.fa"
// Alien sequences
alienseq="$db/alienTrimmerPF8contaminants.fasta"
// Filtering database
filterRef=["danio":"$db/danio_rerio.fna", "human":"$db/homo_sapiens.fna", "mosquito":"$db/anopheles_stephensi.fna", "mouse":"$db/mus_musculus.fna", "phi":"$db/NC_001422.fna"]
// Findley
// http://www.mothur.org/w/images/2/20/Findley_ITS_database.zip
findley="$db/ITSdb.findley.fasta"
// Greengenes
// ftp://greengenes.microbio.me/greengenes_release/gg_13_5/
greengenes="$db/gg_13_5.fasta"
greengenes_taxonomy="$db/gg_13_5_taxonomy.txt"
// Silva
// http://www.arb-silva.de/no_cache/download/archive/release_123/Exports/
silva="$db/SILVA_128_SSURef_Nr99_tax_silva.fasta"
silvalsu="$db/SILVA_128_LSURef_tax_silva.fasta"
underhill="$db/THFv1.3.sequence.fasta"
underhill_taxonomy="$db/THFv1.3.tsv"
unite="$db/sh_general_release_dynamic_s_20.11.2016.fasta"



// PROGRAMS
// Path to the programs
bin='/usr/local/bin'
// AlienTrimmer
alientrimmer="java -jar $bin/AlienTrimmer_0.4.0/src/AlienTrimmer.jar"
// Biom
biom="biom"
// Blastn
blastn="$bin/ncbi-blast-2.5.0+/bin/blastn"
// BMGE ftp://ftp.pasteur.fr/pub/gensoft/projects/BMGE/
BMGE="java -jar $bin/BMGE-1.12/BMGE.jar"
// Bowtie2
bowtie2="$bin/bowtie2-2.2.9/bowtie2"
// Extract fasta
extract_fasta="$bin/extract_fasta/extract_fasta.py"
// Extract result
extract_result="$bin/extract_result/extract_result.py"
// Fastq2fasta
fastq2fasta="$bin/fastq2fasta/fastq2fasta.py"
// Fastqc
fastqc="$bin/FastQC/fastqc"
// Fasttree
FastTreeMP="$bin/FastTree-2.1.9/FastTree"
// FLASH
flash="$bin/FLASH-1.2.11/flash"
// mafft
mafft="$bin/mafft-linux64/mafft.bat"
// get_taxonomy
get_taxonomy="$bin/get_taxonomy/get_taxonomy.py"
// IQ-TREE
iqtree="$bin/iqtree-omp-1.5.1-Linux/bin/iqtree-omp"
// rename_otu
rename_otu="$bin/rename_otu/rename_otu.py"
// rdp classifier
rdp_classifier="java -jar $bin/rdp_classifier_2.12/dist/classifier.jar"
// swarm
swarm="$bin/swarm_bin/bin/swarm"
// swarm2vsearch
swarm2vsearch="$bin/swarm2vsearch/swarm2vsearch.py"
// vsearch
vsearch="$bin/vsearch_bin/bin/vsearch"




















