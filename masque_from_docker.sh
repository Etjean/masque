#!/bin/bash
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html
# ------------------------------------------------------------------
# Author: Amine Ghozlane (amine.ghozlane@pasteur.fr)
# Title:  16S-18S-23S-28S-ITS pipeline
# Description: De novo 16S-18S-23S-28S-ITS pipeline assignation
# ------------------------------------------------------------------

function say_parameters {
    # echo blue
    echo -e "\e[34m# $1\e[0m" >&2
}

function say {
    # echo green
    echo -e "\033[1;32m* $1\033[0m" >&2
}

function error {
    # echo red
    echo -e  "\e[31m* $1\e[0m" >&2
}

function check_log {
    # check if log file is not empty
    if [ -s $1 ]
    then
        error "$1 is not empty !"
        exit 1
    fi
}

function check_integer {
    if [[ $1 != [0-9]* ]]
    then
        error "\"$1\" is not an integer value"
        exit 1
    fi
}

function check_file {
    # Check if result is well produced
    if [ ! -f $1 ] && [ ! -s $1 ]
    then
        error "File \"$1\" does not exist or is empty !"
        exit 1
    fi
}

function check_dir {
    # Check if directory doesnt exist
    if [ ! -d $1 ]
    then
        mkdir $1
        if [ ! -d $1 ]
        then 
            error "The program cannot create the directory \"$1\""
            exit 1
        fi
    fi
}

function check_name {
    if [ "$1" = "" ]
    then
        error "Error failed to parse the sample name:"
        echo "sample name=$1"
        exit 1
    fi
}

function check_soft {
 if type $1 > /dev/null  2>&1;
 then
    echo $1
 else
    echo $2
 fi
}

function get_reverse_file {
    testname=$(echo $1|sed -r "s:R1_([0-9]+).(f*q*):R2_\1.\2:g")
    if [ ! -f "$testname" ] || [ "$testname" == "$1" ] 
    then
        testname=$(echo $1|sed -r "s:R1\.(f):R2.\1:g")
    fi
    if [ ! -f "$testname" ]  || [ "$testname" == "$1" ]
    then
        testname=$(echo $1|sed -r "s:R1\.(f):R2.\1:g")
    fi
    if [ ! -f "$testname" ]  || [ "$testname" == "$1" ]
    then
        testname=$(echo $1|sed -r "s:R1:R2:g")
    fi
    if [ ! -f "$testname" ]  || [ "$testname" == "$1" ]
    then
        error "The program is unable to find the corresdponding R2 pair of $1. Please consider to simplify the pair name (like name_R1.fastq and name_R2.fastq) with no point inside the name. R1 and R2 tag serve to identify pair."
        exit 1
    fi
    echo "$testname" 
}

display_help() {
    if [ "$1" -eq "0" ]
    then
        printf "%-10s $0 %3s %-30s %-35s\n" "16S/18S:" "" "-i </path/to/input/directory/>" "-o </path/to/result/directory/>"
        printf "%-10s $0 %3s %-30s %-35s\n" "23S/28S:" "-l" "-i </path/to/input/directory/>" "-o </path/to/result/directory/>"
        printf "%-10s $0 %3s %-30s %-35s\n" "ITS:" "-f" "-i </path/to/input/directory/>" "-o </path/to/result/directory/>"
        printf "%-10s $0 %3s %-30s %-35s\n" "Amplicon:" "" "-a <amplicon file>" "-o </path/to/result/directory/>"
        echo -e "\e[1;34mAll parameters:\e[m"
        printf "%-25s %-30s\n" "-i" "Provide </path/to/input/directory/>"
        printf "%-25s %-30s\n" "-a" "Provide <amplicon file>"
        printf "%-25s %-30s\n" "-o" "Provide </path/to/result/directory/>"
        printf "%-25s %-30s\n" "-n" "Indicate <project-name> (default: use the name of the input directory or meta)"
        printf "%-25s %-30s\n" "-t" "Number of <thread> (default all cpu will be used)"
        printf "%-25s %-30s\n" "-c" "Contaminant filtering [danio,human,mouse,mosquito,phi] (Default: human,phi)"
        printf "%-25s %-30s\n" "-s" "Perform OTU clustering with swarm"
        printf "%-25s %-30s\n" "-b" "Perform taxonomical annotation with blast (Default vsearch)"
        printf "%-25s %-30s\n" "-l" "Perform taxonomical annotation against LSU databases: Silva/RDP"
        printf "%-25s %-30s\n" "-f" "Perform taxonomical annotation against ITS databases: Unite/Findley/Underhill/RDP"
        printf "%-25s %-30s\n" "--minreadlength" "Minimum read length take in accound in the study (Default 35nt)"
        printf "%-25s %-30s\n" "--minphred" "Qvalue must lie between [0-40] (Default minimum qvalue 20)"
        printf "%-25s %-30s\n" "--minphredperc" "Minimum allowed percentage of correctly called nucleotides [0-100] (Default 80)"
        printf "%-25s %-30s\n" "--NbMismatchMapping" "Maximum number of mismatch when mapping end-to-end against Human genome and Phi174 genome (Default 1 mismatch is accepted)"
        printf "%-25s %-30s\n" "--maxoverlap" "Maximum overlap when paired reads are considered (Default 200)"
        printf "%-25s %-30s\n" "--minoverlap" "Minimum overlap when paired reads are considered (Default 50)"
        printf "%-25s %-30s\n" "--minampliconlength" "Minimum amplicon length (Default 64)"
        printf "%-25s %-30s\n" "--minotusize" "Indicate minimum OTU size (Default 4)"
        printf "%-25s %-30s\n" "--prefixdrep" "Perform prefix dereplication (Default: full length dereplication)"
        printf "%-25s %-30s\n" "--chimeraslayerfiltering" "Use ChimeraSlayer database for chimera filtering (Default Perform a de novo chimera filtering)"
        printf "%-25s %-30s\n" "--otudiffswarm" "Number of difference accepted in an OTU with swarm (Default 1)"
        printf "%-25s %-30s\n" "--evalueTaxAnnot" "Evalue threshold for taxonomical annotation with blast (Default evalue=1E-5)"
        printf "%-25s %-30s\n" "--maxTargetSeqs" "Number of hit per OTU with blast (Default 1)"
        printf "%-25s %-30s\n" "--identityThreshold" "Identity threshold for taxonomical annotation with vsearch (Default 0.75)"
        printf "%-25s %-30s\n" "--conservedPosition" "Percentage of conserved position in the multiple alignment considered for phylogenetic tree (Default 0.8)"
        printf "%-25s %-30s\n" "--accurateTree" "Accurate tree calculation with IQ-TREE instead of FastTree (Default FastTree)"
    else
        display_parameters
    fi
    exit
}

display_parameters() {
    # Display the parameters of the analysis
    say_parameters "Project name [-n]:"                                          
    echo $ProjectName  >&2
    if [ "$input_dir" != "" ]
    then
        say_parameters "Sample input [-i]:"
        echo $input_dir >&2
    elif [ "amplicon" != "" ]
    then
        say_parameters "Amplicon input [-a]:"
        echo $amplicon >&2
    fi
    say_parameters "Result output [-o]:"
    echo $resultDir >&2
    say_parameters "Number of threads [-t]:" >&2
    echo "$NbProc processes will be used" >&2
    say_parameters "Read filtering:" >&2
    echo """Minimum read length [--minreadlength]= $minreadlength
Minimum phred quality [--minphred]= $minphred
Minimum allowed percentage of correctly called nucleotides [--minphredperc]= $minphredperc
Minimum number of mistach for the filtering [--NbMismatchMapping]= $NbMismatchMapping
Filtering databases= ${contaminant[@]}""">&2
    if [ "$paired" -eq "1" ]
    then
        say_parameters "Merge reads parameters:"
        echo """Maxoverlap [--maxoverlap]= $maxoverlap
Minoverlap [--minoverlap]= $minoverlap""" >&2
    fi
    say_parameters "OTU process:" >&2
    if [ "$prefixdrep" -eq "1" ]
    then
        echo "Dereplication is in prefix mode [--prefixdrep]" >&2
    else
        echo "Dereplication is in full length mode" >&2
    fi
    echo """Minimum length of an amplicon [--minampliconlength]= $minampliconlength
Minimum size of an OTU for singleton removal [--minotusize]= $minotusize""" >&2
    if [ "$chimeraslayerfiltering" -eq "1" ]
    then
        echo "Chimera filtering use chimera slayer database for filtering [--chimeraslayerfiltering]" >&2
    else
        echo "Chimera filtering is in de novo mode" >&2
    fi
    if [ "$swarm_clust" -eq "1" ]
    then
        echo """Clustering is performed with swarm [-s]
Number of difference accepted in an OTU with swarm [--otudiffswarm]= $otudiffswarm""">&2
    else
        echo "Clustering is performed with vsearch" >&2
    fi
    if [ "$fungi" -eq "1" ]
    then
        say_parameters "Fungi annotation [-f]" >&2
    elif [ "$lsu" -eq "1" ]
    then
        say_parameters "23S/28S annotation [-l]" >&2
    else
        say_parameters "16S/18S annotation" >&2
    fi
    if [ "$blast_tax" -eq "0" ]
    then
        echo "Identity threshold with vsearch [--identityThreshold]= $identityThreshold" >&2
    else
        echo """E-value with blast [--evalueTaxAnnot]= $evalueTaxAnnot
Maximum number of targets with blast [--maxTargetSeqs]= $maxTargetSeqs""" >&2
    fi
    echo "Conserved position for alignment[--conservedPosition]= $conservedPosition" >&2
    if [ "$accurateTree" -eq "1" ]
    then
        echo "Tree generated in accurate mode with IQ-TREE [--accurateTree]"
    else
        echo "Tree generated in fast mode with FastTree"
    fi
}

function timer()
{
   if [[ $# -eq 0 ]]
   then
         echo $(date '+%s')
   else
      local  stime=$1
      etime=$(date '+%s')
      if [[ -z "$stime" ]]; then stime=$etime; fi
      dt=$((etime - stime))
      ds=$((dt % 60))
      dm=$(((dt / 60) % 60))
      dh=$((dt / 3600))
      printf '%d:%02d:%02d' $dh $dm $ds
  fi
}

SCRIPTPATH=$(dirname "${BASH_SOURCE[0]}")

#############
# Databases #
#############
# ChimeraSlayer reference database
#http://drive5.com/uchime/uchime_download.html
gold="$SCRIPTPATH/databases/gold.fa"
# Alien sequences
alienseq="$SCRIPTPATH/databases/alienTrimmerPF8contaminants.fasta"
# Filtering database
declare -A filterRef
filterRef=(["danio"]="$SCRIPTPATH/databases/danio_rerio.fna" ["human"]="$SCRIPTPATH/databases/homo_sapiens.fna" ["mosquito"]="$SCRIPTPATH/databases/anopheles_stephensi.fna" ["mouse"]="$SCRIPTPATH/databases/mus_musculus.fna" ["phi"]="$SCRIPTPATH/databases/NC_001422.fna")
# Findley
# http://www.mothur.org/w/images/2/20/Findley_ITS_database.zip
findley="$SCRIPTPATH/databases/ITSdb.findley.fasta"
# Greengenes
# ftp://greengenes.microbio.me/greengenes_release/gg_13_5/
greengenes="$SCRIPTPATH/databases/gg_13_5.fasta"
#greengenes="/local/databases/fasta/greengenes.fa"
greengenes_taxonomy="$SCRIPTPATH/databases/gg_13_5_taxonomy.txt"
# RDP
#http://rdp.cme.msu.edu/misc/resources.jsp
# rdp="$SCRIPTPATH/databases/rdp_11_4.fa"
# Silva
#http://www.arb-silva.de/no_cache/download/archive/release_123/Exports/
silva="$SCRIPTPATH/databases/SILVA_128_SSURef_Nr99_tax_silva.fasta"
silvalsu="$SCRIPTPATH/databases/SILVA_128_LSURef_tax_silva.fasta"
underhill="$SCRIPTPATH/databases/THFv1.3.sequence.fasta"
underhill_taxonomy="$SCRIPTPATH/databases/THFv1.3.tsv"
unite="$SCRIPTPATH/databases/sh_general_release_dynamic_s_20.11.2016.fasta"

#######################
# Assembly Parameters #
#######################
accurateTree=0
amplicon=""
blast_tax=0
chimeraslayerfiltering=0
conservedPosition=0.5
contaminant=("human" "phi")
evalueTaxAnnot="1E-5"
fungi=0
identityThreshold=0.75
input_dir=""
lsu=0
maxTargetSeqs=1
maxoverlap=550
minampliconlength=64
minotusize=4
minoverlap=10
minphred=20
minphredperc=80
minreadlength=35
NbMismatchMapping=1
NbProc=$(grep -c ^processor /proc/cpuinfo)
otudiffswarm=1
paired=0
prefixdrep=0
ProjectName=""
swarm_clust=0

############
# Programs #
############
# AlienTrimmer
alientrimmer=$(check_soft "AlienTrimmer" "java -jar $SCRIPTPATH/AlienTrimmer_0.4.0/src/AlienTrimmer.jar")
# Biom
biom="biom"
# Blastn
blastn=$(check_soft "blastn" "$SCRIPTPATH/ncbi-blast-2.5.0+/bin/blastn")
# BMGE ftp://ftp.pasteur.fr/pub/gensoft/projects/BMGE/
BMGE=$(check_soft "BMGE" "java -jar $SCRIPTPATH/BMGE-1.12/BMGE.jar")
# Bowtie2
bowtie2=$(check_soft "bowtie2" "$SCRIPTPATH/bowtie2-2.2.9/bowtie2")
# Extract fasta
extract_fasta="$SCRIPTPATH/extract_fasta/extract_fasta.py"
# Extract result
extract_result="$SCRIPTPATH/extract_result/extract_result.py"
# Fastq2fasta
fastq2fasta="$SCRIPTPATH/fastq2fasta/fastq2fasta.py"
# Fastqc
fastqc=$(check_soft "fastqc" "$SCRIPTPATH/FastQC/fastqc")
# Fasttree
FastTreeMP=$(check_soft "FastTree" "$SCRIPTPATH/FastTree-2.1.9/FastTree")
# FLASH
flash=$(check_soft "flash" "$SCRIPTPATH/FLASH-1.2.11/flash")
# mafft
mafft=$(check_soft "mafft" "$SCRIPTPATH/mafft-linux64/mafft.bat")
# get_taxonomy
get_taxonomy="$SCRIPTPATH/get_taxonomy/get_taxonomy.py"
# IQ-TREE
iqtree=$(check_soft "iqtree-omp" "$SCRIPTPATH/iqtree-omp-1.5.1-Linux/bin/iqtree-omp") 
# otu_tab_size
#otu_tab_size="$SCRIPTPATH/otu_tab_size/otu_tab_size.py"
# rename_otu
rename_otu="$SCRIPTPATH/rename_otu/rename_otu.py"
# rdp classifier
rdp_classifier=$(check_soft "classifier" "java -jar $SCRIPTPATH/rdp_classifier_2.12/dist/classifier.jar")
# swarm
swarm=$(check_soft "swarm" "$SCRIPTPATH/swarm_bin/bin/swarm")
# swarm2vsearch
swarm2vsearch="$SCRIPTPATH/swarm2vsearch/swarm2vsearch.py"
# uc2otutab
#uc2otutab="$SCRIPTPATH/usearch_python_scripts/uc2otutab.py"
# usearch
#usearch="$SCRIPTPATH/usearch8.1.1756_i86linux32"
#usearch -makeudb_utax 16s_ref.fa -output 16s_ref.udb -report 16s_report.txt
# vsearch
#vsearch="$SCRIPTPATH/vsearch-1.4.1-linux-x86_64/bin/vsearch"
vsearch=$(check_soft "vsearch" "$SCRIPTPATH/vsearch_bin/bin/vsearch")

########
# Main #
########
# Execute getopt on the arguments passed to this program, identified by the special character $@
PARSED_OPTIONS=$(getopt -n "$0"  -o hi:o:r:t:a:sblfn:c: --long "help,input_dir:,output:,thread:,minampliconlength:,maxoverlap:,maxTargetSeqs:,minotusize:,minoverlap:,minphred:,minphredperc:,minreadlength:,identityThreshold:,evalueTaxAnnot:,NbMismatchMapping:,amplicon:,swarm,blast,fungi,name:,prefixdrep,chimeraslayerfiltering,conservedPosition:,accurateTree,contaminant:"  -- "$@")

#Check arguments
if [ $# -eq 0 ]
then
    display_help 0
fi

#Bad arguments, something has gone wrong with the getopt command.
if [ $? -ne 0 ];
then
    display_help 1
    exit 1
fi

# A little magic, necessary when using getopt.
eval set -- "$PARSED_OPTIONS"

# Get Cmd line arguments depending on options
while true;
do
  case "$1" in
    -h|--help)
        display_help 0
        shift;;
    -i|--input_dir)
        input_dir=$2
        shift 2;;
    -o|--output)
        resultDir=$2
        readsDir=$resultDir/reads/
        logDir=$resultDir/log/
        errorlogDir=$resultDir/error_log/
        check_dir $resultDir
        check_dir $logDir
        check_dir $errorlogDir
        check_dir $readsDir
        shift 2;;
    -t|--thread)
        check_integer $2
        NbProc=$2
        shift 2;;
    -b|--blast)
        blast_tax=1
        shift;;
    -a|amplicon)
        check_file $2
        amplicon=$2
        shift 2;;
    -s|swarm)
        swarm_clust=1
        shift ;;
    -l|lsu)
        lsu=1
        shift ;;
    -f|fungi)
        fungi=1
        shift ;;
    -n|--name)
        ProjectName=$2
        shift 2;;
    -c|--contaminant)
        contaminant=(${2//,/ })
        shift 2;;
    --chimeraslayerfiltering)
        chimeraslayerfiltering=1
        shift ;;
    --minampliconlength)
        check_integer $2
        minampliconlength=$2
        shift 2;;
    --maxoverlap)
        check_integer $2
        maxoverlap=$2
        shift 2;;
    --minotusize)
        check_integer $2
        minotusize=$2
        shift 2;;
    --minoverlap)
        check_integer $2
        minoverlap=$2
        shift 2;;
    --minphred)
        check_integer $2
        minphred=$2
        shift 2;;
    --minphredperc)
        check_integer $2
        minphredperc=$2
        shift 2;;
    --minreadlength)
        check_integer $2
        minreadlength=$2
        shift 2;;
    --otudiffswarm)
        check_integer $2
        otudiffswarm=$2
        shift 2;;
    --prefixdrep)
        prefixdrep=1
        shift ;;
    --identityThreshold)
        identityThreshold=$2
        shift 2;;
    --evalueTaxAnnot)
        evalueTaxAnnot=$2
        shift 2;;
    --NbMismatchMapping)
        check_integer $2
        NbMismatchMapping=$2
        shift 2;;
    --maxTargetSeqs)
        check_integer $2
        maxTargetSeqs=$2
        shift 2;;
    --conservedPosition)
        conservedPosition=$2
        shift 2;;
    --accurateTree)
        accurateTree=1
        shift ;;
    --)
      shift
      break;;
  esac
done


if [ "$resultDir" = "" ]
then
    error "Please indicate the output directory."
    exit 1
fi

if [ -d "$input_dir" ]
then
    if [ "$ProjectName" = "" ]
    then
        dir_name=$(basename "$input_dir")
        if [ "$dir_name" == "." ]
        then
            ProjectName="meta"
        else
            ProjectName=$(basename "$input_dir")
        fi
    fi
    amplicon="${resultDir}/${ProjectName}_extendedFrags.fasta"
elif [ -f "$amplicon" ]
then
    if [ "$ProjectName" = "" ]
    then
        ProjectName=$(basename $(dirname "$amplicon"))
    fi
else
    error "Error no input dir and no amplicon file given. Please given me one or the other."
    exit 1
fi


# display parameters
display_parameters

# Start timer
say "Start analysis"
wall_time=$(timer)



say "Start working on reads"
all_start_time=$(timer)

if [ -d "$input_dir" ]
then
    list_product_fa=""

    nb_samples=$(ls $input_dir/*R1*.{fastq,fq,fastq.gz,fq.gz} -1  2>/dev/null |wc -l)
    num_sample=0
    if [ "$nb_samples" -eq "0" ]
    then
        nb_samples=$(ls $input_dir/*.{fastq,fq,fastq.gz,fq.gz} -1  2>/dev/null |wc -l)
        for input in $(ls $input_dir/*.{fastq,fq,fastq.gz,fq.gz}  2>/dev/null )
        do
            let "num_sample=$num_sample+1"
            # Get the sample name
            filename=$(basename "$input"|sed "s:.gz::g")
            SampleName="${filename%.*}"
            check_name $SampleName
            list_product_fa+="${resultDir}/reads/${SampleName}_alien_filt.fasta "
            # Triming
            if [ -f "$input" ] && [ ! -f "${readsDir}/${SampleName}_alien.fastq" ] && [ ! -f "${readsDir}/${SampleName}_alien_filt.fastq" ]
            then
                say "$num_sample/$nb_samples - Triming reads with Alientrimmer"
                start_time=$(timer)
                filename=$(basename "$input")
                extension=".${filename##*.}"
                if [ "$extension" == ".gz" ]
                then
                    gunzip -c $input > ${readsDir}/${SampleName}_tmp.fastq
                    input="${readsDir}/${SampleName}_tmp.fastq"
                fi
                $alientrimmer -i $input -o ${readsDir}/${SampleName}_alien.fastq -c $alienseq -l $minreadlength -p $minphredperc -q $minphred > ${logDir}/log_alientrimmer_${SampleName}.txt 2> ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
                check_file ${readsDir}/${SampleName}_alien.fastq
                check_log ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
                rm -f ${readsDir}/${SampleName}_tmp.fastq
                say "$num_sample/$nb_samples - Elapsed time to trim with Alientrimmer: $(timer $start_time)"
            fi
            # Filtering reads against contaminant db
            let "essai=0";
            if [ ! -f "${readsDir}/${SampleName}_alien_filt.fastq" ]
            then
                let "totalFilter=${#contaminant[@]}-1"
                for db in ${contaminant[@]}
                do
                    # Set to lowercase
                    db=$(echo "${db,,}")
                    if [ "${filterRef[$db]}" == "" ]
                    then
                        error "$db does not belong to the list of possible contaminant [danio,human,mouse,mosquito,phi]"
                        echo "${filterRef[$db]}"
                        exit 1
                    else
                        let "num=$essai-1";
                        if [ ! -f "${readsDir}/${SampleName}_${contaminant[${essai}]}_${essai}.fastq" ] && [ -f "${readsDir}/${SampleName}_${contaminant[${num}]}_${num}.fastq" ] && [ ! -f "${readsDir}/${SampleName}_${contaminant[${totalFilter}]}_${totalFilter}.fastq" ]
                        then
                                say "$num_sample/$nb_samples - Filter reads against $db"
                                start_time=$(timer)
                                # Next mapping
                                $bowtie2  -q -N $NbMismatchMapping -p $NbProc -x ${filterRef[$db]} -U ${readsDir}/${SampleName}_${contaminant[${num}]}_${num}.fastq -S /dev/null --un ${readsDir}/${SampleName}_${contaminant[${essai}]}_${essai}.fastq -t --end-to-end --very-fast  > ${logDir}/log_mapping_${SampleName}_${contaminant[${essai}]}_${essai}.txt 2>&1
                                check_file ${readsDir}/${SampleName}_${contaminant[${essai}]}_${essai}.fastq
                                # Remove old file
                                rm -f ${readsDir}/${SampleName}_${num}.fastq
                                say "$num_sample/$nb_samples - Elapsed time to filter reads with $db reference: $(timer $start_time)"
                        elif [ -f "${readsDir}/${SampleName}_alien.fastq" ] && [ "$essai" -eq "0" ] && [ ! -f "${readsDir}/${SampleName}_${contaminant[${essai}]}_${essai}.fastq" ] && [ ! -f "${readsDir}/${SampleName}_${contaminant[${totalFilter}]}_${totalFilter}.fastq" ]
                        then
                                say "$num_sample/$nb_samples - Filter reads against $db"
                                start_time=$(timer)
                                # First mapping
                                $bowtie2 -q -N $NbMismatchMapping -p $NbProc -x ${filterRef[$db]}  -U ${readsDir}/${SampleName}_alien.fastq  -S /dev/null --un ${readsDir}/${SampleName}_${contaminant[${essai}]}_${essai}.fastq -t --end-to-end --very-fast  > ${logDir}/log_mapping_${SampleName}_${contaminant[${essai}]}_${essai}.txt 2>&1
                                check_file ${readsDir}/${SampleName}_${contaminant[${essai}]}_${essai}.fastq
                                rm -f ${readsDir}/${SampleName}_alien.fastq
                                say "$num_sample/$nb_samples - Elapsed time to filter reads with $db reference: $(timer $start_time)"
                        fi
                        let "essai=$essai+1";
                    fi
                done
                mv ${readsDir}/${SampleName}_${contaminant[${totalFilter}]}_${totalFilter}.fastq ${readsDir}/${SampleName}_alien_filt.fastq
            fi
            # Quality control
            if [ -f "${readsDir}/${SampleName}_alien_filt.fastq" ] && [ ! -f "${readsDir}/${SampleName}_alien_filt_fastqc.html" ]
            then
                say "$num_sample/$nb_samples - Quality control with Fastqc"
                start_time=$(timer)
                $fastqc ${readsDir}/${SampleName}_alien_filt.fastq --nogroup -q 2> ${errorlogDir}/error_log_fastqc_${SampleName}.txt
                check_file ${readsDir}/${SampleName}_alien_filt_fastqc.html
                check_log ${errorlogDir}/error_log_fastqc_${SampleName}.txt
                say "$num_sample/$nb_samples - Elapsed time with Fastqc: $(timer $start_time)"
            fi
            # Convert to fasta with the right name
            if [ -f "${readsDir}/${SampleName}_alien_filt.fastq" ] && [ ! -f "${readsDir}/${SampleName}_alien_filt.fasta" ]
            then
                say "$num_sample/$nb_samples - Convert fastq to fasta with fastq2fasta"
                start_time=$(timer)
                $fastq2fasta -i ${readsDir}/${SampleName}_alien_filt.fastq -o ${readsDir}/${SampleName}_alien_filt.fasta -s ${SampleName}  2> ${errorlogDir}/error_log_fastq2fasta_${SampleName}.txt
                check_file ${readsDir}/${SampleName}_alien_filt.fasta
                check_log ${errorlogDir}/error_log_fastq2fasta_${SampleName}.txt
                say "$num_sample/$nb_samples - Elapsed time with fastq2fasta: $(timer $start_time)"
            fi
        done
    else
        paired=1
        for r1_file in $(ls $input_dir/*R1*.{fastq,fq,fastq.gz,fq.gz}  2>/dev/null )
        do
            let "num_sample=$num_sample+1"
            input1=$r1_file
            input2=$(get_reverse_file $input1)
            #input2=$(echo $r1_file|sed "s:R1:R2:g")
            check_file $input1
            check_file $input2
            # Get the sample name
            filename=$(basename "$input1")
            #SampleName=$(echo "${filename%.*}" |sed "s:_L001:@:g"|cut -f 1 -d"@")
            SampleName=$(echo "${filename%.*}" |sed "s:_R1:@:g"|cut -f 1 -d"@")
            check_name $SampleName
            list_product_fa+="${resultDir}/reads/${SampleName}_extendedFrags.fasta "
            # Trimming
            if [ -f "$input1" ] && [ -f "$input2" ] && [ ! -f "${readsDir}/${SampleName}_alien_f.fastq" ] && [ ! -f "${readsDir}/${SampleName}_alien_f_filt.fastq" ]
            then
                say "$num_sample/$nb_samples - Triming reads with Alientrimmer"
                start_time=$(timer)
                filename=$(basename "$input1")
                extension=".${filename##*.}"
                if [ "$extension" == ".gz" ]
                then
                    gunzip -c $input1 > ${readsDir}/${SampleName}_R1_tmp.fastq
                    gunzip -c $input2 > ${readsDir}/${SampleName}_R2_tmp.fastq
                    input1="${readsDir}/${SampleName}_R1_tmp.fastq"
                    input2="${readsDir}/${SampleName}_R2_tmp.fastq"
                fi
                $alientrimmer -if $input1 -ir $input2 -of ${readsDir}/${SampleName}_alien_f.fastq -or ${readsDir}/${SampleName}_alien_r.fastq -os ${readsDir}/${SampleName}_alien_s.fastq -c $alienseq -l $minreadlength -p $minphredperc -q $minphred > ${logDir}/log_alientrimmer_${SampleName}.txt  2> ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
                check_file ${readsDir}/${SampleName}_alien_f.fastq
                check_file ${readsDir}/${SampleName}_alien_r.fastq
                check_log ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
                rm -f ${readsDir}/${SampleName}_R1_tmp.fastq ${readsDir}/${SampleName}_R2_tmp.fastq
                say "$num_sample/$nb_samples - Elapsed time with Alientrimmer: $(timer $start_time)"
            fi
            # Filtering reads against contaminant db
            let "essai=0";
            let "totalFilter=${#contaminant[@]}-1"
            if [ ! -d "${readsDir}/filter_${#filterRef[@]}" ] && [ ! -f "${readsDir}/${SampleName}_alien_f_filt.fastq" ]
            then
                for db in ${contaminant[@]}
                do
                    # Set to lowercase
                    db=$(echo "${db,,}")
                    if [ "${filterRef[$db]}" == "" ]
                    then
                        error "$db does not belong to the list of possible contaminant [danio,human,mouse,mosquito,phi]"
                        exit 1
                    else
                        #TODO to finish
                        let "num=$essai-1";
                        if [ ! -d "${readsDir}/filter_${contaminant[${essai}]}_${essai}" ] && [ -d "${readsDir}/filter_${contaminant[${num}]}_${num}" ] || [ "$essai" -ne "0" ] 
                        then
                            say "$num_sample/$nb_samples - Filter reads against $db"
                            start_time=$(timer)
                            mkdir ${readsDir}/filter_${contaminant[${essai}]}_${essai}
                            # Next mapping
                            $bowtie2  -q -N $NbMismatchMapping -p $NbProc -x ${filterRef[$db]} -1 ${readsDir}/filter_${contaminant[${num}]}_${num}/un-conc-mate.1  -2 ${readsDir}/filter_${contaminant[${num}]}_${num}/un-conc-mate.2 -S /dev/null --un-conc ${readsDir}/filter_${contaminant[${essai}]}_${essai}/ -t --very-fast  > ${logDir}/log_mapping_${SampleName}_${contaminant[${essai}]}_${essai}.txt 2>&1
                            check_file ${readsDir}/filter_${contaminant[${essai}]}_${essai}/un-conc-mate.1
                            # Remove old file
                            rm -rf ${readsDir}/filter_${contaminant[${num}]}_${num}
                            say "$num_sample/$nb_samples - Elapsed time to filter reads with $db reference: $(timer $start_time)"
                        elif [ -f "${readsDir}/${SampleName}_alien_f.fastq" ] && [ -f "${readsDir}/${SampleName}_alien_r.fastq" ]  && [ "$essai" -eq "0" ] && [ ! -d "${readsDir}/filter_${contaminant[${num}]}_${essai}" ]
                        then
                            say "$num_sample/$nb_samples - Filter reads against $db"
                            start_time=$(timer)
                            mkdir  ${readsDir}/filter_${contaminant[${essai}]}_${essai}
                            # First mapping
                            $bowtie2 -q -N $NbMismatchMapping -p $NbProc -x ${filterRef[$db]}  -1 ${readsDir}/${SampleName}_alien_f.fastq -2 ${readsDir}/${SampleName}_alien_r.fastq -S /dev/null --un-conc ${readsDir}/filter_${contaminant[${essai}]}_${essai} -t --very-fast > ${logDir}/log_mapping_${SampleName}_${contaminant[${essai}]}_${essai}.txt 2>&1
                            check_file ${readsDir}/filter_${contaminant[${essai}]}_${essai}/un-conc-mate.1
                            rm -f ${readsDir}/${SampleName}_alien_f.fastq ${readsDir}/${SampleName}_alien_r.fastq ${readsDir}/${SampleName}_alien_s.fastq
                            say "$num_sample/$nb_samples - Elapsed time to filter reads with $db reference: $(timer $start_time)"
                        fi
                        let "essai=$essai+1";
                    fi
                done
                mv ${readsDir}/filter_${contaminant[${totalFilter}]}_${totalFilter}/un-conc-mate.1 ${readsDir}/${SampleName}_alien_f_filt.fastq
                mv ${readsDir}/filter_${contaminant[${totalFilter}]}_${totalFilter}/un-conc-mate.2 ${readsDir}/${SampleName}_alien_r_filt.fastq
                rmdir ${readsDir}/filter_${contaminant[${totalFilter}]}_${totalFilter}
            fi
            # Merging reads
            if [ -f "${readsDir}/${SampleName}_alien_f_filt.fastq" ] && [ -f "${readsDir}/${SampleName}_alien_r_filt.fastq" ] && [ ! -f "${readsDir}/${SampleName}.extendedFrags.fastq" ]
            then
                say "$num_sample/$nb_samples - Merging paired reads with FLASH"
                #say "$num_sample/$nb_samples - Merging paired reads with vsearch"
                start_time=$(timer)
                $flash ${readsDir}/${SampleName}_alien_f_filt.fastq ${readsDir}/${SampleName}_alien_r_filt.fastq -M $maxoverlap -m $minoverlap -d $readsDir/ -o $SampleName -t $NbProc  > ${logDir}/log_flash_${SampleName}.txt
                #$vsearch  --fastq_mergepairs ${readsDir}/${SampleName}_alien_f_filt.fastq --reverse ${readsDir}/${SampleName}_alien_r_filt.fastq --fastqout ${readsDir}/${SampleName}.extendedFrags.fastq --fastq_minovlen $minoverlap --threads $NbProc
                # --label_suffix " ;barcodelabel=${SampleName}" --fastaout ${readsDir}/${SampleName}_extendedFrags_tmp.fasta
                #cut -f 1,3 -d " " ${readsDir}/${SampleName}_extendedFrags_tmp.fasta | sed "s: ::g" >${readsDir}/${SampleName}_extendedFrags.fasta
                check_file ${readsDir}/${SampleName}.extendedFrags.fastq
                #check_file ${readsDir}/${SampleName}_extendedFrags.fasta
                #rm  ${readsDir}/${SampleName}_extendedFrags_tmp.fasta
                #say "$num_sample/$nb_samples - Elapsed time with vsearch: $(timer $start_time)"
                say "$num_sample/$nb_samples - Elapsed time with FLASH: $(timer $start_time)"
            fi
            # Quality control
            if [ -f "${readsDir}/${SampleName}.extendedFrags.fastq" ] && [ ! -f "${readsDir}/${SampleName}.extendedFrags_fastqc.html" ] 
            then
                say "$num_sample/$nb_samples - Quality control with Fastqc"
                start_time=$(timer)
                $fastqc ${readsDir}/${SampleName}.extendedFrags.fastq --nogroup -q 2> ${errorlogDir}/error_log_fastqc_${SampleName}.txt
                #$vsearch --fastq_stats ${readsDir}/${SampleName}.extendedFrags.fastq --log ${readsDir}/${SampleName}.log
                check_file ${readsDir}/${SampleName}.extendedFrags_fastqc.html
                check_log ${errorlogDir}/error_log_fastqc_${SampleName}.txt
                say "$num_sample/$nb_samples - Elapsed time with Fastqc: $(timer $start_time)"
            fi
            # Convert to fasta with the right name
            if [ -f "${readsDir}/${SampleName}.extendedFrags.fastq" ] && [ ! -f "${readsDir}/${SampleName}_extendedFrags.fasta" ]
            then
                say "$num_sample/$nb_samples - Convert fastq to fasta with fastq2fasta"
                start_time=$(timer)
                $fastq2fasta -i ${readsDir}/${SampleName}.extendedFrags.fastq -o ${readsDir}/${SampleName}_extendedFrags.fasta -s ${SampleName}  2> ${errorlogDir}/error_log_fastq2fasta_${SampleName}.txt
                check_file ${readsDir}/${SampleName}_extendedFrags.fasta
                check_log ${errorlogDir}/error_log_fastq2fasta_${SampleName}.txt
                say "$num_sample/$nb_samples - Elapsed time with fastq2fasta: $(timer $start_time)"
            fi
        done
    fi
    say "Elapsed time with read processing: $(timer $all_start_time)"
fi
# Combine all files
if [ ! -f "$amplicon" ]
then
    say "Combine fasta files"
    start_time=$(timer)
    cat $list_product_fa > $amplicon #${resultDir}/${ProjectName}_extendedFrags.fasta
    #check_file ${resultDir}/${ProjectName}_extendedFrags.fasta
    check_file $amplicon
    say "Elapsed time to combine fasta files: $(timer $start_time)"
fi

#if [ -f "${resultDir}/${ProjectName}_extendedFrags.fasta" ] && [ ! -f ${resultDir}/${ProjectName}_reads_vs_rdp.txt ]
#then
#    say "Classify reads with rdp"
#    start_time=$(timer)
#    $rdp_classifier classify  -q ${resultDir}/${ProjectName}_extendedFrags.fasta -o  ${resultDir}/${ProjectName}_reads_vs_rdp.txt
#    check_file ${resultDir}/${ProjectName}_reads_vs_rdp.txt
#    say "Elapsed time to rdp: $(timer $start_time)"
#fi
#[ -f "${resultDir}/${ProjectName}_extendedFrags.fasta" ]
if [ -f "$amplicon" ] && [ ! -f "${resultDir}/${ProjectName}_drep.fasta" ]
then
     say "Dereplication"
     start_time=$(timer)
     #$usearch -derep_fulllength ${resultDir}/${ProjectName}.extendedFrags.fasta -fastaout ${resultDir}/${ProjectName}_drep.fasta -sizeout 
     # -minseqlength 64
     #${resultDir}/${ProjectName}_extendedFrags.fasta
     if [ "$prefixdrep" -eq "1" ]
     then
        $vsearch --derep_prefix $amplicon -output ${resultDir}/${ProjectName}_drep.fasta -sizeout -minseqlength $minampliconlength --strand both
     else
        $vsearch --derep_fulllength $amplicon -output ${resultDir}/${ProjectName}_drep.fasta -sizeout -minseqlength $minampliconlength --strand both
     fi
     check_file ${resultDir}/${ProjectName}_drep.fasta
     say "Elapsed time to dereplicate: $(timer $start_time)"
fi

if [ -f "${resultDir}/${ProjectName}_drep.fasta" ] && [ ! -f "${resultDir}/${ProjectName}_sorted.fasta" ]
then
     say "Abundance sort and discard singletons"
     start_time=$(timer)
     #$usearch -sortbysize ${resultDir}/${ProjectName}_drep.fasta -fastaout ${resultDir}/${ProjectName}_sorted.fasta -minsize 4
     $vsearch -sortbysize ${resultDir}/${ProjectName}_drep.fasta -output ${resultDir}/${ProjectName}_sorted.fasta  -minsize $minotusize
 > ${logDir}/log_search_sort_${ProjectName}.txt 2>&1
     check_file ${resultDir}/${ProjectName}_sorted.fasta
     say "Elapsed time to sort: $(timer $start_time)"
fi

if [ -f "${resultDir}/${ProjectName}_sorted.fasta" ] &&  [ ! -f "${resultDir}/${ProjectName}_nochim.fasta" ]
then
     say "Chimera filtering using reference database"
     start_time=$(timer)
     #$usearch -uchime_ref ${resultDir}/${ProjectName}_otu.fasta -db $gold -strand plus -nonchimeras ${resultDir}/${ProjectName}_otu_nochim.fasta
     if [ "$chimeraslayerfiltering" -eq "0" ]
     then
         $vsearch --uchime_denovo ${resultDir}/${ProjectName}_sorted.fasta --strand both --nonchimeras ${resultDir}/${ProjectName}_nochim.fasta --chimeras ${resultDir}/${ProjectName}_chim.fasta
     else
        $vsearch --uchime_ref ${resultDir}/${ProjectName}_sorted.fasta --db $gold --strand both --nonchimeras ${resultDir}/${ProjectName}_nochim.fasta --chimeras ${resultDir}/${ProjectName}_chim.fasta
     fi
     check_file ${resultDir}/${ProjectName}_nochim.fasta 
     say "Elapsed time to filter chimera: $(timer $start_time)"
fi

#[ ! -f "${resultDir}/${ProjectName}_swarm_representant.fasta" ]
if [ -f "${resultDir}/${ProjectName}_nochim.fasta" ] && [ ! -f "${resultDir}/${ProjectName}_otu.fasta" ] && [ "$swarm_clust" -eq 0 ]
then
     say "OTU clustering with vsearch"
     start_time=$(timer)
     #$usearch -cluster_otus ${resultDir}/${ProjectName}_sorted.fasta -otus ${resultDir}/${ProjectName}_otu.fasta -uparseout ${resultDir}/${ProjectName}_uparse.txt -relabel OTU_ -sizein #-sizeout 
     # --relabel OTU_
     $vsearch --cluster_size ${resultDir}/${ProjectName}_nochim.fasta --id 0.97 --centroids ${resultDir}/${ProjectName}_otu_compl.fasta --sizein --strand both #--sizeout
     python $rename_otu -i ${resultDir}/${ProjectName}_otu_compl.fasta -o ${resultDir}/${ProjectName}_otu.fasta
     check_file ${resultDir}/${ProjectName}_otu.fasta
     say "Elapsed time to OTU clustering with vsearch: $(timer $start_time)"
fi

if [ -f "${resultDir}/${ProjectName}_nochim.fasta" ] && [ ! -f "${resultDir}/${ProjectName}_otu_compl.fasta" ] && [ "$swarm_clust" -eq "1" ]
then
    say "OTU clustering with swarm"
    start_time=$(timer)
    $swarm -t $NbProc -f -z -w ${resultDir}/${ProjectName}_otu_compl.fasta  -o ${resultDir}/${ProjectName}_swarm_clustering.txt -s ${resultDir}/${ProjectName}_swarm_stats.txt -u ${resultDir}/${ProjectName}_swarm_uclust.txt ${resultDir}/${ProjectName}_nochim.fasta
    check_file ${resultDir}/${ProjectName}_otu_compl.fasta 
    say "Elapsed time to OTU clustering with swarm: $(timer $start_time)"
fi

if [ -f "${resultDir}/${ProjectName}_otu_compl.fasta" ] && [ ! -f "${resultDir}/${ProjectName}_otu.fasta" ]  && [ "$swarm_clust" -eq "1" ]
then
     say "Extract OTU clustering with swarm2vsearch"
     start_time=$(timer)
     python $swarm2vsearch -i ${resultDir}/${ProjectName}_otu_compl.fasta   -c ${resultDir}/${ProjectName}_swarm_clustering.txt -o ${resultDir}/${ProjectName}_otu.fasta -oc ${resultDir}/${ProjectName}_otu_swarm_clustering.txt -u ${resultDir}/${ProjectName}_swarm_uclust.txt -ou ${resultDir}/${ProjectName}_otu_swarm_uclust.txt
     check_file ${resultDir}/${ProjectName}_otu.fasta
     say "Elapsed time with swarm2vsearch: $(timer $start_time)"
fi

#[ -f "${resultDir}/${ProjectName}_extendedFrags.fasta" ]
#if [ -f "${resultDir}/${ProjectName}_otu.fasta" ] && [ -f "$amplicon" ] &&  [ ! -f "${resultDir}/${ProjectName}_map.txt" ]
if [ -f "${resultDir}/${ProjectName}_otu.fasta" ] && [ -f "$amplicon" ] && [ ! -f "${resultDir}/${ProjectName}_otu_table.tsv" ]
then
    say "Map reads back to OTUs"
    start_time=$(timer)
    #$usearch -usearch_global ${resultDir}/${SampleName}_extendedFrags.fasta -db ${resultDir}/${SampleName}_otu_nochim.fasta -strand plus -id 0.97 -uc ${resultDir}/${SampleName}_map.txt
    #${resultDir}/${ProjectName}_extendedFrags.fasta
    #$vsearch -usearch_global $amplicon -db ${resultDir}/${ProjectName}_otu.fasta --strand both --id 0.97 -uc ${resultDir}/${ProjectName}_map.txt 
    $vsearch -usearch_global $amplicon -db ${resultDir}/${ProjectName}_otu.fasta --strand both --id 0.97 --otutabout ${resultDir}/${ProjectName}_otu_table.tsv --biomout ${resultDir}/${ProjectName}_count.biom
    #check_file ${resultDir}/${ProjectName}_map.txt
    check_file ${resultDir}/${ProjectName}_otu_table.tsv
    check_file ${resultDir}/${ProjectName}_count.biom
    say "Elapsed time to map reads: $(timer $start_time)"
fi


#if [ -f "${resultDir}/${ProjectName}_map.txt" ] && [ ! -f "${resultDir}/${ProjectName}_otu_table.tsv" ]
#then
#    say "Build OTUs table"
#    start_time=$(timer)
#    python $uc2otutab ${resultDir}/${ProjectName}_map.txt > ${resultDir}/${ProjectName}_otu_table.tsv
#    check_file ${resultDir}/${ProjectName}_otu_table.tsv
#    say "Elapsed time to build OTUs table: $(timer $start_time)"
#fi

#if [ -f "${resultDir}/${ProjectName}_otu_table.tsv" ] && [ ! -f "${resultDir}/${ProjectName}_count.biom" ]
#then
#    say "Convert to biom format"
#    start_time=$(timer)
#    $biom convert -i ${resultDir}/${ProjectName}_otu_table.tsv -o ${resultDir}/${ProjectName}_count.biom --table-type="OTU table" --to-json
#    check_file ${resultDir}/${ProjectName}_count.biom
#    say "Elapsed time to convert the count table to biom format: $(timer $start_time)"
#fi

#if [ -f "${resultDir}/${ProjectName}_otu_table.tsv" ] && [ -f "${resultDir}/${ProjectName}_otu.fasta" ]  && [ ! -f "${resultDir}/${ProjectName}_otu_table_wgl.tsv" ]
#then
#    say "Build OTUs table for gene length normalization"
#    start_time=$(timer)
#    python $otu_tab_size -i ${resultDir}/${ProjectName}_otu_table.tsv -g ${resultDir}/${ProjectName}_otu.fasta -o ${resultDir}/${ProjectName}_otu_table_wgl.tsv
#    check_file ${resultDir}/${ProjectName}_otu_table_wgl.tsv
#    say "Elapsed time to build OTUs table wgl: $(timer $start_time)"
#fi

if [ -f "${resultDir}/${ProjectName}_otu.fasta" ]
then
    if [ ! -f "${resultDir}/${ProjectName}_vs_rdp.tsv" ]
    then
        say "Assign taxonomy with rdp_classifier"
        start_time=$(timer)
        #$usearch -utax ${resultDir}/${ProjectName}_otu.fasta -db $rdp -strand both -taxconfs rdp_16s_short.tc -utaxout ${resultDir}/${ProjectName}_otu_tax_rdp.tsv -utax_cutoff 0.8
        #$vsearch  --usearch_global ${resultDir}/${ProjectName}_otu.fasta --db $rdp --id 0.9 --blast6out ${resultDir}/${ProjectName}_vs_rdp.tsv
        $rdp_classifier classify  -q ${resultDir}/${ProjectName}_otu.fasta -o  ${resultDir}/${ProjectName}_vs_rdp.tsv
        check_file ${resultDir}/${ProjectName}_vs_rdp.tsv
        #python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_silva.tsv -d $rdp -dtype rdp -o ${resultDir}/${ProjectName}_vs_rdp_annotation.tsv
        say "Elapsed time with rdp_classifier: $(timer $start_time)"
    fi

    # SILVA
    if [ ! -f "${resultDir}/${ProjectName}_vs_silva_id_${identityThreshold}.tsv" ] && [ "$blast_tax" -eq "0" ]  && [ "$fungi" -eq "0" ]
    then
        say "Assign taxonomy against silva with vsearch"
        start_time=$(timer)
        #$usearch -utax ${resultDir}/${ProjectName}_otu.fasta -db $silva -strand both -taxconfs silva_16s_short.tc -utaxout ${resultDir}/${ProjectName}_otu_tax_silva.tsv -utax_cutoff 0.8
        if [ "$lsu" -eq "1" ]
        then
            $vsearch --usearch_global ${resultDir}/${ProjectName}_otu.fasta --db $silvalsu --id $identityThreshold --blast6out ${resultDir}/${ProjectName}_vs_silva_id_${identityThreshold}.tsv --strand both
        else
            $vsearch --usearch_global ${resultDir}/${ProjectName}_otu.fasta --db $silva --id $identityThreshold --blast6out ${resultDir}/${ProjectName}_vs_silva_id_${identityThreshold}.tsv --strand both
        fi
        #check_file ${resultDir}/${ProjectName}_vs_silva_id_${identityThreshold}.tsv
        say "Elapsed time with vsearch: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_vs_silva_id_${identityThreshold}.tsv" ] && [ ! -f "${resultDir}/${ProjectName}_vs_silva_annotation_id_${identityThreshold}.tsv" ]
    then
        say "Extract vsearch - silva annotation with get_taxonomy"
        start_time=$(timer)
        if [ "$lsu" -eq "1" ]
        then
            python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_silva_id_${identityThreshold}.tsv -u ${resultDir}/${ProjectName}_otu.fasta -d $silvalsu -o ${resultDir}/${ProjectName}_vs_silva_annotation_id_${identityThreshold}.tsv -ob ${resultDir}/${ProjectName}_vs_silva_annotation_id_${identityThreshold}.biomtsv
        else
            python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_silva_id_${identityThreshold}.tsv -u ${resultDir}/${ProjectName}_otu.fasta -d $silva -o ${resultDir}/${ProjectName}_vs_silva_annotation_id_${identityThreshold}.tsv -ob ${resultDir}/${ProjectName}_vs_silva_annotation_id_${identityThreshold}.biomtsv
        fi
        #check_file ${resultDir}/${ProjectName}_vs_silva_annotation_id_${identityThreshold}.tsv
        say "Elapsed time with get_taxonomy: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_count.biom" ] && [ -f "${resultDir}/${ProjectName}_vs_silva_annotation_id_${identityThreshold}.biomtsv" ] && [ ! -f "${resultDir}/${ProjectName}_silva_id_${identityThreshold}.biom" ]
    then
        say "Build vsearch-silva biom"
        start_time=$(timer)
        $biom add-metadata -i ${resultDir}/${ProjectName}_count.biom -o ${resultDir}/${ProjectName}_silva_id_${identityThreshold}.biom --observation-metadata-fp ${resultDir}/${ProjectName}_vs_silva_annotation_id_${identityThreshold}.biomtsv --observation-header id,taxonomy --sc-separated taxonomy --output-as-json
        check_file ${resultDir}/${ProjectName}_silva_id_${identityThreshold}.biom 
        say "Elapsed time to build vsearch-silva biom: $(timer $start_time)"
    fi
    if [ ! -f "${resultDir}/${ProjectName}_vs_silva_eval_${evalueTaxAnnot}.tsv" ] && [ "$blast_tax" -eq "1" ]  && [ "$fungi" -eq "0" ]
    then
        say "Assign taxonomy against silva with blast"
        start_time=$(timer)
        if [ "$lsu" -eq "1" ]
        then
            $blastn -query ${resultDir}/${ProjectName}_otu.fasta -db $silvalsu -evalue $evalueTaxAnnot -num_threads $NbProc -out ${resultDir}/${ProjectName}_vs_silva_eval_${evalueTaxAnnot}.tsv -max_target_seqs $maxTargetSeqs -task megablast -outfmt "6 qseqid sseqid  pident qcovs evalue" -use_index true 
        else
            $blastn -query ${resultDir}/${ProjectName}_otu.fasta -db $silva -evalue $evalueTaxAnnot -num_threads $NbProc -out ${resultDir}/${ProjectName}_vs_silva_eval_${evalueTaxAnnot}.tsv -max_target_seqs $maxTargetSeqs -task megablast -outfmt "6 qseqid sseqid  pident qcovs evalue" -use_index true 
        fi
        #check_file ${resultDir}/${ProjectName}_vs_silva_eval_${evalueTaxAnnot}.tsv
        say "Elapsed time with blast: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_vs_silva_eval_${evalueTaxAnnot}.tsv" ] && [ ! -f "${resultDir}/${ProjectName}_vs_silva_annotation_eval_${evalueTaxAnnot}.tsv" ]
    then
        say "Extract silva annotation with get_taxonomy"
        start_time=$(timer)
        if [ "$lsu" -eq "1" ]
        then
            python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_silva_eval_${evalueTaxAnnot}.tsv -d $silvalsu -u ${resultDir}/${ProjectName}_otu.fasta -o ${resultDir}/${ProjectName}_vs_silva_annotation_eval_${evalueTaxAnnot}.tsv -ob ${resultDir}/${ProjectName}_vs_silva_annotation_eval_${evalueTaxAnnot}.biomtsv
        else
            python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_silva_eval_${evalueTaxAnnot}.tsv -d $silva -u ${resultDir}/${ProjectName}_otu.fasta -o ${resultDir}/${ProjectName}_vs_silva_annotation_eval_${evalueTaxAnnot}.tsv -ob ${resultDir}/${ProjectName}_vs_silva_annotation_eval_${evalueTaxAnnot}.biomtsv
        fi
        #check_file ${resultDir}/${ProjectName}_vs_silva_annotation_eval_${evalueTaxAnnot}.tsv
        say "Elapsed time with get_taxonomy: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_count.biom" ] && [ -f "${resultDir}/${ProjectName}_vs_silva_annotation_eval_${evalueTaxAnnot}.biomtsv" ] && [ ! -f "${resultDir}/${ProjectName}_silva_eval_${evalueTaxAnnot}.biom" ]
    then
        say "Build blast-silva biom"
        start_time=$(timer)
        $biom add-metadata -i ${resultDir}/${ProjectName}_count.biom -o ${resultDir}/${ProjectName}_silva_eval_${evalueTaxAnnot}.biom --observation-metadata-fp ${resultDir}/${ProjectName}_vs_silva_annotation_eval_${evalueTaxAnnot}.biomtsv --observation-header id,taxonomy --sc-separated taxonomy --output-as-json
        check_file ${resultDir}/${ProjectName}_silva_eval_${evalueTaxAnnot}.biom 
        say "Elapsed time to build blast-silva biom: $(timer $start_time)"
    fi
    # Greengenes
    if [ ! -f "${resultDir}/${ProjectName}_vs_greengenes_id_${identityThreshold}.tsv" ] && [ "$blast_tax" -eq "0" ] && [ "$fungi" -eq "0" ] && [ "$lsu" -eq "0" ]
    then
        say "Assign taxonomy against greengenes with vsearch"
        start_time=$(timer)
        $vsearch --usearch_global ${resultDir}/${ProjectName}_otu.fasta --db $greengenes --id $identityThreshold --blast6out ${resultDir}/${ProjectName}_vs_greengenes_id_${identityThreshold}.tsv --strand both
        #check_file ${resultDir}/${ProjectName}_vs_greengenes_id_${identityThreshold}.tsv 
        say "Elapsed time with vsearch: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_vs_greengenes_id_${identityThreshold}.tsv" ] && [ ! -f "${resultDir}/${ProjectName}_vs_greengenes_annotation_id_${identityThreshold}.tsv" ]
    then
        say "Extract vsearch - greengenes annotation with get_taxonomy"
        start_time=$(timer)
        python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_greengenes_id_${identityThreshold}.tsv -d $greengenes -o ${resultDir}/${ProjectName}_vs_greengenes_annotation_id_${identityThreshold}.tsv -dtype greengenes -t $greengenes_taxonomy -ob ${resultDir}/${ProjectName}_vs_greengenes_annotation_id_${identityThreshold}.biomtsv -u ${resultDir}/${ProjectName}_otu.fasta
        #check_file ${resultDir}/${ProjectName}_vs_greengenes_annotation_id_${identityThreshold}.tsv
        say "Elapsed time with vsearch: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_count.biom" ] && [ -f "${resultDir}/${ProjectName}_vs_greengenes_annotation_id_${identityThreshold}.biomtsv" ] && [ ! -f "${resultDir}/${ProjectName}_greengenes_id_${identityThreshold}.biom" ]
    then
        say "Build vsearch-greengenes biom"
        start_time=$(timer)
        $biom add-metadata -i ${resultDir}/${ProjectName}_count.biom -o ${resultDir}/${ProjectName}_greengenes_id_${identityThreshold}.biom --observation-metadata-fp ${resultDir}/${ProjectName}_vs_greengenes_annotation_id_${identityThreshold}.biomtsv --observation-header id,taxonomy --sc-separated taxonomy --output-as-json
        check_file ${resultDir}/${ProjectName}_greengenes_id_${identityThreshold}.biom 
        say "Elapsed time to build vsearch-greengenes biom: $(timer $start_time)"
    fi
    if [ ! -f "${resultDir}/${ProjectName}_vs_greengenes_eval_${evalueTaxAnnot}.tsv" ] && [ "$blast_tax" -eq "1" ] && [ "$fungi" -eq "0" ] && [ "$lsu" -eq "0" ]
    then
        say "Assign taxonomy against greengenes with blast"
        start_time=$(timer)
        $blastn -query ${resultDir}/${ProjectName}_otu.fasta -db $greengenes -evalue $evalueTaxAnnot -num_threads $NbProc -out ${resultDir}/${ProjectName}_vs_greengenes_eval_${evalueTaxAnnot}.tsv -max_target_seqs $maxTargetSeqs -task megablast -outfmt "6 qseqid sseqid  pident qcovs evalue" -use_index true
        #check_file ${resultDir}/${ProjectName}_vs_greengenes_eval_${evalueTaxAnnot}.tsv
        say "Elapsed time with blast: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_vs_greengenes_eval_${evalueTaxAnnot}.tsv" ] && [ ! -f "${resultDir}/${ProjectName}_vs_greengenes_annotation_eval_${evalueTaxAnnot}.tsv" ]
    then
        say "Extract greengenes annotation with get_taxonomy"
        start_time=$(timer)
        python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_greengenes_eval_${evalueTaxAnnot}.tsv -d $greengenes -u ${resultDir}/${ProjectName}_otu.fasta -o ${resultDir}/${ProjectName}_vs_greengenes_annotation_eval_${evalueTaxAnnot}.tsv -dtype greengenes -t $greengenes_taxonomy -ob ${resultDir}/${ProjectName}_vs_greengenes_annotation_eval_${evalueTaxAnnot}.biomtsv
        #check_file ${resultDir}/${ProjectName}_vs_greengenes_annotation_eval_${evalueTaxAnnot}.tsv
        say "Elapsed time with get_taxonomy: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_count.biom" ] && [ -f "${resultDir}/${ProjectName}_vs_greengenes_annotation_eval_${evalueTaxAnnot}.biomtsv" ] && [ ! -f "${resultDir}/${ProjectName}_greengenes_eval_${evalueTaxAnnot}.biom" ]
    then
        say "Build blast-greengenes biom"
        start_time=$(timer)
        $biom add-metadata -i ${resultDir}/${ProjectName}_count.biom -o ${resultDir}/${ProjectName}_greengenes_eval_${evalueTaxAnnot}.biom --observation-metadata-fp ${resultDir}/${ProjectName}_vs_greengenes_annotation_eval_${evalueTaxAnnot}.biomtsv --observation-header id,taxonomy --sc-separated taxonomy --output-as-json
        check_file ${resultDir}/${ProjectName}_greengenes_eval_${evalueTaxAnnot}.biom 
        say "Elapsed time to build blast-greengenes biom: $(timer $start_time)"
    fi
    if [ ! -f "${resultDir}/${ProjectName}_vs_findley_id_${identityThreshold}.tsv" ] && [ "$blast_tax" -eq "0" ] && [ "$fungi" -eq "1" ]
    then
        say "Assign taxonomy against findley with vsearch"
        start_time=$(timer)
        $vsearch --usearch_global ${resultDir}/${ProjectName}_otu.fasta --db $findley --id $identityThreshold --blast6out ${resultDir}/${ProjectName}_vs_findley_id_${identityThreshold}.tsv --strand both
        #check_file ${resultDir}/${ProjectName}_vs_findley_id_${identityThreshold}.tsv
        say "Elapsed time with vsearch: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_vs_findley_id_${identityThreshold}.tsv" ] && [ ! -f "${resultDir}/${ProjectName}_vs_findley_annotation_id_${identityThreshold}.tsv" ]
    then
        say "Extract vsearch - findley annotation with get_taxonomy"
        start_time=$(timer)
        python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_findley_id_${identityThreshold}.tsv -d $findley  -u ${resultDir}/${ProjectName}_otu.fasta -o ${resultDir}/${ProjectName}_vs_findley_annotation_id_${identityThreshold}.tsv -ob ${resultDir}/${ProjectName}_vs_findley_annotation_id_${identityThreshold}.biomtsv  -dtype findley
        #check_file ${resultDir}/${ProjectName}_vs_findley_annotation_id_${identityThreshold}.tsv
        say "Elapsed time with vsearch: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_count.biom" ] && [ -f "${resultDir}/${ProjectName}_vs_findley_annotation_id_${identityThreshold}.biomtsv" ] && [ ! -f "${resultDir}/${ProjectName}_findley_id_${identityThreshold}.biom" ]
    then
        say "Build vsearch-findley biom"
        start_time=$(timer)
        $biom add-metadata -i ${resultDir}/${ProjectName}_count.biom -o ${resultDir}/${ProjectName}_findley_id_${identityThreshold}.biom --observation-metadata-fp ${resultDir}/${ProjectName}_vs_findley_annotation_id_${identityThreshold}.biomtsv --observation-header id,taxonomy --sc-separated taxonomy --output-as-json
        check_file ${resultDir}/${ProjectName}_findley_id_${identityThreshold}.biom 
        say "Elapsed time to build vsearch-findley biom: $(timer $start_time)"
    fi
    if [ ! -f "${resultDir}/${ProjectName}_vs_findley_eval_${evalueTaxAnnot}.tsv" ] && [ "$blast_tax" -eq "1" ] && [ "$fungi" -eq "1" ]
    then
        say "Assign taxonomy against findley with blast"
        start_time=$(timer)
        $blastn -query ${resultDir}/${ProjectName}_otu.fasta -db $findley -evalue $evalueTaxAnnot -num_threads $NbProc -out ${resultDir}/${ProjectName}_vs_findley_eval_${evalueTaxAnnot}.tsv -max_target_seqs $maxTargetSeqs -task megablast -outfmt "6 qseqid sseqid  pident qcovs evalue" -use_index true
        #check_file ${resultDir}/${ProjectName}_vs_findley_eval_${evalueTaxAnnot}.tsv
        say "Elapsed time with blast: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_vs_findley_eval_${evalueTaxAnnot}.tsv" ] && [ ! -f "${resultDir}/${ProjectName}_vs_findley_annotation_eval_${evalueTaxAnnot}.tsv" ]
    then
        say "Extract findley annotation with get_taxonomy"
        start_time=$(timer)
        python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_findley_eval_${evalueTaxAnnot}.tsv -d $findley -u  ${resultDir}/${ProjectName}_otu.fasta  -o ${resultDir}/${ProjectName}_vs_findley_annotation_eval_${evalueTaxAnnot}.tsv -ob ${resultDir}/${ProjectName}_vs_findley_annotation_eval_${evalueTaxAnnot}.biomtsv -dtype findley
        #check_file ${resultDir}/${ProjectName}_vs_findley_annotation_eval_${evalueTaxAnnot}.tsv
        say "Elapsed time with get_taxonomy: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_count.biom" ] && [ -f "${resultDir}/${ProjectName}_vs_findley_annotation_eval_${evalueTaxAnnot}.biomtsv" ] && [ ! -f "${resultDir}/${ProjectName}_findley_eval_${evalueTaxAnnot}.biom" ]
    then
        say "Build blast-findley biom"
        start_time=$(timer)
        $biom add-metadata -i ${resultDir}/${ProjectName}_count.biom -o ${resultDir}/${ProjectName}_findley_eval_${evalueTaxAnnot}.biom --observation-metadata-fp ${resultDir}/${ProjectName}_vs_findley_annotation_eval_${evalueTaxAnnot}.biomtsv --observation-header id,taxonomy --sc-separated taxonomy --output-as-json
        check_file ${resultDir}/${ProjectName}_findley_eval_${evalueTaxAnnot}.biom 
        say "Elapsed time to build blast-findley biom: $(timer $start_time)"
    fi
    # UNITE
    if [ ! -f "${resultDir}/${ProjectName}_vs_unite_id_${identityThreshold}.tsv" ] && [ "$blast_tax" -eq "0" ] && [ "$fungi" -eq "1" ]
    then
        say "Assign taxonomy against unite with vsearch"
        start_time=$(timer)
        $vsearch --usearch_global ${resultDir}/${ProjectName}_otu.fasta --db $unite --id $identityThreshold --blast6out ${resultDir}/${ProjectName}_vs_unite_id_${identityThreshold}.tsv --strand both
        #check_file ${resultDir}/${ProjectName}_vs_unite_id_${identityThreshold}.tsv
        say "Elapsed time with vsearch: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_vs_unite_id_${identityThreshold}.tsv" ] && [ ! -f "${resultDir}/${ProjectName}_vs_unite_annotation_id_${identityThreshold}.tsv" ]
    then
        say "Extract vsearch - unite annotation with get_taxonomy"
        start_time=$(timer)
        python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_unite_id_${identityThreshold}.tsv -d $unite -u  ${resultDir}/${ProjectName}_otu.fasta  -o ${resultDir}/${ProjectName}_vs_unite_annotation_id_${identityThreshold}.tsv -ob ${resultDir}/${ProjectName}_vs_unite_annotation_id_${identityThreshold}.biomtsv -dtype unite
        #check_file ${resultDir}/${ProjectName}_vs_unite_annotation_id_${identityThreshold}.tsv
        say "Elapsed time with vsearch: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_count.biom" ] && [ -f "${resultDir}/${ProjectName}_vs_unite_annotation_id_${identityThreshold}.biomtsv" ] && [ ! -f "${resultDir}/${ProjectName}_unite_id_${identityThreshold}.biom" ]
    then
        say "Build vsearch-unite biom"
        start_time=$(timer)
        $biom add-metadata -i ${resultDir}/${ProjectName}_count.biom -o ${resultDir}/${ProjectName}_unite_id_${identityThreshold}.biom --observation-metadata-fp ${resultDir}/${ProjectName}_vs_unite_annotation_id_${identityThreshold}.biomtsv --observation-header id,taxonomy --sc-separated taxonomy --output-as-json
        check_file ${resultDir}/${ProjectName}_unite_id_${identityThreshold}.biom 
        say "Elapsed time to build vsearch-unite biom: $(timer $start_time)"
    fi
    if [ ! -f "${resultDir}/${ProjectName}_vs_unite_eval_${evalueTaxAnnot}.tsv" ] && [ "$blast_tax" -eq "1" ] && [ "$fungi" -eq "1" ]
    then
         say "Assign taxonomy against unite with blast"
         start_time=$(timer)
         $blastn -query ${resultDir}/${ProjectName}_otu.fasta -db $unite -evalue $evalueTaxAnnot -num_threads $NbProc -out ${resultDir}/${ProjectName}_vs_unite_eval_${evalueTaxAnnot}.tsv -max_target_seqs $maxTargetSeqs -task megablast -outfmt "6 qseqid sseqid  pident qcovs evalue" -use_index true
         #check_file ${resultDir}/${ProjectName}_vs_unite_eval_${evalueTaxAnnot}.tsv
         say "Elapsed time with blast: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_vs_unite_eval_${evalueTaxAnnot}.tsv" ] && [ ! -f "${resultDir}/${ProjectName}_vs_unite_annotation_eval_${evalueTaxAnnot}.tsv" ]
    then
         say "Extract unite annotation with get_taxonomy"
         start_time=$(timer)
         python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_unite_eval_${evalueTaxAnnot}.tsv -d $unite -u  ${resultDir}/${ProjectName}_otu.fasta  -o ${resultDir}/${ProjectName}_vs_unite_annotation_eval_${evalueTaxAnnot}.tsv -ob ${resultDir}/${ProjectName}_vs_unite_annotation_eval_${evalueTaxAnnot}.biomtsv -dtype unite
         #check_file ${resultDir}/${ProjectName}_vs_unite_annotation_eval_${evalueTaxAnnot}.tsv
         say "Elapsed time with get_taxonomy: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_count.biom" ] && [ -f "${resultDir}/${ProjectName}_vs_unite_annotation_eval_${evalueTaxAnnot}.biomtsv" ] && [ ! -f "${resultDir}/${ProjectName}_unite_eval_${evalueTaxAnnot}.biom" ]
    then
         say "Build blast-unite biom"
         start_time=$(timer)
         $biom add-metadata -i ${resultDir}/${ProjectName}_count.biom -o ${resultDir}/${ProjectName}_unite_eval_${evalueTaxAnnot}.biom --observation-metadata-fp ${resultDir}/${ProjectName}_vs_unite_annotation_eval_${evalueTaxAnnot}.biomtsv --observation-header id,taxonomy --sc-separated taxonomy --output-as-json
         check_file ${resultDir}/${ProjectName}_unite_eval_${evalueTaxAnnot}.biom 
         say "Elapsed time to build blast-unite biom: $(timer $start_time)"
    fi
    #Underhill
    if [ ! -f "${resultDir}/${ProjectName}_vs_underhill_id_${identityThreshold}.tsv" ] && [ "$blast_tax" -eq "0" ] && [ "$fungi" -eq "1" ]
    then
        say "Assign taxonomy against underhill with vsearch"
        start_time=$(timer)
        $vsearch --usearch_global ${resultDir}/${ProjectName}_otu.fasta --db $underhill --id $identityThreshold --blast6out ${resultDir}/${ProjectName}_vs_underhill_id_${identityThreshold}.tsv --strand both
        #check_file ${resultDir}/${ProjectName}_vs_underhill_id_${identityThreshold}.tsv
        say "Elapsed time with vsearch: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_vs_underhill_id_${identityThreshold}.tsv" ] && [ ! -f "${resultDir}/${ProjectName}_vs_underhill_annotation_id_${identityThreshold}.tsv" ]
    then
        say "Extract vsearch - underhill annotation with get_taxonomy"
        start_time=$(timer)
        python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_underhill_id_${identityThreshold}.tsv -d $underhill -u  ${resultDir}/${ProjectName}_otu.fasta -t $underhill_taxonomy   -o ${resultDir}/${ProjectName}_vs_underhill_annotation_id_${identityThreshold}.tsv -ob ${resultDir}/${ProjectName}_vs_underhill_annotation_id_${identityThreshold}.biomtsv -dtype underhill
        #check_file ${resultDir}/${ProjectName}_vs_underhill_annotation_id_${identityThreshold}.tsv
        say "Elapsed time with vsearch: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_count.biom" ] && [ -f "${resultDir}/${ProjectName}_vs_underhill_annotation_id_${identityThreshold}.biomtsv" ] && [ ! -f "${resultDir}/${ProjectName}_underhill_id_${identityThreshold}.biom" ]
    then
        say "Build vsearch-underhill biom"
        start_time=$(timer)
        $biom add-metadata -i ${resultDir}/${ProjectName}_count.biom -o ${resultDir}/${ProjectName}_underhill_id_${identityThreshold}.biom --observation-metadata-fp ${resultDir}/${ProjectName}_vs_underhill_annotation_id_${identityThreshold}.biomtsv --observation-header id,taxonomy --sc-separated taxonomy --output-as-json
        check_file ${resultDir}/${ProjectName}_underhill_id_${identityThreshold}.biom 
        say "Elapsed time to build vsearch-underhill biom: $(timer $start_time)"
    fi
    if [ ! -f "${resultDir}/${ProjectName}_vs_underhill_eval_${evalueTaxAnnot}.tsv" ] && [ "$blast_tax" -eq "1" ] && [ "$fungi" -eq "1" ]
    then
         say "Assign taxonomy against underhill with blast"
         start_time=$(timer)
         $blastn -query ${resultDir}/${ProjectName}_otu.fasta -db $underhill -evalue $evalueTaxAnnot -num_threads $NbProc -out ${resultDir}/${ProjectName}_vs_underhill_eval_${evalueTaxAnnot}.tsv -max_target_seqs $maxTargetSeqs -task megablast -outfmt "6 qseqid sseqid  pident qcovs evalue" -use_index true
         #check_file ${resultDir}/${ProjectName}_vs_underhill_eval_${evalueTaxAnnot}.tsv
         say "Elapsed time with blast: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_vs_underhill_eval_${evalueTaxAnnot}.tsv" ] && [ ! -f "${resultDir}/${ProjectName}_vs_underhill_annotation_eval_${evalueTaxAnnot}.tsv" ]
    then
         say "Extract underhill annotation with get_taxonomy"
         start_time=$(timer)
         python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_underhill_eval_${evalueTaxAnnot}.tsv -d $underhill -u  ${resultDir}/${ProjectName}_otu.fasta -t $underhill_taxonomy  -o ${resultDir}/${ProjectName}_vs_underhill_annotation_eval_${evalueTaxAnnot}.tsv -ob ${resultDir}/${ProjectName}_vs_underhill_annotation_eval_${evalueTaxAnnot}.biomtsv -dtype underhill
         #check_file ${resultDir}/${ProjectName}_vs_underhill_annotation_eval_${evalueTaxAnnot}.tsv
         say "Elapsed time with get_taxonomy: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_count.biom" ] && [ -f "${resultDir}/${ProjectName}_vs_underhill_annotation_eval_${evalueTaxAnnot}.biomtsv" ] && [ ! -f "${resultDir}/${ProjectName}_underhill_eval_${evalueTaxAnnot}.biom" ]
    then
         say "Build blast-underhill biom"
         start_time=$(timer)
         $biom add-metadata -i ${resultDir}/${ProjectName}_count.biom -o ${resultDir}/${ProjectName}_underhill_eval_${evalueTaxAnnot}.biom --observation-metadata-fp ${resultDir}/${ProjectName}_vs_underhill_annotation_eval_${evalueTaxAnnot}.biomtsv --observation-header id,taxonomy --sc-separated taxonomy --output-as-json
         check_file ${resultDir}/${ProjectName}_underhill_eval_${evalueTaxAnnot}.biom 
         say "Elapsed time to build blast-underhill biom: $(timer $start_time)"
    fi
    ##
    # Phylogenetic analysis
    ##
    for annotation in $(ls  ${resultDir}/${ProjectName}_*_annotation_*.tsv ${resultDir}/${ProjectName}_vs_rdp.tsv)
    do
        # Ugly
        soft=$(echo $annotation |sed -r "s:.*vs_(.+)_annotation_.*:\1:g" )
        if [ "$annotation" == "${resultDir}/${ProjectName}_vs_rdp.tsv" ]
        then
            soft="rdp"
        fi
        if [ "$soft" == "" ]
        then
            error "Unable to recognize the software for: $annotation"
            exit 1
        fi
        # Select annotated OTU
        if [ -f "${resultDir}/${ProjectName}_otu.fasta" ] && [ ! -f "${resultDir}/${ProjectName}_otu_${soft}.fasta" ]
        then
            say "Extract OTU annotated with $soft"
            start_time=$(timer)
            python $extract_fasta -d ${resultDir}/${ProjectName}_otu.fasta -i $annotation -o ${resultDir}/${ProjectName}_otu_${soft}.fasta
            check_file ${resultDir}/${ProjectName}_otu_${soft}.fasta
            say "Elapsed time with extract_fasta: $(timer $start_time)"
        fi
        # Alignment
        if [ -f "${resultDir}/${ProjectName}_otu_${soft}.fasta" ] && [ ! -f "${resultDir}/${ProjectName}_otu_${soft}.ali" ]
        then
            say "Align OTU annotated with $soft with mafft"
            start_time=$(timer)
            $mafft --adjustdirectionaccurately --thread $NbProc --genafpair --maxiterate 1000 --ep 0  ${resultDir}/${ProjectName}_otu_${soft}.fasta > ${resultDir}/${ProjectName}_otu_${soft}.ali 2> ${logDir}/log_mafft_${ProjectName}_${soft}.txt
            check_file ${resultDir}/${ProjectName}_otu_${soft}.ali
            sed "s:_R_::g" ${resultDir}/${ProjectName}_otu_${soft}.ali -i
            say "Elapsed time with mafft: $(timer $start_time)"
        fi

        # BMGE
        if [ -f "${resultDir}/${ProjectName}_otu_${soft}.ali" ] && [ ! -f "${resultDir}/${ProjectName}_otu_${soft}_bmge.ali" ]
        then
            say "Cut the alignment with BMGE for $soft annotation"
            start_time=$(timer)
            $BMGE -i ${resultDir}/${ProjectName}_otu_${soft}.ali -t DNA -m ID -h 1 -g $conservedPosition -w 1 -b 1 -of ${resultDir}/${ProjectName}_otu_${soft}_bmge.ali
            check_file ${resultDir}/${ProjectName}_otu_${soft}_bmge.ali
            say "Elapsed time with BMGE: $(timer $start_time)"
        fi

        # Phylogeny
        #if [ -f "${resultDir}/${ProjectName}_otu_${soft}_bmge.ali" ] && [ ! -f "${resultDir}/${ProjectName}_otu_${soft}_bmge.tree" ]
        if [ -f "${resultDir}/${ProjectName}_otu_${soft}_bmge.ali" ] && [ ! -f "${resultDir}/${ProjectName}_otu_${soft}_bmge.ali.treefile" ]
        then
            start_time=$(timer)
            if [ "$accurateTree" == "1" ]
            then  
                say "Compute tree with IQ-TREE for $soft annotation"
                $iqtree -m GTR+I+G4  -nt $NbProc -s ${resultDir}/${ProjectName}_otu_${soft}_bmge.ali > ${logDir}/log_iqtree_${soft}.txt 
            else
                say "Compute tree with FastTree $soft annotation"
                $FastTreeMP -nt ${resultDir}/${ProjectName}_otu_${soft}_bmge.ali > ${resultDir}/${ProjectName}_otu_${soft}_bmge.ali.treefile 2> ${logDir}/log_fasttree_${soft}.txt
            fi
            check_file ${resultDir}/${ProjectName}_otu_${soft}_bmge.ali.treefile
            #$FastTreeMP -nt ${resultDir}/${ProjectName}_otu_${soft}_bmge.ali > ${resultDir}/${ProjectName}_otu_${soft}_bmge.tree
            #check_file ${resultDir}/${ProjectName}_otu_${soft}_bmge.tree
            say "Elapsed time with tree building: $(timer $start_time)"
        fi
    done
fi

# Extract result
if [ -d "$input_dir" ] && [ ! -f "${resultDir}/${ProjectName}_build_process.tsv" ]
then
    say "Extract MASQUE process results with extract_result"
    start_time=$(timer)
    if [ "$paired" -eq "1" ]
    then
        python $extract_result -d ${resultDir}/ -r $input_dir/ -p -o1 ${resultDir}/${ProjectName}_build_process.tsv -o2 ${resultDir}/${ProjectName}_annotation_process.tsv
    elif [ "$paired" -eq "0" ]
    then
        python $extract_result -d ${resultDir}/ -r $input_dir/ -o1 ${resultDir}/${ProjectName}_build_process.tsv -o2 ${resultDir}/${ProjectName}_annotation_process.tsv
#     elif [ ! -d "$input_dir" ] && [ -f "$amplicon" ]
#     then
#         $extract_result -a $input_dir/ -p -o1 ${resultDir}/${ProjectName}_build_process.tsv -o2 ${resultDir}/${ProjectName}_annotation_process.tsv
    fi
    check_file ${resultDir}/${ProjectName}_build_process.tsv
    say "Elapsed time with extract_result: $(timer $start_time)"
fi

say "MASQUE analysis is done. Elapsed time: $(timer $wall_time)"
