#!/bin/bash

help () {
        echo "error_correction.sh [OPTIONS]"
        echo "          -h              print the help page"
        echo "          -t INT  Number of threads (default: 1)"
        echo "          -p PATH path of 'path.conf' file"
        echo "          -o PATH path of outDir"
        echo "          -n INT  Number of iteration (default: 1)"
        echo "          -i PATH path of 'input.txt' file"
        echo "          -s PATH path of 'Bin' "
        exit;
}

while getopts "t:p:o:n:i:s:h" opt
do
        case $opt in
                t) threads=$OPTARG
                        ;;
                p) path=$OPTARG
                        ;;
                o) outDir=$OPTARG
                        ;;
                n) iterNum=$OPTARG
                        ;;
                i) input=$OPTARG
                        ;;
                s) bin=$OPTARG
                        ;;
                h) help ;;
                ?) help ;;
        esac
done

shift $(( $OPTIND - 1))
file=$1


##############################################################################
### TOOL PATH CONFIGURATION ###
##############################################################################

source $path

fasplit=$kent/faSplit

##############################################################################
### INPUT FILE ###
##############################################################################

parameter=$input

INPUT_FASTA=''
P1=()
P2=()

while IFS='' read -a line || [[ -n "$line" ]]
do
        read -r name value <<< "$line"

        if [ "$name" == "fasta" ]; then
                INPUT_FASTA="$value"
        fi

        if [ "$name" == "p1" ]; then
                P1+=( "$value" )
        fi

        if [ "$name" == "p2" ]; then
                P2+=( "$value" )
        fi
done < $parameter

##############################################################################
### OUTPUT PATH ###
##############################################################################


outDir=$outDir
OUTPUT_BWA=$outDir/bwa/
OUTPUT_MERGEBAM=$outDir/mergeBAM/
OUTPUT_SORTBAM=$outDir/sortCOORD/
OUTPUT_MARKDUP=$outDir/markDUP/
OUTPUT_VARIANT=$outDir/variants/
OUTPUT_PILON=$outDir/pilon/
OUTPUT_FASPLIT=$outDir/pilon_byName/
OUTPUT_FINAL=$outDir/finalEdit/

##############################################################################
### MAKE OUTPUT FOLDERS ###
##############################################################################

# PREPROCESSING INPUT DATA
mkdir -p $outDir
mkdir -p $OUTPUT_BWA
mkdir -p $OUTPUT_MERGEBAM
mkdir -p $OUTPUT_SORTBAM
mkdir -p $OUTPUT_MARKDUP
mkdir -p $OUTPUT_VARIANT
mkdir -p $OUTPUT_PILON
mkdir -p $OUTPUT_FASPLIT
mkdir -p $OUTPUT_FINAL

##############################################################################
### PREPROCESSING THE FASTA FILE ###
##############################################################################

echo "    ($iterNum-1) Indexing the meta-assembly for read mapping"

$samtools faidx $INPUT_FASTA 2> $OUTPUT_BWA/log
$bwa index -a bwtsw $INPUT_FASTA 2> $OUTPUT_BWA/log
java -jar $picard CreateSequenceDictionary R=$INPUT_FASTA O=$(dirname $INPUT_FASTA)/$(basename $INPUT_FASTA .fa).dict 2> $OUTPUT_BWA/log

##############################################################################
### BWA MAPPING ###
##############################################################################

echo "    ($iterNum-2) Mapping the reads to meta-assembly"

for idx in ${!P1[@]};do
        fq1=${P1[idx]}
        fq2=${P2[idx]}
        res=$(basename $INPUT_FASTA .fa).$idx
        MERGEINPUT=$MERGEINPUT' INPUT='$OUTPUT_BWA/$res.bam
        $bwa mem $INPUT_FASTA -t $threads -R "@RG\tID:FLOWCELL1.LANE1\tPL:ILLUMINA\tLB:IMAP\tSM:IMAP" $fq1 $fq2 1> $OUTPUT_BWA/$res.sam 2> $OUTPUT_BWA/log
        $samtools view -h -@ $threads -bS $OUTPUT_BWA/$res.sam 1> $OUTPUT_BWA/$res.bam 2>> $OUTPUT_BWA/log
        rm $OUTPUT_BWA/$res.sam
done


##############################################################################
### BAM FILE MERGING ###
##############################################################################

echo "    ($iterNum-3) Merging the bam file(s)"

if [ ${#P1[@]} == 1 ];
then
        in=$res
        res=$(basename $INPUT_FASTA .fa).bamMerge.bam
        mv $OUTPUT_BWA/$in.bam $OUTPUT_MERGEBAM/$res

else
        res=$(basename $INPUT_FASTA .fa).bamMerge.bam
        java -jar $picard MergeSamFiles $MERGEINPUT O=$OUTPUT_MERGEBAM/$res 2> $OUTPUT_MERGEBAM/log
fi


##############################################################################
### SORTING BAM FILES AS COORDINATE ###
##############################################################################

echo "    ($iterNum-4) Sorting the bam file(s) as coordination"

in=$res
res=$(basename $in .bam).sorted.bam
java -jar $picard SortSam I=$OUTPUT_MERGEBAM/$in O=$OUTPUT_SORTBAM/$res SORT_ORDER=coordinate 2> $OUTPUT_SORTBAM/log

##############################################################################
### MARKING THE DUPLICATE ###
##############################################################################

echo "    ($iterNum-5) Marking the duplicate among bam file(s)"

in=$res
res=$(basename $in .bam).dedup.bam
java -jar $picard MarkDuplicates I=$OUTPUT_SORTBAM/$in O=$OUTPUT_MARKDUP/$res METRICS_FILE=$OUTPUT_MARKDUP/$(basename $res .bam).metrics.txt 2> $OUTPUT_MARKDUP/log
java -jar $picard BuildBamIndex I=$OUTPUT_MARKDUP/$res 2> $OUTPUT_MARKDUP/log

##############################################################################
### PILON  ###
##############################################################################

echo "    ($iterNum-6) Error correction"

PILON_BAM_RESULT=$OUTPUT_PILON/bam_result/
PILON_FASTA=$OUTPUT_PILON/fasta/
PILON_PILON=$OUTPUT_PILON/pilon/

mkdir -p $PILON_BAM_RESULT
mkdir -p $PILON_FASTA
mkdir -p $PILON_PILON

$fasplit byname $INPUT_FASTA $PILON_FASTA
cd $PILON_FASTA
ls *.fa -lv | awk '{ split($9,arr,"."); print arr[1] }' > $OUTPUT_PILON/seq_list

perl $bin/batch_split_bam_per_contigs.pl $OUTPUT_PILON/seq_list $OUTPUT_MARKDUP/$res $PILON_BAM_RESULT $threads 2> $OUTPUT_PILON/log
perl $bin/batch_pilon_per_scaffolds.pl $PILON_BAM_RESULT $PILON_FASTA $PILON_PILON $threads 2> $OUTPUT_PILON/log

cd $PILON_PILON
ls *.fasta -lv | awk 'index($9,"sscf") !=0 { print $9 }'  > $OUTPUT_PILON/scf_name
ls *.fasta -lv | awk 'index($9,"sscf") ==0 { print $9 }'  >> $OUTPUT_PILON/scf_name


while read LINE; do cat $PILON_PILON/$LINE; done < $OUTPUT_PILON/scf_name > $OUTPUT_PILON/$(basename $INPUT_FASTA .fa).pilon.fasta 2> $OUTPUT_PILON/log

in=$(basename $INPUT_FASTA .fa).pilon.fasta
res=$(basename $in .fasta).n50.fa
java -jar $picard NormalizeFasta INPUT=$OUTPUT_PILON/$in OUTPUT=$OUTPUT_PILON/$res LINE_LENGTH=50 2> $OUTPUT_PILON/log

##############################################################################
### START VARIANT CALLING INPUT IS PILON RESULT FILE  ###
##############################################################################

echo "    ($iterNum-7) Indexing the corrected meta-assembly for read mapping"

INPUT_FASTA=$OUTPUT_PILON/$res

$samtools faidx $INPUT_FASTA 2> $OUTPUT_BWA/log
$bwa index -a bwtsw $INPUT_FASTA 2> $OUTPUT_BWA/log
java -jar $picard CreateSequenceDictionary R=$INPUT_FASTA O=$(dirname $INPUT_FASTA)/$(basename $INPUT_FASTA .fa).dict 2> $OUTPUT_BWA/log

# $INPUT_FASTA is sigma.spades.edit.pilon.n50.fa

##############################################################################
### BWA ###
##############################################################################

echo "    ($iterNum-8) Mapping the reads to final assembly"

rm $OUTPUT_BWA/*
MERGEINPUT=''
for idx in ${!P1[@]};do
        fq1=${P1[idx]}
        fq2=${P2[idx]}
        res=$(basename $INPUT_FASTA .fa).$idx
        MERGEINPUT=$MERGEINPUT' INPUT='$OUTPUT_BWA/$res.bam
        $bwa mem $INPUT_FASTA -t $threads -R "@RG\tID:FLOWCELL1.LANE1\tPL:ILLUMINA\tLB:IMAP\tSM:IMAP" $fq1 $fq2 1> $OUTPUT_BWA/$res.sam 2> $OUTPUT_BWA/log
        $samtools view -h -@ $threads -bS $OUTPUT_BWA/$res.sam 1> $OUTPUT_BWA/$res.bam 2>> $OUTPUT_BWA/loeds
        rm $OUTPUT_BWA/$res.sam
done

##############################################################################
### BAM FILE MERGING ###
##############################################################################

echo "    ($iterNum-9) Merging the bam file(s)"

rm $OUTPUT_MERGEBAM/*

if [ ${#P1[@]} == 1 ];
then
        in=$res
        res=$(basename $INPUT_FASTA .fa).bamMerge.bam
        mv $OUTPUT_BWA/$in.bam $OUTPUT_MERGEBAM/$res

else
        res=$(basename $INPUT_FASTA .fa).bamMerge.bam
        java -jar $picard MergeSamFiles $MERGEINPUT O=$OUTPUT_MERGEBAM/$res 2> $OUTPUT_MERGEBAM/log
fi


##############################################################################
### SORTING BAM FILES AS COORDINATE ###
##############################################################################

echo "    ($iterNum-10) Sorting the bam file(s) as corrdination"

rm $OUTPUT_SORTBAM/*
in=$res
res=$(basename $in .bam).sorted.bam
java -jar $picard SortSam I=$OUTPUT_MERGEBAM/$in O=$OUTPUT_SORTBAM/$res SORT_ORDER=coordinate 2> $OUTPUT_SORTBAM/log

##############################################################################
### MARKING THE DUPLICATE ###
##############################################################################

echo "    ($iterNum-11) marking the duplication of final assembly"

rm $OUTPUT_MARKDUP/*
in=$res
res=$(basename $in .bam).dedup.bam
java -jar $picard MarkDuplicates I=$OUTPUT_SORTBAM/$in O=$OUTPUT_MARKDUP/$res METRICS_FILE=$OUTPUT_MARKDUP/$(basename $res .bam).metrics.txt 2> $OUTPUT_MARKDUP/log
java -jar $picard BuildBamIndex I=$OUTPUT_MARKDUP/$res 2> $OUTPUT_MARKDUP/log

##############################################################################
### REALIGNING THE TARGET CREATOR ###
##############################################################################

echo "    ($iterNum-12) Realigning the target creator"

in=$res
resList=$(basename $in .dedup.bam).targetintervals.list
java -jar $gatk -T RealignerTargetCreator -R $INPUT_FASTA -I $OUTPUT_MARKDUP/$in -o $OUTPUT_VARIANT/$resList 2> $OUTPUT_VARIANT/log

##############################################################################
### INDEL REALIGNER ###
##############################################################################

echo "    ($iterNum-13) Realigning indel of final assembly"

res=$(basename $in .dedup.bam).realigned.bam
java -jar $gatk -T IndelRealigner -R $INPUT_FASTA -I $OUTPUT_MARKDUP/$in -targetIntervals $OUTPUT_VARIANT/$resList -o $OUTPUT_VARIANT/$res 2> $OUTPUT_VARIANT/log

##############################################################################
### HAPLOTYPE CALLING ###
##############################################################################

echo "    ($iterNum-14) Haplotype calling of final assembly"

in=$res
res=$(basename $in .realigned.bam).vcf
java -jar $gatk -T HaplotypeCaller -R $INPUT_FASTA -I $OUTPUT_VARIANT/$in -o $OUTPUT_VARIANT/$res 2> $OUTPUT_VARIANT/log

##############################################################################
### FASTA FILE SPLIT BY NAME ###
##############################################################################

$fasplit byname $INPUT_FASTA $OUTPUT_FASPLIT/

##############################################################################
### SNP EDITING BY VCF FILE ###
##############################################################################

echo "    ($iterNum-15) Editing the SNP of final assembly"

vcf=$res
res=$(basename $INPUT_FASTA .fa)
$bin/snp_edit $OUTPUT_VARIANT/$vcf $OUTPUT_FASPLIT $OUTPUT_FINAL/$res 2> $OUTPUT_FINAL/log

##############################################################################
### NORMALIZE THE FASTA FILE FINAL OUTPUT ###
##############################################################################

echo "    ($iterNum-16) Normalizing the final assembly"

in=$res.edit.fa
res=$(basename $in .fa).norm.fa
java -jar $picard NormalizeFasta INPUT=$OUTPUT_FINAL/$in OUTPUT=$OUTPUT_FINAL/$res LINE_LENGTH=50 2> $OUTPUT_FINAL/log
cp $OUTPUT_FINAL/$res $outDir/meta_assembled_scaf.fa

rm -rf $PILON_BAM_RESULT
rm -rf $PILON_FASTA
rm -rf $PILON_PILON

exit 0
