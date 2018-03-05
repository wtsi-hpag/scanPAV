#!/bin/bash 

folder=`pwd`
scanPAV_version=$folder/../src/

cpus=30
aligner=bwa
score=550  
output=$aligner\_$score


# Assembly locations:
assembly1=$folder/yeast_data/S288c_with_indels_and_snps.fasta
assembly2=$folder/yeast_data/S288C_reference_sequence_R64-1-1_20110203.fsa
outname=$(basename $assembly1 .fasta)


rm -rf $output
mkdir -p $output 
cd $output

ln -sf $assembly1 .
ln -sf $assembly2 .

echo; echo Looking for Presence PAVs | tee log.txt ; echo
$scanPAV_version/scanPAV -nodes $cpus -score $score -align $aligner $(basename $assembly1) $(basename $assembly2)  pavs_present_in_$outname | tee -a  log.txt
echo; echo Looking for Absence PAVs | tee -a log.txt ; echo
$scanPAV_version/scanPAV -nodes $cpus -score $score -align $aligner $(basename $assembly2) $(basename $assembly1)  pavs_absent_in_$outname  | tee -a  log.txt

echo; echo; echo "Pipeline output can be found in " bwa_550/log.txt
echo;echo Presence/Absence PAVs are in the $output folder:
echo " " $output/pavs_present_in_$outname.fasta : pavs present in $(basename $assembly1) but absent in $(basename $assembly2)
echo "   Details of presence PAVs are listed in the vcf file:" $output/pavs_present_in_$outname 
echo;echo " " $output/pavs_absent_in_$outname.fasta : pavs absent in $(basename $assembly1) but present in $(basename $assembly2)
echo "   Details of absence PAVs are listed in the vcf file:" $output/pavs_absent_in_$outname 
echo
