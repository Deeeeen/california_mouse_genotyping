#!/bin/bash

dir_path=/peromyscus
#### !!! where you keep .fastq sequence data, .fna reference genome, .csv metadata
#### also all the output goes to this directory
#### !!! change the path to all softwares to where you keep them

################################### Demux #####################################
ncpu=23
#### !!! This should be set according to num of node and ppn / CUP you requested.
#### Make sure to leave one cpu unused in case one process uses up all resources
#### and other process can not be processed.

mkdir ${dir_path}/demux
#### Make a directory to keep all demux results

START=$(date +%s)

java -Xmx40G -XX:+AggressiveOpts -XX:+AggressiveHeap \
     -jar fgbio-1.1.0.jar DemuxFastqs \
     --inputs ${dir_path}/Riptide09_S5_L003_R1_001.fastq.gz \
             ${dir_path}/Riptide09_S5_L003_R2_001.fastq.gz \
     --metadata ${dir_path}/2020-03-16-Flowcell_Sample-Barcode_list_Riptide09_Peromyscus.csv \
     --read-structures 8B12M+T 8M+T \
     --output-type=Fastq \
     --threads $ncpu \
     --output ${dir_path}/demux/ \
     --metrics ${dir_path}/demux/peromyscus_riptide_demux_barcode_metrics.txt
#### -Xmx40G memory request for JVM

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Demux time elapsed: $(( $END - $START )) seconds"
#### START and END keep track of how much time were used for thie process.
#### The while loop makes sure all steps in this process are done before
#### entering next process.

########################### Reference preparation #############################
START=$(date +%s)
#### Generate the bwa index 
bwa index ${dir_path}/GCA_007827085.2_ASM782708v2_genomic.fa

#### Generate the fasta file index 
samtools faidx ${dir_path}/GCA_007827085.2_ASM782708v2_genomic.fa

#### Generate the sequence dictionary 
java -Xmx40G -XX:+AggressiveOpts -XX:+AggressiveHeap \
     -jar picard.jar CreateSequenceDictionary \
        REFERENCE=${dir_path}/GCA_007827085.2_ASM782708v2_genomic.fa \
        OUTPUT=${dir_path}/GCA_007827085.2_ASM782708v2_genomic.dict

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Reference preparation time elapsed: $(( $END - $START )) seconds"

###################### Check quality of sequence reads ########################
#### FastQC reports
ncpu=23
mkdir ${dir_path}/qc
mkdir ${dir_path}/qc/fastqc_demux

START=$(date +%s)

fs=$(ls ${dir_path}/demux/*.fastq.gz)
cnt=0
for f in $fs
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   echo -e "\n-----run FastQC on ${cnt}-th file: $f-----"
   fastqc $f --outdir=${dir_path}/qc/fastqc_demux/ &
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "FastQC time elapsed: $(( $END - $START )) seconds"

#### MultiQC reports
mkdir ${dir_path}/qc/multiqc_demux

START=$(date +%s)

multiqc ${dir_path}/qc/fastqc_demux/*.zip -o ${dir_path}/qc/multiqc_demux

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "MultiQC time elapsed: $(( $END - $START )) seconds"

############################# Cutadapt trimming ###############################
ncpu=23

START=$(date +%s)

fs_out=()
fs_in=$(ls ${dir_path}/demux/*.fastq.gz)
mkdir ${dir_path}/trimmed
cnt=0

for f in $fs_in
do
   (( cnt += 1 ))
   temp=$(echo ${f} | rev | cut -d'/' -f 1 | rev)
   fs_out[cnt]=$(cut -d '_' -f 1 <<< $temp)
done

fs_out=$(for f in ${fs_out[@]}; do echo $f; done | sort -u)
cnt=0
for f in $fs_out
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   echo -e "\n-----run Cutadapt on ${cnt}-th file: ${dir_path}/demux/${f} > ${dir_path}/trimmed/${f}-trimmed-----"
   cutadapt -q 20 \
         --minimum-length 50 -pair-filter=both \
         -o ${dir_path}/trimmed/${f}-trimmed_R1.fastq.gz \
         -p ${dir_path}/trimmed/${f}-trimmed_R2.fastq.gz \
         ${dir_path}/demux/${f}_R1.fastq.gz ${dir_path}/demux/${f}_R2.fastq.gz &
done
#### The filters for cutadapt is ignore base quality less than 20
####                             keep the minimum read leagth as 50
while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Cutadapt time elapsed: $(( $END - $START )) seconds"

################## Check quality of trimmed sequence reads ####################
#### FastQC reports
ncpu=23
mkdir ${dir_path}/qc/fastqc_trimmed

START=$(date +%s)

fs=$(ls ${dir_path}/trimmed/*.fastq.gz)
cnt=0
for f in $fs
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   echo -e "\n-----run ${cnt}-th file: $f-----"
   fastqc $f --outdir=${dir_path}/qc/fastqc_trimmed/ &
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "FastQC time elapsed: $(( $END - $START )) seconds"

#### MultiQC reports
mkdir ${dir_path}/qc/multiqc_trimmed

START=$(date +%s)

multiqc ${dir_path}/qc/fastqc_trimmed/*.zip -o ${dir_path}/qc/multiqc_trimmed

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "MultiQC time elapsed: $(( $END - $START )) seconds"