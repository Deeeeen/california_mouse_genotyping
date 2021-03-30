#!/bin/bash

dir_path=/peromyscus
#### !!! where you keep .fastq sequence data, .fna reference genome, .csv metadata
#### also all the output goes to this directory

################# Map sequence reads against reference genome #################
#### Map to reference
ncpu=12
#### !!! This should be set according to num of node and ppn / CUP you requested.

mkdir ${dir_path}/sams

START=$(date +%s)

fs_out=()
fs_in=$(ls ${dir_path}/trimmed/*.fastq.gz)
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
   rg=($(cut -d '-' -f 1 <<< $f) $(cut -d '-' -f 2 <<< $f) $(cut -d '-' -f 3 <<< $f))
   echo -e "\n-----run ${cnt}-th file: ${dir_path}/trimmed/${f} > ${dir_path}/sams/${f}.sam-----"

#### !!! You may want to chage the RG read group based on your fastq files.
   bwa mem -aM -t 2\
     -R "@RG\tID:A00953.86.HC2FMDSXY.3\tLB:Riptide09\tPL:ILLUMINA\tSM:${rg[1]}\tPU:HC2FMDSXY.3.${rg[2]}" \
     ${dir_path}/GCA_007827085.2_ASM782708v2_genomic.fa ${dir_path}/trimmed/${f}_R1.fastq.gz \
     ${dir_path}/trimmed/${f}_R2.fastq.gz > ${dir_path}/sams/${f}.sam &
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

END=$(date +%s)
echo "BWA map to reference Time elapsed: $(( $END - $START )) seconds"

#### Convert a SAM file to a BAM file and sort BAM file by coordinates
ncpu=23
mkdir ${dir_path}/bams
START=$(date +%s)

fs_in=$(ls ${dir_path}/sams/*.sam)
fs_out=()
cnt=0
for f in $fs_in
do
   (( cnt += 1 ))
   temp=$(echo ${f} | rev | cut -d '/' -f 1 | rev)
   fs_out[cnt]=$(echo ${temp} | rev | cut -d '.' -f 2 | rev)
done

cnt=0
for f in $fs_out
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   echo -e "\n-----run ${cnt}-th file: ${dir_path}/sams/${f}.sam > ${dir_path}/bams/${f}.bam-----"
   samtools view -h -b \
      -t ${dir_path}/GCA_007827085.2_ASM782708v2_genomic.fa \
      -o ${dir_path}/bams/${f}.bam ${dir_path}/sams/${f}.sam
   echo -e "\n-----run ${cnt}-th file: ${dir_path}/bams/${f}.bam > ${dir_path}/bams/${f}_sorted.bam-----"
   samtools sort -m 30G \
      -o ${dir_path}/bams/${f}_sorted.bam ${dir_path}/bams/${f}.bam
   rm ${dir_path}/bams/${f}.bam
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

END=$(date +%s)
echo "SAM -> BAM -> sorted BAM time elapsed: $(( $END - $START )) seconds"

################### Mark PCR Duplicates #######################
ncpu=23
START=$(date +%s)

fs_in=$(ls ${dir_path}/bams/*_sorted.bam)
cnt=0
for f in $fs_in
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   fn=$(echo ${f} | rev | cut -d '.' -f 2 | rev)
   echo -e "\n-----run ${cnt}-th file: $f > ${fn}_mkDup.bam-----"
   java -Xmx20G -XX:+AggressiveOpts -XX:+AggressiveHeap\
      -jar picard.jar MarkDuplicates \
      INPUT=$f \
      REMOVE_DUPLICATES=false \
      ASSUME_SORTED=true \
      METRICS_FILE=${dir_path}/bams/${fn}_mkDup_metrics.txt \
      OUTPUT=${fn}_mkDup.bam &
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

END=$(date +%s)
echo "Mark PCR duplicates time elapsed: $(( $END - $START )) seconds"

#### Indexing sorted alignment for fast random access
ncpu=23
START=$(date +%s)

fs_in=$(ls ${dir_path}/bams/*_mkDup.bam)
cnt=0
for f in $fs_in
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   fn=$(echo ${f} | rev | cut -d '.' -f 2 | rev)
   echo -e "\n----- run ${cnt}-th file: $f > ${fn}.bai-----"
   samtools index $f ${fn}.bai
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Index alignments time elapsed: $(( $END - $START )) seconds"

############################### variant calling ###############################
mkdir ${dir_path}/vcfs
ncpu=23
START=$(date +%s)

fs_in=$(ls ${dir_path}/bams/p.cal*_mkDup.bam)
#### we drop the unmatched reads during the demux process here
echo -e "----- 96samples_samtools.vcf -----"

## -C50 -q 25 -Q 30 make things worse?
bcftools mpileup -q 20 -Q 20 \
   --threads $ncpu \
   -f ${dir_path}/GCA_007827085.2_ASM782708v2_genomic.fa \
   -O z \
   $fs_in | \
bcftools call --no-version -m -v \
   --threads $ncpu \
   -O z \
   -o ${dir_path}/vcfs/96samples_samtools.vcf.gz

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

END=$(date +%s)
echo "Variant calling time elapsed: $(( $END - $START )) seconds"