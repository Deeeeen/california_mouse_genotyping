![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# genotyping
## Overview
This folder contains the code for genotyping steps.

## Contents
**[peromyscus_phase1.sh](peromyscus_phase1.sh)**: This script contains the code for the first phase of genotyping. Including the following steps:  
- Reference Genome Preparation 
- Demultiplexing
- Quality Control
- Quality Trimming

**[peromyscus_phase2.sh](peromyscus_phase2.sh)**: This script contains the code for the second phase of genotyping. Including the following steps:  
- Map to Reference Genome 
- Convert SAM to BAM
- Sort BAM Files
- Mark Duplicates
- Index BAM Files
- Variant Calling

## Prerequisites
### software used
Please install the following software before running the script. The version information here is the ones used on our original analysis. You may want to upgrade them to the newest versions when you use them.
- Bcftools 1.9
- BWA 0.7.12
- Cutadapt 1.9.1
- FastQC 0.11.8
- Fgbio 1.1.0
- MultiQC 1.8
- Picard 1.141
- Samtools 1.3  

### Before running the scripts
Please read the comments inside the script with extra attention to the comments start with #### !!! You may need to make some modifications inside the script according to your computer/cluster setup.

### To run the scripts
1. Change the permission of the script
```
chmod u+x peromyscus_phase1.sh
```
2. Run the script
```
./peromyscus_phase1.sh
```  