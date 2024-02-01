\### **Step 4 - Preprocessing/cleanup of loaded data**

\### The reconstruction of this scripts mainly focused on quality control of the

\### original metagenome data, improvement of contaminated sequences, and

\### normalization of the depth of data loss. This rigorous quality control procedure

\### ensures the reliability of the subsequent analysis of the microbial data. 

 

\### This script employed the KneadData tool (version 0.7.5),

\### available at http://huttenhower.sph.harvard.edu/kneaddata, 

\### to ensure the quality and integrity of the microbial data under analysis.

\### Trimmomatic points to trimmomatic software, used for trimming the sequences.

\### Trimmomatic software (version 0.39) points to trimmomatic software, 

\### used for trimming the sequences, employing the following parameters: 

\### SLIDINGWINDOW:4:20, MINLEN:50, LEADING:3, and TRAILING:3. 

 

\### This script runs Kneaddata for the provided fastq files to align them to the

\### Human Genome to eliminate any host-associated reads that could act

\### as potential contaminants using the bowtie2 software (version 2.2.3) . 

 **4.1 keneaddata **

````shell
#/usr/bin/bash

# Setup Variables

# wd points to the folder containing all the fastq files 

wd=/SVanalysis/raw_paired_end_fastq.file

 

# wd_rst points to the folder where results of the run are stored 

wd_rst=/SVanalysis/4.1keneaddata

 

# job_pre is a prefix added to intermediate files

job_pre=SVs_4.1kendadata

 

# keanddata_db points to the folder containing the human reference genome

keanddata_db=/SV/db/kneaddata_database/human/

 

# Trimmomatic points to trimmomatic software, used for trimming the sequences

Trimmomatic=/SV/software/trimmomatic-0.39-1

 

# input_list points to a file containing all the fastq files 

input_list=$wd/gapspace_fatsq.list

 

# create necessary directories if not already present 

mkdir -p $out/$job_pre.spbs

mkdir -p $out/$job_pre.rst

 

# iterate through all the fastq files provided in the list 

cat $input_list | while read line

do 

# store the individual fastq filename in sample_id

 sample_id=$line

# Create necessary directories if not already present 

 mkdir -p "$out/$job_pre.spbs/01kneaddata_pbs"

# create job script 

 oo="$out/$job_pre.spbs/01kneaddata_pbs/$job_pre.$sample_id.sh"

# write job script 

 echo "#!/bin/bash" > $oo

 echo "#PBS -N Kneaddata" >> $oo

 echo "#PBS -l nodes=1:ppn=20" >> $oo

 echo "#PBS -j n" >> $oo

 echo "#PBS -q workq" >> $oo

 echo "#PBS -e step4.1.${sample_id}.CRC_SV_import.err" >> $oo

 echo "#PBS -o step4.1.${sample_id}.CRC_SV_import.out" >> $oo

# Create directories for the results 

 echo "mkdir -p $out/$job_pre.rst/01kneaddata" >> $oo

 echo "mkdir -p $out/$job_pre.rst/01kneaddata/$sample_id" >> $oo

#>> step 4.1: Kneaddata <<

 echo "echo 2 kneaddata started on \`date\`" >> $oo

 echo "source activate kneaddata" >> $oo

# Run Kneaddata 

 echo "kneaddata -i $wd/${sample_id}_1.fastq -i $wd/${sample_id}_2.fastq -db $keanddata_db --reorder -t 20 -p 20 -o $out/$job_pre.rst/01kneaddata/$sample_id --log $out/$job_pre.rst/01kneaddata/$sample_id/$sample_id.log" >> $oo 

done
````

 **4.2 Seqtk **

\### Seqtk is a fast and lightweight tool for processing sequences 

\### in the FASTA or FASTQ format. It seamlessly parses both FASTA 

\### and FASTQ files which can also be optionally compressed by gzip.

\### Then, we conducted a stochastic sampling of 10 million reads from 

\### paired-end fastq files originating from a common sample based on 

\### seqtk tool (version 1.3-r117-dirty) 28, employing a uniformly distributed

\### random number generator initialized with a shared seed value.

```shell
#/usr/bin/bash

# Setup Variables
# wd points to the folder containing all the fastq files after quality control of kneaddata 

wd=/SVanalysis/4.1keneaddata/SVs_4.1kendadata.rst

# wd_rst points to the folder where results of the run are stored

out=/SVanalysis/4.2seqtk

# job_pre is a prefix added to intermediate files

job_pre=SVs_4.2seqtk

input_list=$wd/gapspace_fatsq.list

job_pre=SVs_4.2seqtk


# Create Directories

mkdir -p $out/$job_pre.spbs/03seqtk_pbs

mkdir -p $out/$job_pre.rst/03seqtk

 
#>> step 4.2: Seqtk <<

# Performing Seqtk Sample

cat $input_list | while read line

do

    sample_id=$line

    oo="$out/$job_pre.spbs/03seqtk_pbs/$job_pre.seqtk.sh"

    echo "seqtk sample -s100 $wd/$sample_id/${sample_id}_1_kneaddata_paired_1.fastq 10000000 > $out/$job_pre.rst/03seqtk/${sample_id}_1_kneaddata_seqtk_paired_1.fastq" >> $oo

    echo "seqtk sample -s100 $wd/$sample_id/${sample_id}_1_kneaddata_paired_2.fastq 10000000 > $out/$job_pre.rst/03seqtk/${sample_id}_1_kneaddata_seqtk_paired_2.fastq" >> $oo

done