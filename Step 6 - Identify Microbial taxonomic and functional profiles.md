**### Step 6 - Identify Microbial taxonomic and functional profiles**

 

\### In this step, we conducted taxonomic and functional profiling of 

\### bacterial organisms using Kraken2 (version 2.1.2) and 

\### The HMP Unified Metabolic Analysis Network 3.0 (HUMAnN 3.0).

 

\### Taxonomic Classification: Normalized metagenomic reads obtained

\### from the step4 were subjected to taxonomic classification using Kraken2. 

\### This step assigned taxonomic labels to each read, identifying the 

\### microbial species and genera they belong to.

 

\### Taxonomic Abundance Estimation:To estimate taxonomic abundance

\### with high accuracy, particularly at the species and genus levels, 

\### we employed Bracken (version 2.6.2). Bracken built upon the 

\### taxonomic tree provided by Kraken2 and generated abundance 

\### estimates for each taxonomic category. The resulting read counts

\### for each microbial species were transformed into relative abundance values.

 

\###

\###

 

\```shell

\#/usr/bin

\# wd points to the folder containing all the fastq files after quality control and normalization

wd=/SVanalysis/4.2seqtk/SVs_4.2seqtk

 

\# wd_rst points to the folder where results of the run are stored 

out=/SVanalysis/6.1kraken

 

\# job_pre is a prefix added to intermediate files

job_pre=SVs_6.1kraken

 

input_list=$wd/gapspace_fatsq.list

job_pre=SVs_6.1kraken

 

\# Create Directories

mkdir -p $wd/$job_pre.spbs/06.1kraken_pbs

mkdir -p $wd/$job_pre.rst/06.1kraken

 

\# Kraken analysis

cat $input_list | while read line

do

​     sample_id=$line

​     oo="$wd/$job_pre.spbs/06.1kraken_pbs/$job_pre.$sample_id.step6_kraken.sh"

​     

​     echo "#!/bin/bash" > $oo

​     echo "#PBS -N kraken.${sample_id}" >> $oo

​     echo "#PBS -l nodes=1:ppn=20" >> $oo

​     echo "#PBS -j n" >> $oo

​     echo "#PBS -q workq" >> $oo

​     echo "#PBS -e step6.1_kraken_${sample_id}.import.err" >> $oo

​     echo "#PBS -o step6.1_kraken_${sample_id}.import.out" >> $oo

 

 

 

​     \#>> step 6.1: kraken <<

​     \# Run kraken2 and assign labels to reads with the constructed database

​     echo "echo Step 6.1- kraken started on \`date\`" >> $oo

​     echo "kraken2 --db /public/apps/anaconda2/bin/kraken2_database --threads 20 --report $wd/$job_pre.rst/06.1kraken/${sample_id}.kraken.report --output $wd/$job_pre.rst/06.1kraken/${sample_id}.kraken.output --paired $sample_id/${sample_id}_1_kneaddata_paired_1.fastq $sample_id/${sample_id}_1_kneaddata_paired_2.fastq" >> $oo

​     echo "echo Step 6.1 - kraken finished on \`date\`" >> $oo

done

 

\# Braken analysis

\# wd points to the folder containing all the report files after kraken

wd=/SVanalysis/4.2seqtk/SVs_6.1kraken

 

\# wd_rst points to the folder where results of the run are stored 

out=/SVanalysis/6.2braken

 

\# job_pre is a prefix added to intermediate files

job_pre=SVs_6.2braken

 

input_list=$wd/gapspace_fatsq.list

job_pre=SVs_6.2braken

 

\# Create Directories

mkdir -p $wd/$job_pre.spbs/6.2braken_pbs

mkdir -p $wd/$job_pre.rst/6.2braken

 

\# Loop through the input list of sample IDs to execute braken

cat $input_list | while read line

do

​     sample_id=$line

​     oo="$wd/$job_pre.spbs/08braken_pbs/$job_pre.$sample_id.step7_braken.sh"

​     

​     echo "#!/bin/bash" > $oo

​     echo "#PBS -N kraken.${sample_id}" >> $oo

​     echo "#PBS -l nodes=1:ppn=20" >> $oo

​     echo "#PBS -j n" >> $oo

​     echo "#PBS -q workq" >> $oo

​     echo "#PBS -e step6_kraken_${sample_id}.import.err" >> $oo

​     echo "#PBS -o step6_kraken_${sample_id}.import.out" >> $oo

 

​     \#>> step 6.2: braken <<

​     

​     echo "echo Step 6.2 - braken started on \`date\`" >> $oo

​     echo "/public/home/cuiying/03Project/03project_public4/assis_CRC/fq/Bracken/bracken -d /public/apps/anaconda2/bin/kraken2_database -i $wd/$job_pre.rst/06.1kraken/${sample_id}.kraken.report -o $wd/$job_pre.rst/07.1braken/${sample_id}.bracken -w $wd/$job_pre.rst/07.1braken/${sample_id}.bracken.report -r 150 -l S" >> $oo

​     echo "/public/home/cuiying/03Project/03project_public4/assis_CRC/fq/KrakenTools/kreport2mpa.py --display-header -r $wd/$job_pre.rst/07.1braken/${sample_id}.bracken.report -o $wd/$job_pre.rst/07.1braken/${sample_id}.bracken.new.report" >> $oo

​     echo "echo Step 6.2 - braken finished on \`date\`" >> $oo

done

\```

 

\### We employed the HMP Unified Metabolic Analysis Network 3.0 

\### (HUMAnN 3.0) for the analysis of functional annotations. To begin, 

\### the clean paired-end sequencing data were merged into a single 

\### fastq file, streamlining the data for subsequent analysis.

\#

\### Functional Annotation:The HUMAnN 3.0 toolkit was then utilized 

\### to annotate the complete functionalities of genes and proteins. 

\### This annotation process was facilitated by leveraging both Bowtie 

\### and DIAMOND (version 31).

\#

\### Database Source:For gene and protein annotation, we employed 

\### the UniRef90 database (version 2021.03) as the primary source 

\### of reference. This database was used to map the sequences 

\### to known functional elements.

\#

\### Quantification:Genes and pathways were quantified using 

\### units of RPKs (reads per kilobase). Additionally, the tool provided 

\### the flexibility to normalize the data to either relative abundance 

\### or copies per million (CPM) units, depending on the specific

\### analysis requirements.



#!/bin/bash

\### Iterate through sample IDs in 'id.txt'

for sample_id in $(cat id.txt); do

​    perl concat_paired_end.pl -p 4 --no_R_match -o ./paired_fastq/ cat_reads rst/03seqtk/${sample_id}*.fastq
done

#!/bin/bash

\### Iterate through sample IDs in 'id.txt'

for sample_id in $(cat id.txt); do

​    humann --input ./paired_fastq/${sample_id}.fastq --output ./humann_out/${sample_id} \
​    --metaphlan-options "--bowtie2db /metaphlan_databases/metaphlan_databases/" --threads 20
done





 