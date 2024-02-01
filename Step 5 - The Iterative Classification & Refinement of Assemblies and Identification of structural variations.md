**### Step 5 - The Iterative Classification & Refinement of Assemblies and Identification of structural variations**

 

**5.1 The Iterative Classification & Refinement of Assemblies **

\### The identification of structural variations in this part of the protocol 

\### primarily relies on the approach taken in the article of David Zeevi et al., 2019. 

\### Importantly, the sequencing depth may vary to some extent between 

\### datasets, and the successful detection of structural variations depends in part

\### on the sequencing depth of the data. Therefore, the sequencing depth should

\### generally not be less than 10 billion. The following examples correspond 

\### to the best choices of SGVFinder and parameters for analysis used in 

\### *David Zeevi et al., 2019* . 

\### In this step, preprocessed metagenomic sequencing reads are mapped 

\### to the reference database comprising 3,953 representative microbial genomes 

\### of proGenomes database which were common commensal bacteria in the gut.

\### To ensure accuracy in read assignments, the Iterative Classification and 

\### Refinement of Assemblies (ICRA) framework will employed subsequent 

\### to completion of preliminary mapping by GEM Mapper. ICRA involved the

\### reassignment of ambiguous reads with multiple alignments based on 

\### their respective mapping quality and genomic coverage.

 

##### ### ICRA

\```shell

\#Settings for ICRA generally

\# Setup Variables

\#/usr/bin

\# wd points to the folder containing all the fastq files after quality control of kneaddata 

wd=/SVanalysis/4.2seqtk/SVs_4.2seqtk

\# wd_rst points to the folder where results of the run are stored 

out=/SVanalysis/5.1ICRA

\# job_pre is a prefix added to intermediate files

job_pre=SVs_5.1ICRA

input_list=$wd/gapspace_fatsq.list

job_pre=SVs_5.1ICRA

 

 

\# Create directories to store PBS jobs and ICRA result files

mkdir -p $wd/$job_pre.spbs/02icra_pbs

mkdir -p $wd/$job_pre.rst/02icra

 

cat $input_list | while read line

do 

  \# Get sample name

​     sample_id=$line

​     oo="$wd/$job_pre.spbs/$job_pre.$sample_id.step3_ICRA.sh"

​     \# Set PBS job name and parameters

​     echo "#!/bin/bash" > $oo

​     echo "#PBS -N ICRA.${sample_id}" >> $oo

​     echo "#PBS -q workq" >> $oo

​     echo "#PBS -e step5.1_CRC_SV_${sample_id}.import.err" >> $oo

echo "#PBS -o step5.1_CRC_SV_${sample_id}.import.out" >> $oo

​     echo "#PBS -l nodes=1:ppn=40" >> $oo

​     echo "#PBS -j n" >> $oo

 

​     echo "#PBS -o step3_CRC_SV_${sample_id}.import.out" >> $oo

 

​     \#>> step 5.1: icra <<

​     \# Execute ICRA

​     echo "source activate python2.7.9" >> $oo

​     echo "echo step 5.1 - icra started on \`date\`" >> $oo

​     echo "python ICRA_cmd.py $wd/$job_pre.rst/02icra $wd/$job_pre.rst/01kneaddata/$sample_id/${sample_id}_1_kneaddata_paired --pe --use_theta --debug" >> $oo

​     echo "echo step 5.1 - icra finished on \`date\`" >> $oo

done

\```

 

**5.2 Identification of structural variations **

\### Based on the reassigned metagenomic reads,the microbial SVs 

\### of metagenomic samples using SGVFinder with default parameters 

\### will be futher identified. SGVFinder was devised to split the reference 

\### genomes into genomic bins and then examining the coverage of 

\### genomic bins across all samples to identify highly variable 

\### genomic segments and detect SVs. SGVFinder incorporates 

\### a detection framework capable of identifying two distinct classes 

\### of structural variations (SVs): deletion SVs (dSVs) and variable SVs (vSVs). 

 

\### The methodology involves evaluating the deletion percentage of 

\### a genomic segment across the population. Based on specific thresholds,

\### different analytical approaches are employed. When the deletion 

\### percentage is below 25%, the standardized coverage is computed for 

\### the vSV. For deletion percentages between 25% and 75%, only the 

\### presence or absence status of the genomic segment (dSV) is 

\### considered. Regions with deletion percentages exceeding 75% are

\### excluded from the analysis. This systematic approach allows for 

\### precise SV classification and analysis.

 

\### In our datasets, we detected SVs using default parameters. All bacterial 

\### species subjected to structural variation (SV) calling exhibited a 

\### minimum presence of 5% across the entire sample set.

 

**### Perfile**

\```shell

\#/usr/bin

\# wd points to the folder containing all the GEM files after reassigned by ICRA

wd=/SVanalysis/5.1ICRA/SVs_5.1ICRA

 

\# wd_rst points to the folder where results of the run are stored 

out=/SVanalysis/5.2Perfile

 

\# job_pre is a prefix added to intermediate files

job_pre=SVs_5.2Perfile

input_list=$wd/gapspace_fatsq.list

job_pre=SVs_5.2Perfile

 

\# Create Directories

mkdir -p $wd/$job_pre.spbs/5.2Perfile_pbs

mkdir -p $wd/$job_pre.rst/5.2Perfile

 

cat $input_list | while read line

do

  \# Read input_list and run commands for each line

​     sample_id=$line

​     oo="$wd/$job_pre.spbs/03perfile_pbs/$job_pre.$sample_id.step5.2_PerFile.sh"

​     

​     \# Set PBS job name and parameters

​     echo "#!/bin/bash" > $oo

​     echo "#PBS -N ${sample_id}.PerFile" >> $oo

​     echo "#PBS -l nodes=1:ppn=10" >> $oo

​     echo "#PBS -j n" >> $oo

​     echo "#PBS -q workq" >> $oo

​     echo "#PBS -e step5.2_CRC_SV_${sample_id}.import.err" >> $oo

​     echo "#PBS -o step5.2_CRC_SV_${sample_id}.import.out" >> $oo

 

​     \#>> step 5.2: PerFile <<

​     \# Execute PerFile

​     echo "source activate python2.7.9" >> $oo

​     echo "echo step 5 - PerFile started on \`date\`" >> $oo

​     echo "python ./SGVFinder/src/SGVF_PerFile_cmd.py $wd/SVs_perSample_CRC-HC.rst/02icra/${sample_id}_1_kneaddata_paired.jsdel $wd/$job_pre.rst/03perfile/${sample_id}.PerFile.jsdel 150 --x_coverage 0.01 --rate_param 10" >> $oo

​     echo "echo step 5 - PerFile finished on \`date\`" >> $oo

done

\```

 

**### SGVFinder**

\```shell

\#/usr/bin

\# wd points to the folder containing all the GEM files after reassigned by ICRA

wd=/SVanalysis/5.1ICRA/SVs_5.2Perfile

 

\# wd_rst points to the folder where results of the run are stored 

out=/SVanalysis/5.3sgvf

 

\# job_pre is a prefix added to intermediate files

job_pre=SVs_5.3sgvf

 

input_list=$wd/gapspace_fatsq.list

job_pre=SVs_5.3sgvf

 

 

\# Create Directories

mkdir -p $wd/${job_pre}.rst/5.3sgvf

mkdir -p $wd/${job_pre}.spbs/5.3sgvf_pbs

 

cat $input_list | while read line

do

​     read -r -a array <<< "$line"

​     sample_id="${array[0]}"

​     glob="${array[1]}"

​     n="${array[2]}"

​     oo="$wd/${job_pre}.spbs/04sgvf_pbs/${job_pre}.${sample_id}.step5_sgvf.sh"

​     

​     \# Create script

​     echo "#!/bin/bash" > $oo

​     echo "#PBS -N ${sample_id}.SVGF" >> $oo

​     echo "#PBS -l nodes=1:ppn=40" >> $oo

​     echo "#PBS -j n" >> $oo

​     echo "#PBS -q workq" >> $oo

​     echo "#PBS -e step5.3_SGVF_${sample_id}_import.err" >> $oo

​     echo "#PBS -o step5.3_SGVF_${sample_id}_step5_import.out" >> $oo

 

​     \#>> step 5.3: SGVF_cmd <<

​     \# Call SGVF_cmd program

​     echo "source activate python2.7.9" >> $oo

​     echo "echo Step 6 - SVGF started on \`date\`" >> $oo

​     echo "mkdir -p $wd/${job_pre}.rst/04sgvf/${sample_id}/${sample_id}.html" >> $oo

​     echo "python ./SGVFinder/src/SGVF_cmd_changed.py \"$glob\" $wd/$job_pre.rst/04sgvf/$sample_id/$sample_id.dsgv.csv $wd/$job_pre.rst/04sgvf/$sample_id/$sample_id.vsgv.csv --min_samp_cutoff $n --x_coverage 0.01 --rate_param 10 --browser_path $wd/$job_pre.rst/04sgvf/$sample_id/$sample_id.html --csv_output --byorig" >> $oo

​     echo "echo Step 5.3 - SVGF finished on \`date\`" >> $oo

done

\```

 

\### **5.3 transform the df file format to CSV file**

\```python

import pandas as pd

dsgv_df = pd.read_pickle("dsgv.df")

dsgv_df.to_csv("dsgv.csv")

vsgv_df = pd.read_pickle("vsgv.df")

vsgv_df.to_csv("vsgv.csv")

\```

 

 