**### Step 9 - Analysis of Species Abundance Differences**

\### Table processing, which processes the downloaded tables into 

\### specific formats that need to be used for the analysis of the following 

\### scripts in lefse. - c 1 represents the use of the first row as the 

\### grouping category (main grouping)- S 2 represents the subcategory 

\### (secondary grouping) in the second row of the table, and - u 3 

\### represents the sample name in the third row. If the - s and - u 

\### parameters do not need to be given a parameter value of -1 

\### (cannot be left blank)- The o parameter represents a total abundance 

\### of 1000000 for species at the same classification level in the sample.

 

\### Prepare input files

\### sample.txt

\### S.taxName.count.tsv, an abundance file containing species classification information

 

\```shell

lefse_format_input.py hmp_aerobiosis_small.txt hmp_aerobiosis_small.in -c 1 -s 2 -u 3 -o 1000000  

 

\### Data preparation - data section

awk -F "\t" -v 'OFS=\t' '{$1=$NF; $NF=""; {print $0}}' S.taxName.count.tsv | \

sed 1d | sed 's/; /|/g' > S.lefse.tmp

 

\### Data preparation - class section

head -n 1 S.lefse.tmp | sed 's/\t/\n/g' | sed '1d' | \

while read sp ;do awk '$2=="'$sp'" {printf "\t"$1 }' sample.txt ;done > S.lefse.header

 

\### Merge

cat S.lefse.header S.lefse.tmp > S.lefse

 

\### Format conversion

format_input.py S.lefse \ # Input

S.lefse.in \ # Output

-c 1 \ # Class in this row

-s -1 \ # Subclass in this row

-u 2 \ # Subject in this row

-o 1000000 # Normalize to 1M

 

\### Run LEfSe

run_lefse.py S.lefse.in \ # Input

S.lefse.out \ # Output

-l 2 # LDA threshold

 

\### LDA plot

plot_res.py S.lefse.out \ # Input file

S.lefse.LDA.pdf \ # Output file

--format pdf # Output format

 

\### Phylogenetic tree plot

plot_cladogram.py S.lefse.out \ # Input file

S.lefse.cladogram.pdf \ # Output file

--format pdf \ # Output format

--labeled_start_lev 1

\```