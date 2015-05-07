# peak-tool
1st commit v0.3

Tool in c++ to parse the gencode database file and a list of ChIP-Seq peaks in bed file and report detailed promoter / intergenic occupancy.

a) Tested with MACS bed file (more accurate if summit data is used)
b) You can combine two bed files with bedtools and use the output with the co-ChIP version of this program
c) If your system has < 24GB RAM please make sure you have ~25% of your RAM allocated as swap with ZRAM

Use ZRAM on a system with 16GB RAM or use a system with >= 24 GB RAM to avoid swapping and unnecessary frustrations

Make sure you have g++ installed on your system

usage:

gunzip gencode.v19.annotation.gtf.gz 

make ./peak_tool 

./peak_tool ./test.bed > test.genes

Output: The output is in the following format

chrnum peak_mid peak_name peak_score peak_location gene_name strand isomer_name isomer_coding_type tss dist_tss

You can safely ignore rest of the columns which are used with co-chip peak files
