# peak-tool
Tool in c++ to parse the gencode database file and a list of ChIP-Seq peaks in bed file and report detailed promoter / intergenic occupancy 

1st commit v0.3, RAM USAGE 24GB / use zRAM on a system with 16GB RAM or use a system with atleast 24GB to avoid swapping and unnecessary frustations

Make sure you have g++ installed on your system

usage:

gunzip gencode.v19.annotation.gtf.gz 
make ./peak_tool 
./test.bed > test.genes

Output: The output is in the following format chrnum peak_mid peak_name peak_score peak_location gene_name strand isomer_name isomer_coding_type tss dist_tss

You can safely ignore rest of the columns which are used with co-chip peak files
