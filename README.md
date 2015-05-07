# peak-tool
1st commit v0.3

Tool in c++ to parse the gencode database file and a list of ChIP-Seq peaks in bed file and report detailed promoter / intergenic occupancy.

a) Tested with MACS bed file (Output is more accurate if summit data is used)

b) Simplistic assumption of Promoter regions (-2kb - +0.5kb) (can be changed in the program)

b) You can combine two bed files with bedtools and use the output with the co-ChIP version of this program

c) Works perfectly on a Macbook Pro Retina laptop with 16GB RAM and MAC OS X > (10.9.x).
   The OS nicely compresses the index data-structures in memory.

Linux specific 

d) Use ZRAM on a Linux system with 16GB RAM or use a system with >= 20 GB RAM to avoid swapping and unnecessary frustrations

e) If your Linux system has < 20GB RAM please make sure you have ~30% of your RAM allocated as swap with ZRAM

edit the following file in Ubuntu Linux
vim /etc/init/zram-config.conf
> mem=$(((totalmem * 30/100 / ${NRDEVICES}) * 1024))


f) Make sure you have g++  / XCode (OS X) installed on your system

Usage:

gunzip gencode.v19.annotation.gtf.gz 

make ./peak_tool 

./peak_tool ./test.bed > test.genes

Output: The output is in the following format

chrnum peak_mid peak_name peak_score peak_location gene_name strand isomer_name isomer_coding_type tss dist_tss

You can safely ignore rest of the columns which are used with co-chip peak files
