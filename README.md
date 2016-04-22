Comments / Criticisms / Feature requests: https://www.biostars.org/p/157346/

#peak-tool-multi 
Tool to annotate human hg19 or mouse mm10 aligned ChIP-Seq peak files. This tool also takes multiple ChIP-Seq peak files from different experiments and finds neighbors of the primary peak file and annotate it.

20GB RAM required on Linux 

16GB RAM required on MAC OS X >= 10.9

This is a tool in c++ to parse the gencode annotation database file and list(s) of ChIP-Seq peaks in bed file format.
Report detailed promoter / gene-body / intergenic / enhancer occupancy in human or mouse. 
Multi peak option reports neighboring peaks in order to elicit co-acting transcription factors. For e.g. With this feature you can correlate your peaks with ENCODE ChIP-Seq data or a related ChIP-Seq data.

Compiling::

>gunzip gencode.v19.annotation.gtf.gz 

>gunzip enhancers.bed.gz 

>make  

Running::

Example1:
./peak_tool_multi ./test.bed > test.genes.txt
-------------------------------------------------
Example Input in BED format

chr1	1778750	1778751	MACS_peak_51	102.12

chr1	1933483	1933484	MACS_peak_57	93.87

chr1	3446145	3446146	MACS_peak_91	85.75

chr1	4003155	4003156	MACS_peak_104	58.09

chr1	5787471	5787472	MACS_peak_121	1325.16

chr1	6473142	6473143	MACS_peak_138	988.16

chr1	7259083	7259084	MACS_peak_154	60.32

chr1	8031408	8031409	MACS_peak_181	750.45

chr1	8319346	8319347	MACS_peak_195	3100.00

chr1	8360245	8360246	MACS_peak_198	261.50

The output is in the following format
-------------------------------------------------
chrnum peak_mid peak_name peak_score peak_location gene_name strand isoform isoform_coding_type tss dist_tss
-------------------------------------------------
chr1	1778750	MACS_peak_51	102.12	INTRON	GNB1	-	GNB1-001	protein_coding	1822495	43745

chr1	1933483	MACS_peak_57	93.87	INTRON	C1orf222	-	C1orf222-007	retained_intron	1935276	1793

chr1	3446145	MACS_peak_91	85.75	INTRON	MEGF6	-	MEGF6-001	protein_coding	3448012	1867

chr1	4003155	MACS_peak_104	58.09	ENHANCER	.	.	.	.	.	.

chr1	5787471	MACS_peak_121	1325.16	INTERGENIC	.	.	.	.	.	.

chr1	6473142	MACS_peak_138	988.16	EXON	HES2	-	HES2-002	protein_coding	6484730	11588

chr1	7259083	MACS_peak_154	60.32	INTRON	CAMTA1	+	CAMTA1-001	protein_coding	6845384	413699

chr1	8031408	MACS_peak_181	750.45	EXON	PARK7	+	PARK7-004	protein_coding	8014351	17057

chr1	8319346	MACS_peak_195	3100	ENHANCER	.	.	.	.	.	.
-------------------------------------------------
Example2:
./peak_tool_multi EXP1.bed EXP2.bed EXP3.bed EXP4.bed > EXP1_EXP2_EXP3_EXP4.genes.txt
--------------------------------------------------
Output of multi-peak feature
-------------------------------------------------
chr1	714304	MACS_peak_2	158.84	INTERGENIC	.	.	.	.	.	.	714017	MACS_peak_1	68.78	-287	.	714039	MACS_peak_1	170.93	-265	.	714023	MACS_peak_1	245.28	-281	.

chr1	769360	MACS_peak_8	421.42	INTERGENIC	.	.	.	.	.	.	769283	MACS_peak_3	55.43	-77	.	769273	MACS_peak_5	154.26	-87	.	769292	MACS_peak_3	71.24	-68	.

chr1	840155	MACS_peak_10	96.66	INTERGENIC	.	.	.	.	.	.	840097	MACS_peak_6	165.5	-58	.	840026	MACS_peak_11	57.72	-129	.	840075	MACS_peak_8	134.82	-80	.

chr1	840738	MACS_peak_11	137.93	ENHANCER	.	.	.	.	.	.	840097	MACS_peak_6	165.5	-641	.	840026	MACS_peak_11	57.72	-712	.	840075	MACS_peak_8	134.82	-663	.

chr1	911680	MACS_peak_15	241.99	PROMOTER	C1orf170	-	C1orf170-002	retained_intron	912021	341	911740	MACS_peak_20	271.3	60	281	911707	MACS_peak_31	285.65	27	314	911709	MACS_peak_23	251.1	29	312

chr1	994706	MACS_peak_18	95.46	ENHANCER	.	.	.	.	.	.	994613	MACS_peak_35	66.76	-93	.	995290	MACS_peak_52	57.49	584	.	994655	MACS_peak_37	61.84	-51	.

chr1	1003263	MACS_peak_19	71.56	INTERGENIC	.	.	.	.	.	.	1003034	MACS_peak_37	252.84	-229	.	1003088	MACS_peak_54	281.96	-175	.	1003078	MACS_peak_39	643.19	-185	.

chr1	1003982	MACS_peak_20	104.28	ENHANCER	.	.	.	.	.	.	1003034	MACS_peak_37	252.84	-948	.	1003088	MACS_peak_54	281.96	-894	.	1003078	MACS_peak_39	643.19	-904	.

chr1	1098334	MACS_peak_23	149.14	INTERGENIC	.	.	.	.	.	.	1098987	MACS_peak_50	185.58	653	.	1098482	MACS_peak_70	130.48	148	.	1098535	MACS_peak_51	138.52	201	.
-------------------------------------------------
Notes::

a) Make sure you have g++ / clang installed on your system 

b) Tested with MACS bed file (Peak summit is preferable)

c) Simplistic assumption of Promoter regions (-2kb - +0.5kb) (can be changed in the program)

d) To check multiple peak files (1k region, can be changed in the program) and elicit neighboring peaks use peak_tool_multi in the multi-ChIP directory.

e) Works perfectly on a Macbook Pro Retina laptop with 16GB RAM and MAC OS X > (10.9.x).
   The OS efficiently compresses the indexes data-structure in memory.

Linux specific notes

f) Use ZRAM on a Linux system with 16GB RAM or use a system with >= 20 GB RAM to avoid swapping and unnecessary frustrations
http://askubuntu.com/questions/174579/how-do-i-use-zram

g) If your Linux system has < 20GB RAM please make sure you have ~30%-50% of your RAM allocated as swap with ZRAM


edit the following zram config file in Ubuntu Linux
vim /etc/init/zram-config.conf
> mem=$(((totalmem * 30/100 / ${NRDEVICES}) * 1024))



