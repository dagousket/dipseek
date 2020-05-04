# Dipseek

This repository conatains an R script built to detect valleys in histone signal. Currently in development.

```
Usage: detect_valley.R [options]


Options:
	-f CHARACTER, --file=CHARACTER
		dataset file name

	-o CHARACTER, --out=CHARACTER
		output file name [default= valleys.out]

	-l INTEGER, --read_length=INTEGER
		Average read extended length used for coverage computing

	-n LOGICAL, --testrun=LOGICAL
		If set as TRUE, will sample 10 peaks for test run [default= FALSE]

	-c CHARACTER, --chromosome=CHARACTER
		Regex of the chromosome you want to keep for the analysis

	-h, --help
		Show this help message and exit
```

The usual pipeline I am currently testing goes through the following steps :

1. Get BAM for H3K27ac ChIP data and if possible input.

2. Call peaks with MACS with both the Broad and Narrow options (keep in mind the read extension used)

3. Extract Broad peaks overlapping at least one narrow peak with minimum threshold score (e.g. 5% FDR).

4. Use Deeptools to generate a Bedgraph file of H3K27ac signal, with no normalisation (gives the raw read count) and the highest resolution (1bp). Make sure you use the proper read extension as used in MACS2.

5. Intersect the Bedgraph file with the selected Broad regions while keeping the assigned peak name.

6. Run Dip first script, merge result files without header and run Dip second script.

7. Extract broad peaks and valleys not overlapping TSS and/or TSS upstream region (1.5kb for human , 500bp for drosophila). If needed, select on FDR threshold the valleys (by default I don't). Extract Fasta sequences from the peaks. Mask repeat sequence without subtracting them with bedtools maskfasta.

8. Run oligo-analysis within peak motifs in RSAT. Test both Dip and Broad without control and Dip vs Broad.



