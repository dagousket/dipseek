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

