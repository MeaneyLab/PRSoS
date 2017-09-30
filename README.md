# PRS-on-SPARK
PRS-on-SPARK (PRSoS) generates polygenic risk scores for large genotype data, including imputed genotype dosages.



## Installation

To clone the repository, use :
```
git clone https://github.com/seriousNel/PRS-on-SPARK.git
```

## Software requirements

The notebooks and scripts require the following to run :
+ spark-2.0.0 +
+ Pyhon 2.7

Instructions for installing Apache Spark on Linux can be found [here](https://www.santoshsrinivas.com/installing-apache-spark-on-ubuntu-16-04/)

The prerequisite to install spark are :
+ java8 (oracle)
+ scala 
+ sbt

Some extra libraries are required for regression and plotting. To install them, first make sure pip is installed on your computer, then type :
```
cd PRS-on-SPARK
pip install -r requirements.txt
```

## What this pipeline does :
+ Calculate PRS from a genotype file (in .gen or .vcf format) and a GWAS file
+ Correct the strand alignment discrepancies between genotype and GWAS data. 

## What this pipeline cannot do :
+ Perform quality control of genotype data

## Default format :
### GWAS
By default, the GWAS should have the same format as that of a GWAS file obtained from Psychiatric Genomics Consortium (PGC). 

|     snpid|   pval|    or| a1| a2|    CEUaf|
|----------|-------|------|---|---|---------|
| rs3131972| 0.2032| 1.047|  A|  G|  0.16055|
| rs3131969|0.08597| 1.067|  A|  G| 0.133028|
| rs3131967|0.06683| 1.077|  T|  C|        .|
| rs1048488| 0.2808|0.9617|  T|  C| 0.836449|
|rs12562034| 0.8489|0.9931|  A|  G|0.0925926|



You can change your GWAS to the same format, or use optional parameter flags to let the script know about the format you are using. More details below.

### .gen file
from [www.shapeit.fr](http://www.shapeit.fr/pages/m02_formats/gensample.html) :
A .gen file is a SPACE delimited file. Each line corresponds to a single SNP. The first 5 columns are:
Chromosome number [integer]
SNP ID [string]
SNP physical position (bp) [integer]
First allele [string]
Second allele [string]

### .vcf file 
This is a default format for the genotype data returned from [Sanger Institute](https://imputation.sanger.ac.uk/). 
Details about the format can be found [here](http://samtools.github.io/hts-specs/VCFv4.2.pdf). 
PRSoS uses posterior genotype probabilities (FORMAT/GP) as genotyping data.

### Output file
By default, the output format of PRS results and SNP logs is csv. 

## Running command-line script PRS_run.py
### Parameters



A description of the parameters for the script can be obtained by typing : 
```
python PRS_run.py --help
```
Command-line type flags are used to specify how the scores are calculated. 

To run the script, use ```spark-submit```. You can add other parameters for spark before the script if desired. 

```
spark-submit PRS_run.py 
```
Followed by three positional parameters (mandatory) :
```
  GENO                  Name of the Genotype files, can be a name or path, or name patterns with '*'
  GWAS                  Name of the GWAS file, can be a name or path.
  OUTPUT                The path and name for the output file.
```
Followed by some optional parameters. By default, the pipeline assumes the following : 

* To specify the format of genotype file :
```
  --filetype VCF  
```
The type of genotype file used as input, choose between VCF and GEN. Default is VCF.

* To specify the p-value thresholds of the score :
```
  --thresholds 0.5 0.2 0.1 0.05, 0.01 0.001 0.0001
```
Enter one or more float numbers separated by space. Default is 0.5 0.2 0.1 0.05 0.01 0.001 0.0001.

Alternatively you can specify a sequence of thresholds to use :

```
  --threshold_seq 0.1 0.5 0.01
```
After the flag, the first number is the starting point of the sequence, the second is the end point of the sequence, the third number denotes the step size. The above example would yield the sequence 0.1,0.11,0.12.....0.49,0.5. Note the interval is inclusive of the endpoints.
  
  --GWAS_no_header      Adding this parameter signals that there is no header 
                        in the GWAS data input. The default is to assume that 
                        GWAS has column names.
                        
  --log_or              Adding this parameter tells the script to log (base=10) 
                        the effect sizes provided in the GWAS summary-level 
                        data. For example, this would be applied to odds ratios 
                        to get log odds or the beta values of logistic 
                        regression.
                        
  --no_check_ref        Adding this option tells the script to not check 
                        reference allele when determining genoypte calls. 
                        Default is checking.
                        
  --no_maf              The pipeline calculates the allele frequencies in the 
                        genotype population by default, which is used to help 
                        retain as many ambiguous SNPs as possible by comparing 
                        the allele frequencies in the GWAS to make the best 
                        estimate of the strand alignment. Use this flag to 
                        disable this feature and all ambiguous SNPs that would 
                        have been used for PRS caculation will be discarded.
                        
  --snp_log             Add this flag to record the SNPs that are used at each 
                        threshold. It will also report whether the a1 or a2 
                        allele in the genotype data was used as reference for 
                        the risk effect, indicated as 'keep' or 'flip'. Any SNPs
                        that meet the p-value threshold criteria but has allele 
                        names that do not match the allele names in the GWAS 
                        description are indicated in the 'discard' column. This 
                        record will be saved to a file with the name specified 
                        in the Output flag, with .snplog as file extension.

  --sample_file SAMPLE_FILE
                        Path and name of the file that contain the sample
                        labels. It is assumed that the sample labels are
                        already in the same order as in the genotype file.
                        
  --sample_delim SAMPLE_DELIM
                        Delimiter of the sample file. Default is comma. Set 
                        quotation marks around the delimiter when applied.




### Examples :
To calculate PRS from a series of .vcf files, while checking the allele alignment between the genotype and the GWAS, and log transform risk effect, using p-value thresholds of 0.2, 0.1, 0.05:
```
spark-submit PRS_run.py "VCF_number*.vcf" pgc.mdd.clump.txt output.csv --sample_file samplefile.csv --sample_file_id 0 --log_or --thresholds  0.2 0.1 0.05
```
To calculate PRS from a series of .gen files, without checking allele alignments, using a GWAS with no header, and not transform the risk effect, using p-value thresholds of 0.2, 0.1, 0.05:

```
spark-submit PRS_run.py "VCF_number*.vcf" pgc.mdd.clump.txt output.csv --filetype GEN --sample_file samplefile.csv --sample_file_id 0 --no_check_ref --GWAS_no_header --thresholds 0.2 0.1 0.05
```




### Full list of parameters when type `python PRS_run.py --help` :
```
positional arguments:
  GENO                  Name of the genotype files, can be a name or path, or
                        name patterns with wildcard character.
  GWAS                  Name of the GWAS file, can be a name or path.
  OUTPUT                The path and name stem for the output files. One name
                        will be used for the score output, the snp log 
                        (optional), and the regression output. This is similar 
                        to the --out flag in Plink.

optional arguments:
  -h, --help            Show this help message and exit.
  -v, --version         Show program's version number and exit.
  --gwas_id GWAS_ID     Column number in your GWAS that contains SNP ID, with
                        first column being 0. Default is 0.
  --gwas_p GWAS_P       Column number in your GWAS that contains p-value, with
                        first column being 0. Default is 1.
  --gwas_or GWAS_OR     Column number in your GWAS that contains odds-
                        ratio/beta, with first column being 0. Default is 2.
  --gwas_a1 GWAS_A1     Column number in your GWAS that contains allele A1,
                        with first column being 0. Default is 3.
  --gwas_a2 GWAS_A2     Column number in your GWAS that contains allele A2,
                        with first column being 0. Default is 4.
  --gwas_a1f GWAS_A1F   Column number in your GWAS that contains frequency of
                        A1, with first column being 0. Default is 5.
  --filetype {GEN,VCF}  The type of genotype file used as input, choose
                        between VCF and GEN. Default is VCF.
  --thresholds THRESHOLD1 THRESHOLD2 THRESHOLD3 ...
                        The p-value thresholds that filters which SNPs are
                        used from the GWAS. Specifying the p-values simply by
                        input one after another. Default is 0.5, 0.2, 0.1,
                        0.05, 0.01, 0.001, 0.0001.
  --threshold_seq LOWERBOUND UPPERBOUND STEPSIZE
                        Defines a sequence that contains all the p-value
                        thresholds that filters which SNPs are used from the
                        GWAS. Input is three numbers separated by space: lower
                        bound, upper bound, step size. Default is None (i.e., 
                        not used). Defining a sequence automatically overwrites 
                        the threshold list defined under --thresholds.
  --GWAS_delim "GWAS_DELIM"
                        Delimiter of the GWAS file. Default is tab-delimited. 
                        Set quotation marks around the delimiter when applied.
  --GWAS_no_header      Adding this parameter signals that there is no header 
                        in the GWAS data input. The default is to assume that 
                        GWAS has column names.
  --log_or              Adding this parameter tells the script to log (natural 
                        base) the effect sizes provided in the GWAS 
                        summary-level data. For example, this would be applied
                        to odds ratios to get log odds or the beta values of 
                        logistic regression.
  --no_check_ref        Adding this option tells the script to not check 
                        reference allele when determining genotype calls. 
                        Default is checking.
  --app_name APP_NAME   Give your spark application a name. Default is PRS.
  --sample_file SAMPLE_FILE
                        Path and name of the file that contain the sample
                        labels. It is assumed that the sample labels are
                        already in the same order as in the genotype file.
  --sample_file_delim "SAMPLE_DELIM"
                        Delimiter of the sample file. Default is comma. Set 
                        quotation marks around the delimiter when applied.
  --sample_file_ID ID_COLUMN1 ID_COLUMN2 ID_COLUMN3 ...
                        Specify which columns in the sample file are used as
                        labels. Can use one integer to specify one column, or
                        multiple integers to specify multiple columns, with 
                        first column being 0. Default is 0.
  --sample_file_skip SAMPLE_SKIP
                        Specify how many header lines to ignore in the sample 
                        file. Default is 1, which assumes that the sample files 
                        has column names and the labels start on the second 
                        line.
  --no_maf              The pipeline calculates the allele frequencies in the 
                        genotype population by default, which is used to help 
                        retain as many ambiguous SNPs as possible by comparing 
                        the allele frequencies in the GWAS to make the best 
                        estimate of the strand alignment. Use this flag to 
                        disable this feature and all ambiguous SNPs that would 
                        have been used for PRS caculation will be discarded.
  --snp_log             Add this flag to record the SNPs that are used at each 
                        threshold. It will also report whether the A1 or A2 
                        allele in the genotype data was used as reference for 
                        the risk effect, indicated as 'keep' or 'flip'. Any SNPs
                        that meet the p-value threshold criteria but has allele 
                        names that do not match the allele names in the GWAS 
                        description are indicated in the 'discard' column. This 
                        record will be saved to a file with the name specified 
                        in the Output flag, with .snplog as file extension.
  --check_dup           Add this flag if you want to check for and discard
                        SNPs that are duplicated, which will take extra time.
                        By default, the script will assume there is no
                        duplicate SNPs.
  --pheno_file PHENO_FILE
                        Specify the path to the data file for the phenotype.
                        It is assumed that the phenotype data is organized in
                        the same order as the samples in the genotype file.
  --pheno_columns PHENO_COLUMN1 PHENO_COLUMN2 PHENO_COLUMN3 ...
                        Specify which columns that the phenotype data is in
                        the provided phenotype data file. Multiple column
                        numbers can be specified to conduct regression with
                        multiple phenotypes. Default is the first column.
  --pheno_delim "PHENO_DELIM"
                        Specify the delimiter for the phenotype data file.
                        Default is comma. Set quotation marks around the 
                        delimiter when applied.
  --pheno_no_header     Specify whether the phenotype has a header row
  --covar_columns COVAR_COLUMN1 COVAR_COLUMN2 COVAR_COLUMN3 ...
                        Specify which columns that the phenotype data is in
                        the provided phenotype data file. Multiple column
                        numbers can be specified to conduct regression with
                        multiple phenotypes. Default is the first column.

``` 
