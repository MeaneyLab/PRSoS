# PRS-on-SPARK
PRS-on-SPARK (PRSoS) generates polygenic risk scores (PRS) for large genotype data, including imputed genotype posterior probabilites. 
It can use multiple cores to increase processing efficiency (i.e., reduce processing time). 
PRSoS is compatible with Linux, Mac OS, and Windows. It runs using [Apache Spark](https://spark.apache.org) and [Python](https://www.python.org/download/releases/2.7.2/).

## Contact Information
Lawrence M. Chen: lawrence.m.chen@mail.mcgill.ca


## Installation

To clone the repository, use:
```
git clone https://github.com/MeaneyLab/PRSoS.git
```

## Software requirements and installation instructions

The notebooks and scripts require the following to run:
+ Spark-2.0.0 + (Apache)
+ Python 2.7 (not Python 3.0)

The prerequisite to install Spark are:
+ Java 8 (Oracle)
+ Scala 
+ sbt

The required Python modules are:
+ matplotlib
+ pandas
+ numpy
+ statsmodels

### Linux

Instructions for installing Apache Spark on Linux can be found [here](https://www.santoshsrinivas.com/installing-apache-spark-on-ubuntu-16-04/).

The Python modules can be installed independently. To install them all at once, first make sure pip is installed on your computer, then run:
```
cd PRSoS
pip install -r requirements.txt
```

### Mac OS

Requires [Homebrew](https://brew.sh) to perform the required installations.

1. Install code-select by typing the following in Terminal:
    ```
    xcode-select --install
    ```

2. Install Scala:
    ```
    brew install scala
    ```

3. Install Spark:
    ```
    brew install apache-spark
    export SPARK_HOME=/usr/local/Cellar/apache-spark/2.1.1/libexec
    export PYTHONPATH=/usr/local/Cellar/apache-spark/2.1.1/libexec/python/:$PYTHONPATH
    ```

    (Note: The version of Spark should be changed to be the one you install.)

4. Install Python:
    ```
    brew install python
    pip install psutil
    ```

### Windows (Win10)

1. Download the installation file for Windows from the following link: https://www.continuum.io/downloads

2. Install Anaconda Python by double clicking the downloaded file and use the default options.

3. Download and install some prerequisite applications:
    + scala: http://www.scala-lang.org/download/
    + java7 sdk (if neccessary): http://www.oracle.com/technetwork/java/javase/downloads/index.html

4. Download winutils.exe from HortonWorks repo or Steve Loughran's GitHub (https://github.com/steveloughran/winutils) and store it in a bin directory under a created Hadoop home directory (e.g., C:\Users\hadoop\bin).

5. Afterwards, open cmd and run the following:
    ```
    setx SPARK_HOME "C:\Users\spark-2.1.1-bin-hadoop2.7"
    setx PATH " %SPARK_HOME%\bin;%PATH%"
    setx HADOOP_HOME "C:\Users\hadoop"
    ```
    
6. To install PySpark, run: 
    ```
    conda install -c conda-forge pyspark=2.1.1
    ```

    (If "conda" command is not found, run the above command line in Anaconda Prompt--you can find this through the Taskbar search menu)
    
7. To install the Python modules, you can use conda or pip.

Alternative instructions for installing Spark can be found [here](http://www.informit.com/articles/article.aspx?p=2755929&seqNum=3).

## What this pipeline does
+ Match the strand alignment between genotype and GWAS data (by default), and then
+ Calculate PRS from the genotype data set (in .gen or .vcf file format), weighted by the associated effect size (e.g., log odds ratio), and subsetted by the selected associated p-value thresholds in the GWAS

## What this pipeline cannot do
+ Perform quality control of genotype data (depending on the file format, quality control can be performed using [QCTOOL](http://www.well.ox.ac.uk/~gav/qctool_v2/), [PLINK](https://www.cog-genomics.org/plink/1.9/), [VCFtools](http://vcftools.sourceforge.net), etc.)
+ Perform linkage disequilibrium clumping or pruning (linkage disequilibrium clumping or pruning can can be performed using [PLINK](https://www.cog-genomics.org/plink/1.9/). However, note that dosage data (e.g., posterior probabilities) are not compatible with PLINK file formats)

## Default format
### GWAS

The GWAS data input requires the following columns:

snpid = rs ID of the SNP

pval = p-value of the association in the GWAS data

or = odds ratio in the GWAS data, or it can be log odds ratio or beta (if it is the odds ratio, add `--log_or` flag to use log odds ratio)

a1 = reference allele for the odds ratio

a2 = alternative allele


|     snpid|   pval|    or| a1| a2|   a1freq|
|----------|-------|------|---|---|---------|
| rs3131972| 0.2032| 1.047|  A|  G|  0.16055|
| rs3131969|0.08597| 1.067|  A|  G| 0.133028|
| rs3131967|0.06683| 1.077|  T|  C|        .|
| rs1048488| 0.2808|0.9617|  T|  C| 0.836449|
|rs12562034| 0.8489|0.9931|  A|  G|0.0925926|

Allele frequency of a1 (a1freq) is optional (if provided in the GWAS data), but it is necessary if you would like to include strand-ambiguous SNPs using the allele frequency information to inform strand alignment.

You can change your GWAS to the same format, or use optional parameter flags to let the script know about the format you are using. Header names are optional. More details below.

### Genotype 
.gen or .vcf file formats can be used for the genotype data input.

#### .gen file
A full description can be found on [www.shapeit.fr](http://www.shapeit.fr/pages/m02_formats/gensample.html). A .gen file is a space-delimited file with each line corresponding to a single SNP. The first five columns are:

Chromosome number [integer]

SNP ID [string]

SNP physical position (bp) [integer]

First allele [string]

Second allele [string]

Starting from the sixth column every three columns indicate the genotype (or posterior genotype probabilities in the case of using imputed genotype data) of one subject sample. 
The first of the three columns indicates the likelihood that the sample carries homozygous first allele, 
the middle column indicates the likelihood that the sample carries heterozygous alleles, 
and the last column indicates the likelihood that the sample carries homozygous second allele.

#### .vcf file 
This is a default format for the genotype data returned from [Sanger Institute](https://imputation.sanger.ac.uk/). 
Details about the format can be found [here](http://samtools.github.io/hts-specs/VCFv4.2.pdf). 
PRSoS uses posterior genotype probabilities (see FORMAT = GP) as genotype input.

### Output file
The output format of PRS results and SNP logs is comma-delimited. 

#### .score

A .score file contains the sample ID, the PRS and the SNP count for each p-value threshold applied.

#### .snplog

A .snplog file contains the SNP list and their effect allele for each PRS p-value threshold applied. It also contains the list of shared SNPs between the GWAS and genotype data that are discarded because of strand-ambiguity, mismatching reference alleles, or SNP duplication.


## Running command-line script PRS_run.py
### Spark-submit command

Use ```spark-submit``` to run the PRS calculations script. You can add other Spark parameters before the script name to control how Spark operates (e.g., --master local[K] to use K number of local worker threads). More Spark-submit options are found [here](http://spark.apache.org/docs/latest/submitting-applications.html).

```
spark-submit PRS_run.py
```

Followed by three positional parameters (mandatory):

```
  GENO                  Name of the genotype files, can be a name or path, or name patterns with '*'
  GWAS                  Name of the GWAS file, can be a name or path.
  OUTPUT                The path and name for the output file.
```

Followed by some optional parameters (full list further below). By default, the pipeline assumes the following: 

* To specify the format of genotype file:
```
  --filetype VCF  
```

The type of genotype file used as input, choose between VCF and GEN. Default is VCF.

* To specify the p-value thresholds of the score:
```
  --thresholds 0.5 0.2 0.1 0.05 0.01 0.001 0.0001
```

Enter one or more float numbers separated by space. Default is 0.5 0.2 0.1 0.05 0.01 0.001 0.0001.

Alternatively you can specify a sequence of thresholds to use:

```
  --threshold_seq 0.1 0.5 0.01
```

After the flag, the first number is the starting point of the sequence, the second is the end point of the sequence, the third number denotes the step size. 
The above example would yield the sequence 0.1, 0.11 ,0.12, ... 0.49, 0.5. Note that the interval is inclusive of the endpoints.

#### Examples
To calculate PRS using the provided test sample files and generate the SNP log output:

```
spark-submit PRS_run.py test_sample.gen gwas.clump.txt test_output --sample test_sample.sample --filetype GEN --snp_log
```

To calculate PRS from a series of .vcf files, while checking the allele alignment between the genotype and the GWAS, and log transform risk effect, using p-value thresholds of 0.2, 0.1, 0.05:

```
spark-submit PRS_run.py "VCF_number*.vcf" gwas.clump.txt output.csv --sample samplefile.csv --sample_id 0 --log_or --thresholds 0.2 0.1 0.05
```

To calculate PRS from a series of .gen files, without checking allele alignments, using a GWAS with no header, and not transform the risk effect, using p-value thresholds of 0.2, 0.1, 0.05:

```
spark-submit PRS_run.py "GEN_number*.gen" gwas.clump.txt output.csv --filetype GEN --sample samplefile.csv --sample_id 0 --no_check_ref --gwas_no_header --thresholds 0.2 0.1 0.05
```

### Parameters

A description of the parameters for the script can be obtained by typing: 
```
python PRS_run.py --help
```
Command-line type flags are used to specify how the scores are calculated. 

#### Full list of parameters when type `python PRS_run.py --help`
```
Positional arguments:
  GENO                  Name of the genotype files, can be a name or path, or
                        name patterns with wildcard character.
  GWAS                  Name of the GWAS file, can be a name or path.
  OUTPUT                The path and name stem for the output files. One name
                        will be used for the score output, the snp log 
                        (optional), and the regression output. This is similar 
                        to the --out flag in PLINK.

Optional arguments:
  -h, --help            Show this help message and exit.
  -v, --version         Show program's version number and exit.
  --app_name APP_NAME   Give your spark application a name. Default is PRS.
  --filetype {GEN,VCF}  The type of genotype file used as input, choose
                        between VCF and GEN. Default is VCF.
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
  --gwas_delim "GWAS_DELIM"
                        Delimiter of the GWAS file. Default is tab-delimited. 
                        Set quotation marks around the delimiter when applied.
  --gwas_no_header      Adding this parameter signals that there is no header 
                        in the GWAS data input. The default is to assume that 
                        GWAS has column names.
  --sample SAMPLE_FILE
                        Path and name of the file that contain the sample
                        labels. It is assumed that the sample labels are
                        already in the same order as in the genotype file.
  --sample_delim "SAMPLE_DELIM"
                        Delimiter of the sample file. Default is comma. Set 
                        quotation marks around the delimiter when applied.
  --sample_id ID_COLUMN1 ID_COLUMN2 ID_COLUMN3 ...
                        Specify which columns in the sample file are used as
                        labels. Can use one integer to specify one column, or
                        multiple integers to specify multiple columns, with 
                        first column being 0. Default is 0.
  --sample_skip_header SAMPLE_HEADER_LINES
                        Specify how many header lines to ignore in the sample 
                        file. Default is 2, which assumes that the sample labels
                        start on the third line.
  --pheno PHENO_FILE
                        Specify the path to the data file for the phenotype.
                        It is assumed that the phenotype data is organized in
                        the same order as the samples in the genotype file.
  --pheno_delim "PHENO_DELIM"
                        Specify the delimiter for the phenotype data file.
                        Default is comma. Set quotation marks around the 
                        delimiter when applied.
  --pheno_no_header     Specify whether the phenotype has a header row.
  --pheno_columns PHENO_COLUMN1 PHENO_COLUMN2 PHENO_COLUMN3 ...
                        Column number(s) in the phenotype file that contain the 
                        phenotype data. Multiple column numbers can be specified
                        to conduct regression with multiple phenotypes, with 
                        first column being 0. Default is 0.
  --covar_columns COVAR_COLUMN1 COVAR_COLUMN2 COVAR_COLUMN3 ...
                        Column number(s) in the phenotype file that contain the 
                        covariate data. Multiple column numbers can be specified
                        to conduct regression with multiple covariates, with 
                        first column being 0. No column number is set as 
                        default.
  --thresholds THRESHOLD1 THRESHOLD2 THRESHOLD3 ...
                        The p-value thresholds that filters which SNPs are
                        used from the GWAS. Specifying the p-values simply by
                        input one after another. Default is 0.5, 0.2, 0.1,
                        0.05, 0.01, 0.001, 0.0001. Note that the interval is 
                        inclusive of the endpoints.
  --threshold_seq LOWERBOUND UPPERBOUND STEPSIZE
                        Defines a sequence that contains all the p-value
                        thresholds that filters which SNPs are used from the
                        GWAS. Input is three numbers separated by space: lower
                        bound, upper bound, step size. Default is None (i.e., 
                        not used). Defining a sequence automatically overwrites 
                        the threshold list defined under --thresholds.
  --log_or              Adding this parameter tells the script to log (natural 
                        base) the effect sizes provided in the GWAS 
                        summary-level data. For example, this would be applied
                        to odds ratios to get log odds or the beta values of 
                        logistic regression.
  --no_a1f              The pipeline calculates the allele frequencies in the 
                        genotype population by default, which is used to help 
                        retain as many ambiguous SNPs as possible by comparing 
                        the allele frequencies in the GWAS to make the best 
                        estimate of the strand alignment. Use this flag to 
                        disable this feature and all ambiguous SNPs that would 
                        have been used for PRS caculation will be discarded.
  --no_check_ref        Adding this option disables checking reference alleles 
                        between GENO and GWAS when determining genotype calls. 
                        When this is used, first allele column in GENO and GWAS 
                        will be used for scoring.
  --check_dup           Add this flag if you want to check for and discard
                        SNPs that are duplicated, which will take extra time.
                        By default, the script will assume there is no
                        duplicate SNPs.
  --snp_log             Add this flag to record the SNPs that are used at each 
                        threshold. It will also report whether the A1 or A2 
                        allele in the genotype data was used as reference for 
                        the risk effect. Any SNPs that meet the p-value 
                        threshold criteria but has allele names that do not 
                        match the allele names in the GWAS description are 
                        indicated in the 'discard' column. This record will be 
                        saved to a file with the name specified in the OUTPUT 
                        flag, with .snplog as file extension.

``` 

## Notes

+ Allele names in the allele columns of the GENO and GWAS files should match. If they do not match, they will be discarded (see the SNP log "--snp_log").
+ Allele names should be written using nucleotide initials with capital letters.
