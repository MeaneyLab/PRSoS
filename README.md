# PRS-on-Spark
PRS-on-Spark (PRSoS) generates polygenic risk scores (PRS) for large genotype data, including imputed genotype posterior probabilites. 
It can use multiple cores to increase processing efficiency (i.e., reduce processing time). 
PRSoS is compatible with Linux, Mac OS, and Windows. It runs using [Apache Spark](https://spark.apache.org) and [Python](https://www.python.org/download/releases/2.7.2/).

## Contact Information
Lawrence M. Chen: lawrence.m.chen@mail.mcgill.ca


## Installation

To clone the repository, use:
```
git clone https://github.com/MeaneyLab/PRSoS.git
```

## Software requirements and detailed installation instructions

There are many ways to configure Spark (local execution, standalone cluster, or Yarn modes), each different on Windows, MacOS, and Linux. 
For simplicity it is now easier to use a quick [PyPi (pip)](https://pypi.org/project/pip/) that is common on any platform where Python is installed, and use its bundled version of spark-submit for client runs on the local file system. 
If you favour configuring Spark to run as a service with [Hadoop Distributed File System (HDFS)](https://searchdatamanagement.techtarget.com/definition/Hadoop-Distributed-File-System-HDFS) and workers / executors, 
there are many resources for this available online.

PRSoS scripts require the following to run:
+ [Spark-2.0.0 + (Apache)](https://spark.apache.org/downloads.html)
+ Python 2.7 (**not Python 3.0,** and preferably in [Anaconda](https://www.anaconda.com/what-is-anaconda/), which is a data science distribution of Python)
+ Java 8 (preferably [Oracle](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html))

The required Python modules are:
+ matplotlib
+ pandas
+ numpy
+ statsmodels
+ pyspark

Below are platform-specific installation guides:
+ [Windows](#win10)
+ [Mac OS](#mac)
+ [Linux](#linux)

### <a name="win10"></a>Windows (Win10)

1. Download and install Oracle JDK8 from [here](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html).
2. Download and install Anaconda2 v5.10 Python from [here](https://repo.anaconda.com/archive/Anaconda2-5.1.0-Windows-x86_64.exe). Use default options.
3. To install Spark, open an elevated command prompt and issue a pip or conda install:
    ```
    pip install pyspark
    ```
    or
    ```
    conda install -c conda-forge pyspark=2.2.1
    ```
    
4. Once installed, try issuing the following command. 
   ```
   spark-submit --version
   ```
   If you get nothing, reopen your elevated command prompt and try again.

5. Then git clone (or [download zipfile](https://github.com/MeaneyLab/PRSoS/archive/master.zip) and unzip), access the folder, and install required modules:
   ```
   git clone https://github.com/MeaneyLab/PRSoS.git
   cd PRSoS*
   pip install –r requirements.txt
   ```
   
6. Finally, perform this PRSoS test run to check if the installation was successful:
   ```
   spark-submit --master local[*] PRSoS.py examples/example.vcf examples/gwasfile.txt test_output
   ```

Note 1: While it is not necessary to use Hadoop and its HDFS filesystem, running Spark in Windows will give some messages such as “**Failed to locate the winutils binary**” when it fails to find %HADOOP_HOME and %HADOOP_HOME%\bin in the Windows environment. 
[Here](http://alvincjin.blogspot.ca/2016/08/setup-hadoophome-in-windows.html) is how to resolve this:
 + Download [winutils.exe] (https://github.com/steveloughran/winutils) and store it in a Hadoop bin directory (e.g., C:\Users\hadoop\bin).
 + Then open an elevated command prompt and run the following: 

    ```
    setx HADOOP_HOME "C:\Users\hadoop"
    setx PATH "%HADOOP_HOME%\bin;%PATH%"
    ```

Note 2: Hadoop and Spark versions should be compatible, e.g., use hadoop-2.7.1 for Spark 2 ([here is the direct link to winutils.exe binary](https://github.com/steveloughran/winutils/blob/master/hadoop-2.7.1/bin/winutils.exe)).

Alternative instructions for installing Spark can be found [here](http://www.informit.com/articles/article.aspx?p=2755929&seqNum=3).

### <a name="mac"></a>Mac OS

[Homebrew](https://brew.sh) package manager in MacOS can be used to install everything necessary for PRSoS. 
For this reason, Spark is available both via "brew" as well as PythonPI's "pip install". 
After installing Homebrew, and brew installing Java8 and Anaconda2, we will use the pip install method to get Spark.

If you already have Homebrew, remember to run "brew doctor" before anything to fix possible issues. 
Copying this into a MacOS Terminal will install Homebrew:
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Then you may run [this PRSoS preparation script](http://www.meaney.lab.mcgill.ca/PRSoS-prep-bash.sh) with step by step installation for PRSoS. 
This installation script is compatible with MacOS and Linux. 
Execute the installation script from the Terminal using:
```
bash PRSoS-prep-bash.sh
```

<details><summary>The installation script essentially runs the following setup. Click here to expand/collapse the code.</summary>
<pre><code>
## Pre-requisite installations and preparations
cd ~/
xcode-select --install # this installs Xcode if it's not already installed
brew tap caskroom/cask
brew cask install java8 # this installs Java 8
brew install python@2 # this installs Python 2
echo `export PATH="/usr/local/opt/python2/libexec/bin:$PATH"` >> ~/.bash_profile
source ~/.bash_profile
#pip2 install --upgrade pip setuptools --user
python --version # this should indicate that you have Python 2
pip install pyspark # you can append "--user" to install psypark to user's home folder
spark-submit --version # this tests that "spark-submit" command is working; it should display Spark profile
## We will now install PRSoS
wget https://github.com/MeaneyLab/PRSoS/archive/master.zip # "git clone https://github.com/MeaneyLab/PRSoS.git" also works if you have git installed
unzip master.zip # this is not necessary if you used the "git" command
cd ~/PRSoS*
pip install -r requirements.txt # you can append "--user" to install to user’s home folder
## We will now submit a test run on all available cores
spark-submit --master local[*] PRSoS.py examples/example.vcf examples/gwasfile.txt test_output

</code></pre>
</details>

### <a name="linux"></a>Linux

You may run [this PRSoS preparation script](http://www.meaney.lab.mcgill.ca/PRSoS-prep-bash.sh) with step by step installation for PRSoS. 
This installation script is compatible with MacOS and Linux. 
Execute the installation script from the Terminal using:
```
bash PRSoS-prep-bash.sh
```

For the classical method of installing Spark locally using pre-built binaries, or how to install with "conda" look [here](http://sigdelta.com/blog/how-to-install-pyspark-locally/).

<details><summary>The joy of Spark installs in Linux is the automation. A quick install can be done by running the following in Terminal. Click here to expand/collapse the code.</summary>
<pre><code>
# ** Anaconda2 install
cd ~/
wget https://repo.anaconda.com/archive/Anaconda2-5.1.0-Linux-x86_64.sh
bash Anaconda2-5.1.0-Linux-x86_64.sh -b -p $HOME/anaconda2
echo -e '##**Added for PRSoS anaconda python**##\nexport PATH="$HOME/anaconda2/bin:$PATH" # comment out beginning of this line with "#" to disable Anaconda and revert to old python' >> ~/.bashrc
source ~/.bashrc
python --version #This should now say Anaconda python
pip --version #needed to install pyspark and dependency modules
# ** Spark install
pip install pyspark --user # you can omit "--user" if you have root / sudo (i.e., system administrator) privileges
spark-submit --version # this tests that "spark-submit" command is working; it should display Spark profile with Scala and Java versions
# **git clone PRSoS, make sure you have installed git !
wget https://github.com/MeaneyLab/PRSoS/archive/master.zip # "git clone https://github.com/MeaneyLab/PRSoS.git" also works if you have git installed
unzip master.zip # this is not necessary if you used the "git" command
cd PRSoS*
pip install -r requirements.txt --user # omit --user if you have root / sudo priveleges
cd ~/PRSoS* ## now do a PRSoS test run using spark-submit
sleep 5
spark-submit --master local[*] PRSoS.py examples/example.vcf examples/gwasfile.txt test_output

</code></pre>
</details>
<br>

At this point you should see a successful PRSoS spark-submit.

<details><summary>A SparkUI service is viewable in your browser at [http://127.0.0.1:4040](http://127.0.0.1:4040/) for the duration of the run. 
Most errors at this point are due to missing java ($JAVA_HOME) or incompatible versions. 
If you have not yet installed Oracle JDK 8, then you can do this quickly by pasting the following. Click here to expand/collapse the code.</summary>
<pre><code>
java -version
cd ~/
wget --no-check-certificate -c --header "Cookie: oraclelicense=accept-securebackup-cookie" http://download.oracle.com/otn-pub/java/jdk/8u171-b11/512cd62ec5174c3487ac17c61aaa89e8/jdk-8u171-linux-x64.tar.gz
tar xvzf jdk-8u171-linux-x64.tar.gz
echo "#**Added for PRSoS Java**" >> ~/.bashrc
echo "export JAVA_HOME=~/jdk1.8.0_171" >> ~/.bashrc
echo 'export PATH="$JAVA_HOME/bin/:$PATH"' >> ~/.bashrc
source ~/.bashrc
java -version # this should return java version "1.8.0_171"

</code></pre>
</details>
<br>

For further help on installing Apache Spark in Linux as an administrator, look [here](https://www.santoshsrinivas.com/installing-apache-spark-on-ubuntu-16-04/).

## What this pipeline does
+ Match the strand alignment between genotype and GWAS data (by default), and then
+ Calculate PRS from the genotype data set (in .gen or .vcf file format), weighted by the associated effect size (e.g., log odds ratio), and subsetted by the selected associated p-value thresholds in the GWAS

## What this pipeline cannot do
+ Perform quality control of genotype data (depending on the file format, quality control can be performed using [QCTOOL](http://www.well.ox.ac.uk/~gav/qctool_v2/), [PLINK](https://www.cog-genomics.org/plink/1.9/), [VCFtools](http://vcftools.sourceforge.net), etc.)
+ Perform linkage disequilibrium clumping or pruning (linkage disequilibrium clumping or pruning can can be performed using [PLINK](https://www.cog-genomics.org/plink/1.9/). However, note that dosage data (e.g., posterior probabilities) are not compatible with PLINK file formats)

## Guideline for linkage disequilibrium clumping
PRSoS does not perform linkage disequilibrium (LD) corrections as it can be performed using various methods and may or may not be necessary for your PRS. 
[PLINK's LD clumping algorithm](https://www.cog-genomics.org/plink/1.9/postproc#clump) is currently the most popular approach. 
The following guide (compatible with MacOS and Linux Terminal.app) can be used to prepare the required GWAS file with LD clump:
```
# for the present example, change the current working directory to the PRSoS/examples/ directory
cd ~/PRSoS/examples # the path may be different depending on where you located the PRSoS directory

# make binary PLINK files from GEN files for clumping
plink --data example --keep-allele-order --make-bed --out example_gen2plink 

# run clumping
plink \
--bfile example_gen2plink \
--clump gwasfile.txt \
--clump-p1 1 --clump-p2 1 --clump-kb 500 --clump-r2 0.2 `# set your LD clumping parameters` \
--clump-snp-field "SNP" --clump-field "P" `# this is not necessary if the SNP ID and p-value column headers are already labeled "SNP" and "P", respectively; you can change these field names in the command to match your GWAS file if they are different` \
--out example_gwas_clump_output

# create GWAS file with only index SNPs from the clumping output
awk 'NR==FNR {FILE1[$3]=$0; next} ($1 in FILE1) {print $0}' `# "($1 in FILE1)" assumes the SNP ID column is first column in the GWAS file; change "$1" to the SNP ID column number if it is different` \
example_gwas_clump_output.clumped gwasfile.txt > gwasfile_example_clumped.txt
```

## Default format
### GWAS

The GWAS data input requires the following columns:

snpid = rs ID of the SNP

pval = p-value of the association in the GWAS data

or = odds ratio in the GWAS data, or it can be log odds ratio or beta (if it is the odds ratio, add `--log_or` flag to use log odds ratio)

a1 = reference/effect allele (to which the odds ratio belong)

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

### Output files
The output format of PRS results and SNP logs is comma-delimited. 

#### .score.csv

A .score.csv file contains the sample ID, the PRS and the SNP count for each p-value threshold applied.

E.g.,

|     ID1|     ID2|     SNP_count_0.01|     PRS_0.01|     SNP_count_0.05|     PRS_0.05|
|--------|--------|-------------------|-------------|-------------------|-------------|
| Family1|   Child|                 13|  -0.02797601|                101|  -0.00647242|
| Family1|  Mother|                 13|  -0.02933779|                101|  -0.00509516|
| Family1|  Father|                 13|  -0.02380963|                101|  -0.00641128|
| Family2|   Child|                 13|  -0.03056985|                101|  -0.00888138|
| Family2|  Father|                 13|  -0.03028115|                101|  -0.00732302|
| Family3|   Child|                 13|  -0.02817242|                101|  -0.00796898|

The columns are as follows:

ID1, ID2, ID..., etc. = sample IDs provided from the sample file (see --sample_id description below)

SNP_count_0.01, SNP_count_0.05, SNP_count_..., etc. = the number of SNPs that contribute to the PRS_0.01, PRS_0.05, PRS_..., etc.

PRS_0.01, PRS_0.05, PRS_..., etc. = the polygenic risk score at the p-value threshold 0.01, 0.05, ..., etc. (see --thresholds description below)

#### .snplog.csv

A .snplog.csv file contains the SNP list and their effect allele for each PRS p-value threshold applied. 
It also contains the list of shared SNPs between the GWAS and genotype data that are discarded because of strand-ambiguity, mismatching reference alleles, or SNP duplication.

E.g.,

|     PRS_0.01|     PRS_0.01_flag|     PRS_0.05|     PRS_0.05_flag|     Discard|
|-------------|------------------|-------------|------------------|------------|
|   rs10017793|                A1|   rs10017793|                A1|   rs1276324|
|   rs11157862|                A1|   rs10501467|                A2|   rs2278342|
|   rs17329328|                A2|   rs10514162|                A2|   rs9919144|
|   rs17772344|                A1|   rs11157862|                A1|    rs997850|
|    rs1994908|                A2|   rs17329328|                A2|            |

The columns are as follows:

PRS_0.01, PRS_0.05, PRS_..., etc. = list of SNPs that contribute the indicated PRS p-value thresholds (see --thresholds description below)

PRS_0.01_flag, PRS_0.05_flag, PRS_..._flag, etc. = the allele in the genotype file that was matched with the effect allele in the GWAS data and weighted by the log odds or beta

Discard = list of SNPs that are discarded because of strand-ambiguity, mismatching reference alleles, or SNP duplication

## Running command-line script PRSoS.py
### Spark-submit command

Use ```spark-submit``` to run the PRS calculations script. 
You can add other Spark parameters before the script name to control how Spark operates (e.g., --master local[K] to use K number of local worker threads). 
More Spark-submit options are found [here](http://spark.apache.org/docs/latest/submitting-applications.html).

```
spark-submit PRSoS.py
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
spark-submit PRSoS.py test_sample.gen gwas.clump.txt test_output --sample test_sample.sample --filetype GEN --snp_log
```

To calculate PRS from a series of .vcf files, while checking the allele alignment between the genotype and the GWAS, and log transform risk effect, using p-value thresholds of 0.2, 0.1, 0.05:

```
spark-submit PRSoS.py "VCF_number*.vcf" gwas.clump.txt output.csv --sample samplefile.csv --sample_id 0 --log_or --thresholds 0.2 0.1 0.05
```

To calculate PRS from a series of .gen files, without checking allele alignments, using a GWAS with no header, and not transform the risk effect, using p-value thresholds of 0.2, 0.1, 0.05:

```
spark-submit PRSoS.py "GEN_number*.gen" gwas.clump.txt output.csv --filetype GEN --sample samplefile.csv --sample_id 0 --no_check_ref --gwas_no_header --thresholds 0.2 0.1 0.05
```

### Parameters

A description of the parameters for the script can be obtained by typing: 
```
python PRSoS.py --help
```
Command-line type flags are used to specify how the scores are calculated. 

#### Full list of parameters when type `python PRSoS.py --help`
```
Positional arguments:
  GENO                  Name of the genotype files, can be a name or path, or
                        name patterns with wildcard character.
  GWAS                  Name of the GWAS file, can be a name or path.
  OUTPUT                The path and name stem for the output files. One name
                        will be used for the score output and the snp log 
                        (optional). This is similar 
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
                        Delimiter of the GWAS file. Default is tab. Set 
                        quotation marks around the delimiter when applied.
  --gwas_no_header      Adding this parameter signals that there is no header 
                        in the GWAS data input. The default is to assume that 
                        GWAS has column names.
  --sample SAMPLE_FILE
                        Path and name of the file that contain the sample
                        IDs. It is assumed that the sample IDs are
                        already in the same order as in the genotype file.
  --sample_delim "SAMPLE_DELIM"
                        Delimiter of the sample file. Default is space. Set 
                        quotation marks around the delimiter when applied.
  --sample_id ID_COLUMN1 ID_COLUMN2 ID_COLUMN3 ...
                        Specify which columns in the sample file are used as
                        IDs. Can use one integer to specify one column, or
                        multiple integers to specify multiple columns, with 
                        first column being 0. Default is 0.
                                Tip: 
                                    This can be used to add other variables from
                                    the sample file to the PRS output (e.g., 
                                    gender or family membership)
  --sample_skip_header SAMPLE_HEADER_LINES
                        Specify how many header lines to ignore in the sample 
                        file. Default is 2, which assumes that the sample IDs
                        start on the third line.
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
                        flag, with .snplog.csv as file extension.

``` 

## Notes

+ Allele names in the allele columns of the GENO and GWAS files should match. If they do not match, they will be discarded (see the SNP log "--snp_log").
+ Allele names should be written using nucleotide initials with capital letters.
