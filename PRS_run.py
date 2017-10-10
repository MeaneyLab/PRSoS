 # -*- coding: utf-8 -*-
from __future__ import division


from operator import add
from math import log
import csv
import pickle
import sys
from collections import Counter
import re
import glob, os

import ntpath
import functools
import itertools

import time
import argparse


# Packages for the regression
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plottings import *
import statsmodels.formula.api as smf
import logging

# Takes a line in the genotype file and return the frequency of A1 allele
def getA1f(geno):
    AA=geno[0::3]
    AB=geno[1::3]
    AA2=[x*2 for x in AA]
    A1count=map(add, AA2, AB)
    A1F=sum(A1count)/(float(len(AA2))*2)
    return A1F

def getCall(line):
    AA=line[0::3]
    AB=line[1::3]
    BB=line[2::3]
    calls=[x+y+z for x,y,z in zip(AA,AB,BB)]
    calls=[int(round(x)) for x in calls]
    return calls

# Takes an array of triplets and convert to number of A1 alleles
def makeGenotype(line):
    AA=line[0::3]
    AB= line[1::3]
    AA2=[x*2 for x in AA]
    genotype=[AA2[i]+AB[i] for i in range(len(AA2))]
    return  genotype

# Takes an array of triplets and convert to number of A1 alleles, while checking for reference alignments
def makeGenotypeCheckRef(line, checkMap, toDF=False):
    rsid=line[0]
    gen=line[1]
    if rsid in checkMap:
        if checkMap[rsid]=="A1":    # applies 2 1 0 weights on the genotype triplet codes when the SNP is tagged with "A1" (i.e., A1 in genotype data is same or complementary to A1 in GWAS)
            AA=gen[0::3]
            AB=gen[1::3]
            AA2=[x*2 for x in AA]
            genotype=list(map(add, AA2, AB))

        elif checkMap[rsid]=="A2":    # applies 0 1 2 weights on the genotype triplet codes when the SNP is tagged with "A2" (i.e., A2 in genotype data is same or complementary to A1 in GWAS)
            AA=gen[2::3]
            AB=gen[1::3]
            AA2=[x*2 for x in AA]
            genotype=list(map(add, AA2, AB))


    #else:
        # fail-safe, in case a SNP does not have a flag but still exist in the genotype file
        #print("SNP {} was not accounted for in the alignment checking step, discarding this SNP".format(rsid))
        #genotype=[0.0]*(len(gen)/3)
    if genotype:
        if toDF:  # future work, in case want to use SQL methods on Spark DataFrame
            return [rsid]+genotype
        else:
            return (rsid, genotype)


def filterGWASByP(GWASRdd, pcolumn,  pHigh, oddscolumn,idcolumn, pLow=0, logOdds=False):
    GWAS_Pfiltered=GWASRdd.filter(lambda line: (float(eval(line[pcolumn]))<=pHigh) and (float(line[pcolumn])>=pLow))
    if logOdds:

        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],log(float(line[oddscolumn]))))
    else:

        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],  float(line[oddscolumn])))
    GWASoddspair=GWAS_Odds.collectAsMap()
    return GWASoddspair


def filterGWASByP_DF(GWASdf, pcolumn,  pHigh, oddscolumn,idcolumn, pLow=0, logOdds=False):
    GWAS_Pfiltered=GWASdf.rdd.filter(lambda line: (float(line[pcolumn])<=pHigh) and (float(line[pcolumn])>=pLow))
    if logOdds:
        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],log(float(line[oddscolumn]))))

    else:

        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],  float(line[oddscolumn])))  # make tuples of the format (snpid, effect size)
    GWASoddspair=GWAS_Odds.collectAsMap()   # make a python dictionary of the format { snpid: effect_size }
    return GWASoddspair


# Checking for reference allele alignment (i.e., is the GWAS A1 the same as A1 in the genotype)
def checkAlignmentDF(dataframe, bpMap):
    snpid=dataframe[0]
    genoA1=dataframe[1].strip()
    genoA2=dataframe[2].strip()
    genoA1F=dataframe[3]
    # no number 4 because the 5th column has the snpid from the GWAS
    gwasA1=dataframe[5].strip()
    gwasA2=dataframe[6].strip()
    gwasA1Fstring=dataframe[7].strip()
    if re.compile(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$').match(gwasA1Fstring):
        gwasA1F=float(gwasA1Fstring)    # keep gwasA1Fstring as gwasA1F if it is a number
    else:
        gwasA1F="NO"    # otherwise gwasA1F is indicated as "NO" because it is not provided

    if genoA1 in bpMap:
        # discard if A1 and A2 in GENO or GWAS are the same
        if genoA1==genoA2 or gwasA1==gwasA2:
            flag='discard'
        # for ambiguous SNPs
        elif genoA2==bpMap[genoA1] and (genoA1==bpMap[gwasA1] or genoA1==bpMap[gwasA2]) and (genoA2==bpMap[gwasA2] or genoA2==bpMap[gwasA1]):
            if gwasA1F=="NO":
                flag="discard"    # discard ambiguous SNP is A1 frequency is not provided
            else:
                gwasA1F=float(gwasA1F)
                genoA1Fdiff=genoA1F*10-5
                gwasA1Fdiff=float(gwasA1F)*10-5
                # if the allele frequency in the GWAS or genotype is between 0.4 and 0.6 (or have an absolute A1Fdiff's < 1), it will be discarded because we cannot be confident about the strand alignment due to genotype differences in the cohorts
                if abs(genoA1Fdiff)<1 or abs(gwasA1Fdiff)<1:
                    flag="discard"
                else:
                    if genoA1Fdiff * gwasA1Fdiff>0:    # if A1Fdiff's are both positives or both negatives, it implies the genotypes in GWAS and genotype data are aligned on the same strand
                        flag="A1"
                    else:
                        flag="A2"
        # for non-ambiguous SNPs
        elif (genoA2==gwasA1 and genoA1==gwasA2) or (genoA2==bpMap[gwasA1] and genoA1==bpMap[gwasA2]):
            flag="A2"

        elif (genoA1==gwasA1 and genoA2==gwasA2) or (genoA1==bpMap[gwasA1] and genoA2==bpMap[gwasA2]):
            flag="A1"

        else:
            flag="discard"    # discard if one allele name does not match (for example, if geno has (A,G) and GWAS has (G,T))
    else :
        flag="discard"    # discard everything else
        #logger.warn("Invalid Genotypes for SNP {}; genotype Alleles:{},{}; GWAS alleles: {},{}".format(snpid, genoA1, genoA2, gwasA1, gwasA2))

    return (snpid,flag)

# Check reference allele alignment without using A1F in the GWAS. Ambiguous SNPS are automatically discarded
def checkAlignmentDFnoMAF(dataframe, bpMap):
    snpid=dataframe[0]
    genoA1=dataframe[1].strip()
    genoA2=dataframe[2].strip()
    gwasA1=dataframe[4].strip()
    gwasA2=dataframe[5].strip()

    if genoA1 in bpMap:
        if genoA2==bpMap[genoA1] or genoA1==genoA2:  # discard ambiguous case and the case where A1 = A2 in the genotype data
            flag="discard"
        elif (genoA2==gwasA1 and genoA1==gwasA2) or (genoA2==bpMap[gwasA1] and genoA1==bpMap[gwasA2]):    # for non-ambiguous SNPs
            flag="A2"
        elif (genoA1==gwasA1 and genoA2==gwasA2) or (genoA1==bpMap[gwasA1] and genoA2==bpMap[gwasA2]):    # for non-ambiguous SNPs
            flag="A1"
        else:
            flag='discard'

    else:
        flag="discard"
        #print("Invalid Genotypes for SNP {}; Genotype alleles:{},{}; GWAS alleles: {},{}".format(snpid, genoA1, genoA2, gwasA1, gwasA2))

    return (snpid,flag)


def getSampleNames(sampleFileName, sampleDelim, sampleIDCol, skip=0):
    labels=[]
    with open(sampleFileName, "r") as f:
        subjList=[item.split(sampleDelim) for item in f.read().splitlines()]
        counter=1
        for i in sampleIDCol:
            subjNames=[x[i] for x in subjList[skip::]]
            subjNames=[name.strip('"') for name in subjNames]
            column=["Label-"+str(counter)]+subjNames
            labels.append(column)
            counter+=1
    return labels

# Remove duplicates from flag list
def rmDup(checkList):
    checkCount=Counter([x[0] for x in checkList])
    nodupMap={snp:flag for snp, flag in checkList if checkCount[snp]==1}
    return nodupMap

# Output each PRS for each sample in the form of [sample, *scores]
def writePRS(prsResults, outputFile, logger, samplenames=None, dialect=None):
    scorefile=outputFile+".score"
    onescore=list(prsResults.values())[0][1]
    samplesize=len(onescore)
    if not samplenames:
        logger.warn("No sample names provided, generating sample names")
        samplenames=[["Label"]+["Sample"+str(i+1) for i in range(samplesize)]]
    labelnumber=len(samplenames[0])-1
    logger.info("Collected {} sample labels".format(labelnumber))
    if labelnumber==samplesize:
        outputdata=samplenames
        pvaluelist=sorted(list(prsResults.keys()))
        for pvalue in pvaluelist:
            snpcounts=[int(x) for x in prsResults[pvalue][0]]
            outputdata.append(["SNP_count_{}".format(pvalue)]+snpcounts)
            outputdata.append(["PRS_{}".format(pvalue)]+prsResults[pvalue][1])

        try:
            with open(scorefile, "w") as f:
                csvwriter=csv.writer(f, dialect=dialect)
                for row in zip(*outputdata):
                    csvwriter.writerow(row)
                logger.info("Successfully wrote scores to "+ scorefile)
        except:
            e = sys.exc_info()[0]
            logger.warn( "Error: %s" % e )
            logger.warn("Data output was unsuccessful.")
            logger.warn("All is not lost, final results saved as binary format in file 'PRSOutput.pk'")
            with open(os.path.dirname(scorefile)+"/PRSOutput.pk", "wb") as f:
                pickle.dump(outputdata, f)
    else:
        logger.warn("Unequal number of labels extracted from sample sheet and number of samples detected in the genotype data, saving results to PRSOutput.pk")
        outputdata=False
        with open(os.path.join(os.path.dirname(outputFile),"PRSOutput.pk"), "wb") as f:
            pickle.dump({"Samples": samplenames, "results": prsResults}, f)

    return outputdata

# Output SNP log
def writeSNPlog(snpidmap, outputFile, logger, flagMap=None, dialect=None):
    snplogfile=outputFile+".snplog"
    outputdata=[]
    maxT=max(snpidmap.keys())
    maxLen=len(snpidmap[maxT])
    for pvalue in sorted(list(snpidmap.keys())):

        sortedlist=sorted(snpidmap[pvalue])
        outputdata.append(["PRS_"+str(pvalue)]+sortedlist+[""]*(maxLen-len(sortedlist)))

        # get the flag for each SNP in the SNP log
        if flagMap:
            flaglist=["PRS_"+str(pvalue)+"_flag"]+[flagMap[snpid] for snpid in sortedlist]+[""]*(maxLen-len(sortedlist))
            outputdata.append(flaglist)
    if flagMap:
        discardlist=[snp for snp in flagMap.keys() if flagMap[snp]=="discard"]
        outputdata.append(["Discard"]+discardlist+[""]*(maxLen-len(discardlist)))

    try:
        with open(snplogfile, "w") as f:
            csvwriter=csv.writer(f, dialect=dialect)
            for row in zip(*outputdata):
                csvwriter.writerow(row)
            logger.info("Successfully wrote log to "+ snplogfile)
    except:
        e = sys.exc_info()[0]
        logger.warn( "Error: %s" % e )
        logger.warn("SNP log output was unsuccessful.")
        logger.warn("All is not lost, logs are saved as binary format in file 'SNPlog.pk'")
        with open(os.path.join(os.path.dirname(snplogfile),"SNPlog.pk"), "wb") as f:
            pickle.dump(outputdata, f)
    return outputdata


def regression(scoreMap,phenoFile, phenoDelim, phenoColumns, phenoNoHeader, covarColumns, outputName, logger):
    samplesize=len(scoreMap[list(scoreMap.keys())[0]][1])
    outputFile="{}.regression".format(outputName)
    if phenoNoHeader:
        header=None
    else:
        header=0
    pvalueList=sorted(list(scoreMap.keys()))
    prsData=pd.DataFrame({pvalue:scoreMap[pvalue][1] for pvalue in pvalueList})
    prsData.columns=["prs{}".format(re.sub("\.","_",str(x))) for x in prsData.columns.tolist()]
    thresholdNames=prsData.columns.tolist()

    phenodata=pd.read_table(phenoFile, header=header, sep=phenoDelim)
    phenodata.columns=[re.sub("\.","_",str(x)) for x in phenodata.columns]


    assert samplesize==phenodata.shape[0], "Unequal sample size in pheno file and PRS data"

    covariateNames=phenodata.columns[covarColumns].tolist()
    covar=phenodata[covariateNames]
    phenotypes=[]
    thresholds=pvalueList
    r2All=[]
    pAll=[]
    with open(outputFile, "w") as f:
        f.write("\n")

    for columnNumber in phenoColumns:

        pheno=phenodata.iloc[:,columnNumber]
        phenoName=phenodata.columns[columnNumber]
        logger.info("Regression with phenotype {}".format(phenoName))
        phenotypes.append(phenoName)
        regressdata=pd.concat([pheno, prsData, covar], axis=1)
        regressdataClean=regressdata.dropna(axis=0)
        logger.info("After removing rows with missing data, {} sample removed, {} samples remain".format(regressdata.shape[0]-regressdataClean.shape[0], regressdataClean.shape[0]))

        plist=[]
        r2list=[]

        for pvalue in thresholdNames:
            formula="{} ~ {}+1".format(phenoName, "+".join([pvalue]+covariateNames))
            lm = smf.ols(formula, data=regressdataClean).fit()
            plist.append(lm.pvalues[1])
            r2list.append(lm.rsquared)

            summary=lm.summary2()
            oldindex=summary.tables[1].index.tolist()
            oldindex[1]=re.sub("_", ".", oldindex[1])
            summary.tables[1].index=oldindex
            summary.title="{} : {} + {}".format(summary.title, phenoName, oldindex[1])
            with open(outputFile, "a") as f:
              f.write("\n")
              f.write(summary.as_text())
              f.write("\n")
            logger.info("Regression finished using {}. Summary written to {}".format(re.sub("_", ".", pvalue),outputFile))
        pAll.append(plist)
        r2All.append(r2list)

    logger.info("All regression(s) finished")
    return phenotypes, thresholds, r2All, pAll


if __name__=="__main__":
    import pyspark
    from pyspark.sql import SparkSession
    from pyspark.sql import Row
    from pyspark.sql.functions import udf
    from pyspark.sql.types import *
    from pyspark import SparkConf, SparkContext
    # ATTN: python index starts at 0, so if you want to specify the second column, use 1
    # Define column number for contents in GWAS
    parser = argparse.ArgumentParser(description='PARAMETERS', version='1.7')
    # Mandatory positional arguments
    parser.add_argument("GENO", action="store", help="Name of the genotype files, can be a name or path, or name patterns with wildcard character.")
    parser.add_argument("GWAS", action="store", help="Name of the GWAS file, can be a name or path.")
    parser.add_argument("OUTPUT", action="store", help="The path and name stem for the output files. One name will be used for the score output, the snp log (optional), and the regression output. This is similar to the --out flag in PLINK.")

    # Optional arguments
    parser.add_argument("--gwas_id", action="store", default=0, dest="gwas_id",type=int, help="Column number in your GWAS that contains SNP ID, with first column being 0. Default is 0.")
    parser.add_argument("--gwas_p", action="store", default=1, dest="gwas_p", type=int, help="Column number in your GWAS that contains p-value, with first column being 0. Default is 1.")
    parser.add_argument("--gwas_or", action="store", default=2, dest="gwas_or", type=int, help="Column number in your GWAS that contains odds-ratio/beta, with first column being 0. Default is 2.")
    parser.add_argument("--gwas_a1", action="store", default=3, dest="gwas_a1", type=int, help="Column number in your GWAS that contains allele A1, with first column being 0. Default is 3.")
    parser.add_argument("--gwas_a2", action="store", default=4, dest="gwas_a2", type=int, help="Column number in your GWAS that contains allele A2, with first column being 0. Default is 4.")
    parser.add_argument("--gwas_a1f", action="store", default=5, dest="gwas_a1f", type=int, help="Column number in your GWAS that contains frequency of A1, with first column being 0. Default is 5.")
    parser.add_argument("--filetype", action="store",default="VCF", dest="filetype", help="The type of genotype file used as input, choose between VCF and GEN. Default is VCF.", choices=set(["VCF", "GEN"]))

    parser.add_argument("--thresholds", action="store", default=[0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001], dest="thresholds", help="The p-value thresholds that controls which SNPs are used from the GWAS. Specifying the p-values simply by input one after another. Default is 0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001.", nargs="+", type=float)

    parser.add_argument("--threshold_seq", action="store", default=None, dest="threshold_seq", help="As an alternative to --thresholds, use this flag to define a sequence that contains all the p-value thresholds that filters which SNPs are used from the GWAS. Input is three numbers separated by space: lower bound, upper bound, step size. Default is None (i.e., not used). This flag automatically overwrites the threshold list defined under --thresholds.", nargs="+", type=float)

    parser.add_argument("--GWAS_delim", action="store", default="\t", dest="GWAS_delim", help="Delimiter of the GWAS file. Default delimiter is tab. Set quotation marks around the delimiter when applied.")

    parser.add_argument("--GWAS_no_header", action="store_false", default=True, dest="GWAS_header", help="Add this flag to signal no header in the GWAS data input. The default is to assume that GWAS has column names.")

    parser.add_argument("--log_or", action="store_true", default=False, dest="log_or", help="Add this flag to log (natural base) the effect sizes provided in the GWAS summary-level data. For example, this would be applied to odds ratios to get log odds or the beta values of logistic regression.")

    parser.add_argument("--no_check_ref", action="store_false", default=True, dest="check_ref", help="Adding this option disables checking reference alleles between GENO and GWAS when determining genotype calls. When this is used, first allele column in GENO and GWAS will be used for scoring.")

    parser.add_argument("--app_name", action="store", default="PRS", dest="app_name", help="Give your spark application a name. Default is PRS.")

    parser.add_argument("--sample_file", action="store", dest="sample_file", default="NOSAMPLE",help="Path and name of the file that contain the sample labels. It is assumed that the sample labels are already in the same order as in the genotype file.")

    parser.add_argument("--sample_file_delim", action="store", default=" ", dest="sample_delim", help="Delimiter of the sample file. Default delimiter is space. Set quotation marks around the delimiter when applied.")

    parser.add_argument("--sample_file_ID", action="store", default=[0], type=int, nargs="+", dest="sample_file_ID", help="Specify which columns in the sample file are used as labels. Can use one integer to specify one column or multiple integers to specify multiple columns, with first column being 0. Default is 0.")

    parser.add_argument("--sample_file_skip", action="store", default=2, type=int, dest="sample_skip", help="Specify how many lines to skip in the sample file, i.e., which row do the labels start. Default is 2, which assumes that the labels start on the third line")

    parser.add_argument("--no_maf", action="store_false", default=True, dest="use_maf", help="The pipeline calculates the allele frequencies in the genotype population by default, which is used to help retain as many ambiguous SNPs as possible by comparing the allele frequencies in the GWAS to make the best estimate of the strand alignment. Use this flag to disable this feature and all ambiguous SNPs that would have been used for PRS caculation will be discarded.")

    parser.add_argument("--snp_log", action="store_true", default=False, dest="snp_log", help="Add this flag to record the SNPs that are used at each threshold. It will also report whether the A1 or A2 allele in the genotype data was used as reference for the risk effect. Any SNPs that meet the p-value threshold criteria but has allele names that do not match the allele names in the GWAS description are indicated in the 'discard' column. This record will be saved to a file with the name specified in the OUTPUT flag, with .snplog as file extension.")

    parser.add_argument("--check_dup", action="store_true", default=False, dest="checkdup", help="Add this flag to discard SNPs that are duplicated, which will take extra time. By default, the script will assume there are no duplicated SNPs.")


    parser.add_argument("--pheno_file", action="store", default=None, dest="pheno_file", help="Specify the path to the data file for the phenotype. It is assumed that the phenotype data is organized in the same order as the samples in the genotype file.")

    parser.add_argument("--pheno_columns", action="store", default=[0], type=int, nargs="+", dest="pheno_columns", help="Specify which columns that the phenotype data is in the provided phenotype data file. Multiple column numbers can be specified to conduct regression with multiple phenotypes. Default is the first column.")

    parser.add_argument("--pheno_delim", action="store", default=",", dest="pheno_delim", help="Specify the delimiter for the phenotype data file. Default delimiter is comma. Set quotation marks around the delimiter when applied.")

    parser.add_argument("--pheno_no_header", action="store_true", default=False, dest="pheno_no_header", help="Specify whether the phenotype has a header row.")

    parser.add_argument("--covar_columns", action="store", default=[], type=int, nargs="+", dest="covar_columns", help="Specify which columns that the phenotype data is in the provided phenotype data file. Multiple column numbers can be specified to conduct regression with multiple phenotypes. Default is the first column.")


    results=parser.parse_args()

    outputPath=results.OUTPUT

    """ configure logging control """

    logger = logging.getLogger(results.app_name)
    logger.setLevel(logging.DEBUG)
    # Create file handler which logs even debug messages
    fh = logging.FileHandler(outputPath+".log")
    fh.setLevel(logging.DEBUG)
    # Create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # Create formatter and add it to the handlers
    formatter1 = logging.Formatter('%(asctime)s %(levelname)s : %(message)s')
    formatter2 = logging.Formatter('%(asctime)s %(levelname)s : %(message)s')

    ch.setFormatter(formatter1)
    fh.setFormatter(formatter2)
    # Add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    # Type of files, VCF or GEN
    filetype=results.filetype

    # Log file
    snp_log=results.snp_log
    
    # Setting parameters
    gwas_id=results.gwas_id    # column of SNP ID
    gwas_p=results.gwas_p      # column of P value
    gwas_or=results.gwas_or    # column of odds ratio
    gwas_a1=results.gwas_a1    # column of A1 in the GWAS
    gwas_a2=results.gwas_a2    # column of A2 in the GWAS
    gwas_a1f=results.gwas_a1f  # column index of A1 frequency in the GWAS

    # Define column number for contents in genfile
    if filetype.lower()=="vcf":
        geno_id= 2 # column number with rsID
        geno_start=9 # column number of the 1st genotype, in the raw vcf files, after separated by the delimiter of choice
        geno_a1 = 3 # column number that contains the reference allele
        GENO_delim= "\t"
    elif filetype.lower()=="gen":
        geno_id = 1
        geno_start=5
        geno_a1=3
        GENO_delim= " "

    step=0.01
    # List of thresholds:
    thresholds=results.thresholds
    threshold_seq=results.threshold_seq
    if threshold_seq is not None:
        if len(threshold_seq)==3:
            lower=min(threshold_seq[0:2])
            upper=max(threshold_seq[0:2])
            step=threshold_seq[2]
            thresholds=np.arange(lower, upper+step, step).tolist()
        else:
            raise("Invalid input for threshold sequence parameters")
            logger.error("Invalid input for threshold sequence parameters")

    # GWAS file information
    gwasFiles=results.GWAS
    GWAS_has_header=results.GWAS_header
    GWAS_delim=results.GWAS_delim

    # Programme parameter
    log_or=results.log_or        # specify whether you want to log your odds ratios
    check_ref=results.check_ref  # if you know that there are mismatch between the top strand in the genotypes and that of the GWAS, set True. Not checking the reference allele will improve the speed
    use_maf=results.use_maf      # whether to use MAF to check reference allele

    # Sample file path and name
    sampleFilePath=results.sample_file    # include the full/relative path and name of the sample file
    sampleFileDelim=results.sample_delim  # sample File Delimiter
    sampleFileID=results.sample_file_ID   # which column in the sample file has the ID
    sample_skip=results.sample_skip       # how many lines to skip so that the sample names can be matched to the genotypes 1-to-1, taking into account the header of the sample file
    
    # Output file information
    outputPath=results.OUTPUT

    # Specify whether to check for duplicate SNPs
    checkDup=results.checkdup


    # Get the name of the genotype files
    genoFileNamePattern=results.GENO
    if "file:/" in genoFileNamePattern:
        genoFileNamePattern=re.sub("file://", "", genoFileNamePattern)


    # Get the whole list of the file names
    genoFileNames=glob.glob(genoFileNamePattern)

    # Parameter for phenotype regression
    pheno_file=results.pheno_file
    pheno_columns=results.pheno_columns
    pheno_delim=results.pheno_delim
    pheno_no_header=results.pheno_no_header
    covar_columns=results.covar_columns


    ''' ####  start spark   ##### '''
    ## Start timing
    totalstart=time.time()
    ## Start spark context
    APP_NAME=results.app_name

    spark=SparkSession.builder.appName(APP_NAME).getOrCreate()

    # If using spark < 2.0.0, use the pyspark module to make Spark context
    #conf=pyspark.SparkConf().setAppName(APP_NAME).set()#.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")

    sc  = spark.sparkContext

    #sc = spark.sparkContext
    sc.setLogLevel("WARN")
    log4jLogger = sc._jvm.org.apache.log4j
    LOGGER = log4jLogger.LogManager.getLogger(__name__)
    logger.info("Start reading files")
    logger.info("Using these genotype files: ")

    for filename in genoFileNames[:min(24, len(genoFileNames))]:
        logger.info(filename)
    if len(genoFileNames)>23:
        logger.info("and more...")

    logger.info("Total of {} genotype files".format(str(len(genoFileNames))))
    # 1. Load files

    # Read the raw data
    genodata=sc.textFile(genoFileNamePattern)
    #logger.info("Using the GWAS file: {}".format(ntpath.basename(gwasFiles)))
    logger.info("Using the GWAS file: {}".format(gwasFiles))
    gwastable=spark.read.option("header",GWAS_has_header).option("delimiter",GWAS_delim).csv(gwasFiles).cache()
    logger.info("Showing top 5 rows of GWAS file")
    gwastable.show(5)

    logger.info("System recognizes the following information in the GWAS based on user inputs:")
    logger.info("SNP ID : Column {}".format(gwas_id))
    logger.info("P-values : Column {}".format(gwas_p))
    logger.info("Effect size : Column {}".format(gwas_or))
    logger.info("Allele A1 : Column {}".format(gwas_a1))
    logger.info("Allele A2 : Column {}".format(gwas_a2))
    if use_maf:
        logger.info("A1 allele frequencies : Column {}".format(gwas_a1f))

    # 1.1 Filter GWAS and prepare odds ratio

    # Filter the genotype to contain only the SNPs less than the maximum p value threshold in the GWAS
    maxThreshold=max(thresholds)  # maximum p value
    gwasOddsMapMax=filterGWASByP_DF(GWASdf=gwastable, pcolumn=gwas_p, idcolumn=gwas_id, oddscolumn=gwas_or, pHigh=maxThreshold, logOdds=log_or)
    gwasOddsMapMaxCA=sc.broadcast(gwasOddsMapMax).value  # Broadcast the map

    # ### 2. Initial processing
    # At this step, the genotypes are already filtered to keep only the ones in 'gwasOddsMapMax'
    bpMap={"A":"T", "T":"A", "C":"G", "G":"C"}
    tic=time.time()
    if filetype.lower()=="vcf":
        logger.info("Genotype data format : VCF ")

        # [chrom, bp, snpid, A1, A2, *genotype]
        genointermediate=genodata.filter(lambda line: ("#" not in line)).map(lambda line: line.split(GENO_delim)).filter(lambda line: line[geno_id] in gwasOddsMapMaxCA).map(lambda line: line[0:5]+[chunk.strip('"').split(":")[3] for chunk in line[geno_start::]]).map(lambda line: line[0:5]+[triplet.split(",") for triplet in line[5::]])

        ## (snpid, [genotypes])
        genotable=genointermediate.map(lambda line: (line[geno_id], list(itertools.chain.from_iterable(line[5::])))).mapValues(lambda geno: [float(x) for x in geno])

        if check_ref:
            if use_maf:

                logger.info("Determining strand alignment using MAF")
                genoA1f=genointermediate.map(lambda line: (line[geno_id], (line[geno_a1], line[geno_a1+1]), [float(x) for x in list(itertools.chain.from_iterable(line[5::]))])).map(lambda line: (line[0], line[1][0], line[1][1], getA1f(line[2]))).toDF(["Snpid_geno", "GenoA1", "GenoA2", "GenoA1f"])

                # 'GwasA1F' means the allele of the A1 frequency in the GWAS
                gwasA1f=gwastable.rdd.map(lambda line:(line[gwas_id], line[gwas_a1], line[gwas_a2], line[gwas_a1f])).toDF(["Snpid_gwas", "GwasA1", "GwasA2", "GwasA1F"])

                # checktable = [ geno_snpid, genoA1, genoA2, genoA1f, gwas_snpid, gwasA1, gwasA2, gwasA1f]
                checktable=genoA1f.join(gwasA1f, genoA1f["Snpid_geno"]==gwasA1f["Snpid_gwas"], "inner").cache()
                if checkDup:
                    logger.info("Searching and removing duplicated SNPs")
                    flagList = checktable.rdd.map(lambda line: checkAlignmentDF(line, bpMap)).collect()  #  (snpid, flag)
                    flagMap = rmDup(flagList)
                else:
                    flagMap = checktable.rdd.map(lambda line: checkAlignmentDF(line, bpMap)).collectAsMap()

            else:
                logger.info("Determining strand alignment without using MAF; SNPs with alleles that are reverse complements will be discarded")
                genoalleles=genointermediate.map(lambda line: (line[geno_id], (line[geno_a1], line[geno_a1+1]), [float(x) for x in list(itertools.chain.from_iterable(line[5::]))])).map(lambda line: (line[0], line[1][0], line[1][1])).toDF(["Snpid_geno", "GenoA1", "GenoA2"])

                gwasalleles=gwastable.rdd.map(lambda line:(line[gwas_id], line[gwas_a1], line[gwas_a2])).toDF(["Snpid_gwas", "GwasA1", "GwasA2"])

                checktable=genoalleles.join(gwasalleles, genoalleles["Snpid_geno"]==gwasalleles["Snpid_gwas"], "inner").cache()

                if checkDup:
                    logger.info("Searching and removing duplicated SNPs")
                    flagList = checktable.rdd.map(lambda line: checkAlignmentDFnoMAF(line, bpMap)).collect()
                    flagMap = rmDup(flagList)
                else:
                    # no need to check the duplicates if the data is preprocessed
                    flagMap = checktable.rdd.map(lambda line: checkAlignmentDFnoMAF(line, bpMap)).collectAsMap()

            logger.info("Generating genotype dosage while taking into account difference in strand alignment")
            flagMap=sc.broadcast(flagMap).value
            genotypeMax=genotable.filter(lambda line: line[0] in flagMap and flagMap[line[0]]!="discard").map(lambda line: makeGenotypeCheckRef(line, checkMap=flagMap)).cache()

        else:
            logger.info("Generating genotype dosage without checking reference allele alignments")
            genotypeMax=genotable.mapValues(lambda line: makeGenotype(line)).cache()
            flagMap=False
            if checkDup:
                logger.info("Searching and removing duplicated SNPs")
                genotypeCount=genotypeMax.map(lambda line: (line[0], 1)).reduceByKey(lambda a,b: a+b).filter(lambda line: line[1]==1).collectAsMap()
                genotypeMax=genotypeMax.filter(lambda line: line[0] in genotypeCount)

    elif filetype.lower() == "gen":
        logger.info("Genotype data format : GEN")
        genotable=genodata.map(lambda line: line.split(GENO_delim)).filter(lambda line: line[geno_id] in gwasOddsMapMaxCA).map(lambda line: (line[geno_id], line[geno_start::])).mapValues(lambda geno: [float(call) for call in geno])
        if check_ref:
            if use_maf:
                logger.info("Determining strand alignment using MAF")
                genoA1f=genodata.map(lambda line: line.split(GENO_delim)).map(lambda line: (line[geno_id], line[geno_a1], line[geno_a1+1], getA1f([float(x) for x in line[geno_start::]]))).toDF(["Snpid_geno", "GenoA1", "GenoA2", "GenoA1f"])
                gwasA1f=gwastable.rdd.map(lambda line:(line[gwas_id], line[gwas_a1], line[gwas_a2], line[gwas_a1f])).toDF(["Snpid_gwas", "GwasA1", "GwasA2", "GwasA1f" ])
                checktable=genoA1f.join(gwasA1f, genoA1f["Snpid_geno"]==gwasA1f["Snpid_gwas"], "inner").cache()
                if checkDup:
                    logger.info("Searching and removing duplicated SNPs")
                    flagList = checktable.rdd.map(lambda line: checkAlignmentDF(line, bpMap)).collect()
                    flagMap = rmDup(flagList)
                else:
                    flagMap = checktable.rdd.map(lambda line: checkAlignmentDF(line, bpMap)).collectAsMap()
            else:
                logger.info("Determining strand alignment without using MAF; SNPs with alleles that are reverse complements will be discarded")
                genoalleles=genodata.map(lambda line: line.split(GENO_delim)).map(lambda line: (line[geno_id], line[geno_a1], line[geno_a1+1])).toDF(["Snpid_geno", "GenoA1", "GenoA2"])
                gwasalleles=gwastable.rdd.map(lambda line:(line[gwas_id], line[gwas_a1], line[gwas_a2])).toDF(["Snpid_gwas", "GwasA1", "GwasA2"])
                checktable=genoalleles.join(gwasalleles, genoalleles["Snpid_geno"]==gwasalleles["Snpid_gwas"], "inner").cache()

                if checkDup:
                    logger.info("Searching and removing duplicated SNPs")
                    flagList = checktable.rdd.map(lambda line: checkAlignmentDFnoMAF(line, bpMap)).collect()
                    flagMap = rmDup(flagList)
                else:
                    flagMap = checktable.rdd.map(lambda line: checkAlignmentDFnoMAF(line, bpMap)).collectAsMap()

            logger.info("Generating genotype dosage while taking into account difference in strand alignment")
            flagMap=sc.broadcast(flagMap).value
            genotypeMax=genotable.filter(lambda line: line[0] in flagMap and flagMap[line[0]]!="discard" ).map(lambda line: makeGenotypeCheckRef(line, checkMap=flagMap)).cache()

        else:
            logger.info("Generating genotype dosage without checking allele alignments")
            genotypeMax=genotable.mapValues(lambda line: makeGenotype(line)).cache()
            flagMap=False
            if checkDup:
                logger.info("Searching and removing duplicated SNPs")
                genotypeCount=genotypeMax.map(lambda line: (line[0], 1)).reduceByKey(lambda a,b: a+b).filter(lambda line: line[1]==1).collectAsMap()
                genotypeMax=genotypeMax.filter(lambda line: line[0] in genotypeCount)

    logger.info("Dosage generated in {:.1f} seconds".format(time.time()-tic) )
    samplesize=int(len(genotypeMax.first()[1]))
    logger.info("Detected {} samples in genotype data" .format(str(samplesize)))

    #genoa1f.map(lambda line:"\t".join([line[0], "\t".join(line[1]), str(line[2])])).saveAsTextFile("../name_maf")

    # Calculate PRS at the specified thresholds
    if flagMap:
      genocalltable=genotable.filter(lambda line: line[0] in flagMap and flagMap[line[0]]!="discard" ).mapValues(lambda geno: getCall(geno)).cache()
    else:
      genocalltable=genotable.mapValues(lambda geno: getCall(geno))

    assert len(genocalltable.first()[1])==samplesize, "Bug found, size of genotype and call table differ"



    def calcPRSFromGeno(genotypeRDD, oddsMap, samplenum,calltable=False):
        assert calltable, "**** No call table found *****"
        totalcount=calltable.map(lambda line: line[1]).reduce(lambda a,b: map(add, a, b))
        multiplied=genotypeRDD.map(lambda line:[call * oddsMap[line[0]] for call in line[1]])
        PRS=multiplied.reduce(lambda a,b: map(add, a, b))
        normalizedPRS=[x/count if count != 0 else x for x, count in zip(PRS, totalcount)]
        return (totalcount,normalizedPRS)

    def calcAll(genotypeRDD, gwasRDD, thresholdlist, logsnp, samplenum,calltableRDD=False):
        logger.info("Started calculating PRS at each threshold")
        prsMap={}
        thresholdNoMaxSorted=sorted(thresholdlist, reverse=True)
        thresholdmax=max(thresholdlist)
        idlog={}
        start=time.time()
        for threshold in thresholdNoMaxSorted:
            tic=time.time()
            gwasFilteredBC=sc.broadcast(filterGWASByP_DF(GWASdf=gwasRDD, pcolumn=gwas_p, idcolumn=gwas_id, oddscolumn=gwas_or, pHigh=threshold, logOdds=log_or))
            #gwasFiltered=spark.sql("SELECT snpid, gwas_or_float FROM gwastable WHERE gwas_p_float < {:f}".format(threshold)
            logger.info("Filtered GWAS at threshold of {}. Time spent : {:.1f} seconds".format(str(threshold), time.time()-tic))
            checkpoint=time.time()
            filteredgenotype=genotypeRDD.filter(lambda line: line[0] in gwasFilteredBC.value)
            assert calltableRDD, "Error, calltable must be provided"
            filteredcalltable=calltableRDD.filter(lambda line: line[0] in gwasFilteredBC.value )
            if not filteredgenotype.isEmpty():
              #assert filteredcalltable.count()==filteredgenotype.count(), "Error, call table have different size from genotype"
              if logsnp:
                idlog[threshold]=filteredgenotype.map(lambda line:line[0]).collect()

              prsMap[threshold]=calcPRSFromGeno(filteredgenotype, gwasFilteredBC.value,samplenum=samplenum, calltable=filteredcalltable)

              logger.info("Finished calculating PRS at threshold of {}. Time spent : {:.1f} seconds".format(str(threshold), time.time()-checkpoint))

            else:
              logger.warn("No SNPs left at threshold {}" .format(threshold))
        return prsMap, idlog

    prsDict, snpids=calcAll(genotypeMax,gwastable, thresholds, logsnp=snp_log, samplenum=samplesize,calltableRDD=genocalltable)

    # Log which SNPs are used in PRS
    if snp_log:
        if flagMap:
            logoutput=writeSNPlog(snpids, outputPath, logger, flagMap=flagMap)
        else:
            logoutput=writeSNPlog(snpids, outputPath, logger)

    # Generate labels for samples
    #if filetype.lower()=="vcf":
        #subjNames=genodata.filter(lambda line: "#CHROM" in line).map(lambda line: line.split(GENO_delim)[9::]).collect()[0]
        #output=writePRS(prsDict,  outputPath, samplenames=subjNames)

    if sampleFilePath!="NOSAMPLE":
        # get sample name from the provided sample file
        subjNames=getSampleNames(sampleFilePath,sampleFileDelim,sampleFileID, skip=sample_skip)
        
        logger.info("This sample file was used to provide sample ID : {}".format(sampleFilePath))
        output=writePRS(prsDict,  outputPath, logger, samplenames=subjNames)
    else:
        output=writePRS(prsDict,  outputPath,logger=logger, samplenames=None)


    if pheno_file is not None:
        phenotypes, thresholds, r2All, pAll=regression(prsDict,pheno_file, pheno_delim, pheno_columns, pheno_no_header, covarColumns=covar_columns, outputName=outputPath, logger=logger)

        r_square_plots(phenotypes,r2All,pAll, thresholds, outputName=outputPath, width = 3,bar_width = step)

    sc.stop()
    seconds=time.time()-totalstart
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    logger.info("Total calculation time : {:d} hrs {:02d} min {:02d} sec".format(int(h), int(m), int(round(s))))

# Remove "spark-warehouse" if it exists and is empty
try:
    os.rmdir(os.path.dirname(os.path.expanduser(outputPath))+"/spark-warehouse")
except OSError:
    pass
