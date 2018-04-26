## CALCULATE PRS

### using vcf file to create nine polygenic scores at sequential p-value thresholds 0.1, 0.15, 0.2, ..., 0.5

spark-submit ../PRS_run.py \
example.vcf \
gwasfile.txt \
example_vcf \
--sample example.sample \
--threshold_seq 0.1 0.5 0.05 \
--snp_log

### using gen file to create PRS while removing ambiguous SNPs

spark-submit ../PRS_run.py \
example.gen \
gwasfile.txt \
example_gen \
--sample example.sample \
--filetype GEN \
--no_a1f \
--snp_log
