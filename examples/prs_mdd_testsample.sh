## CALCULATE PRS SCORE

### using vcf file
spark-submit PRS_run.py \
test_sample.vcf \
GWAS_MDD_beta_clump_withaf.txt \
test_sample_vcf_GWAS_MDD_beta_clump_withaf \
--sample_file test_sample.samples \
--sample_file_skip 2 \
--sample_file_delim " " \
--filetype VCF \
--thresholds 1 0.5 0.4 0.3 0.2 0.1 0.05 0.01 0.001 0.0001 \
--GWAS_delim "\t" \
--snp_log

### using gen file
spark-submit PRS_run.py \
test_sample.gen \
GWAS_MDD_beta_clump_withaf.txt \
test_sample_gen_GWAS_MDD_beta_clump_withaf \
--sample_file test_sample.samples \
--sample_file_skip 2 \
--sample_file_delim " " \
--filetype GEN \
--thresholds 1 0.5 0.4 0.3 0.2 0.1 0.05 0.01 0.001 0.0001 \
--GWAS_delim "\t" \
--snp_log
