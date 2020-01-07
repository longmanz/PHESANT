# 1. Run PHESANT on the cloud using neale_lab_parsed_and_restricted_to_QCed_samples.tsv as the phenotype file
See parsing_raw_pheno_initiation_script.

# 2. Run 01_run_summarise_phenotypes_initial_PHESANT_run.r 
Get the summary file for the PHESANT run on the cloud.

# 3. Run 02_run_all_sexes_to_male_female_exomes.r
This restricts the PHESANT run to males and females and removes phenotypes with too few examples. Note that this creates IRNT - here a new renormalised version for each sex. To get the pheno on the same scale as both_sexes, you need to restrict those phenotypes to men/women.
Note that this may segfault if there is not sufficient RAM.

# 4. Run 03_run_summarise_phenotypes_cloud.r. 
This then generates summary files for the sex-specific phenotype files that have been parsed by PHESANT, and restricts the summary files to the sex-specific.

# 5. Run 04_run_combine_PHESANT_output.r
Combine the categorical (1-n_chunks), the IRNT cts, and the raw cts variables.
