# First ensure that the PHESANT output is moved to a bucket on the cloud
# e.g. 	gsutil cp multi_ancestry_jan_2020*gz gs://ukb-diverse-pops/Phenotypes/Everyone/PHESANT_output/jan2020/
# 		gsutil cp multi_ancestry_jan_2020*log gs://ukb-diverse-pops/Phenotypes/Everyone/PHESANT_output/jan2020/

# Navigate to PHESANT/post_PHESANT_multi_ancestry_pipeline and run this code.
# 01_run_summarise_multi_ancestry_phenotypes_initial_PHESANT_run.r

n_chunks <- 10
# Updated version :
variable_info <- "../variable-info/outcome_info_final_multi_ancestry_jan2020.tsv"
coding_info_file <- "../variable-info/data-coding-ordinal-info-nov2019-update.txt"

# Copy down from the cloud.
cloud_filename_root <- "gs://ukb-diverse-pops/Phenotypes/Everyone/PHESANT_output/jan2020/multi_ancestry_jan_2020."
filename_root <- "../../multi_ancestry_jan_2020."

intermediate_output_location <- "gs://ukb-diverse-pops/Phenotypes/Everyone/PHESANT_intermediate_output/jan_2020/"
plot_location <- "gs://ukb-diverse-pops/Phenotypes/Everyone/plots/jan_2020/"

source("01_run_summarise_multi_ancestry_phenotypes_initial_PHESANT_run.r")

# 02_run_all_sexes_to_male_female_multi_ancestry

sex_df_path <- "../../neale_lab_parsed_QC_Oct2019_sex_sampleID.tsv"

n_cts_min <- 5000
n_cat_ordered <- 5000
n_cat_min <- 100

QCed_io_name <- '../../neale_lab_parsed_QC_Oct2019'
QCed_io_name_cloud <- 'gs://ukb-diverse-pops/Phenotypes/Everyone/neale_lab_parsed_QC_Oct2019.tsv'

source("02_run_all_sexes_to_male_female_multi_ancestry.r")

only_males_file <- "../should_only_be_in_males.tsv"
only_females_file <- "../should_only_be_in_females.tsv"

source("03_run_summarise_multi_ancestry_phenotypes_cloud.r")

date <- 'January_2020_plus_pharma_and_updated_codings'
final_output <- '../../phesant_output_multi_ancestry_'
final_output_location <- "gs://ukb-diverse-pops/Phenotypes/Everyone/PHESANT_final_output/"

source("04_run_combine_PHESANT_multi_ancestry_output.r")