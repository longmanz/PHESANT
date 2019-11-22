library(data.table)
source("~/Repositories/PHESANT/summarise_phenotypes.r")

# Latest.
filename_root <- "ukb11214_final_january_reference_QC_more_phenos_and_corrected."
# July run
# Make sure I'm in the right subdirectory.
filename_root <- "ukb11214_final_july_reference_QC_more_phenos_and_corrected."

phenotype_info_file <- "variable-info/outcome_info_final_round2.tsv"
phenotype_info_file <- "../variable-info/outcome_info_final_round2.tsv"

coding_info_file <- "variable-info/data-coding-ordinal-info.txt"
coding_info_file <- "../variable-info/data-coding-ordinal-info.txt"

for (i in 1:4)
{
	# Both sexes
	hist_filename <- paste0(filename_root, i, "_hist")
	pheno_summary <- paste0(filename_root, i, "_phenosummary.tsv")

	filename <- paste0(filename_root, i)
	tsv_filename <- paste(filename, ".tsv", sep="")
	log_filename <- paste0(filename_root, i)
	log_file <- paste(log_filename, ".log", sep="")
	tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
	names(tsv_data)[1] <- "userId"

	outcome_info <- read.table(phenotype_info_file, sep='\t', quote="", comment.char="", header=TRUE)
	print(paste0("both sexes ", i))
	summary_11214 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=4)
	write.table(summary_11214, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
}