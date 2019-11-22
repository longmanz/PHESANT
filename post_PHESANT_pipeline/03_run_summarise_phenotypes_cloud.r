library(data.table)
source("~/Repositories/PHESANT/summarise_phenotypes.r")

root_file <- "ukb11214_final_january_reference_QC_more_phenos_and_corrected"
root_file <- "ukb11214_final_july_reference_QC_more_phenos_and_corrected"

phenotype_info_file <- "variable-info/outcome_info_final_round2.tsv"
phenotype_info_file <- "../variable-info/outcome_info_final_round2.tsv"

coding_info_file <- "variable-info/data-coding-ordinal-info.txt"
coding_info_file <- "../variable-info/data-coding-ordinal-info.txt"
out_root_file <- "neale_lab_parsed_and_restricted_to_QCed_samples"

only_males_file <- "../should_only_be_in_males.tsv"
only_females_file <- "../should_only_be_in_females.tsv"

for (i in 1:4)
{
	# Both sexes
	hist_filename <- paste0(out_root_file, "_cat_variables_both_sexes.", i, "_hist")
	pheno_summary <- paste0(out_root_file, "_cat_variables_both_sexes.", i, "_phenosummary.tsv")

	filename <- paste0(out_root_file, "_cat_variables_both_sexes.", i)
	tsv_filename <- paste(filename, ".tsv", sep="")
	log_filename <- paste0(root_file, ".", i)
	log_file <- paste(log_filename, ".log", sep="")
	tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
	names(tsv_data)[1] <- "userId"

	outcome_info <- read.table(phenotype_info_file, sep='\t', quote="", comment.char="", header=TRUE)
	print(paste0("both sexes ", i))
	summary_11214 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
	write.table(summary_11214, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

	# Males
	hist_filename <- paste0(out_root_file, "_cat_variables_males.", i, "_hist")
	pheno_summary <- paste0(out_root_file, "_cat_variables_males.", i, "_phenosummary.tsv")

	filename <- paste0(out_root_file, "_cat_variables_males.", i)
	tsv_filename <- paste(filename, ".tsv", sep="")
	log_filename <- paste0(root_file, ".", i)
	log_file <- paste(log_filename, ".log", sep="")
	tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
	names(tsv_data)[1] <- "userId"

	outcome_info <- read.table(phenotype_info_file, sep='\t', quote="", comment.char="", header=TRUE)
	print(paste0("males ", i))
	summary_11214 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
	write.table(summary_11214, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

	# Females
	hist_filename <- paste0(out_root_file, "_cat_variables_females.", i, "_hist")
	pheno_summary <- paste0(out_root_file, "_cat_variables_females.", i, "_phenosummary.tsv")

	filename <- paste0(out_root_file, "_cat_variables_females.", i)
	tsv_filename <- paste(filename, ".tsv", sep="")
	log_filename <- paste0(root_file, ".", i)
	log_file <- paste(log_filename, ".log", sep="")
	tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
	names(tsv_data)[1] <- "userId"

	outcome_info <- read.table(phenotype_info_file, sep='\t', quote="", comment.char="", header=TRUE)
	print(paste0("females ", i))
	summary_11214 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
	write.table(summary_11214, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
}

create_cts_summary_file_from_get_hist_and_notes_output <- function(get_hists_and_notes_output, only_males, only_females, males_females_or_both='both')
{
	get_hists_and_notes_output <- data.frame(get_hists_and_notes_output, stringsAsFactors = FALSE)
	get_hists_and_notes_output$variable.type <- get_hists_and_notes_output$PHESANT.notes
	get_hists_and_notes_output$variable.type[grep("IRNT", get_hists_and_notes_output$PHESANT.notes)] <- "CONTINOUS IRNT"
	
	if(males_females_or_both == "both") {
		get_hists_and_notes_output <- get_hists_and_notes_output[!(rownames(get_hists_and_notes_output) %in% c(only_males$FullFieldID, only_females$FullFieldID)),]
	} else if(males_females_or_both == "males") {
		get_hists_and_notes_output <- get_hists_and_notes_output[!(rownames(get_hists_and_notes_output) %in% only_females$FullFieldID),]
	} else {
		get_hists_and_notes_output <- get_hists_and_notes_output[!(rownames(get_hists_and_notes_output) %in% only_males$FullFieldID),]
	}

	get_hists_and_notes_output_cp <- get_hists_and_notes_output
	get_hists_and_notes_output_cp$variable.type[grep("IRNT", get_hists_and_notes_output$PHESANT.notes)] <- "CONTINOUS RAW"
	get_hists_and_notes_output_cp$PHESANT.notes <- gsub(" IRNT \\|\\|", "", get_hists_and_notes_output_cp$PHESANT.notes)
	rownames(get_hists_and_notes_output_cp) <- paste0(rownames(get_hists_and_notes_output_cp), "_raw")
	rownames(get_hists_and_notes_output) <- paste0(rownames(get_hists_and_notes_output), "_irnt")
	get_hists_and_notes_output <- rbind(get_hists_and_notes_output, get_hists_and_notes_output_cp)
	return(get_hists_and_notes_output)
}

# Now, remove variables that are in male only or female only from both sexes and the opposite sex phenotype summary file.
only_males <- fread(only_males_file, header=FALSE, sep='\t')
only_females <- fread(only_females_file, header=FALSE, sep='\t')
names(only_males) <- c("FullFieldID_rowname","FullFieldID", "FieldID", "Field", "SubFieldID", "SubField",
	"N.non.missing.males", "N.missing.males", "N.controls.males", "N.cases.males",
	 "N.non.missing", "N.missing", "N.controls", "N.cases")
names(only_females) <- c("FullFieldID_rowname", "FullFieldID", "FieldID", "Field", "SubFieldID", "SubField",
	"N.non.missing.females", "N.missing.females", "N.controls.females", "N.cases.females",
	 "N.non.missing", "N.missing", "N.controls", "N.cases")

# We also need to do the cts phenotypes - new format now because of how I parsed them.
system(paste0("awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' ", root_file, ".*.log > ", root_file, "_all.log"))

# Both sexes
hist_filename <- paste0(out_root_file, "_cts_irnt_both_sexes_hist")
pheno_summary <- paste0(out_root_file, "_cts_irnt_both_sexes_phenosummary.tsv")
summary_11214_file_sex_specific <- paste0(out_root_file, "_cts_both_sexes_phenosummary.tsv")

filename <- paste0(out_root_file, "_cts_irnt")
tsv_filename <- paste(filename, ".tsv", sep="")
log_filename <- paste0(root_file, "_all")
log_file <- paste(log_filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
names(tsv_data)[1] <- "userId"

summary_11214 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
write.table(summary_11214, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
final_summary_11214 <- create_cts_summary_file_from_get_hist_and_notes_output(summary_11214, only_males, only_females)
write.table(final_summary_11214, file=summary_11214_file_sex_specific, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

# Males
hist_filename <- paste0(out_root_file, "_cts_irnt_males_hist")
pheno_summary <- paste0(out_root_file, "_cts_irnt_males_phenosummary.tsv")
summary_11214_file_sex_specific <- paste0(out_root_file, "_cts_males_phenosummary.tsv")

filename <- paste0(out_root_file, "_cts_irnt_males")
tsv_filename <- paste(filename, ".tsv", sep="")
log_filename <- paste0(root_file, "_all")
log_file <- paste(log_filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
names(tsv_data)[1] <- "userId"

summary_11214 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
write.table(summary_11214, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
final_summary_11214 <- create_cts_summary_file_from_get_hist_and_notes_output(summary_11214, only_males, only_females, "males")
write.table(final_summary_11214, file=summary_11214_file_sex_specific, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

# Females
hist_filename <- paste0(out_root_file, "_cts_irnt_females_hist")
pheno_summary <- paste0(out_root_file, "_cts_irnt_females_phenosummary.tsv")
summary_11214_file_sex_specific <- paste0(out_root_file, "_cts_females_phenosummary.tsv")

filename <- paste0(out_root_file, "_cts_irnt_females")
tsv_filename <- paste(filename, ".tsv", sep="")
log_filename <- paste0(root_file, "_all")
log_file <- paste(log_filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
names(tsv_data)[1] <- "userId"

summary_11214 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
write.table(summary_11214, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
final_summary_11214 <- create_cts_summary_file_from_get_hist_and_notes_output(summary_11214, only_males, only_females, "females")
write.table(final_summary_11214, file=summary_11214_file_sex_specific, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

for (i in 1:4)
{
	# Both sexes.
	pheno_summary_file <- paste0(out_root_file, "_cat_variables_both_sexes.", i, "_phenosummary.tsv")
	pheno_summary_codings_included_file <- paste0(out_root_file, "_cat_variables_both_sexes_phesant_recodings.", i, "_phenosummary.tsv")
	pheno_summary_codings_included_file_sex_specific <- paste0(out_root_file, "_cat_variables_both_sexes_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")
	pheno_summary_codings_included <- include_PHESANT_reassignment_names(pheno_summary_file, outcome_info)
	write.table(pheno_summary_codings_included, file=pheno_summary_codings_included_file,
		col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

	pheno_summary_codings_included <- pheno_summary_codings_included[!(rownames(pheno_summary_codings_included) %in% c(only_males$FullFieldID, only_females$FullFieldID)),]
	pheno_summary_codings_included$variable.type <- pheno_summary_codings_included$PHESANT.notes
	pheno_summary_codings_included$variable.type[grep("CAT-ORD", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-ORDERED"
	pheno_summary_codings_included$variable.type[grep("CAT-SINGLE-BINARY", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-UNORDERED"
	pheno_summary_codings_included$variable.type[grep("CAT-MUL-BINARY-VAR", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-UNORDERED"
	pheno_summary_codings_included$variable.type[grep("two bins and treat as binary", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-UNORDERED"

	write.table(pheno_summary_codings_included, file=pheno_summary_codings_included_file_sex_specific,
		col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
	# Males.
	pheno_summary_file <- paste0(out_root_file, "_cat_variables_males.", i, "_phenosummary.tsv")
	pheno_summary_codings_included_file <- paste0(out_root_file, "_cat_variables_males_phesant_recodings.", i, "_phenosummary.tsv")
	pheno_summary_codings_included_file_sex_specific <- paste0(out_root_file, "_cat_variables_males_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")

	pheno_summary_codings_included <- include_PHESANT_reassignment_names(pheno_summary_file, outcome_info)
	write.table(pheno_summary_codings_included, file=pheno_summary_codings_included_file,
		col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
	pheno_summary_codings_included <- pheno_summary_codings_included[!(rownames(pheno_summary_codings_included) %in% only_females$FullFieldID),]
	pheno_summary_codings_included$variable.type <- pheno_summary_codings_included$PHESANT.notes
	pheno_summary_codings_included$variable.type[grep("CAT-ORD", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-ORDERED"
	pheno_summary_codings_included$variable.type[grep("CAT-SINGLE-BINARY", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-UNORDERED"
	pheno_summary_codings_included$variable.type[grep("CAT-MUL-BINARY-VAR", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-UNORDERED"
	pheno_summary_codings_included$variable.type[grep("two bins and treat as binary", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-UNORDERED"
	
	write.table(pheno_summary_codings_included, file=pheno_summary_codings_included_file_sex_specific,
		col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

	# Females.
	pheno_summary_file <- paste0(out_root_file, "_cat_variables_females.", i, "_phenosummary.tsv")
	pheno_summary_codings_included_file <- paste0(out_root_file, "_cat_variables_females_phesant_recodings.", i, "_phenosummary.tsv")
	pheno_summary_codings_included_file_sex_specific <- paste0(out_root_file, "_cat_variables_females_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")

	pheno_summary_codings_included <- include_PHESANT_reassignment_names(pheno_summary_file, outcome_info)
	write.table(pheno_summary_codings_included, file=pheno_summary_codings_included_file,
		col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
	pheno_summary_codings_included <- pheno_summary_codings_included[!(rownames(pheno_summary_codings_included) %in% only_males$FullFieldID),]
	pheno_summary_codings_included$variable.type <- pheno_summary_codings_included$PHESANT.notes
	pheno_summary_codings_included$variable.type[grep("CAT-ORD", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-ORDERED"
	pheno_summary_codings_included$variable.type[grep("CAT-SINGLE-BINARY", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-UNORDERED"
	pheno_summary_codings_included$variable.type[grep("CAT-MUL-BINARY-VAR", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-UNORDERED"
	pheno_summary_codings_included$variable.type[grep("two bins and treat as binary", pheno_summary_codings_included$PHESANT.notes)] <- "CAT-UNORDERED"
	write.table(pheno_summary_codings_included, file=pheno_summary_codings_included_file_sex_specific,
		col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

}

# for (i in 1:4)
# {
# 	filename <- paste0(root_file, ".", i)
# 	tsv_filename <- paste(filename, ".tsv", sep="")
# 	log_file <- paste(filename, ".log", sep="")
# 	tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')

# 	outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final_round2.tsv",
# 						   sep='\t', quote="", comment.char="", header=TRUE)

# 	categorical_11214 <-  get_barplot_numbers(tsv_data, log_file, outcome_info, codings_tables)

# 	if(i == 1) {
# 		categorical_full <- categorical_11214
# 	} else {
# 		categorical_full <- rbind(categorical_full, categorical_11214)
# 	}

# 	print(dim(categorical_full))
# }

# write.table(categorical_full, file="~/Repositories/PHESANT/variable-info/PHESANT_categoricals.tsv",
# 	col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
