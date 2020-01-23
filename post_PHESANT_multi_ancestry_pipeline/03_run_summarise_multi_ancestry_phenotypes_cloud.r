library(data.table)
library(hyperSpec)

source("../summarise_phenotypes.r")
source("functions_for_run_all_sexes_male_female.r")

# Read in all the data table codings
to_read <- paste("../WAS/codings/", dir("../WAS/codings"), sep="")
codings_tables <- list()
filename_root <- gsub("\\.$", "", filename_root)

for(i in to_read) {
	name <- gsub("../WAS/codings/coding(.*).tsv", "\\1", i)
	codings_tables[[name]] <- read.table(i, header=TRUE, sep='\t', quote="", stringsAsFactors=FALSE)
}

for (i in 1:n_chunks)
{
	# Both sexes
	hist_filename <- paste0(QCed_io_name, "_cat_variables_both_sexes.", i, "_hist")
	pheno_summary <- paste0(QCed_io_name, "_cat_variables_both_sexes.", i, "_phenosummary.tsv")

	filename <- paste0(QCed_io_name, "_cat_variables_both_sexes.", i)
	tsv_filename <- paste(filename, ".tsv", sep="")
	log_filename <- paste0(filename_root, ".", i)
	log_file <- paste(log_filename, ".log", sep="")
	
	# If the file doesn't exists, then in that chunk, no phenotypes made it through PHESANT.
	if(!file.exists(tsv_filename)) next

	tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
	names(tsv_data)[1] <- "userId"

	outcome_info <- read.table(variable_info, sep='\t', quote="", comment.char="", header=TRUE)
	cat(paste0("both sexes ", i, '\n'))
	summary <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
	write.table(summary, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

	# Males
	hist_filename <- paste0(QCed_io_name, "_cat_variables_males.", i, "_hist")
	pheno_summary <- paste0(QCed_io_name, "_cat_variables_males.", i, "_phenosummary.tsv")

	filename <- paste0(QCed_io_name, "_cat_variables_males.", i)
	tsv_filename <- paste(filename, ".tsv", sep="")
	log_filename <- paste0(filename_root, ".", i)
	log_file <- paste(log_filename, ".log", sep="")

	if(file.exists(tsv_filename))
	{
		tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
		names(tsv_data)[1] <- "userId"

		outcome_info <- read.table(variable_info, sep='\t', quote="", comment.char="", header=TRUE)
		cat(paste0("males ", i, '\n'))
		summary <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
		write.table(summary, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
	}

	# Females
	hist_filename <- paste0(QCed_io_name, "_cat_variables_females.", i, "_hist")
	pheno_summary <- paste0(QCed_io_name, "_cat_variables_females.", i, "_phenosummary.tsv")

	filename <- paste0(QCed_io_name, "_cat_variables_females.", i)
	tsv_filename <- paste(filename, ".tsv", sep="")
	log_filename <- paste0(filename_root, ".", i)
	log_file <- paste(log_filename, ".log", sep="")

	if(file.exists(tsv_filename))
	{
		tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
		names(tsv_data)[1] <- "userId"

		outcome_info <- read.table(variable_info, sep='\t', quote="", comment.char="", header=TRUE)
		cat(paste0("females ", i, '\n'))
		summary <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
		write.table(summary, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
	}
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
system(paste0("awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' ", filename_root, ".*.log > ", filename_root, "_all.log"))

# Both sexes
hist_filename <- paste0(QCed_io_name, "_cts_irnt_both_sexes_hist")
pheno_summary <- paste0(QCed_io_name, "_cts_irnt_both_sexes_phenosummary.tsv")
summary_file_sex_specific <- paste0(QCed_io_name, "_cts_both_sexes_phenosummary.tsv")

filename <- paste0(QCed_io_name, "_cts_irnt")
tsv_filename <- paste(filename, ".tsv", sep="")
log_filename <- paste0(filename_root, "_all")
log_file <- paste(log_filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
names(tsv_data)[1] <- "userId"

outcome_info <- read.table(variable_info, sep='\t', quote="", comment.char="", header=TRUE)
summary <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
write.table(summary, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
final_summary <- create_cts_summary_file_from_get_hist_and_notes_output(summary, only_males, only_females)
write.table(final_summary, file=summary_file_sex_specific, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

# Males
hist_filename <- paste0(QCed_io_name, "_cts_irnt_males_hist")
pheno_summary <- paste0(QCed_io_name, "_cts_irnt_males_phenosummary.tsv")
summary_file_sex_specific <- paste0(QCed_io_name, "_cts_males_phenosummary.tsv")

filename <- paste0(QCed_io_name, "_cts_irnt_males")
tsv_filename <- paste(filename, ".tsv", sep="")
log_filename <- paste0(filename_root, "_all")
log_file <- paste(log_filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
names(tsv_data)[1] <- "userId"

summary <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
write.table(summary, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
final_summary <- create_cts_summary_file_from_get_hist_and_notes_output(summary, only_males, only_females, "males")
write.table(final_summary, file=summary_file_sex_specific, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

# Females
hist_filename <- paste0(QCed_io_name, "_cts_irnt_females_hist")
pheno_summary <- paste0(QCed_io_name, "_cts_irnt_females_phenosummary.tsv")
summary_file_sex_specific <- paste0(QCed_io_name, "_cts_females_phenosummary.tsv")

filename <- paste0(QCed_io_name, "_cts_irnt_females")
tsv_filename <- paste(filename, ".tsv", sep="")
log_filename <- paste0(filename_root, "_all")
log_file <- paste(log_filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
names(tsv_data)[1] <- "userId"

summary <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=2)
write.table(summary, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
final_summary <- create_cts_summary_file_from_get_hist_and_notes_output(summary, only_males, only_females, "females")
write.table(final_summary, file=summary_file_sex_specific, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

for (i in 1:n_chunks)
{
	# Both sexes.
	pheno_summary_file <- paste0(QCed_io_name, "_cat_variables_both_sexes.", i, "_phenosummary.tsv")

	if(file.exists(pheno_summary_file)) {
		n_lines <- strsplit(system(paste0("wc -l ", pheno_summary_file), intern=TRUE), split=" ")[[1]]
		n_lines <- as.integer(n_lines[length(n_lines)-1])
		if(n_lines > 1)
		{
			pheno_summary_codings_included_file <- paste0(QCed_io_name, "_cat_variables_both_sexes_phesant_recodings.", i, "_phenosummary.tsv")
			pheno_summary_codings_included_file_sex_specific <- paste0(QCed_io_name, "_cat_variables_both_sexes_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")
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
		}
	}

	# Males.
	pheno_summary_file <- paste0(QCed_io_name, "_cat_variables_males.", i, "_phenosummary.tsv")

	if(file.exists(pheno_summary_file)) {
		n_lines <- strsplit(system(paste0("wc -l ", pheno_summary_file), intern=TRUE), split=" ")[[1]]
		n_lines <- as.integer(n_lines[length(n_lines)-1])
		if(n_lines > 1)
		{
			pheno_summary_codings_included_file <- paste0(QCed_io_name, "_cat_variables_males_phesant_recodings.", i, "_phenosummary.tsv")
			pheno_summary_codings_included_file_sex_specific <- paste0(QCed_io_name, "_cat_variables_males_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")

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
		}
	}

	# Females.
	pheno_summary_file <- paste0(QCed_io_name, "_cat_variables_females.", i, "_phenosummary.tsv")

	if(file.exists(pheno_summary_file)) {
		n_lines <- strsplit(system(paste0("wc -l ", pheno_summary_file), intern=TRUE), split=" ")[[1]]
		n_lines <- as.integer(n_lines[length(n_lines)-1])
		if(n_lines > 1)
		{
			pheno_summary_codings_included_file <- paste0(QCed_io_name, "_cat_variables_females_phesant_recodings.", i, "_phenosummary.tsv")
			pheno_summary_codings_included_file_sex_specific <- paste0(QCed_io_name, "_cat_variables_females_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")

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
	}

}
