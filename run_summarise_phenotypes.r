source("summarise_phenotypes.r")

# Running Verneri's application through.
hist_filename <- "~/results/ukb1859_hist"
pheno_summary <- "~/results/ukb1859_phenosummary.tsv"

filename <- "~/results/ukb7127_output"
tsv_filename <- paste(filename, ".tsv", sep="")
log_file <- paste(filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')

qc_data <- read.table("~/results/ukb1859_qc.tsv", header=TRUE)

outcome_info <- read.table("PHESANT/variable-info/outcome_info_final.tsv",
					   sep='\t', quote="", comment.char="", header=TRUE)

samples_for_removal <- as.character(read.table("~/results/w1859_20170726_participantwithdrawallist.csv")$V1)

notes_for_manny_1859 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal)

write.table(notes_for_manny_1859, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

# Running Joel's application through.
hist_filename <- "~/results/ukb1189_hist"
pheno_summary <- "~/results/ukb1189_phenosummary.tsv"

filename <- "~/results/ukb1189_output"
tsv_filename <- paste(filename, ".tsv", sep="")
log_file <- paste(filename, ".log", sep="")
tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')

qc_data <- read.table("~/results/ukb1189_qc.tsv", header=TRUE)

outcome_info <- read.table("PHESANT/variable-info/outcome_info_final.tsv",
					   sep='\t', quote="", comment.char="", header=TRUE)

samples_for_removal <- as.character(read.table("~/results/w1189_20170726_participantwithdrawallist.csv")$V1)

notes_for_manny_1189 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal)

write.table(notes_for_manny_1189, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

pheno_summary <- "~/results/ukb1859_phenosummary.tsv"
pheno_summary_PHESANT <- "~/results/ukb1859_phenosummary_final.tsv"
notes_for_manny_1859_PHESANT_codings_included <- include_PHESANT_reassignment_names(pheno_summary, outcome_info)
write.table(notes_for_manny_1859_PHESANT_codings_included, file=pheno_summary_PHESANT, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

pheno_summary <- "~/results/ukb1189_phenosummary.tsv"
pheno_summary_PHESANT <- "~/results/ukb1189_phenosummary_final.tsv"
notes_for_manny_1189_PHESANT_codings_included <- include_PHESANT_reassignment_names(pheno_summary, outcome_info)
write.table(notes_for_manny_1189_PHESANT_codings_included, file=pheno_summary_PHESANT, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

# Run the ICD10 codes in Joel's application through.
hist_filename <- "~/results/ukb1189_icd10_hist"
pheno_summary <- "~/results/ukb1189_icd10_phenosummary.tsv"

filename <- "~/results/ukb1189_icd10_flags.tsv"
log_file <- FALSE
tsv_data <- read.table(filename, header=TRUE, sep='\t')
names(tsv_data)[1] <- "userId"
names(tsv_data)[-1] <- paste("X41202", names(tsv_data)[-1], sep="_")

qc_data <- read.table("~/results/ukb1189_qc.tsv", header=TRUE)

outcome_info <- read.table("PHESANT/variable-info/outcome_info_final.tsv",
					   sep='\t', quote="", comment.char="", header=TRUE)

samples_for_removal <- as.character(read.table("~/results/w1189_20170726_participantwithdrawallist.csv")$V1)

notes_for_manny_1189_icd10 <- get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal)
write.table(notes_for_manny_1189_icd10, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
