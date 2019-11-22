# Quick and dirty code to clean up the PHESANT data when restricting to men and women

library(data.table)

root_file <- "ukb11214_final_january_reference_QC_more_phenos_and_corrected."
root_file <- "ukb11214_final_july_reference_QC_more_phenos_and_corrected."
# root_file_old <- "ukb11214_final_QC_more_phenos_and_corrected."
root_file_old <- "new_and_corrected_PHESANT_phenos/ukb11214_rw_restricted_to_QC_more_phenos_and_corrected."

# Need the phenotype summary file as well.
# Use run_summarise_phenotypes_cloud.r if these have not yet been created.

# These are just used to pull out the males and females.
# root_file_males <- "ukb11214_final_QC_males_more_phenos_and_corrected."
# root_file_females <- "ukb11214_final_QC_females_more_phenos_and_corrected."
root_file_males <- "new_and_corrected_PHESANT_phenos/ukb11214_rw_restricted_to_QC_males_more_phenos_and_corrected."
root_file_females <- "new_and_corrected_PHESANT_phenos/ukb11214_rw_restricted_to_QC_females_more_phenos_and_corrected."
# This is for checking I got the right answer when in 'summarise_phenotypes_output_july copy'
root_file_males <- "../new_and_corrected_PHESANT_phenos/ukb11214_rw_restricted_to_QC_males_more_phenos_and_corrected."
root_file_females <- "../new_and_corrected_PHESANT_phenos/ukb11214_rw_restricted_to_QC_females_more_phenos_and_corrected."

phenotype_info_file <- "variable-info/outcome_info_final_round2.tsv"
phenotype_info_file <- "../variable-info/outcome_info_final_round2.tsv"

coding_info_file <- "variable-info/data-coding-ordinal-info.txt"
coding_info_file <- "../variable-info/data-coding-ordinal-info.txt"

minimum_bin <- function(df_column, minimum=100) {
	
	if(any(df_column=="", na.rm=TRUE)){
		df_column <- df_column[-which(df_column=="")]
		if(length(df_column) == 0) return(FALSE)
	}
	
	if(sum(is.na(df_column)) == length(df_column)) return(FALSE)
	result <- table(df_column)

	if(length(result) == 1) {
		# This is the case when everyone is in a single bin.
		return(FALSE)
	} else {
		# Otherwise, ask the size of the smallest bin.
		return(min(result) > minimum)
	}
}

keep_cat_ordered <- function(df_column, minimum=5000) {
	
	if(any(df_column=="", na.rm=TRUE)){
		df_column <- df_column[-which(df_column=="")]
		if(length(df_column) == 0) return(FALSE)
	}

	if(sum(!is.na(df_column)) >= 5000) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

get_raw_cts <- function(cts_variable, file_header, file="neale_lab_parsed_and_restricted_to_QCed_samples.tsv", single=TRUE)
{	
	# Get matches in the header file.
	match <- grep(paste0('x', cts_variable, '_'), header)

	if(single) {
		if(length(match)==1) {
			return(match)
		} else {
			return(c())
		}
	} else {
		if(length(match)==1) {
			return(c())
		} else {
			return(match)
		}
	}
}

reassign_all_values <- function(cts_variables, phenofile) {
	where <- phenofile$FieldID %in% cts_variables
	changes <- cbind(phenofile$FieldID[where], phenofile$DATA_CODING[where])
	changes <- changes[!is.na(changes[,2]),]
	# create a list of variables for each coding
	codings <- list()
	for(i in names(table(changes[,2]))) {
		codings[[i]] <- changes[which(changes[,2]==strtoi(i)),1]
	}
	return(codings)
}

get_reassignment <- function(reassignment) {
	matrix(as.integer(unlist(strsplit(strsplit(reassignment, split='\\|')[[1]], split='='))), ncol=2, byrow=TRUE)
}

make_the_changes <- function(reassignment_matrix, cts_variables, data_frame)
{
	for(i in 1:nrow(reassignment_matrix)) {
		matches <- data_frame[cts_variables] == reassignment_matrix[i,1]
		if(length(matches) > 0)
			data_frame[cts_variables][matches] <- reassignment_matrix[i,2]
	}
	return(data_frame)
}

change_values <- function(codings, data_frame, coding_info_file) {
	# Need to read in the encoding and determine the changes to make
	reassignments <- fread(coding_info_file, sep=',', header=TRUE, data.table=FALSE)
	reassignments <- reassignments[reassignments$dataCode %in% names(codings),]
	reassignment_matrices <- sapply(reassignments$reassignments, get_reassignment)
	names(reassignment_matrices) <- reassignments$dataCode
	codings <- codings[order(names(codings))]
	reassignment_matrices <- reassignment_matrices[order(names(reassignment_matrices))]

	for(i in 1:length(reassignment_matrices)) {
		data_frame <- make_the_changes(reassignment_matrices[[i]], as.character(codings[[i]]), data_frame)
	}

	return(data_frame)
}

average_over_cts_multi <- function(cts_variable, data_frame)
{
	where <- grep(cts_variable, names(data_frame))
	column <- rowMeans(data_frame[,where], na.rm=TRUE)
	column[is.nan(column)] <- NA
	return(column)
}

irnt <- function(cts_variable) {
    set.seed(1234) # This is the same as was used by PHESANT - for checking.
    n_cts <- length(which(!is.na(cts_variable)))
    quantile_cts <- (rank(cts_variable, na.last = "keep", ties.method = "random") - 0.5) / n_cts
    # use the above to check, but also use frank for the real thing
    cts_IRNT <- qnorm(quantile_cts)	
    return(cts_IRNT)
}

look_for_logical <- function(column) {
	return("TRUE" %in% column | "FALSE" %in% column)
}

# Let's read in the phenotype information file
phenofile <- fread(phenotype_info_file, sep='\t', header=TRUE)

# Want to pull out variables that end up as cts variables in the both_sex PHESANT file.
# First, let's get the header.
file <- "neale_lab_parsed_and_restricted_to_QCed_samples.tsv"
header <- strsplit(system(paste0("head -n 1 ", file), intern=TRUE), split='\t')[[1]]

single_cts_columns <- c()
multi_cts_columns <- c()

for(i in 1:4) {
	pheno_summary <- paste0(root_file, i, "_phenosummary.tsv")
	cts_variables <- system(paste("grep IRNT", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)

	# These columns can be written to a new file as is (no need to perform any averaging)
	single_cts_columns <- c(single_cts_columns, unlist(lapply(cts_variables, get_raw_cts, header, single=TRUE)))
	multi_cts_columns <- c(multi_cts_columns, unlist(lapply(cts_variables, get_raw_cts, header, single=FALSE)))
}

# awk out the single cts columns and write them to a file.
# Include the userId
single_cts_columns <- c(1, single_cts_columns)
multi_cts_columns <- c(1, multi_cts_columns)
outfile_single <- "neale_lab_parsed_and_restricted_to_QCed_samples_cts_single.tsv"
system(paste0("awk -F $'\t' -v OFS=$'\t' '{print ", paste0("$", single_cts_columns, collapse=","), "}' ", file, " > ", outfile_single))
outfile_multi <- "neale_lab_parsed_and_restricted_to_QCed_samples_cts_multi.tsv"
system(paste0("awk -F $'\t' -v OFS=$'\t' '{print ", paste0("$", multi_cts_columns, collapse=","), "}' ", file, " > ", outfile_multi))

# Now, create the average columns for the other cts variables...
# Get the names of the variables, and then grep for them and take the average.
cts_multi <- fread(outfile_multi, header=TRUE, sep='\t', data.table=FALSE)
cts_multi_to_average_over <- unique(gsub("_.*", "_", names(cts_multi)[-1]))
codings <- reassign_all_values(substr(cts_multi_to_average_over, 2, nchar(cts_multi_to_average_over)-1), phenofile)
# Need extra step for cts_multi
new_codings <- list()
for(i in 1:length(codings)) {
	for(j in 1:length(codings[[i]])) {
		if(j ==1) {
			new_codings[[i]] <- names(cts_multi)[grep(paste0("x", codings[[i]][j]), names(cts_multi))]
		} else {
			# Changed this c(new_codings[[i]] from c(codings[[i]]...
			new_codings[[i]] <- c(new_codings[[i]], names(cts_multi)[grep(paste0("x", codings[[i]][j]), names(cts_multi))])
		}
	}
}

names(new_codings) <- names(codings)
codings <- new_codings

cts_multi <- change_values(codings, cts_multi, coding_info_file)
cts_multi <- data.frame(userId=cts_multi$userId, sapply(cts_multi_to_average_over, average_over_cts_multi, cts_multi))
names(cts_multi) <- c("userId", gsub("x(.*)_.*", "\\1", cts_multi_to_average_over))

# Now, write these, along with userID, to disk.
cts_columns <- fread(outfile_single, header=TRUE, sep='\t', data.table=FALSE)
names(cts_columns) <- gsub("x([^_]*)_.*", "\\1", names(cts_columns))

# Make the alterations before the transform
codings <- reassign_all_values(names(cts_columns), phenofile)
cts_columns <- change_values(codings, cts_columns, coding_info_file)

cts_columns <- merge(cts_multi, cts_columns, by='userId')

# Pull out the make and female userIds for out application.
males <- system(paste0("awk 'BEGIN { FS=\"\t\" } { print $1 }' ", root_file_males, 1, ".tsv"), intern=TRUE)
females <- system(paste0("awk 'BEGIN { FS=\"\t\" } { print $1 }' ", root_file_females, 1, ".tsv"), intern=TRUE)
males <- males[-1]
females <- females[-1]

# As a check, I want to IRNT the cts variables in all sexes, this should be the same as before.
# Not anymore - as we've made some changes to the encodings and PHESANT!
# The male and female IRNT should be different.
cts_output <- data.frame(sapply(cts_columns, irnt))
cts_output_males <- data.frame(sapply(cts_columns[cts_columns$userId %in% males,], irnt))
cts_output_females <- data.frame(sapply(cts_columns[cts_columns$userId %in% females,], irnt))

names(cts_output) <- names(cts_columns)
names(cts_output_males) <- names(cts_columns)
names(cts_output_females) <- names(cts_columns)

# Finally, need to remove all the cts variables that have less than 5000 instances of non-NAs.
cat(paste("started with", ncol(cts_output), "cts.\n"))
to_keep_all_sexes_cts <- apply(cts_output, 2, keep_cat_ordered)
cat(paste(sum(to_keep_all_sexes_cts), "cts remain for both sexes.\n"))
to_keep_all_sexes_cts[1] <- TRUE
to_keep_males_cts <- apply(cts_output_males, 2, keep_cat_ordered)
cat(paste(sum(to_keep_males_cts), "cts remain for males.\n"))
to_keep_males_cts[1] <- TRUE
to_keep_females_cts <- apply(cts_output_females, 2, keep_cat_ordered)
cat(paste(sum(to_keep_females_cts), "cts remain for females.\n"))
to_keep_females_cts[1] <- TRUE

fwrite(cts_output[,to_keep_all_sexes_cts], file="neale_lab_parsed_and_restricted_to_QCed_samples_cts_irnt.tsv", sep='\t')
fwrite(cts_output_males[, to_keep_males_cts], file="neale_lab_parsed_and_restricted_to_QCed_samples_cts_irnt_males.tsv" , sep='\t')
fwrite(cts_output_females[, to_keep_females_cts], file="neale_lab_parsed_and_restricted_to_QCed_samples_cts_irnt_females.tsv" , sep='\t')

# Now create the raw cts files.
cts_output_raw <- cts_columns
cts_output_males_raw <- cts_columns[cts_columns$userId %in% males,]
cts_output_females_raw <- cts_columns[cts_columns$userId %in% females,]

# Finally, need to remove all the cts variables that have less than 5000 instances of non-NAs.
cat(paste("started with", ncol(cts_output), "cts.\n"))
to_keep_all_sexes_cts_raw <- apply(cts_output_raw, 2, keep_cat_ordered)
cat(paste(sum(to_keep_all_sexes_cts_raw), "cts remain for both sexes.\n"))
to_keep_all_sexes_cts_raw[1] <- TRUE
to_keep_males_cts_raw <- apply(cts_output_males_raw, 2, keep_cat_ordered)
cat(paste(sum(to_keep_males_cts_raw), "cts remain for males.\n"))
to_keep_males_cts_raw[1] <- TRUE
to_keep_females_cts_raw <- apply(cts_output_females_raw, 2, keep_cat_ordered)
cat(paste(sum(to_keep_females_cts_raw), "cts remain for females.\n"))
to_keep_females_cts_raw[1] <- TRUE

fwrite(cts_output_raw[,to_keep_all_sexes_cts_raw], file="neale_lab_parsed_and_restricted_to_QCed_samples_cts_raw.tsv", sep='\t')
fwrite(cts_output_males_raw[, to_keep_males_cts_raw], file="neale_lab_parsed_and_restricted_to_QCed_samples_cts_raw_males.tsv" , sep='\t')
fwrite(cts_output_females_raw[, to_keep_females_cts_raw], file="neale_lab_parsed_and_restricted_to_QCed_samples_cts_raw_females.tsv" , sep='\t')

# Perform the check
for(i in 1:4) {
	all_sexes <- fread(paste0(root_file, i, ".tsv"), header=TRUE, sep='\t', data.table=FALSE)
	to_look <- names(all_sexes)[names(all_sexes) %in% names(cts_output)]
	all_sexes_check <- all_sexes[,names(all_sexes) %in% names(cts_output)]
	for(j in 2:ncol(all_sexes_check)){
		damn <- max(abs(cts_output[,to_look[j]] - all_sexes_check[,j]),na.rm=TRUE)
		if(damn > 1e-14) {
			print(paste(to_look[j], names(all_sexes_check)[j]))
			print(damn)
		}
	}
}


for(i in 1:4)
{
	all_sexes <- fread(paste0(root_file, i, ".tsv"), header=TRUE, sep='\t', data.table=FALSE)
	log_file_both_sex <- paste0(root_file, i, ".log")
	pheno_summary <- paste0(root_file, i, "_phenosummary.tsv")
	cts_variables <- system(paste("grep IRNT", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)

	# Now, look in the PHESANT file for IRNT - and exclude those from the male and females 
	all_sexes <- all_sexes[,-which(names(all_sexes) %in% c(cts_variables, "age", "sex"))]

	PHESANT_males <- all_sexes[all_sexes$userId %in% males,]
	PHESANT_females <- all_sexes[all_sexes$userId %in% females,]

	# Check to make sure that all remaining variables are accounted for - they should be mapped to either Binary, or ordered categorical.
	cat_ordered <- system(paste("grep CAT-ORD", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)
	cat_unordered <- system(paste("grep CAT-SINGLE-BINARY", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)
	cat_mul_unordered <- system(paste("grep CAT-MUL-BINARY-VAR", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)
	edge_cases_cat_unordered <- system(paste("grep 'Combine .* two bins and treat as binary'", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)

	# Check that the grepping for cat-unordered, cat-ordered, and cat-mul-unordered are disjoint.
	if(length(c(cat_unordered, cat_ordered, cat_mul_unordered))!= length(unique(c(cat_unordered, cat_ordered, cat_mul_unordered))))
		cat("Error: non-unique greps for CAT-ORD, CAT-SINGLE-BINARY and CAT-MUL-BINARY-VAR")

	cat_unordered <- c(cat_unordered, cat_mul_unordered, edge_cases_cat_unordered)

	cat_unordered_all_sexes <- all_sexes[,which(names(all_sexes) %in% cat_unordered)]
	cat_unordered_males <- PHESANT_males[,which(names(PHESANT_males) %in% cat_unordered)]
	cat_unordered_females <- PHESANT_females[,which(names(PHESANT_females) %in% cat_unordered)]

	cat_ordered_all_sexes <- all_sexes[,which(names(all_sexes) %in% cat_ordered)]
	cat_ordered_males <- PHESANT_males[,which(names(PHESANT_males) %in% cat_ordered)]
	cat_ordered_females <- PHESANT_females[,which(names(PHESANT_females) %in% cat_ordered)]

	cat(paste("started with", ncol(cat_unordered_all_sexes), "categorical unordered.\n"))
	to_keep_all_sexes <- apply(cat_unordered_all_sexes, 2, minimum_bin)
	cat(paste(sum(to_keep_all_sexes), "categorical unordered remain for both sexes.\n"))
	to_keep_male <-  apply(cat_unordered_males, 2, minimum_bin)
	cat(paste(sum(to_keep_male), "categorical unordered remain for males.\n"))
	to_keep_female <- apply(cat_unordered_females, 2, minimum_bin)
	cat(paste(sum(to_keep_female), "categorical unordered remain for females.\n"))

	# Finally, need to remove ordered variables that are categorical ordered and have less that 5000 (the default in PHESANT) non-NA values.
	cat(paste("started with", ncol(cat_ordered_all_sexes), "categorical ordered.\n"))
	to_keep_all_sexes_ordered <- apply(cat_ordered_all_sexes, 2, keep_cat_ordered)
	cat(paste(sum(to_keep_all_sexes_ordered), "categorical ordered remain for both sexes.\n"))
	to_keep_males_ordered <- apply(cat_ordered_males, 2, keep_cat_ordered)
	cat(paste(sum(to_keep_males_ordered), "categorical ordered remain for males.\n"))
	to_keep_females_ordered <- apply(cat_ordered_females, 2, keep_cat_ordered)
	cat(paste(sum(to_keep_females_ordered), "categorical ordered remain for females.\n"))

	fwrite(cbind(all_sexes$userId, cat_ordered_all_sexes[, to_keep_all_sexes_ordered], cat_unordered_all_sexes[, to_keep_all_sexes]),
		sep='\t', file=paste0('neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_both_sexes.', i, '.tsv'))
	fwrite(cbind(PHESANT_males$userId, cat_ordered_males[, to_keep_males_ordered], cat_unordered_males[, to_keep_male]),
		sep='\t', file=paste0('neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_males.', i, '.tsv'))
	fwrite(cbind(PHESANT_females$userId, cat_ordered_females[, to_keep_females_ordered], cat_unordered_females[, to_keep_female]),
		sep='\t', file=paste0('neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_females.', i, '.tsv'))

	all_sexes_check <- all_sexes[,-c(1,which(names(all_sexes) %in% c(cat_unordered, cat_ordered, cat_mul_unordered)))]
	if(ncol(all_sexes_check)!=0)
		cat("Error: there are variables that have not been accounted for!")
}

# The code below was to determine which subset of the 
# phenotypes should be run again, based on the previous run.

# cts_file <- fread("neale_lab_parsed_and_restricted_to_QCed_samples_cts_irnt.tsv", sep='\t', header=TRUE, data.table=FALSE)

# Determine which variables were not removed, that were removed previously.
# old_vars_males <- c()
# old_vars_females <- c()
# old_vars_both_sexes <- c()
# new_vars_males <- c()
# new_vars_females <- c()
# new_vars_both_sexes <- c()

# for(i in 1:4) {
# 	old_vars_males <- c(old_vars_males, gsub("\"", "", strsplit(system(paste0("head -n 1 ", paste0(root_file_males, i, ".tsv")),
# 		intern=TRUE), split="\t")[[1]][-1]))
# 	old_vars_females <- c(old_vars_females, gsub("\"", "", strsplit(system(paste0("head -n 1 ", paste0(root_file_females, i, ".tsv")),
# 		intern=TRUE), split="\t")[[1]][-1]))
# 	new_vars_males <- c(new_vars_males, strsplit(system(paste0("head -n 1 ", 'neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_males.', i, '.tsv'),
# 		intern=TRUE), split="\t")[[1]][-1])
# 	new_vars_females <- c(new_vars_females, strsplit(system(paste0("head -n 1 ", 'neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_females.', i, '.tsv'),
# 		intern=TRUE), split="\t")[[1]][-1])

# 	old_vars_both_sexes <- c(old_vars_both_sexes, gsub("\"", "", strsplit(system(paste0("head -n 1 ", paste0(root_file_old, i, ".tsv")), 
# 		intern=TRUE), split="\t")[[1]][-1]))
# 	new_vars_both_sexes <- c(new_vars_both_sexes, strsplit(system(paste0("head -n 1 ", 'neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_both_sexes.', i, '.tsv'),
# 		intern=TRUE), split="\t")[[1]][-1])
# }

# # New variables to definitely include (they're new).
# to_include_both_sexes <- setdiff(new_vars_both_sexes, c(old_vars_both_sexes, names(cts_file)))
# to_include_males <- setdiff(new_vars_males, c(old_vars_males, names(cts_file)))
# to_include_females <- setdiff(new_vars_females, c(old_vars_females, names(cts_file)))

# # Old variables to remove (they got deleted - likely too few in CAT-MULTIPLE).
# to_remove_males <- setdiff(old_vars_males, c(new_vars_males, names(cts_file), "age", "sex"))
# to_remove_females <- setdiff(old_vars_females, c(new_vars_females, names(cts_file), "age", "sex"))
# to_remove_both_sexes <- setdiff(old_vars_both_sexes, c(new_vars_both_sexes, names(cts_file), "age", "sex"))

# changed_variables_males <- c()
# changed_variables_females <- c()
# changed_variables_both_sexes <- c()

# for(i in 1:4) {
# 	new_file <- fread(file=paste0('neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_both_sexes.', i, '.tsv'), sep='\t', header=TRUE, data.table=FALSE)
# 	old_file <- fread(paste0(root_file_old, i, ".tsv"), header=TRUE, sep='\t', data.table=FALSE)
# 	old_file <- old_file[,!(names(old_file) %in% c("age", "sex"))]

# 	is_logical <- sapply(old_file, look_for_logical)
# 	old_file[is_logical] <- sapply(old_file[is_logical], as.logical)
# 	old_file[!is_logical] <- sapply(old_file[!is_logical], as.numeric)

# 	index_both <- intersect(names(new_file), names(old_file))
# 	differences <- sapply(abs(new_file[index_both] - old_file[index_both]), max, na.rm=TRUE)
# 	NA_differences <- unique(which(is.na(new_file[index_both]) != is.na(old_file[index_both]), arr.ind=TRUE)[,2])
# 	changed_variables_both_sexes <- c(changed_variables_both_sexes, index_both[NA_differences], index_both[differences > 1e-10])
# 	print(changed_variables_both_sexes)
# 	# Next, get the differences between the male and female specific files.
# 	new_file <- fread(file=paste0('neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_males.', i, '.tsv'), sep='\t', header=TRUE, data.table=FALSE)
# 	old_file <- fread(paste0(root_file_males, i, ".tsv"), header=TRUE, sep='\t', data.table=FALSE)
	
# 	is_logical <- sapply(old_file, look_for_logical)
# 	old_file[is_logical] <- sapply(old_file[is_logical], as.logical)
# 	old_file[!is_logical] <- sapply(old_file[!is_logical], as.numeric)

# 	index_both <- intersect(names(new_file), names(old_file))
# 	differences <- sapply(abs(new_file[index_both] - old_file[index_both]), max, na.rm=TRUE)
# 	NA_differences <- unique(which(is.na(new_file[index_both]) != is.na(old_file[index_both]), arr.ind=TRUE)[,2])
# 	changed_variables_males <- c(changed_variables_males, index_both[NA_differences], index_both[differences > 1e-10])
# 	print(changed_variables_males)
# 	# Next, get the differences between the male and female specific files.
# 	new_file <- fread(file=paste0('neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_females.', i, '.tsv'), sep='\t', header=TRUE, data.table=FALSE)
# 	old_file <- fread(paste0(root_file_females, i, ".tsv"), header=TRUE, sep='\t', data.table=FALSE)
	
# 	is_logical <- sapply(old_file, look_for_logical)
# 	old_file[is_logical] <- sapply(old_file[is_logical], as.logical)
# 	old_file[!is_logical] <- sapply(old_file[!is_logical], as.numeric)

# 	index_both <- intersect(names(new_file), names(old_file))
# 	differences <- sapply(abs(new_file[index_both] - old_file[index_both]), max, na.rm=TRUE)
# 	NA_differences <- unique(which(is.na(new_file[index_both]) != is.na(old_file[index_both]), arr.ind=TRUE)[,2])
# 	changed_variables_females <- c(changed_variables_females, index_both[NA_differences], index_both[differences > 1e-10])
# 	print(changed_variables_females)
# }

# # Finally, want to check what has changed for the cts variables.
# cts_males_file <- fread("neale_lab_parsed_and_restricted_to_QCed_samples_cts_irnt_males.tsv", sep='\t', header=TRUE, data.table=FALSE)
# cts_females_file <- fread("neale_lab_parsed_and_restricted_to_QCed_samples_cts_irnt_females.tsv", sep='\t', header=TRUE, data.table=FALSE)

# # Things that are now continuous that were previously categorical.
# to_include_males_cts <- c(to_include_males, setdiff(names(cts_males_file), old_vars_males))
# to_include_females_cts <- c(to_include_females, setdiff(names(cts_females_file), old_vars_females))

# for(i in 1:4) {
# 	# Next, get the differences between the male and female specific files.
# 	old_file <- fread(paste0(root_file_old, i, ".tsv"), header=TRUE, sep='\t', data.table=FALSE)
# 	index_both <- intersect(names(cts_file), names(old_file))
# 	old_file <- old_file[index_both]
# 	is_logical <- sapply(old_file, look_for_logical)
# 	old_file[is_logical] <- sapply(old_file[is_logical], as.logical)
# 	old_file[!is_logical] <- sapply(old_file[!is_logical], as.numeric)

# 	cts_tmp <- cts_file[index_both]
# 	differences <- sapply(abs(cts_tmp - old_file), max, na.rm=TRUE)
# 	NA_differences <- unique(which(is.na(cts_tmp) != is.na(old_file), arr.ind=TRUE)[,2])
# 	changed_variables_both_sexes <- c(changed_variables_both_sexes, c(names(cts_tmp)[NA_differences], names(cts_tmp)[differences > 1e-10]))

# 	# Next, get the differences between the male and female specific files.
# 	old_file <- fread(paste0(root_file_males, i, ".tsv"), header=TRUE, sep='\t', data.table=FALSE)
# 	index_both <- intersect(names(cts_males_file), names(old_file))
# 	old_file <- old_file[index_both]
# 	is_logical <- sapply(old_file, look_for_logical)
# 	old_file[is_logical] <- sapply(old_file[is_logical], as.logical)
# 	old_file[!is_logical] <- sapply(old_file[!is_logical], as.numeric)

# 	cts_males_tmp <- cts_males_file[index_both]
# 	differences <- sapply(abs(cts_males_tmp - old_file), max, na.rm=TRUE)
# 	NA_differences <- unique(which(is.na(cts_males_tmp) != is.na(old_file), arr.ind=TRUE)[,2])
# 	changed_variables_males <- c(changed_variables_males, c(names(cts_males_tmp)[NA_differences], names(cts_males_tmp)[differences > 1e-10]))

# 	# Next, get the differences between the male and female specific files.
# 	old_file <- fread(paste0(root_file_females, i, ".tsv"), header=TRUE, sep='\t', data.table=FALSE)
# 	index_both <- intersect(names(cts_females_file), names(old_file))
# 	old_file <- old_file[index_both]
# 	is_logical <- sapply(old_file, look_for_logical)
# 	old_file[is_logical] <- sapply(old_file[is_logical], as.logical)
# 	old_file[!is_logical] <- sapply(old_file[!is_logical], as.numeric)

# 	cts_females_tmp <- cts_females_file[index_both]
# 	differences <- sapply(abs(cts_females_tmp - old_file), max, na.rm=TRUE)
# 	NA_differences <- unique(which(is.na(cts_females_tmp) != is.na(old_file), arr.ind=TRUE)[,2])
# 	changed_variables_females <- c(changed_variables_females, c(names(cts_females_tmp)[NA_differences], names(cts_females_tmp)[differences > 1e-10]))

# }

# print(changed_variables_females)
# print(changed_variables_males)
# print(changed_variables_both_sexes)

# changed_variables_males <- changed_variables_males[-which(changed_variables_males == "userId")]
# changed_variables_females <- changed_variables_females[-which(changed_variables_females == "userId")]
# changed_variables_both_sexes <- changed_variables_both_sexes[-which(changed_variables_both_sexes == "userId")]

# to_include_males_cts <- to_include_males_cts[-which(to_include_males_cts == "userId")]
# to_include_females_cts <- to_include_females_cts[-which(to_include_females_cts == "userId")]

# fwrite(list(changed_variables_both_sexes), file='change_variables_both_sexes.tsv', sep='\t')
# fwrite(list(changed_variables_males), file='change_variables_males.tsv',sep='\t')
# fwrite(list(changed_variables_females), file='change_variables_females.tsv', sep='\t')

# fwrite(list(to_include_both_sexes),file='to_include_both_sexes.tsv', sep='\t')
# fwrite(list(to_include_males_cts), file='to_include_males.tsv', sep='\t')
# fwrite(list(to_include_females_cts), file='to_include_females.tsv', sep='\t')

# fwrite(list(to_remove_both_sexes), file='to_remove_both_sexes.tsv', sep='\t')
# fwrite(list(to_remove_males), file='to_remove_males.tsv', sep='\t')
# fwrite(list(to_remove_females), file='to_remove_females.tsv', sep='\t')

