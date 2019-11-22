library(hyperSpec)

trim <- function (x) {
	x <- gsub("^\\s+|\\s+$", "", x)
	x <- gsub("^\\\"|\"$", "", x)
	return(x)
}

remove_excess_whitespace <- function(x) x <- gsub("\\s+", " ", x)

remove_excess_quotes <- function(x) x <- gsub("\"+", "\"", x)

get_subtype <- function(x) paste(x[2:length(x)],collapse=" ")

get_barplot_numbers <- function(tsv_data, log_file, outcome_info, codings_tables, start_column=4)
{	
	good_samples <- tsv_data$userId
	tsv_data <- tsv_data[tsv_data$userId %in% good_samples,]

	print(dim(tsv_data))

	if (ncol(tsv_data) > (start_column-1)) {
		# Create character matrix 'notes' that we will write to disk.
		notes <- matrix(nrow=(ncol(tsv_data) - start_column+1), ncol=4)
		colnames(notes) <- c("Field", "Max.Category", "Min.Category", "Dist.Category")
		# Rownames are the FieldIDs
		rownames(notes) <- colnames(tsv_data)[start_column:ncol(tsv_data)]
		samples <- nrow(tsv_data)
		
		k <- 1
		for(i in colnames(tsv_data)[start_column:ncol(tsv_data)]){
			
			type <- class(tsv_data[i][,1])

			if (type=="logical" | type=="integer")
			{
				# Get the variable name.
				var <- strsplit(i, split="_")[[1]][1]
				# Remove the X.
				var <- substr(var, 2, nchar(var))
				# Subvariable, if it exists.
				subvar <- strsplit(i, split="_")[[1]][2]
				where <- which(outcome_info$FieldID == var)

				i_name <- paste(trim(outcome_info$Field[where]),"-", gsub("X", "", i))
				coding <- as.character(outcome_info$Coding[where])

				if(nchar(coding)>0 && !is.na(subvar)) {

					if(is.null(codings_tables[coding][[1]]))
						print(paste("Error: Data coding table for coding:", coding, "not found!"))

					where_coding <- which(codings_tables[coding][[1]]$coding == subvar)

					i_subname <- ifelse(
						length(where_coding) > 0, 
						codings_tables[coding][[1]]$meaning[where_coding],
						"PHESANT recoding")

					notes[k,1] <- paste(trim(outcome_info$Field[where]),": ",i_subname, sep="")
				} else {
					# The first column is the Field name.
					notes[k,1] <- trim(outcome_info$Field[where])
				}

				y <- table(tsv_data[i][,1])
				min_y <- min(y)
				max_y <- max(y)
				dist_y <- paste(y, collapse="|")

				notes[k,2] <- min_y
				notes[k,3] <- max_y
				notes[k,4] <- dist_y

			} else {
				print("Error: not one of the specified types!")
				print(type)
			}

			k <- k+1
		}

		rownames(notes) <- substr(rownames(notes), 2, nchar(rownames(notes)))
		return(notes[-Reduce(intersect, list(which(is.na(notes[,2])), which(is.na(notes[,3])), which(is.na(notes[,4])))),])
	}
}

get_hists_and_notes <- function(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data=FALSE,
	samples_for_removal=c(), samples_for_inclusion=FALSE, check=TRUE, start_column=4)
{	
	if(samples_for_inclusion == FALSE) {
		# First, let's restrict to the samples that we want to parse.
		# - in.white.British.ancestry.subset==1
		# - used.in.pca.calculation==1
		# - excess.relatives==0
		# - putative.sex.chromosome.aneuploidy==0
		# ...this should leave you with 337208 samples

		where_good_samples <- which(qc_data$in.white.British.ancestry.subset==1 &
							 	     qc_data$used.in.pca.calculation == 1 &
							 	     qc_data$excess.relatives == 0 & 
							 	     qc_data$putative.sex.chromosome.aneuploidy == 0)
		# Check
		if(length(where_good_samples) != 337208) {
			return("failed")
		}

		# Check that after removal of redacted samples, we have 337205.
		where_redacted_samples <- which(qc_data$eid %in% c("-1", "-2", "-3"))
		if(length(where_redacted_samples) != 3) {
			return("failed - couldn't find all 3 redacted samples")
		}

		# Check after final removal of the last 6 samples, we have 337199.
		where_samples_for_removal <- which(qc_data$eid %in% samples_for_removal)

		where_good_samples <- setdiff(where_good_samples, where_redacted_samples)
		where_good_samples <- setdiff(where_good_samples, where_samples_for_removal)

		if(length(where_good_samples) != 337199) {
			return("failed - not the correct number of samples after removal of redacted samples and individuals who removed consent")
		}

		good_samples <- qc_data$eid[where_good_samples]

	} else {
		if(samples_for_inclusion == TRUE) {
			# Include all the samples.
			good_samples <- tsv_data$userId
		} else {
			good_samples <- samples_for_inclusion
		}
	}

	tsv_data <- tsv_data[tsv_data$userId %in% good_samples,, drop=FALSE]

	if(check == TRUE) {
		if(dim(tsv_data)[1] != 337199) {
			return("failed - not the correct number of samples after removal of redacted samples and individuals who removed consent")
		}
	}

	if (ncol(tsv_data) > (start_column-1)) {
		pdf(file=paste(hist_filename,".pdf",sep=""), width=8, height=5)
		par(oma=c(4,0,0,0))
		# Create character matrix 'notes' that we will write to disk and pass to Manny.
		notes <- matrix(nrow=(ncol(tsv_data)-start_column+1),ncol=9)
		colnames(notes) <- c("Field", "N.non.missing", "N.missing", "N.cases", "N.controls",
			"Notes", "PHESANT.notes", "PHESANT.reassignments", "warning.for.case.control")
		# Rownames are the FieldIDs
		rownames(notes) <- colnames(tsv_data)[start_column:ncol(tsv_data)]
		samples <- nrow(tsv_data)
		k <- 1
		for(i in colnames(tsv_data)[start_column:ncol(tsv_data)]){
			type <- class(tsv_data[i][,1])
			if(type == "numeric") {
				print("numeric")
				# Get the variable name.
				var <- substr(i,2,nchar(i))
				where <- which(outcome_info$FieldID == var)
				i_name <- paste(trim(outcome_info$Field[where]),"-", gsub("X", "", i))
				# The first column is the Field name (not the ID).
				notes[k,1] <- trim(outcome_info$Field[where])
				# The sixth column is the 'Notes' field in variable-info file:
				notes[k,6] <- remove_excess_whitespace(as.character(trim(outcome_info$Notes[where])))
				notes[k,6] <- remove_excess_quotes(notes[k,6])
				
				if (log_file != FALSE){
					# The seventh column is information about how the data is parsed using PHESANT:

					matching_line <- trim(grep(paste('^', var, '_', sep=""), readLines(log_file), value=TRUE))
					notes[k,7] <- matching_line

					# The eighth column is usually empty, unless there are PHESANT reassignments, in which
					# case we detail those reassignments here:
					if (length(grep("reassignments", matching_line)) > 0) {
						notes[k,8] <- trim(gsub("^.*(reassignments: .*?)\\|\\|.*", "\\1", matching_line))
					}
				}

				# Create a simple histogram of this continous data.
				main_cex = 1
				if(nchar(i_name) > 75) {
					main_cex = 75/nchar(i_name)
				}
				hist(tsv_data[i][,1], main=i_name, col="grey", breaks=100, xlab="value")

			} else if (type=="logical" | type=="integer") {
				
				print("logical or integer")
				# Get the variable name.
				var <- strsplit(i, split="_")[[1]][1]
				# Remove the X.
				var <- substr(var, 2, nchar(var))
				# Subvariable, if it exists.
				subvar <- strsplit(i, split="_")[[1]][2]
				where <- which(outcome_info$FieldID == var)
				print(where)
				
				# Add the notes:
				# The sixth column is the 'notes' field in variable-info file:
				notes[k,6] <- remove_excess_whitespace(as.character(trim(outcome_info$Notes[where])))
				notes[k,6] <- remove_excess_quotes(notes[k,6])
				
				if (log_file != FALSE) {
					# The seventh column is information about how the data is parsed using PHESANT:
					matching_line <- trim(grep(paste('^', var, '_', sep=""), readLines(log_file), value=TRUE))
					
					# Do some cleaning up.
					if (length(grep("^.*CAT-MULTIPLE", matching_line)) > 0) {
						new_matching_line <- paste(gsub(paste("^(.*?)\\|.*(CAT-MUL-BINARY-VAR", subvar, ".*?)CAT.*"), "\\1 \\|\\| \\2", matching_line))
						if(matching_line == new_matching_line){
							new_matching_line <- paste(gsub(paste("^(.*?)\\|.*(CAT-MUL-BINARY-VAR", subvar, ".*)"), "\\1 \\|\\| \\2", matching_line))
						}
						notes[k,7] <- trim(new_matching_line)
					} else if (length(grep("CAT-SINGLE-UNORDERED", matching_line)) > 0) {
						new_matching_line <- gsub(paste("^(.*?)\\|.*(Inc.*?: ", subvar, "\\([0-9]+\\)).*", sep=""),
							paste("\\1 \\|\\| CAT-SINGLE \\|\\| CAT-SINGLE-BINARY-VAR:", subvar, " \\|\\| \\2 \\|\\|"), matching_line)
						notes[k,7] <- trim(new_matching_line)
					} else {
						print(matching_line)
						notes[k,7] <- trim(matching_line)
					}

					# The eighth column is usually empty, unless there are PHESANT reassignments, 
					# in which case we detail those reassignments here:
					if (length(grep("reassignments", matching_line)) > 0) {
						notes[k,8] <- trim(gsub("^.*(reassignments: .*?)\\|\\|.*", "\\1", matching_line))
					}
				}

				i_name <- paste(trim(outcome_info$Field[where]),"-", gsub("X", "", i))
				coding <- as.character(outcome_info$Coding[where])

				if(nchar(coding)>0 && !is.na(subvar)) {

					if(is.null(codings_tables[coding][[1]]))
						print(paste("Error: Data coding table for coding:", coding, "not found!"))

					where_coding <- which(codings_tables[coding][[1]]$coding == subvar)

					i_subname <- ifelse(
						length(where_coding) > 0, 
						codings_tables[coding][[1]]$meaning[where_coding],
						"PHESANT recoding")

					notes[k,1] <- paste(trim(outcome_info$Field[where]),": ",i_subname, sep="")
				} else {
					# The first column is the Field name.
					notes[k,1] <- trim(outcome_info$Field[where])
					i_subname <- ""
				}

				colour <- "grey"
				y <- table(tsv_data[i][,1])
				main <- ifelse(type=="logical", paste(i_name,"\n",i_subname), i_name)
				main_cex = 1
				if( (nchar(i_name) > 75) | (nchar(i_subname) > 75)) {
					main_cex = 75/max(nchar(i_name), nchar(i_subname))
				}
				xx <- barplot(height = y, main=main, ylim=c(0, 1.1*max(y)), cex.main=main_cex)
				text(x = xx, y = y, label = y, pos = 3, cex = 0.8, col = "red")

				# Finally, include the case and control numbers for the logical and integer phenotypes
				# with just two categories.
				if(type=="logical") {
					notes[k,4] <- y[names(y) == "TRUE"]
					notes[k,5] <- y[names(y) == "FALSE"]
				} else {
					if(length(y) == 2) {
						# I assume that 1 encodes a positive here:
						if(all(names(y) == c("0", "1"))) {
							notes[k,4] <- y[names(y) == "1"]
							notes[k,5] <- y[names(y) == "0"]
						} else {
							# All bets are off if the variable is not 0,1, 
							# but I report the numbers in each of the two categories.
							# This should be the larger number encoding the positive, according to the 
							# laws of the table function.
							notes[k,4] <- y[1]
							notes[k,5] <- y[2]
							notes[k,9] <- "YES"
						}
					}
				}

			} else {
				print("Error: not one of the specified types!")
			}
			# The third column is the number of missing data for this phenotype.
			n_miss <- sum(is.na(tsv_data[i][,1]))
			notes[i,3] <- n_miss
			# The second column is the number of non-missing data for this phenotype.
			n_non_miss <- sum(!is.na(tsv_data[i][,1]))
			notes[i,2] <- n_non_miss
			p <- n_miss/samples
			colour <- ifelse(p > 0.95, "red", "darkgreen")
			mtext(side=1, text=paste("Missing = ", n_miss, ", Proportion = ", round(p, 2), sep=""), line=0, outer=TRUE, col=colour)
			if(type=="integer") mtext(side=1, text=i_subname, line=-1, outer=TRUE)
			# print(i_subname)
			k <- k+1
		}
		dev.off()
		# Get rid of the X that are prepended to the start of each variable.
		rownames(notes) <- substr(rownames(notes), 2, nchar(rownames(notes)))
		return(notes)
	}
}

include_PHESANT_reassignment_names <- function(pheno_summary_file, outcome_info)
{	
	pheno_summary <- read.table(pheno_summary_file, sep='\t', comment.char="", quote="", header=TRUE, stringsAsFactors=FALSE) 

	# Find the variables that have had a reassignment.
	where_reassignment <- which(!is.na(pheno_summary[,8]))
	reassignments <- as.character(pheno_summary[where_reassignment,8])
	
	reassignments_before_after_list <- list()
	reassignments <- gsub("reassignments: ", "", reassignments)
	reassignments_list <- strsplit(reassignments, split="\\|")

	if(length(reassignments_list) > 0) {
		for(i in 1:length(reassignments_list)) {
			reassignments_before_after_list[[i]] <- matrix(unlist(strsplit(reassignments_list[[i]], "=")), ncol=2, byrow=TRUE)
		}
	}

	if(length(where_reassignment) > 0) {
		for(i in 1:length(where_reassignment)) {
			if(length(grep("_", rownames(pheno_summary)[where_reassignment[i]])) > 0) {
				for(j in 1:nrow(reassignments_before_after_list[[i]])) {
					if(length(grep(paste("_", reassignments_before_after_list[[i]][j,2], sep=""),
						rownames(pheno_summary)[where_reassignment[i]])) > 0) {
						# We've found this recoded variable in this row.
						# I want to now find what this codes for - so need to find the coding table for this variable
						# from the outcome-info file, then look at what the uncoded version of the phenotype is.
						var <- strsplit(rownames(pheno_summary)[where_reassignment[i]], split="_")[[1]][1]
						where_var <- which(outcome_info$FieldID == var)

						coding <- as.character(outcome_info$Coding[where_var])
						subvar <- reassignments_before_after_list[[i]][j,1]

						if(nchar(coding)>0 && !is.na(subvar)) {

							if(is.null(codings_tables[coding][[1]]))
								print(paste("Error: Data coding table for coding:", coding, "not found!"))

							where_coding <- which(codings_tables[coding][[1]]$coding == subvar)

							i_subname <- ifelse(length(where_coding) > 0, 
								codings_tables[coding][[1]]$meaning[where_coding],
								"PHESANT recoding")
							print(pheno_summary[where_reassignment[i],1])
							print(gsub("PHESANT recoding", i_subname, pheno_summary[where_reassignment[i],1]))
							pheno_summary[where_reassignment[i],1] <- gsub("PHESANT recoding", i_subname, pheno_summary[where_reassignment[i],1])
							
						}				
					}
				}
			}
		}
	}
	return(pheno_summary)
}
