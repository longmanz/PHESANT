trim <- function (x) {
	x <- gsub("^\\s+|\\s+$", "", x)
	x <- gsub("^\\\"|\"$", "", x)
	return(x)
}

remove_excess_whitespace <- function(x) x <- gsub("\\s+", " ", x)

# Read in all the data table codings
to_read <- paste("~/Repositories/PHESANT/WAS/codings/",
	dir("~/Repositories/PHESANT/WAS/codings"), sep="")
codings_tables <- list()

for(i in to_read) {
	name <- gsub("PHESANT/WAS/codings/coding(.*).tsv", "\\1", i)
	codings_tables[[name]] <- read.table(i, header=TRUE, sep='\t', quote="", stringsAsFactors=FALSE)
}

get_subtype <- function(x) paste(x[2:length(x)],collapse=" ")

get_hists_and_notes <- function(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data,
	samples_for_removal, samples_for_inclusion=FALSE, check=TRUE, start_column=4)
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
		good_samples <- samples_for_inclusion
	}

	tsv_data <- tsv_data[tsv_data$userId %in% good_samples,]

	if(check == TRUE) {
		if(dim(tsv_data)[1] != 337199) {
			return("failed - not the correct number of samples after removal of redacted samples and individuals who removed consent")
		}
	}

	if (ncol(tsv_data) > (start_column-1)) {
		pdf(file=paste(hist_filename,".pdf",sep=""), width=5, height=5)
		par(oma=c(4,0,0,0))
		# Create character matrix 'notes' that we will write to disk and pass to Manny.
		notes <- matrix(nrow=(ncol(tsv_data)-3),ncol=8)
		colnames(notes) <- c("Field", "N.non.missing", "N.missing", "N.cases", "N.controls",
			"Notes", "PHESANT.notes", "PHESANT.reassignments")
		# Rownames are the FieldIDs
		rownames(notes) <- colnames(tsv_data)[start_column:ncol(tsv_data)]
		samples <- nrow(tsv_data)
		k <- 1
		for(i in colnames(tsv_data)[start_column:ncol(tsv_data)]){
			type <- class(tsv_data[i][,1])
			if(type == "numeric") {
				var <- substr(i,2,nchar(i))
				where <- which(outcome_info$FieldID == var)
				i_name <- paste(trim(outcome_info$Field[where]),"-",i)
				# The first column is the Field name.
				notes[k,1] <- trim(outcome_info$Field[where])
				# The sixth column is the 'Notes' field in variable-info file:
				notes[k,6] <- remove_excess_whitespace(as.character(trim(outcome_info$Notes[where])))
				
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

				hist(tsv_data[i][,1], main=i_name, col="grey", breaks=100, xlab="value")
			} else if (type=="logical" | type=="integer") {
				var <- strsplit(i, split="_")[[1]][1]
				var <- substr(var, 2, nchar(var))
				subvar <- strsplit(i, split="_")[[1]][2]
				where <- which(outcome_info$FieldID == var)
				
				# Add the notes:
				# The sixth column is the 'notes' field in variable-info file:
				notes[k,6] <- remove_excess_whitespace(as.character(trim(outcome_info$Notes[where])))
				
				if (log_file != FALSE) {
					# The seventh column is information about how the data is parsed using PHESANT:
					matching_line <- trim(grep(paste('^', var, '_', sep=""), readLines(log_file), value=TRUE))
					
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
						notes[k,7] <- trim(matching_line)
					}

					# The eighth column is usually empty, unless there are PHESANT reassignments, in which
					# case we detail those reassignments here:
					if (length(grep("reassignments", matching_line)) > 0) {
						notes[k,8] <- trim(gsub("^.*(reassignments: .*?)\\|\\|.*", "\\1", matching_line))
					}
				}

				i_name <- paste(trim(outcome_info$Field[where]), "-", i)

				i_subname <- as.character(outcome_info$Coding[where])
				if(nchar(i_subname)>0 && !is.na(subvar)) {
					to_parse <- strsplit(i_subname, split="\\|")[[1]]
					if(length(to_parse) > 1) {
						where_subvar <- which(sapply(strsplit(trim(to_parse)," "), "[[",1) == subvar) 
						i_subname <- ifelse(length(where_subvar) > 0,
							get_subtype(strsplit(trim(to_parse[where_subvar])," ")[[1]]),
							"PHESANT recoding")
						# The first column is the Field name - if we've entered this if statement,
						# we need to include the subfield (as the categorical has been split into 
						# loads of booleans).
						notes[k,1] <- paste(trim(outcome_info$Field[where]),": ",i_subname, sep="")
					} else {
						notes[k,1] <- trim(outcome_info$Field[where])
					}
				} else {
					# The first column is the Field name.
					notes[k,1] <- trim(outcome_info$Field[where])
				}
				colour <- "grey"
				if(length(grep("Too many", i_subname))>0) {
					# Then we need to look in the corresponding tsv_data-coding table to get the label.
					coding <- trim(gsub("^.*id=","",i_subname))
					where_coding <- which(codings_tables[coding][[1]]$coding == subvar)
					i_subname <- codings_tables[coding][[1]]$meaning[where_coding]
					# The first column is the Field name - if we've entered this if statement,
					# we need to include the subfield (as the categorical has been split into 
					# loads of booleans) - this is accessed from the associated coding table for the 
					# phenotype.
					notes[k,1] <- paste(trim(outcome_info$Field[where]),": ",i_subname, sep="")
				}
				y <- table(tsv_data[i][,1])
				main <- ifelse(type=="logical", paste(i_name,"\n",i_subname), i_name)
				xx <- barplot(height = y, main=main, ylim=c(0, 1.1*max(y)))
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
			k <- k+1
		}
		dev.off()
		# Get rid of the X that's appended at the start of each variable.
		rownames(notes) <- substr(rownames(notes), 2, nchar(rownames(notes)))
		return(notes)
	}
}

include_PHESANT_reassignment_names <- function(pheno_summary, outcome_info)
{	
	pheno_summary <- read.table(pheno_summary, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE) 

	# Find the variables that have had a reassignment.
	where_reassignment <- which(!is.na(pheno_summary[,8]))
	reassignments <- as.character(pheno_summary[where_reassignment,8])
	
	reassignments_before_after_list <- list()
	reassignments <- gsub("reassignments: ", "", reassignments)
	reassignments_list <- strsplit(reassignments, split="\\|")
	for(i in 1:length(reassignments_list)) {
		reassignments_before_after_list[[i]] <- matrix(unlist(strsplit(reassignments_list[[i]], "=")),ncol=2, byrow=TRUE)
	}

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

					i_subname <- as.character(outcome_info$Coding[where_var])
					subvar <- reassignments_before_after_list[[i]][j,1]

					if(nchar(i_subname)>0 && !is.na(subvar)) {
						to_parse <- strsplit(i_subname, split="\\|")[[1]]
						if(length(to_parse) > 1) {
							where_subvar <- which(sapply(strsplit(trim(to_parse)," "), "[[",1) == subvar) 
							i_subname <- ifelse(length(where_subvar) > 0,
								get_subtype(strsplit(trim(to_parse[where_subvar])," ")[[1]]),
								"PHESANT recoding")
							
							pheno_summary[where_reassignment[i],1] <- gsub("PHESANT recoding", i_subname, pheno_summary[where_reassignment[i],1])
						}
					}
				
					if(length(grep("Too many", i_subname))>0) {
						# Then we need to look in the corresponding tsv_data-coding table to get the label.
						coding <- trim(gsub("^.*id=","",i_subname))
						where_coding <- which(codings_tables[coding][[1]]$coding == reassignments_before_after_list[[i]][j,1])
						i_subname <- codings_tables[coding][[1]]$meaning[where_coding]

						pheno_summary[where_reassignment[i],1] <- paste(pheno_summary[where_reassignment[i],1], i_subname, sep="")
					}
				}
			}
		}
	}
	return(pheno_summary)
}
