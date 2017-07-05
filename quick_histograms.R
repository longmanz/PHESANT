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
	name <- gsub("codings/coding(.*).tsv", "\\1", i)
	codings_tables[[name]] <- read.table(i, header=TRUE, sep='\t', quote="", stringsAsFactors=FALSE)
}

get_subtype <- function(x) paste(x[2:length(x)],collapse=" ")

get_hists_and_notes <- function(hist_filename, tsv_data, log_file, outcome_info, codings_tables)
{
	if (ncol(tsv_data) > 4) {
		pdf(file=paste(hist_filename,".pdf",sep=""), width=5, height=5)
		par(oma=c(4,0,0,0))

		notes <- matrix(nrow=(ncol(tsv_data)-4),ncol=3)
		rownames(notes) <- colnames(tsv_data)[5:ncol(tsv_data)]
		samples <- nrow(tsv_data)
		k <- 1
		for(i in colnames(tsv_data)[5:ncol(tsv_data)]){
			type <- class(tsv_data[i][,1])
			if(type == "numeric") {
				var <- substr(i,2,nchar(i))
				where <- which(outcome_info$FieldID == var)
				i_name <- paste(trim(outcome_info$Field[where]),"-",i)
				# Add the notes:
				notes[k,1] <- remove_excess_whitespace(as.character(trim(outcome_info$Notes[where])))
				matching_line <- trim(grep(paste('^', var, '_', sep=""), readLines(log_file), value=TRUE))
				notes[k,2] <- matching_line

				if (length(grep("reassignments", matching_line)) > 0) {
					notes[k,3] <- trim(gsub("^.*(reassignments: .*?)\\|\\|.*", "\\1", matching_line))
				}

				hist(tsv_data[i][,1], main=i_name, col="grey", breaks=100, xlab="value")
			} else if (type=="logical" | type=="integer") {
				var <- strsplit(i, split="_")[[1]][1]
				var <- substr(var, 2, nchar(var))
				subvar <- strsplit(i, split="_")[[1]][2]
				where <- which(outcome_info$FieldID == var)
				# Add the notes:
				notes[k,1] <- remove_excess_whitespace(as.character(trim(outcome_info$Notes[where])))
				matching_line <- trim(grep(paste('^', var, '_', sep=""), readLines(log_file), value=TRUE))
				if (length(grep("^.*CAT-MULTIPLE", matching_line)) > 0) {
					new_matching_line <- paste(gsub(paste("^(.*?)\\|.*(CAT-MUL-BINARY-VAR", subvar, ".*?)CAT.*"), "\\1 \\|\\| \\2", matching_line))
					if(matching_line == new_matching_line){
						new_matching_line <- paste(gsub(paste("^(.*?)\\|.*(CAT-MUL-BINARY-VAR", subvar, ".*)"), "\\1 \\|\\| \\2", matching_line))
					}
					notes[k,2] <- trim(new_matching_line)
				} else if (length(grep("CAT-SINGLE-UNORDERED", matching_line)) > 0) {
					new_matching_line <- gsub(paste("^(.*?)\\|.*(Inc.*?: ", subvar, "\\([0-9]+\\)).*", sep=""),
						paste("\\1 \\|\\| CAT-SINGLE \\|\\| CAT-SINGLE-BINARY-VAR:", subvar, " \\|\\| \\2 \\|\\|"), matching_line)
					notes[k,2] <- trim(new_matching_line)
				} else {
					notes[k,2] <- trim(matching_line)
				}

				if (length(grep("reassignments", matching_line)) > 0) {
					notes[k,3] <- trim(gsub("^.*(reassignments: .*?)\\|\\|.*", "\\1", matching_line))
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
					}
				}
				colour <- "grey"
				if(length(grep("Too many", i_subname))>0) {
					# Then we need to look in the corresponding tsv_data-coding table to get the label.
					coding <- trim(gsub("^.*id=","",i_subname))
					where_coding <- which(codings_tables[coding][[1]]$coding == subvar)
					i_subname <- codings_tables[coding][[1]]$meaning[where_coding]
				}
				y <- table(tsv_data[i][,1])
				main <- ifelse(type=="logical", paste(i_name,"\n",i_subname), i_name)
				xx <- barplot(height = y, main=main, ylim=c(0, 1.1*max(y)))
				text(x = xx, y = y, label = y, pos = 3, cex = 0.8, col = "red")
			} else {
				print("Error: not one of the specified types!")
			}
			n_miss <- sum(is.na(tsv_data[i][,1]))
			p <- n_miss/samples
			colour <- ifelse(p > 0.95, "red", "darkgreen")
			mtext(side=1, text=paste("Missing = ", n_miss, ", Proportion = ", round(p, 2), sep=""), line=0, outer=TRUE, col=colour)
			if(type=="integer") mtext(side=1, text=i_subname, line=-1, outer=TRUE)
			k <- k+1
		}
		dev.off()
		return(notes)
	}
}

tsv_filename <- "~/results/ukb9403"
hist_filename <- "~/results/ukb9403"

for(i in 1:20) {
	tsv_filename_i <- paste(tsv_filename, ".", i, ".tsv", sep="")
	log_file <- paste(tsv_filename, ".", i, ".log", sep="")
	tsv_data <- read.table(tsv_filename_i, header=TRUE, sep='\t')
	# tsv_data <- as.data.frame(tsv_data)
	
	outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final.tsv",
						   sep='\t', quote="", comment.char="", header=TRUE)

	get_hists_and_notes(paste(hist_filename,".",i,sep=""), tsv_data, log_file, outcome_info, codings_tables)

}
