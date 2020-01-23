library(data.table)

source("../summarise_phenotypes.r")

# Only copy if the file isn't there already
if (!file.exists(paste0(filename_root, "1.tsv"))) {
	system(paste0("gsutil cp ", cloud_filename_root,"*tsv.gz ../../"))
	system(paste0("gsutil cp ", cloud_filename_root,"*log ../../"))
	# Decompress the bgzipped files.
	# system(paste0('for i in {1..', n_chunks, '}; do bgzip -d ', filename_root, '$i.tsv.gz ;done'))
	# Decompress the gzipped files.
	system(paste0('for i in {1..', n_chunks, '}; do gzip -d ', filename_root, '$i.tsv.gz ;done'))
}

# Read in all the data table codings
to_read <- paste("../WAS/codings/", dir("../WAS/codings"), sep="")
codings_tables <- list()

for(i in to_read) {
	name <- gsub("../WAS/codings/coding(.*).tsv", "\\1", i)
	codings_tables[[name]] <- read.table(i, header=TRUE, sep='\t', quote="", stringsAsFactors=FALSE)
}

for (i in 1:n_chunks)
{
	# Both sexes
	hist_filename <- paste0(filename_root, i, "_hist")
	pheno_summary <- paste0(filename_root, i, "_phenosummary.tsv")

	filename <- paste0(filename_root, i)
	tsv_filename <- paste0(filename, ".tsv")
	log_filename <- paste0(filename_root, i)
	log_file <- paste0(log_filename, ".log")
	tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
	names(tsv_data)[1] <- "userId"

	outcome_info <- read.table(variable_info, sep='\t', quote="", comment.char="", header=TRUE)
	print(paste0("both sexes ", i))
	summary_file <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=4)
	write.table(summary_file, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
	# Copy the intermediate summary .tsv files across to the bucket
	system(paste0("gsutil cp ", pheno_summary, " ", intermediate_output_location))
	# Copy the PHESANT run log files across to the bucket
	system(paste0("gsutil cp ", log_file, " ", intermediate_output_location))
	# Copy the PHESANT output files across to the bucket
	system(paste0("gsutil cp ", tsv_filename, " ", intermediate_output_location))
	# Copy the intermediate .pdf files across to the bucket
	system(paste0("gsutil cp ", hist_filename, ".pdf ", plot_location))
}
