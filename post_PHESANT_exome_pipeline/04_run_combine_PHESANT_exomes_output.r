library(data.table)
library(dplyr)

# Move relevant files to the cloud if running out of RAM locally
# for (i in c("both_sexes", "male", "female")) {
# 	for (j in 1:n_chunks) {
# 		system(paste0("gzip ", QCed_io_name, "_cat_variables_", i, ".", j, ".tsv"))
# 	}
# 	system(paste0("gzip ", QCed_io_name, "_cts_irnt.tsv"))
# 	system(paste0("gzip ", QCed_io_name, "_cts_raw.tsv"))

# 	system(paste0("gsutil cp ", QCed_io_name, "_cat_variables_", i, "*gz gs://ukb-diverse-pops/Phenotypes/Everyone/PHESANT_intermediate_output/jan_2020/to_combine/"))
# 	system(paste0("gsutil cp ", QCed_io_name, "_cts_*gz gs://ukb-diverse-pops/Phenotypes/Everyone/PHESANT_intermediate_output/jan_2020/to_combine/"))
# 	system(paste0("gsutil cp ", QCed_io_name, "_cat_variables_", i, "_phesant_recodings_remove_sex_specific", "*_phenosummary.tsv gs://ukb-diverse-pops/Phenotypes/Everyone/PHESANT_intermediate_output/jan_2020/to_combine/"))
# 	system(paste0("gsutil cp ", QCed_io_name, "_cts_", i, "_phenosummary.tsv gs://ukb-diverse-pops/Phenotypes/Everyone/PHESANT_intermediate_output/jan_2020/to_combine/"))
# }
# Now, need to load a VM and move it across. 
# Navigate to the repo/post_PHESANT_multi_ancestry_pipeline/ and run
# gsutil cp gs://ukb-diverse-pops/Phenotypes/Everyone/PHESANT_intermediate_output/jan_2020/to_combine/* ../../
# for(i in c("both_sexes")) {#, "male", "female")) {
# 	for(j in 1:n_chunks) {
# 		system(paste0("gzip -d ", QCed_io_name, "_cat_variables_", i, ".", j, ".tsv.gz"))
# 	}
# }
# system(paste0("gzip -d ", QCed_io_name, "_cts_irnt.tsv.gz"))
# system(paste0("gzip -d ", QCed_io_name, "_cts_raw.tsv.gz"))

# If the file doesn't exists, then in that chunk, no cat phenotypes made it through PHESANT.
pheno_file <- paste0(QCed_io_name, "_cat_variables_both_sexes.1.tsv")
cat('Combining and checking cat variables...\n')
cat('Chunk 1...\n')
if(file.exists(pheno_file))
	pheno <- fread(pheno_file, sep='\t', data.table=FALSE, header=TRUE)


if(n_chunks > 1) {
	for(i in 2:n_chunks)
	{	
		cat(paste0('Chunk ', i, '...\n'))
		pheno_file <- paste0(QCed_io_name, "_cat_variables_both_sexes.", i, ".tsv")
		
		if(!file.exists(pheno_file))
			next
		
		pheno_tmp <- fread(pheno_file, sep='\t', data.table=FALSE, header=TRUE)
		# check 
		cat('Checking:\n')
		cat(paste0(all(pheno_tmp[,1] == pheno[,1]), '\n'))
		pheno <- cbind(pheno, pheno_tmp[,-1, drop=FALSE])
	}
}

# Remove sex specific 
males_only <- fread(only_males_file, sep='\t', data.table=FALSE)$V1
females_only <- fread(only_females_file, sep='\t', data.table=FALSE)$V1

if(exists("pheno"))
	pheno <- pheno[,-which(colnames(pheno) %in% c(males_only, females_only))]

# Next the cts raw
cat('Combining and checking cts raw variables...\n')
if(file.exists(paste0(QCed_io_name, "_cts_raw.tsv")))
{
	pheno_tmp <- fread(paste0(QCed_io_name, "_cts_raw.tsv"), sep='\t', data.table=FALSE, header=TRUE)
	
	if (any(colnames(pheno_tmp) %in% c(males_only, females_only)))
		pheno_tmp <- pheno_tmp[,-which(colnames(pheno_tmp) %in% c(males_only, females_only))]

	if(!exists("pheno")) {
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_raw')
		pheno <- pheno_tmp
	} else {
		cat('Checking:\n')
		cat(paste0(all(order(pheno_tmp[,1] == pheno[,1])), '\n'))
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_raw')
		pheno <- cbind(pheno, pheno_tmp[,-1, drop=FALSE])
	}
}

# Finally, cts IRNT
if(!exists("pheno")) 
	cat("No phenotype output files exist!\n")

cat('Combining and checking cts irnt variables...\n')
if(file.exists(paste0(QCed_io_name, "_cts_irnt.tsv")))
{
	pheno_tmp <- fread(paste0(QCed_io_name, "_cts_irnt.tsv"), sep='\t', data.table=FALSE, header=TRUE)

	if (any(colnames(pheno_tmp) %in% c(males_only, females_only)))
		pheno_tmp <- pheno_tmp[,-which(colnames(pheno_tmp) %in% c(males_only, females_only))]

	if(!exists("pheno")) {
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_irnt')
		pheno <- pheno_tmp
	} else {
		cat('Checking:\n')
		cat(paste0(all(order(pheno_tmp[,1]) == order(pheno[,1])), '\n'))
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_irnt')
		pheno <- cbind(pheno, pheno_tmp[,-1, drop=FALSE])
	}
}

colnames(pheno)[1] <- 'userId'

# Combine the phenotype summary files (these have already been restricted to sex-specific).
cat('Next the summary files...\nCombining and checking cat variables...\n')
cat('Chunk 1...\n')
pheno_summary_file <- paste0(QCed_io_name, "_cat_variables_both_sexes_phesant_recodings_remove_sex_specific.1_phenosummary.tsv")
if(file.exists(pheno_summary_file))
	pheno_summary <- read.table(pheno_summary_file, sep='\t', quote="", comment.char="", header=TRUE, stringsAsFactors=FALSE)

if(n_chunks > 1) {
	for(i in 2:n_chunks)
	{
		cat(paste0('Chunk ', i, '...\n'))
		pheno_summary_file <- paste0(QCed_io_name, "_cat_variables_both_sexes_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")

		if(!file.exists(pheno_summary_file))
			next

		pheno_summary <- rbind(pheno_summary, read.table(pheno_summary_file, sep='\t', quote="", comment.char="", header=TRUE, stringsAsFactors=FALSE))
	}
}

cat('Combining and checking cts variables...\n')
pheno_summary_file <- paste0(QCed_io_name, "_cts_both_sexes_phenosummary.tsv")
if(file.exists(pheno_summary_file)) {
	pheno_summary_cts <- read.table(pheno_summary_file, sep='\t', quote="", comment.char="", header=TRUE, stringsAsFactors=FALSE)
	if(!exists("pheno_summary")) {
		pheno_summary <- pheno_summary_cts
	} else{
		pheno_summary <- rbind(pheno_summary, pheno_summary_cts)	
	}
}

# Now, check what's in this file that isn't in the summary files...
colnames(pheno)[which(!(colnames(pheno) %in% rownames(pheno_summary)))]
rownames(pheno_summary)[which(!(rownames(pheno_summary) %in% colnames(pheno)))]

cat('Writing the results...\n')
fwrite(pheno_summary, file=paste0(final_output, 'combined_both_sexes_no_sex_specific_summary.tsv'), quote=FALSE, sep='\t', row.names=TRUE)
fwrite(pheno, file=paste0(final_output, 'combined_both_sexes_no_sex_specific.tsv'), quote=FALSE, sep='\t')

rm("pheno")
rm("pheno_summary")

# Males
cat('Next for males...\nCombining and checking cat variables...\n')
cat('Chunk 1...\n')
pheno_file <- paste0(QCed_io_name, "_cat_variables_males.1.tsv")

if(file.exists(pheno_file))
	pheno <- fread(pheno_file, sep='\t', data.table=FALSE, header=TRUE)

if(n_chunks > 1) {
	for(i in 2:n_chunks)
	{
		cat(paste0('Chunk ', i, '...\n'))
		pheno_file <- paste0(QCed_io_name, "_cat_variables_males.", i, ".tsv")
		
		if(!file.exists(pheno_file))
			next
		
		pheno_tmp <- fread(pheno_file, sep='\t', data.table=FALSE, header=TRUE)
		# check 
		cat('Checking:\n')
		cat(paste0(all(pheno_tmp[,1] == pheno[,1]), '\n'))
		pheno <- cbind(pheno, pheno_tmp[,-1, drop=FALSE])
	}
}

# Next the cts raw
cat('Combining and checking cts raw variables...\n')
if(file.exists(paste0(QCed_io_name, "_cts_raw_males.tsv")))
{
	pheno_tmp <- fread(paste0(QCed_io_name, "_cts_raw_males.tsv"), sep='\t', data.table=FALSE, header=TRUE)

	if(!exists("pheno")) {
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_raw')
		pheno <- pheno_tmp
	} else {
		cat('Checking:\n')
		cat(paste0(all(order(pheno_tmp[,1]) == order(pheno[,1])), '\n'))
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_raw')
		pheno <- cbind(pheno, pheno_tmp[,-1, drop=FALSE])
	}
}

# Finally, cts IRNT
cat('Combining and checking cts raw variables...\n')
if(file.exists(paste0(paste0(QCed_io_name, "_cts_irnt_males.tsv"))))
{
	pheno_tmp <- fread(paste0(QCed_io_name, "_cts_irnt_males.tsv"), sep='\t', data.table=FALSE, header=TRUE)
	
	if(!exists("pheno")) {
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_irnt')
		pheno <- pheno_tmp
	} else {
		cat('Checking:\n')
		cat(paste0(all(order(pheno_tmp[,1]) == order(pheno[,1])), '\n'))
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_irnt')
		pheno <- cbind(pheno, pheno_tmp[,-1, drop=FALSE])
	}
}

if(exists("pheno"))
	colnames(pheno)[1] <- 'userId'

# Now, check what's in this file that isn't in the summary files...
cat('Next the summary files...\nCombining and checking cat variables...\n')
cat('Chunk 1...\n')
pheno_summary_file <- paste0(QCed_io_name, "_cat_variables_males_phesant_recodings_remove_sex_specific.1_phenosummary.tsv")
if(file.exists(pheno_summary_file))
	pheno_summary <- read.table(pheno_summary_file, sep='\t', quote="", comment.char="", header=TRUE, stringsAsFactors=FALSE)

if(n_chunks > 1) {
	for(i in 2:n_chunks)
	{
		cat(paste0('Chunk ', i, '...\n'))
		pheno_summary_file <- paste0(QCed_io_name, "_cat_variables_males_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")

		if (!file.exists(pheno_summary_file))
			next

		pheno_summary <- rbind(pheno_summary, read.table(pheno_summary_file, sep='\t', quote="", comment.char="", header=TRUE, stringsAsFactors=FALSE))
	}
}

cat('Combining and checking cts variables...\n')
pheno_summary_file <- paste0(QCed_io_name, "_cts_males_phenosummary.tsv")
if(file.exists(pheno_summary_file)) {
	pheno_summary_cts <- read.table(pheno_summary_file, sep='\t', quote="", comment.char="", header=TRUE, stringsAsFactors=FALSE)
	if(!exists("pheno_summary")) {
		pheno_summary <- pheno_summary_cts
	} else{
		pheno_summary <- rbind(pheno_summary, pheno_summary_cts)	
	}
}

colnames(pheno)[which(!(colnames(pheno) %in% rownames(pheno_summary)))]
rownames(pheno_summary)[which(!(rownames(pheno_summary) %in% colnames(pheno)))]

cat('Writing the results...\n')
fwrite(pheno_summary, file=paste0(final_output, 'combined_males_summary.tsv'), quote=FALSE, sep='\t', row.names=TRUE)
fwrite(pheno, file=paste0(final_output, 'combined_males.tsv'), quote=FALSE, sep='\t')

rm("pheno")
rm("pheno_summary")

# Females
cat('Next for females...\nCombining and checking cat variables...\n')
cat('Chunk 1...\n')
pheno_file <- paste0(QCed_io_name, "_cat_variables_females.1.tsv")

if(file.exists(pheno_file))
	pheno <- fread(pheno_file, sep='\t', data.table=FALSE, header=TRUE)

if(n_chunks > 1) {
	for(i in 2:n_chunks)
	{
		cat(paste0('Chunk ', i, '...\n'))
		pheno_file <- paste0(QCed_io_name, "_cat_variables_females.", i, ".tsv")
		
		if(!file.exists(pheno_file))
			next
		
		pheno_tmp <- fread(pheno_file, sep='\t', data.table=FALSE, header=TRUE)
		# check 
		cat('Checking:\n')
		cat(paste0(all(pheno_tmp[,1] == pheno[,1])), '\n')
		pheno <- cbind(pheno, pheno_tmp[,-1, drop=FALSE])
	}
}

# Next the cts raw
cat('Combining and checking cts raw variables...\n')
if(file.exists(paste0(QCed_io_name, "_cts_raw_females.tsv")))
{
	pheno_tmp <- fread(paste0(QCed_io_name, "_cts_raw_females.tsv"), sep='\t', data.table=FALSE, header=TRUE)

	if(!exists("pheno")) {
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_raw')
		pheno <- pheno_tmp
	} else {
		cat('Checking:\n')
		cat(paste0(all(order(pheno_tmp[,1]) == order(pheno[,1])), '\n'))
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_raw')
		pheno <- cbind(pheno, pheno_tmp[,-1, drop=FALSE])
	}
}

# Finally, cts IRNT
cat('Combining and checking cts irnt variables...\n')
if(file.exists(paste0(paste0(QCed_io_name, "_cts_irnt_females.tsv"))))
{
	pheno_tmp <- fread(paste0(QCed_io_name, "_cts_irnt_females.tsv"), sep='\t', data.table=FALSE, header=TRUE)
	
	if(!exists("pheno")) {
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_irnt')
		pheno <- pheno_tmp
	} else {
		cat('Checking:\n')
		cat(paste0(all(order(pheno_tmp[,1]) == order(pheno[,1])), '\n'))
		names(pheno_tmp) <- paste0(names(pheno_tmp), '_irnt')
		pheno <- cbind(pheno, pheno_tmp[,-1, drop=FALSE])
	}
}

if(exists("pheno"))
	colnames(pheno)[1] <- 'userId'

# Now, check what's in this file that isn't in the summary files...
cat('Next the summary files...\nCombining and checking cat variables...\n')
cat('Chunk 1...\n')
pheno_summary_file <- paste0(QCed_io_name, "_cat_variables_females_phesant_recodings_remove_sex_specific.1_phenosummary.tsv")
if(file.exists(pheno_summary_file))
	pheno_summary <- read.table(pheno_summary_file, sep='\t', quote="", comment.char="", header=TRUE, stringsAsFactors=FALSE)

if(n_chunks > 1) {
	for(i in 2:n_chunks)
	{
		cat(paste0('Chunk ', i, '...\n'))
		pheno_summary_file <- paste0(QCed_io_name, "_cat_variables_females_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")

		if (!file.exists(pheno_summary_file))
			next

		pheno_summary <- rbind(pheno_summary, read.table(pheno_summary_file, sep='\t', quote="", comment.char="", header=TRUE, stringsAsFactors=FALSE))
	}
}

cat('Combining and checking cts variables...\n')
pheno_summary_file <- paste0(QCed_io_name, "_cts_females_phenosummary.tsv")
if(file.exists(pheno_summary_file)) {
	pheno_summary_cts <- read.table(pheno_summary_file, sep='\t', quote="", comment.char="", header=TRUE, stringsAsFactors=FALSE)
	if(!exists("pheno_summary")) {
		pheno_summary <- pheno_summary_cts
	} else{
		pheno_summary <- rbind(pheno_summary, pheno_summary_cts)	
	}
}

colnames(pheno)[which(!(colnames(pheno) %in% rownames(pheno_summary)))]
rownames(pheno_summary)[which(!(rownames(pheno_summary) %in% colnames(pheno)))]

cat('Writing the results...\n')
fwrite(pheno_summary, file=paste0(final_output, 'combined_females_summary.tsv'), quote=FALSE, sep='\t', row.names=TRUE)
fwrite(pheno, file=paste0(final_output, 'combined_females.tsv'), quote=FALSE, sep='\t')

rm("pheno")
rm("pheno_summary")

system(paste0("bgzip ", final_output, 'combined_both_sexes_no_sex_specific.tsv'))
system(paste0("bgzip ", final_output, 'combined_males.tsv'))
system(paste0("bgzip ", final_output, 'combined_females.tsv'))

system(paste0("gsutil cp ", final_output, '*combined*gz ', final_output_location, date, "/"))
system(paste0("gsutil cp ", final_output, '*_both_sexes_no_sex_specific_summary.tsv ', final_output_location, date, "/"))
system(paste0("gsutil cp ", final_output, '*_males_summary.tsv ', final_output_location, date, "/"))
system(paste0("gsutil cp ", final_output, '*_females_summary.tsv ', final_output_location, date, "/"))

# system(paste0("gsutil cp ", final_output, '*combined* gs://phenotype_pharma/PHESANT_output/', n_exomes, "/", date, "/"))
