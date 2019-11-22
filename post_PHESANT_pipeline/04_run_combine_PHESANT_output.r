library(data.table)

only_males_file <- "../should_only_be_in_males.tsv"
only_females_file <- "../should_only_be_in_females.tsv"

root_file <- "neale_lab_parsed_and_restricted_to_QCed_samples"

pheno <- fread(paste0(root_file, "_cat_variables_both_sexes.1.tsv"), sep='\t', data.table=FALSE, header=TRUE)

for(i in 2:4) {
	pheno_tmp <- fread(paste0(root_file, "_cat_variables_both_sexes.", i, ".tsv"), sep='\t', data.table=FALSE, header=TRUE)
	# check 
	print(all(pheno_tmp[,1] == pheno[,1]))
	pheno <- cbind(pheno, pheno_tmp[,-1])
}

# Remove sex specific 
males_only <- fread(only_males_file, sep='\t', data.table=FALSE)$V1
females_only <- fread(only_females_file, sep='\t', data.table=FALSE)$V1

pheno <- pheno[,-which(colnames(pheno) %in% c(males_only, females_only))]
# Next the cts IRNT

pheno_tmp <- fread(paste0(root_file, "_cts_irnt.tsv"), sep='\t', data.table=FALSE, header=TRUE)
if (any(colnames(pheno_tmp) %in% c(males_only, females_only)))
	pheno_tmp <- pheno_tmp[,-which(colnames(pheno_tmp) %in% c(males_only, females_only))]
print(all(order(pheno_tmp[,1]) == order(pheno[,1])))
names(pheno_tmp) <- paste0(names(pheno_tmp), '_irnt')
pheno <- cbind(pheno, pheno_tmp[,-1])

# Finally, cts raw
pheno_tmp <- fread(paste0(root_file, "_cts_raw.tsv"), sep='\t', data.table=FALSE, header=TRUE)
if (any(colnames(pheno_tmp) %in% c(males_only, females_only)))
	pheno_tmp <- pheno_tmp[,-which(colnames(pheno_tmp) %in% c(males_only, females_only))]
print(all(pheno_tmp[,1] == pheno[,1]))
names(pheno_tmp) <- paste0(names(pheno_tmp), '_raw')
pheno <- cbind(pheno, pheno_tmp[,-1])

colnames(pheno)[1] <- 'userID'

# Combine the phenotype summary files (these have already been restricted to sex-specific).
pheno_summary_file <- paste0(root_file, "_cat_variables_both_sexes_phesant_recodings_remove_sex_specific.1_phenosummary.tsv")
pheno_summary <- read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE)
for(i in 2:4) {
	pheno_summary_file <- paste0(root_file, "_cat_variables_both_sexes_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")
	pheno_summary <- rbind(pheno_summary, read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE))
}

pheno_summary_file <- paste0(root_file, "_cts_both_sexes_phenosummary.tsv")
pheno_summary_cts <- read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE)
pheno_summary <- rbind(pheno_summary, read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE))

# Now, check what's in this file that isn't in the summary files...
colnames(pheno)[which(!(colnames(pheno) %in% rownames(pheno_summary)))]
rownames(pheno_summary)[which(!(rownames(pheno_summary) %in% colnames(pheno)))]

fwrite(pheno_summary, file='phesant_output_combined_both_sexes_no_sex_specific_summary.tsv', quote=FALSE, sep='\t', row.names=TRUE)
fwrite(pheno, file='phesant_output_combined_both_sexes_no_sex_specific.tsv', quote=FALSE, sep='\t')

# Males
pheno <- fread(paste0(root_file, "_cat_variables_males.1.tsv"), sep='\t', data.table=FALSE, header=TRUE)

for(i in 2:4) {
	pheno_tmp <- fread(paste0(root_file, "_cat_variables_males.", i, ".tsv"), sep='\t', data.table=FALSE, header=TRUE)
	# check 
	print(all(pheno_tmp[,1] == pheno[,1]))
	pheno <- cbind(pheno, pheno_tmp[,-1])
}

# Next the cts IRNT
pheno_tmp <- fread(paste0(root_file, "_cts_irnt_males.tsv"), sep='\t', data.table=FALSE, header=TRUE)
print(all(order(pheno_tmp[,1]) == order(pheno[,1])))
names(pheno_tmp) <- paste0(names(pheno_tmp), '_irnt')
pheno <- cbind(pheno, pheno_tmp[,-1])

# Finally, cts raw
pheno_tmp <- fread(paste0(root_file, "_cts_raw_males.tsv"), sep='\t', data.table=FALSE, header=TRUE)
print(all(pheno_tmp[,1] == pheno[,1]))
names(pheno_tmp) <- paste0(names(pheno_tmp), '_raw')
pheno <- cbind(pheno, pheno_tmp[,-1])

colnames(pheno)[1] <- 'userID'

# Now, check what's in this file that isn't in the summary files...
pheno_summary_file <- paste0(root_file, "_cat_variables_males_phesant_recodings_remove_sex_specific.1_phenosummary.tsv")
pheno_summary <- read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE)
for(i in 2:4) {
	pheno_summary_file <- paste0(root_file, "_cat_variables_males_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")
	pheno_summary <- rbind(pheno_summary, read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE))
}

pheno_summary_file <- paste0(root_file, "_cts_males_phenosummary.tsv")
pheno_summary_cts <- read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE)
pheno_summary <- rbind(pheno_summary, read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE))

colnames(pheno)[which(!(colnames(pheno) %in% rownames(pheno_summary)))]
rownames(pheno_summary)[which(!(rownames(pheno_summary) %in% colnames(pheno)))]

fwrite(pheno_summary, file='phesant_output_combined_males_summary.tsv', quote=FALSE, sep='\t', row.names=TRUE)
fwrite(pheno, file='phesant_output_combined_males.tsv', quote=FALSE, sep='\t')

# Females
pheno <- fread(paste0(root_file, "_cat_variables_females.1.tsv"), sep='\t', data.table=FALSE, header=TRUE)

for(i in 2:4) {
	pheno_tmp <- fread(paste0(root_file, "_cat_variables_females.", i, ".tsv"), sep='\t', data.table=FALSE, header=TRUE)
	# check 
	print(all(pheno_tmp[,1] == pheno[,1]))
	pheno <- cbind(pheno, pheno_tmp[,-1])
}

# Next the cts IRNT
pheno_tmp <- fread(paste0(root_file, "_cts_irnt_females.tsv"), sep='\t', data.table=FALSE, header=TRUE)
print(all(order(pheno_tmp[,1]) == order(pheno[,1])))
names(pheno_tmp) <- paste0(names(pheno_tmp), '_irnt')
pheno <- cbind(pheno, pheno_tmp[,-1])

# Finally, cts raw
pheno_tmp <- fread(paste0(root_file, "_cts_raw_females.tsv"), sep='\t', data.table=FALSE, header=TRUE)
print(all(pheno_tmp[,1] == pheno[,1]))
names(pheno_tmp) <- paste0(names(pheno_tmp), '_raw')
pheno <- cbind(pheno, pheno_tmp[,-1])

colnames(pheno)[1] <- 'userID'

# Now, check what's in this file that isn't in the summary files...
pheno_summary_file <- paste0(root_file, "_cat_variables_females_phesant_recodings_remove_sex_specific.1_phenosummary.tsv")
pheno_summary <- read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE)
for(i in 2:4) {
	pheno_summary_file <- paste0(root_file, "_cat_variables_females_phesant_recodings_remove_sex_specific.", i, "_phenosummary.tsv")
	pheno_summary <- rbind(pheno_summary, read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE))
}

pheno_summary_file <- paste0(root_file, "_cts_females_phenosummary.tsv")
pheno_summary_cts <- read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE)
pheno_summary <- rbind(pheno_summary, read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE))

colnames(pheno)[which(!(colnames(pheno) %in% rownames(pheno_summary)))]
rownames(pheno_summary)[which(!(rownames(pheno_summary) %in% colnames(pheno)))]

fwrite(pheno_summary, file='phesant_output_combined_females_summary.tsv', quote=FALSE, sep='\t', row.names=TRUE)
fwrite(pheno, file='phesant_output_combined_females.tsv', quote=FALSE, sep='\t')
