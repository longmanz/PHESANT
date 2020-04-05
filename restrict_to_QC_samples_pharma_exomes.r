library('data.table')

n_exomes <- '300K'

df <- fread("pharma_exomes_parsed.tsv", header=TRUE, sep='\t', na.strings=c("NA", ""), data.table=FALSE)
df_biomarkers <- fread("pharma_biomarkers_parsed.tsv", header=TRUE, sep='\t', na.strings=c("NA", ""), data.table=FALSE)
df_biomarkers <- df_biomarkers[,c(which(names(df_biomarkers)=="userId"), grep("^x30", names(df_biomarkers)))]

if(n_exomes == '100K') {
	# 100k exomes
	df_restrict <- fread("Project_26041_bridge.csv", header=TRUE, data.table=FALSE)
	# For the 100k we actually ran biomarkers and the rest separately as the biomarkers were added
	# later, but we should merge them before passing to PHESANT
} else if (n_exomes == '200K') {
	# 200k exomes
	df_restrict <- fread("linking_file_200K_withbatch.csv", header=TRUE, data.table=FALSE)
} else if (n_exomes == '300K') {
	# 300k exomes
	df_restrict <- fread("linking_file_300K_withbatch.csv", header=TRUE, data.table=FALSE)
} else {
	cat('number of exomes not recognised!\n'))
}

df_restrict <- df_restrict$eid_26041
df_restrict <- data.frame(userId = df_restrict)

print(nrow(df))
df <- merge(df, df_restrict, by='userId')
df <- merge(df, df_biomarkers, by='userId')

# Check
print(nrow(df_restrict))
print(nrow(df))

fwrite(df, file=paste0("pharma_parsed_and_restricted_to_", n_exomes, "_sample_subset.tsv"), row.names=FALSE, quote=FALSE, sep='\t')
