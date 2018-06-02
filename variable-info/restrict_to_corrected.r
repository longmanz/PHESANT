library(data.table)

round2_variable_info <- "outcome_info_final_round2.tsv"
new_and_fixed_variable_info <- "new_and_fixed_outcome_info_final_round2.tsv"

corrected_phenos <- c("20492",
	"4935", "4946", "4957", "4968", "4979", "4990", "5001", "5012",
	"5556", "5699", "5779", "5790", "5866", "20165","20167", "20169", "20171","20173",
	"20175","20177","20179","20181","20183","20185","20187","20189","20193",
	"20240",
	"20242", "20244", "20245",
	"20526", "20527", "20528", "20529", "20531", "20530")

new_phenotypes <- fread("new_phenotypes_may_2018_excluded.tsv", sep='\t', header=TRUE, data.table=FALSE)

new_phenotypes <- as.character(new_phenotypes$FieldID[which(new_phenotypes$EXCLUDED=="")])
phenotypes_to_run <- unique(c(corrected_phenos, new_phenotypes))

# Read in the latest variable info file and set everything that isn't in this list to NOT-NEW-OR-CORRECTED.
df <- fread(round2_variable_info, header=TRUE, data.table=FALSE)
df$EXCLUDED[!(as.character(df$FieldID) %in% phenotypes_to_run)] <- "NOT-NEW-OR-CORRECTED"

# Write the result to a new file.
fwrite(df, sep='\t', file = new_and_fixed_variable_info)