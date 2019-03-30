library(data.table)
exclusion_name <- "YES-NEALELAB-MARCH-19"
# The old version of the outcome info file
old_variable_info <- "variable-info/outcome_info_final_round2.tsv"
# The latest version of the Data Dictionary Showcase - download from UKB if updating.
# wget http://biobank.ctsu.ox.ac.uk/%7Ebbdatan/Data_Dictionary_Showcase.csv
new_variable_info <- "Data_Dictionary_Showcase.csv"
PHEASANT_variable_info <- "variable-info/outcome_info_PHESANT_main.tsv"

old_df <- fread(old_variable_info, header=TRUE, data.table=FALSE)
new_df <- fread(new_variable_info, header=TRUE, data.table=FALSE)
# We only read the following in to makes sure that the overall format is right, and we've convinced ourselves that 
# we've included/removed the right things.
PH_df <- fread(PHEASANT_variable_info, header=TRUE, data.table=FALSE)
PH_df <- PH_df[,!(names(PH_df) %in%
	c("Path", "Category", "Participants", "Items", "Stability", 
	  "ValueType", "Units", "ItemType", "Strata", "Sexed", "Instances", 
	  "Array", "Coding", "Notes", "Link"))]

# Merge the two data frames.

df <- merge(x=old_df, y=new_df, by="FieldID", all=TRUE)

checking_differences <- function(field_to_check, dataframe=df) {
	x <- paste0(field_to_check, '.x')
	y <- paste0(field_to_check, '.y')

	where <- which(df[x] != df[y])

	return(list(where=where, fields=df[where, c('FieldID', 'Field.x', 'Field.y', x,y)]))
}

# Manually checked through, and things have nicely matched up, with the .y variables being updated versions of the .x variables.
# To merge, we can simply remove all instances of the .x variables, replacing them with the .y variables.

# Before we do that, find out which phenotypes have been added, and write this subset to disk to determine which of these to include.
df_added <- df[is.na(df$Field.x),]
df_removed <- df[is.na(df$Field.y),]

# Now write the FieldID and the Field to disk.
names(df_added)[names(df_added) == 'Field.y'] <- 'Field'
names(df_removed)[names(df_removed) == 'Field.x'] <- 'Field'
df_to_check <- rbind(df_added[,c('FieldID', 'Field')], df_removed[,c('FieldID', 'Field')])
fwrite(df_to_check, sep='\t', file='variable-info/new_phenotypes.tsv')

# We then add an 'EXCLUDED' column and fill it in manually.
df <- df[,-grep('\\.x', names(df))]
names(df) <- gsub('\\.y', '', names(df))

# Now, remove all the column names that are no longer present in the variable info file.
# The columns that remain are the new columns (or columns that didn't change), and the extra columns that PHESANT requires.

# Now, we're done...let's just double check that this matches the type of file expected by PHESANT.
df <- merge(x=PH_df, y=df, by="FieldID", all=TRUE)

checking_differences('TRAIT_OF_INTEREST')
checking_differences('CAT_MULT_INDICATOR_FIELDS')
checking_differences('CAT_SINGLE_TO_CAT_MULT')
checking_differences('DATA_CODING')
# Do we believe that our changes are the right ones - if so, we're all good.

# This portion is to check if we've made any mistakes - just sanity checking what we removed.
# There have been some changes - so change to the latest PHESANT version.
exclude <- checking_differences('EXCLUDED')$fields
# There are two differences that aren't NEALELAB exclusions, so change these back to being included?
exclude[which(exclude$EXCLUDED.y != "YES-NEALELAB"),]
# No, they're the ICD9 codes.

names(df)[which(names(df) == "EXCLUDED.y")] <- "EXCLUDED"

# Now I can do the same as before to remove unwanted columns.
df <- df[,-grep('\\.x', names(df))]
names(df) <- gsub('\\.y', '', names(df))

# Now read in and merge in the manually curated list of new variables to be excluded.
# Need to manually create this file using the output 'variable-info/new_phenotypes.tsv' file: add in a 
# new 'EXCLUDED' column and manually curate.
manual_df <- fread("variable-info/new_phenotypes_march_2019_excluded.tsv", sep='\t', header=TRUE, data.table=FALSE)
df <- merge(df, manual_df, by="FieldID", all=TRUE)

# Double check that there's no overlap in the EXCLUDED - they should be completely disjoint.
if (any(is.na(df$EXCLUDED.y) & is.na(df$EXCLUDED.x)))
	print('ERROR')

df$EXCLUDED.x[which(df$EXCLUDED.y != "")] <- exclusion_name
names(df)[which(names(df) == "EXCLUDED.x")] <- "EXCLUDED"
df <- df[,-grep('\\.y', names(df))]
names(df) <- gsub('\\.x', '', names(df))

# Finally, set 'CODING' equal to the newer 'coding'
df$DATA_CODING <- df$Coding

# Also, need to fill in information for the categorical multiple variables that have been added.

# Write out, and make sure it's tab separated.
fwrite(df, sep='\t', file = "variable-info/outcome_info_final_round3.tsv")
# NOTE: Have to read in with Excel and write to .tsv to get the correct behavious with PHESANT...

