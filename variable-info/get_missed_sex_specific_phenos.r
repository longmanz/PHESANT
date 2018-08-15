get_phenos <- function(tsv_file, cloud=TRUE) {
	if(cloud) {
		header <- strsplit(system(paste("gsutil cat", tsv_file, "| head -n 1"), intern=TRUE),split='\t')[[1]]
		header <- gsub("\"", "", header)[-c(1,2,3)]
	} else {
		header <- strsplit(system(paste("cat", tsv_file, "| head -n 1"), intern=TRUE),split='\t')[[1]]
		header <- gsub("\"", "", header)[-c(1,2,3)]
	}
	return(header)
}

remove_excess_whitespace <- function(x) x <- gsub("\\s+", " ", x)
remove_trailing_whitespace <- function(x) x <- gsub(" $", "", x)

# Read in all the data table codings
get_codings <- function(folder)
{
	to_read <- paste0(folder, "/", dir(folder))
	codings_tables <- list()

	for(i in to_read) {
		name <- gsub(paste0(folder, "/coding(.*).tsv"), "\\1", i)
		codings_tables[[name]] <- read.table(i, header=TRUE, sep='\t', quote="", stringsAsFactors=FALSE)
	}
	return(codings_tables)
}

get_df <- function(PHESANT_IDs, outcome_info)
{
	var_subvar <- sapply(PHESANT_IDs, strsplit, split="_")
	var <- sapply(var_subvar, "[", 1)
	subvar <- sapply(var_subvar, "[", 2)

	var_meaning <- rep("", length(var))
	subvar_meaning <- rep("", length(var))

	for(i in 1:length(var)) {
		if(!is.na(subvar[i])) {
			# Consider the subvar
			coding <- as.character(outcome_info$Coding[which(outcome_info$FieldID == var[i])])
			where_coding <- which(codings_tables[coding][[1]]$coding == subvar[i])
			subvar_meaning[i] <- ifelse(length(where_coding) > 0, 
							codings_tables[coding][[1]]$meaning[where_coding],
							"PHESANT recoding")
		}
		var_meaning[i] <- as.character(outcome_info$Field[which(outcome_info$FieldID == var[i])])
	}

	PHESANT_IDs_df <- data.frame(FullFieldID=PHESANT_IDs, FieldID=var, Field=var_meaning, SubFieldID=subvar, SubField=subvar_meaning)
	return(PHESANT_IDs_df)
}

get_manual_df <- function(specific_fields, df)
{
	where <- rep(0, length(specific_fields))
	for(i in 1:length(specific_fields)) {
		where[i] <- which(remove_trailing_whitespace(specific_fields[i]) ==
			sapply(paste(df$Field, df$SubField), remove_trailing_whitespace))
	}
	return(df[where,])
}

get_incorrect_phenos_in_both_sex_analysis <- function(df_manual, pheno_summary_both_sex, pheno_summary_single_sex,
	label)
{
	where_single_sex <- c()
	where <- c()
	for(i in 1:length(df_manual$FullFieldID)) {
		where_single_sex <- c(where_single_sex, which(rownames(pheno_summary_single_sex) == df_manual$FullFieldID[i]))
		where <- c(where, which(rownames(pheno_summary_both_sex) == df_manual$FullFieldID[i]))
	}

	append(names(df_manual), c(paste0(c('N.non.missing.', 'N.missing.', 'N.controls.', 'N.cases.'), label), 'Check'))

	df_manual[paste0('N.non.missing.', label)] <- pheno_summary_single_sex$N.non.missing[where_single_sex]
	df_manual[paste0('N.missing.', label)] <- pheno_summary_single_sex$N.missing[where_single_sex]
	df_manual[paste0('N.controls.', label)] <- pheno_summary_single_sex$N.controls[where_single_sex]
	df_manual[paste0('N.cases.', label)] <-pheno_summary_single_sex$N.cases[where_single_sex]

	df_manual$N.non.missing <- pheno_summary_both_sex$N.non.missing[where]
	df_manual$N.missing <- pheno_summary_both_sex$N.missing[where]
	df_manual$N.controls <- pheno_summary_both_sex$N.controls[where]
	df_manual$N.cases <- pheno_summary_both_sex$N.cases[where]
	df_manual['Check'] <- pheno_summary_both_sex$Field[where]

	# Now, with these, should be able to look into the both sexes analyses and see if the numbers are the same, or different.
	# include another column to be the number of cases observed in the different log files.
	where_different_controls <- which(df_manual$N.controls!=df_manual[paste0('N.controls.', label)])
	where_different_missing <- which(df_manual[paste0('N.non.missing.', label)]!=df_manual$N.non.missing)
	where_different <- unique(where_different_controls, where_different_missing)
	return(df_manual[where_different,])
}

codings_tables <- get_codings("~/Repositories/PHESANT/WAS/codings")

# Get the variable info file and the data_codings files as in summarise_phenotypes.r.
in_men_not_women <- c()
in_women_not_men <- c()

for(i in 1:4) {
	# Want to compare the male and female phenotypes to the all sexes phenotypes
	male_phenos <- paste0("ukb11214_final_QC_males_more_phenos_and_corrected.", i, ".tsv")
	female_phenos <- paste0("ukb11214_final_QC_females_more_phenos_and_corrected.", i, ".tsv")	
	phenos <- paste0("ukb11214_final_QC_more_phenos_and_corrected.", i, ".tsv")
	male_phenos <- get_phenos(male_phenos, cloud=FALSE)
	female_phenos <- get_phenos(female_phenos, cloud=FALSE)
	phenos <- get_phenos(phenos, cloud=FALSE)
	# Now ask which phenotypes are in one of the sexes, but not the other?
	in_men_not_women <- c(in_men_not_women, setdiff(male_phenos, female_phenos))
	in_women_not_men <- c(in_women_not_men, setdiff(female_phenos, male_phenos))
}

outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final_round2.tsv",
						   sep='\t', quote="", comment.char="", header=TRUE)

df_men_not_women <- get_df(in_men_not_women, outcome_info)
df_women_not_men <- get_df(in_women_not_men, outcome_info)

# Manually curate a sex-specific list from these.
male_specific <- c(
	"Ever had prostate specific antigen (PSA) test",
	"Relative age of first facial hair ",
	"Relative age voice broke ",
	"Hair/balding pattern Pattern 2",
	"Hair/balding pattern Pattern 3",
	"Hair/balding pattern Pattern 4",
	"Number of children fathered",
	"\"Cancer code, self-reported\" prostate cancer",
	"\"Cancer code, self-reported\" testicular cancer",
	"\"Non-cancer illness code, self-reported\" enlarged prostate",
	"\"Non-cancer illness code, self-reported\" prostate problem (not cancer)",
	"\"Non-cancer illness code, self-reported\" testicular problems (not cancer)",
	"\"Non-cancer illness code, self-reported\" bph / benign prostatic hypertrophy",
	"\"Non-cancer illness code, self-reported\" prostatitis",
	"\"Non-cancer illness code, self-reported\" undescended testicle",
	"\"Non-cancer illness code, self-reported\" erectile dysfunction / impotence",
	"Treatment/medication code tamsulosin",
	"Treatment/medication code alfuzosin",
	"Treatment/medication code finasteride",
	"Treatment/medication code testosterone product",
	"Treatment/medication code viagra 100mg tablet",
	"Treatment/medication code dutasteride",
	"Treatment/medication code cardura 1mg tablet",
	"Treatment/medication code cialis 10mg tablet",
	"Treatment/medication code zoladex 3.6mg implant",
	"Treatment/medication code saw palmetto product",
	"Treatment/medication code flomax mr 400micrograms m/r capsule",
	"Treatment/medication code tadalafil",
	"Treatment/medication code cialis 20mg tablet",
	"Treatment/medication code viagra 50mg tablet",
	"Treatment/medication code viagra 25mg tablet",
	"Treatment/medication code testogel 50mg gel 5g sachet",
	"Treatment/medication code xatral 2.5mg tablet",
	"Treatment/medication code sildenafil",
	"Treatment/medication code vardenafil",
	"Underlying (primary) cause of death: ICD10 C61 Malignant neoplasm of prostate"
	)

female_specific <- c(
	"Ever had breast cancer screening / mammogram",
	"Years since last breast cancer screening / mammogram",
	"Ever had cervical smear test",
	"Years since last cervical smear test",
	"Age when periods started (menarche)",
	"Had menopause",
	"Number of live births",
	"Age at first live birth",
	"Age at last live birth",
	"\"Ever had stillbirth, spontaneous miscarriage or termination\"",
	"Ever taken oral contraceptive pill",
	"Age started oral contraceptive pill",
	"Bilateral oophorectomy (both ovaries removed)",
	"Age at menopause (last menstrual period)",
	"Ever had hysterectomy (womb removed)",
	"Length of menstrual cycle",
	"Number of stillbirths",
	"Number of spontaneous miscarriages",
	"Number of pregnancy terminations",
	"Age of primiparous women at birth of child",
	"Gestational diabetes only",
	"\"Cancer code, self-reported\" breast cancer",
	"\"Cancer code, self-reported\" cervical cancer",                                                                                                                      
	"\"Cancer code, self-reported\" cin/pre-cancer cells cervix",                                                                                                          
	"\"Cancer code, self-reported\" uterine/endometrial cancer",                                                                                                           
	"\"Cancer code, self-reported\" ovarian cancer",
	"\"Non-cancer illness code, self-reported\" benign breast lump",
	"\"Non-cancer illness code, self-reported\" cervical intra-epithelial neoplasia (cin) / precancerous cells cervix",
	"\"Non-cancer illness code, self-reported\" uterine fibroids",
	"\"Non-cancer illness code, self-reported\" ovarian cyst or cysts",
	"\"Non-cancer illness code, self-reported\" polycystic ovaries/polycystic ovarian syndrome",
	"\"Non-cancer illness code, self-reported\" endometriosis",
	"\"Non-cancer illness code, self-reported\" dysmenorrhoea / dysmenorrhea",
	"\"Non-cancer illness code, self-reported\" breast disease (not cancer)",                                                                                              
	"\"Non-cancer illness code, self-reported\" cervical polyps",
	"\"Non-cancer illness code, self-reported\" menopausal symptoms / menopause",
	"\"Non-cancer illness code, self-reported\" uterine polyps",                                                                                                           
	"\"Non-cancer illness code, self-reported\" vaginal prolapse/uterine prolapse",                                                                                        
	"\"Non-cancer illness code, self-reported\" gestational hypertension/pre-eclampsia",
	"\"Non-cancer illness code, self-reported\" breast cysts",                                                                                                          
	"\"Non-cancer illness code, self-reported\" gynaecological disorder (not cancer)",
	"\"Non-cancer illness code, self-reported\" abnormal smear (cervix)",
	"\"Non-cancer illness code, self-reported\" post-natal depression",
	"\"Non-cancer illness code, self-reported\" female infertility",
	"\"Non-cancer illness code, self-reported\" menorrhagia (unknown cause)",
	"\"Non-cancer illness code, self-reported\" miscarriage",
	"\"Non-cancer illness code, self-reported\" ectopic pregnancy",
	"\"Non-cancer illness code, self-reported\" breast fibroadenoma",
	"\"Non-cancer illness code, self-reported\" fibrocystic disease",
	"\"Non-cancer illness code, self-reported\" gestational diabetes",
	"Treatment/medication code tamoxifen",
	"Treatment/medication code cerazette 75micrograms tablet",
	"Treatment/medication code femulen tablet",
	"Treatment/medication code estraderm mx 25 patch",
	"Treatment/medication code mirena 52mg intrauterine system",
	"Treatment/medication code anastrozole",
	"Treatment/medication code implanon 68mg subdermal implant",
	"Treatment/medication code prempak 0.625 tablet",
	"Treatment/medication code climaval 1mg tablet",
	"Treatment/medication code depo-provera 50mg/1ml injection",
	"Treatment/medication code micronor tablet",
	"Treatment/medication code oestrogen product",
	"Treatment/medication code mirena 20mcg/24hrs intrauterine system",
	"Treatment/medication code arimidex 1mg tablet",
	"Treatment/medication code ovestin 0.1% vaginal cream",
	"Treatment/medication code elleste duet conti tablet",
	"Treatment/medication code microgynon 30 tablet",
	"Treatment/medication code estradiol product",
	"Treatment/medication code femoston 1/10 tablet",
	"Treatment/medication code norethisterone",
	"Treatment/medication code climagest 1mg tablet",
	"Treatment/medication code premarin 625micrograms tablet",
	"Treatment/medication code premique 0.625mg/5mg tablet",
	"Treatment/medication code evorel 25 patch",
	"Treatment/medication code starflower oil",
	"Treatment/medication code logynon tablet",
	"Treatment/medication code elleste-solo 1mg tablet",
	"Treatment/medication code femseven 50 patch",
	"Treatment/medication code vagifem 25mcg pessary",
	"Treatment/medication code evorel conti patch",
	"Treatment/medication code premique cycle 10mg tablet",
	"Treatment/medication code kliofem tablet",
	"Treatment/medication code mercilon tablet",
	"Treatment/medication code ortho-gynest 500micrograms pessary",
	"Treatment/medication code hormonin tablet",
	"Treatment/medication code noriday tablet",
	"Treatment/medication code kliovance 1mg/0.5mg tablet",
	"Treatment/medication code evista 60mg tablet",
	"Treatment/medication code livial 2.5mg tablet",
	"Treatment/medication code climesse tablet",
	"Treatment/medication code progynova 1mg tablet",
	"Treatment/medication code estriol product",
	"Treatment/medication code cilest tablet",
	"Treatment/medication code tibolone",
	"Treatment/medication code tranexamic acid",
	"Treatment/medication code loestrin 20 tablet",
	"Treatment/medication code raloxifene hydrochloride",
	"Treatment/medication code norgeston tablet",
	"Treatment/medication code indivina 1mg/2.5mg tablet",
	"Treatment/medication code zumenon 1mg tablet",
	"Treatment/medication code exemestane",
	"Treatment/medication code nuvelle tablet",
	"Treatment/medication code progesterone product",
	"Treatment/medication code letrozole",
	"Treatment/medication code menophase tablet",
	"Anaesthetics administered during delivery Spinal anaesthetic",
	"Anaesthetics administered during delivery Epidural or caudal anaesthetic",                                                                                           
	"Anaesthetics administered during delivery Other",             
	"Anaesthetics administered during delivery General anaesthetic",                                                                                                       	
	"Anaesthetics administered during delivery Epidural or caudal, and spinal anaesthetic",                                                                               	
	"Anaesthetics administered post delivery Other",               
	"Anaesthetics administered post delivery General anaesthetic", 
	"Anaesthetics administered post delivery Epidural or caudal anaesthetic",                                                                                              
	"Anaesthetics administered post delivery Spinal anaesthetic",  
	"Delivery methods Emergency caesarean section",                
	"Delivery methods Other forceps, not breech",                  
	"Delivery methods Spontaneous other cephalic",                 
	"Delivery methods Spontaneous vertex",                         
	"Delivery methods Elective caesarean section",                 
	"Delivery methods Low forceps, not breech",                    
	"Delivery methods Other than those specified above",           
	"Delivery methods Ventouse, vacuum extraction",                
	"Hospital episode type Delivery episode",                      
	"Destinations on discharge from hospital (recoded) Transfer to other NHS provider: Obstetrics"
	)

df_manual_men_not_women <- get_manual_df(male_specific, df_men_not_women)
df_manual_women_not_men <- get_manual_df(female_specific, df_women_not_men)

# Next, want to read through the phenotype summary files and add that information.
for(i in 1:4) {
	pheno_summary_file_males <- paste0("ukb11214_final_QC_males_more_phenos_and_corrected_phesant_recodings.", i, "_phenosummary.tsv")
	pheno_summary_file_females <- paste0("ukb11214_final_QC_females_more_phenos_and_corrected_phesant_recodings.", i, "_phenosummary.tsv")
	pheno_summary_file <- paste0("ukb11214_final_QC_more_phenos_and_corrected_phesant_recodings.", i, "_phenosummary.tsv")
	if(i==1) {
		pheno_summary_males <- read.table(pheno_summary_file_males, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE)
		pheno_summary_females <- read.table(pheno_summary_file_females, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE)
		pheno_summary <- read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE)
	} else {
		pheno_summary_males <- rbind(pheno_summary_males, read.table(pheno_summary_file_males, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE))
		pheno_summary_females <- rbind(pheno_summary_females, read.table(pheno_summary_file_females, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE))
		pheno_summary <- rbind(pheno_summary, read.table(pheno_summary_file, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE))
	}
}

df_remove_should_only_be_in_males <- get_incorrect_phenos_in_both_sex_analysis(df_manual_men_not_women, pheno_summary, pheno_summary_males, 'males')
df_remove_should_only_be_in_females <- get_incorrect_phenos_in_both_sex_analysis(df_manual_women_not_men, pheno_summary, pheno_summary_females, 'females')
