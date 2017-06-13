# Loads phenotype and trait of interest data files
# Creates phenotype / trait of interest data frame
# Creates confounder data frame
# Returns an object holding these two data frames
loadData <- function()
{
    # Load phenotype
    cat("Loading phenotypes...\n")
    
    if (length(grep('.tsv', opt$phenofile)) == 1 | length(grep('.tab', opt$phenofile)) == 1) {
        phenotype <- read.table(opt$phenofile, header=1,sep='\t')
    } else if (length(grep('.csv', opt$phenofile)) == 1) {
            phenotype <- read.table(opt$phenofile, header=1,sep=',')
    } else if (length(grep('.Rdata', opt$phenofile)) == 1) {
        load(opt$phenofile)
        if (!exists("phenotype")) {
            stop("Error: phenotype not found in .Rdata file")
        }

    } else {
        stop("Cannot detect the format of the phenotype file")
    }
    if (ncol(phenotype) == 1)
        stop("Number of columns of the phenotype file is 1: is it tab separated?")

    validatePhenotypeInput(phenotype)

    # Load SNPs.
    if (is.null(opt$traitofinterestfile)) {
        cat("Extracting trait of interest from pheno file...\n")
        
        if (opt$traitofinterest %in% colnames(phenotype)) {

            cat("Trait of interest found in pheno file.\n")
            # Set name of trait of interest to geno in phenotype file
            idxTOI <- which(colnames(phenotype) == opt$traitofinterest)
            colnames(phenotype)[idxTOI] <- "geno"		
            datax <- phenotype

            # Remove all rows with no trait of interest
            idxNotEmpty <- which(!is.na(phenotype[,idxTOI]))
            cat(cat("Phenotype file has", nrow(phenotype), "rows with",
                    length(idxNotEmpty), "not NA for trait of interest ("),
                opt$traitofinterest,").", sep="")
            phenotype <- phenotype[idxNotEmpty,]
        } else {
            stop(cat("Trait of interest (", opt$traitofinterest,")",
                     " not found in phenotype file ",opt$phenofile, ". ",
                     "Trait of interest should either be in phenotype file ",
                     "or seperate trait of interest file specified in ",
                     "traitofinterestfile arg.", sep=""), call.=FALSE)
        }

    } else {
        cat("Loading trait of interest file...\n")
        snpScores <- read.table(opt$traitofinterestfile, sep=",", header=1)
        validateTraitInput(snpScores)

        # Keep only the userID and exposure variable
        idx1 <- which(names(snpScores) == opt$userId)
        idx2 <- which(names(snpScores) == opt$traitofinterest)

        # Remove all rows with no trait of interest
        idxNotEmpty <- which(!is.na(snpScores[,idx2]))
        cat("Trait of interest has", nrow(snpScores), "rows with",
            length(idxNotEmpty), "not NA.")
        snpScores <- snpScores[idxNotEmpty,]

        snpScores <- cbind.data.frame(snpScores[,idx1], snpScores[,idx2])
        colnames(snpScores)[1] <- opt$userId
        colnames(snpScores)[2] <- "geno"

        print("Merging trait of interest and phenotype data")
        # Merge to one matrix
        datax <- merge(snpScores, phenotype, by=opt$userId, all=FALSE)
    }

    if (nrow(datax)==0) {
    	stop("No examples with row in both trait of interest and phenotype files", call.=FALSE)
    }

    datax <- fixOddFieldsToCatMul(datax)

    # Variables we adjust for.
    age <- datax[,"x21022_0_0"]
    sex <- datax[,"x31_0_0"]
    confounders <- cbind.data.frame(age,sex)

    # If genetic trait of interest then adjust for genotype chip
    # and also let user choose sensitivity analysis that also adjusts 
    # for top 10 genetic principal components and assessment centre.
    if (opt$genetic == TRUE)
    {
    	genoBatch <- datax[,"x22000_0_0"]

    	# Chip comes from batch field 22000
    	genoChip <- rep.int(NA,nrow(datax))
    	idxForVar <- which(genoBatch<0)
    	genoChip[idxForVar] <- 0
    	idxForVar <- which(genoBatch>=0 & genoBatch<2000)
    	genoChip[idxForVar] <- 1
    	confounders <- cbind.data.frame(confounders, genoChip)

        if (opt$sensitivity==TRUE) {
            genoPCs <- cbind(datax[,"x22009_0_1"], datax[,"x22009_0_2"],
                             datax[,"x22009_0_3"], datax[,"x22009_0_4"],
                             datax[,"x22009_0_5"], datax[,"x22009_0_6"],
                             datax[,"x22009_0_7"], datax[,"x22009_0_8"],
                             datax[,"x22009_0_9"], datax[,"x22009_0_10"])
            assessCenter <- datax[,"x54_0_0"]
            confounders <- cbind.data.frame(confounders, genoPCs, assessCenter)
            cat("Adjusting for age, sex, genotype chip, top 10 genetic principal",
                  "components and assessment centre.\n")
        } else {
            print("Adjusting for age, sex and genotype chip")
        }
    } else {
        # Non genetic trait of interest, then sensitivity adjusts for assessment center
        if (opt$sensitivity==TRUE) {
            assessCenter <- datax[,"x54_0_0"]
            confounders <- cbind.data.frame(confounders, assessCenter)
            cat("Adjusting for age, sex and assessment centre.\n")
        } else {
            cat("Adjusting for age and sex.\n")
        }
    }

    d <- list(datax=datax, confounders=confounders)

    return(d)

}
