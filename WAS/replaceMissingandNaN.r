# Replace negative values with NA as these are assumed to be missing

replaceMissingCodes <- function(pheno)
{
    phenoReplaced <- pheno
    uniqVar <- unique(na.omit(phenoReplaced))

    # Variable values <0 are `missing' codes
    for (u in uniqVar) {
        if (u < 0) {
            idxU <- which(phenoReplaced == u)
            phenoReplaced[idxU] <- NA
        }
    }

    return(phenoReplaced)
}

# Replace NaN and empty values with NA in pheno
replaceNaN <- function(pheno)
{
    if (is.factor(pheno)) {
        phenoReplaced <- pheno
        nanStr <- which(phenoReplaced == "NaN")
        phenoReplaced[nanStr] <- NA 
        emptyx <- which(phenoReplaced == "")
        phenoReplaced[emptyx] <- NA
    } else {
        phenoReplaced <- pheno
        nanx <- which(is.nan(phenoReplaced))
        phenoReplaced[nanx] <- NA
        emptyStr<- which(phenoReplaced == "")   
        phenoReplaced[emptyStr]<- NA  
    }

    return(phenoReplaced)
}

allNAs <- function(pheno_row) {
    all(is.na(pheno_row))
}