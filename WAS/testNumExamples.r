# Remove variable values if less than 10 examples have this value
testNumExamples <- function(pheno, varlogfile)
{
    # Loop through values and remove if has < 10 examples
    uniqVar <- unique(na.omit(pheno))
    
    for (u in uniqVar) {
        withValIdx <- which(pheno==u)
        numWithVal <- length(withValIdx)
        if (numWithVal<10) {
            pheno[withValIdx] <- NA
            cat("Removed ",u ,": ", numWithVal, "<10 examples || ", sep="",
                file=varlogfile, append=TRUE)
        } else {
            cat("Inc(>=10): ", u, "(", numWithVal, ") || ", sep="",
                file=varlogfile, append=TRUE)
        }
    }
    return(pheno)
}
