# Remove variable values if less than 10 examples have this value
testNumExamples <- function(pheno, varlogfile)
{
    # Loop through values and remove if has < opt$mincategorysize examples
    uniqVar <- unique(na.omit(pheno))
    
    for (u in uniqVar) {
        withValIdx <- which(pheno==u)
        numWithVal <- length(withValIdx)
        if (numWithVal < opt$mincategorysize) {
            pheno[withValIdx] <- NA
            cat("Removed ",u ,": ", numWithVal, "<", opt$mincategorysize, " examples || ", sep="",
                file=varlogfile, append=TRUE)
        } else {
            cat("Inc(>=", opt$mincategorysize"): ", u, "(", numWithVal, ") || ", sep="",
                file=varlogfile, append=TRUE)
        }
    }
    return(pheno)
}
