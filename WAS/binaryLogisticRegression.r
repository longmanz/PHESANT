# Performs binary logistic regression on the phenotype stored in this data 
# and stores result in 'results-logistic-binary' results file.

binaryLogisticRegression <- function(varName, varType, thisdata, varlogfile)
{   
    to_save <- thisdata[, phenoStartIdx]
    phenoFactor <- factor(thisdata[, phenoStartIdx])
    facLevels <- levels(phenoFactor)

    # Assert variable has exactly two distinct values
    if (length(facLevels) != 2) {
        cat("BINARY-NOT2LEVELS- (", length(facLevels), ") || ",
            sep="", file=varlogfile, append=TRUE)
        incrementCounter("binary.nottwolevels")
    }

    idxTrue <- length(which(phenoFactor == facLevels[1]))
    idxFalse <- length(which(phenoFactor == facLevels[2]))
    numNotNA <- length(which(!is.na(phenoFactor)))

    if (idxTrue < 10 || idxFalse < 10) {
        cat("BINARY-LOGISTIC-SKIP-10 (", idxTrue, "/", idxFalse, ") || ",
            sep="", file=varlogfile, append=TRUE)
        incrementCounter("binary.10")
        return(NULL)
    } else if (numNotNA < 500) {	
        cat("BINARY-LOGISTIC-SKIP-500 (", numNotNA, ") || ",
            sep="", file=varlogfile, append=TRUE)
        incrementCounter("binary.500")
        return(NULL)
    } else {
        cat("sample ", idxTrue, "/", idxFalse, "(", numNotNA, ") || ",
            sep="", file=varlogfile, append=TRUE)

        tryCatch(
        { 
            incrementCounter("success.binary")
            return(list(phenoFactor, varName))
            # END TRYCATCH
        }, error = function(e) {
            print(paste("ERROR:", varName, gsub("[\r\n]", "", e)) )
            incrementCounter("binary.error")
            return(NULL)
        })
    }
}

