# Tests an unordered categorical phenotype with multinomial regression
# and saves this result in the multinomial logistic results file.
testCategoricalUnordered <- function(varName, varType, thisdata, varlogfile) {

    pheno <- thisdata[,phenoStartIdx:ncol(thisdata)]
    geno <- thisdata[,"geno"]

    idxNA <- which(is.na(pheno))
    numNotNA <- length(which(!is.na(pheno)))
    numRows <- nrow(thisdata)

    if (numNotNA < 500) {
        cat("CATUNORD-SKIP-500 (", numNotNA, ") || ",sep="",
            file=varlogfile, append=TRUE)
        incrementCounter("unordCat.500")
        return(NULL)
    } else {
        # Check there are not too many levels and skip if there are.
        numUnique <- length(unique(na.omit(pheno)))
        if (numUnique > 1000) {
            cat("Too many levels: ", numUnique, " > 1000 || SKIP ", sep="",
                file=varlogfile, append=TRUE)
            incrementCounter("unordCat.cats")
            return(NULL)
        }

        phenoFactor <- chooseReferenceCategory(pheno, varlogfile)
        loop <- levels(phenoFactor)[-1]
        varBinarymat <- as.data.frame(matrix(ncol=length(loop), nrow=numRows))

        j <- 1
        binaryMatColNames <- paste(varName, loop, sep='_')
        binaryMatColNames <- gsub("[[:space:]]", "", binaryMatColNames)

        for (i in loop) {
            idxForVar <- which(phenoFactor == i)
            cat(" CAT-SINGLE-BINARY-VAR: ", i, " || ", sep="",
                file=varlogfile, append=TRUE)
            incrementCounter("catSingle.binary")
            # Make zero vector and set 1s for those with this variable value
            varBinarymat[,j] <- rep.int(FALSE, numRows)
            varBinarymat[idxForVar,j] <- TRUE
            j <- j+1
        }

        varBinarymat[idxNA,] <- NA

        # BEGIN TRYCATCH
        tryCatch(
        {
            incrementCounter("success.unordCat")
            return(list(varBinarymat, binaryMatColNames))
            # END TRYCATCH
        }, error = function(e) {
            cat(paste("ERROR:", varName, gsub("[\r\n]", "", e), sep=" "))
            incrementCounter("unordCat.error")
    	})
    }
}

# Find reference category - category with most number of examples
chooseReferenceCategory <- function(pheno, varlogfile)
{
    uniqVar <- unique(na.omit(pheno))
    phenoFactor <- factor(pheno, ordered=FALSE)
    maxFreq <- 0
    maxFreqVar <- ""

    for (u in uniqVar) {
        withValIdx <- which(pheno==u)
        numWithVal <- length(withValIdx)

        if (numWithVal>maxFreq) {
            maxFreq <- numWithVal
            maxFreqVar <- u
        }
    }
    cat("reference: ", maxFreqVar,"=",maxFreq, " || ", sep="",
        file=varlogfile, append=TRUE)

    # Choose reference (category with largest frequency)
    phenoFactor <- relevel(phenoFactor, ref = paste("", maxFreqVar, sep=""))
    return(phenoFactor)
}
