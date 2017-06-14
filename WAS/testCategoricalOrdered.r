# Performs ordered logistic regression test and saves results in ordered logistic results file
testCategoricalOrdered <- function(varName, varType, thisdata,
    varlogfile, orderStr="")
{
    pheno <- thisdata[,phenoStartIdx:ncol(thisdata)]
    geno <- thisdata[,"geno"]

    cat("CAT-ORD || ", file=varlogfile, append=TRUE)
    incrementCounter("ordCat")

    doCatOrdAssertions(pheno)
    uniqVar <- unique(na.omit(pheno))

    # Log the ordering of categories used
    orderStr <- setOrderString(orderStr, uniqVar)
    cat("order: ", orderStr, " || ", sep="", file=varlogfile, append=TRUE)

    # Check sample size
    numNotNA <- length(which(!is.na(pheno)))
    if (numNotNA < 500) {
        cat("CATORD-SKIP-500 (", numNotNA, ") || ", sep="",
            file=varlogfile, append=TRUE)
        incrementCounter("ordCat.500")
        return(NULL)
    } else {
        # Test this cat ordered variable with ordered logistic regression	
        phenoFactor <- factor(pheno)
        cat("num categories:", length(unique(na.omit(phenoFactor))), "||",
            file=varlogfile, append=TRUE)

        # BEGIN TRYCATCH
        tryCatch(
        {    
            incrementCounter("success.ordCat")
            return(list(phenoFactor, varName))
        # END TRYCATCH
        }, error = function(e) {
            cat("ERROR:", varName, gsub("[\r\n]", "", e))
            incrementCounter("ordCat.error")
        })
    }
}

# Check that the phenotype is valid - that there are more than two categories
# and that these all have at least 10 cases.
# Something has gone wrong if this is the case
doCatOrdAssertions <- function(pheno)
{
    # Assert variable has only one column    
    if (!is.null(dim(pheno))) stop("More than one column for categorical ordered")

    uniqVar <- unique(na.omit(pheno))

    # Assert more than 2 categories
    if (length(uniqVar) <= 1) stop("1 or zero values")
    if (length(uniqVar) == 2) stop("this variable is binary")

    # Assert each value has >= 10 examples
    for (u in uniqVar) {
        withValIdx <- which(pheno == u)
        numWithVal <- length(withValIdx)
        if (numWithVal < 10) stop("value with <10 examples")
    }
}

# If data coding file does not specify an order then we use the default order as in coding defined by Biobank
# and this function just generates a string with this order for logging purposes
setOrderString <- function(orderStr, uniqVar) {

    if (is.na(orderStr) || nchar(orderStr) == 0) {
    	orderStr <- ""
    	# Create order str by appending each value
        uniqVarSorted <- sort(uniqVar)
        first <- 1

        for (i in uniqVarSorted) {
            if (first == 0) orderStr <- paste(orderStr, "|", sep="")
    		if (i >= 0) orderStr <- paste(orderStr, i, sep="")
            first <- 0
        }
    }
    return(orderStr)
}
