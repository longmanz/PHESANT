# Processing integer fields, namely:
# 1) Re-assigning values as specified in the data code information file
# 2) Generate a single value if there are several values (arrays) by taking the mean
# 3) Treating this field as continuous if at least 20 distinct values.
# Otherwise treat as binary or ordered categorical if 2 or more than two values. 

testInteger <- function(varName, varType, thisdata, varlogfile) {
    cat("INTEGER || ", file=varlogfile, append=TRUE)

    pheno <- thisdata[,phenoStartIdx:ncol(thisdata)]

    if (!is.numeric(as.matrix(pheno))) {
        cat("SKIP Integer type but not numeric", sep="",
            file=varlogfile, append=TRUE)
        return(NULL)
    }

    pheno <- reassignValue(pheno, varName, varlogfile)

    # Average if multiple columns
    if (!is.null(dim(pheno))) {
        phenoAvg <- rowMeans(pheno, na.rm=TRUE)
        # If participant only has NA values then NaN is generated so we 
        # convert back to NA
        phenoAvg <- replaceNaN(phenoAvg)
    } else {
        phenoAvg <- pheno
    }

    uniqVar <- unique(na.omit(phenoAvg))

    # If >=20 separate values then treat as continuous
    if (length(uniqVar) >= 20) {
        thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], phenoAvg)
        data_to_add <- testContinuous2(varName, varType, thisdatanew, varlogfile)
        incrementCounter("int.continuous")
        return(data_to_add)
    } else {
        # Remove categories if < 10 examples
        phenoAvg <- testNumExamples(phenoAvg, varlogfile)

        # Binary if 2 distinct values, else ordered categorical
        phenoFactor <- factor(phenoAvg)
        numLevels <- length(levels(phenoFactor))
        if (numLevels <= 1) {
            cat("SKIP (number of levels: ", numLevels, ")", sep="",
                file=varlogfile, append=TRUE)
            incrementCounter("int.onevalue")
            return(NULL)
        } else if (numLevels==2) {
            incrementCounter("int.binary")
            # Binary so logistic regression
            thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], phenoFactor)
            data_to_add <- binaryLogisticRegression(varName, varType, thisdatanew, varlogfile)
            return(data_to_add)
        } else {
            incrementCounter("int.catord")
            cat("3-20 values || ", file=varlogfile, append=TRUE)
            # We don't use equal sized bins just the original integers 
            # (that have >=10 examples) as categories
            thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], phenoFactor)
            # Treat as ordinal categorical
            data_to_add <- testCategoricalOrdered(varName, varType, thisdatanew, varlogfile)
            return(data_to_add)
        }
    }
}
