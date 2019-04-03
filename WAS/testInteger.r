# Processing integer fields, namely:
# 1) Re-assigning values as specified in the data code information file
# 2) Generate a single value if there are several values (arrays) by taking the mean
# 3) Treating this field as continuous if at least 20 distinct values.
# Otherwise treat as binary or ordered categorical if 2 or more than two values. 

testInteger <- function(varName, varType, thisdata, varlogfile) {
    cat("INTEGER || ", file=varlogfile, append=TRUE)
    print("this is actually the current variable")
    print(varName)
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

    # If opt$inttocontcutoff separate values then treat as continuous
    if (length(uniqVar) >= opt$inttocontcutoff) {
        print("convert to cts")
        thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], phenoAvg)
        data_to_add <- testContinuous2(varName, varType, thisdatanew, varlogfile)
        incrementCounter("int.continuous")
        return(data_to_add)
    } else {
        print("don't convert to cts")
        # Remove categories if < opt$mincategorysize examples
        phenoAvg <- testNumExamples(phenoAvg, varlogfile)

        # Binary if 2 distinct values, else ordered categorical
        phenoFactor <- factor(phenoAvg)
        numLevels <- length(levels(phenoFactor))
        print("number of levels")
        print(numLevels)
        if (numLevels <= 1) {
            cat("SKIP (number of levels: ", numLevels, ")", sep="",
                file=varlogfile, append=TRUE)
            incrementCounter("int.onevalue")
            return(NULL)
        } else if (numLevels==2) {
            print("num levels: 2")
            incrementCounter("int.binary")
            # Binary so logistic regression
            cat("INT-BINARY-VAR || ", sep="", file=varlogfile, append=TRUE)
            thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], phenoFactor)
            data_to_add <- binaryLogisticRegression(varName, varType, thisdatanew, varlogfile)
            return(data_to_add)
        } else {
            print("int cat ord")
            incrementCounter("int.catord")
            cat("3-", opt$inttocontcutoff, " values || ", file=varlogfile, append=TRUE)
            # We don't use equal sized bins just the original integers 
            # (that have >=opt$inttocontcutoff examples) as categories
            thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], phenoFactor)
            # Treat as ordinal categorical
            data_to_add <- testCategoricalOrdered(varName, varType, thisdatanew, varlogfile)
            return(data_to_add)
        }
    }
}
