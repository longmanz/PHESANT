# Main function called for continuous fields
testContinuous <- function(varName, varType, thisdata, varlogfile)
{
    cat("CONTINUOUS MAIN || ", file=varlogfile, append=TRUE)	

    pheno <- thisdata[,phenoStartIdx:ncol(thisdata)]
    # reassign values
    pheno <- reassignValue(pheno, varName, varlogfile)
    thisdata[,phenoStartIdx:ncol(thisdata)] <- pheno

    data_to_add <- testContinuous2(varName, varType, thisdata, varlogfile)
    return(data_to_add)
}

# Main code used to process continuous fields, or integer fields that have been reassigned as continuous because they have >20 distinct values.
# This is needed because we have already reassigned values for integer fields, so do this in the function above for continuous fields.
testContinuous2 <- function(varName, varType, thisdata, varlogfile)
{
    cat("CONTINUOUS || ", file=varlogfile, append=TRUE)

    pheno <- thisdata[,phenoStartIdx:ncol(thisdata)]

    if (!is.null(dim(pheno))) {
        phenoAvg <- rowMeans(pheno, na.rm = TRUE)
    } else {
        phenoAvg <- pheno
    }

    ## Recode NaN to NA, which is generated if all cols of pheno are NA for a given person
    idxNan <- which(is.nan(phenoAvg))
    phenoAvg[idxNan] <- NA
    numNotNA <- length(na.omit(phenoAvg))

    ## Check whether > 20% examples with same value
    uniqVar <- unique(na.omit(phenoAvg))
    valid <- TRUE
    for (uniq in uniqVar) {
        numWithValue <- length(which(phenoAvg == uniq))
        if (numWithValue / numNotNA >= 0.2) {
            valid <- FALSE
            break
        }
    }

    if (valid == FALSE) {
        # Treat as ordinal categorical
        cat(">20% IN ONE CATEGORY || ", file=varlogfile, append=TRUE)

        # If >2 unique values then treat as ordered categorical
        numUniqueValues <- length(uniqVar)

        # Straight forward case that there are two (or one) values		
        if (numUniqueValues <= 2) {
            # Treat as binary or skip (binary requires>=10 per category)

            # Remove categories if < 10 examples to see if this should be binary
            # or not, but if ordered categorical then we include all values when 
            # generating this.
            phenoAvgMoreThan10 <- testNumExamples(phenoAvg, varlogfile)

            # Binary if 2 distinct values, else ordered categorical
            phenoFactor <- factor(phenoAvgMoreThan10)
            numLevels <- length(unique(na.omit(phenoAvgMoreThan10)))

            if (numLevels<=1) {
                cat("SKIP (number of levels: ", numLevels, ")",
                    sep="", file=varlogfile, append=TRUE)
                incrementCounter("cont.onevalue")
                return(NULL)
            } else if (numLevels == 2) {
            	# Binary so logistic regression
                incrementCounter("cont.binary")
            	thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], phenoFactor)			
            	data_to_add <- binaryLogisticRegression(varName, varType,
                    thisdatanew, varlogfile)
                return(data_to_add)
            }

        } else {
            ## Try to treat as ordered categorical
            incrementCounter("cont.ordcattry")
            ## Equal sized bins
            phenoBinned <- equalSizedBins(phenoAvg, varlogfile)
            
            # Check number of people in each bin
            bin0Num <- length(which(phenoBinned == 0))
            bin1Num <- length(which(phenoBinned == 1))
            bin2Num <- length(which(phenoBinned == 2))

            if (bin0Num >= 10 & bin1Num >= 10 & bin2Num >= 10) {
                # Successful binning. >=10 examples in each of the 3 bins
                incrementCounter("cont.ordcattry.ordcat")
                thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], phenoBinned)
                data_to_add <- testCategoricalOrdered(varName, varType,
                    thisdatanew, varlogfile)
                return(data_to_add)
            } else {
                # Try to treat as binary because not enough examples in each bin
                if (bin0Num < 10 & bin2Num < 10) {
                    # Skip - not possible to create binary variable because first
                    # and third bins are too small ie. could merge bin1 with bin 2 
                    # but then bin3 still too small etc...
                    cat("SKIP 2 bins are too small || ", file=varlogfile, append=TRUE)
                    incrementCounter("cont.ordcattry.smallbins")
                    return(NULL)
                } else if ((bin0Num < 10 | bin1Num < 10) &
                          (bin0Num + bin1Num) >= 10) {
                    # Combine first and second bin to create binary variable
                    incrementCounter("cont.ordcattry.binsbinary")
                    cat("Combine first two bins and treat as binary || ",
                        file=varlogfile, append=TRUE)
                    phenoBinned[which(phenoBinned==0)] <- 1	
                    # Binary so logistic regression.
                    thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], phenoBinned)
                    data_to_add <- binaryLogisticRegression(varName, varType,
                        thisdatanew, varlogfile)
                    return(data_to_add)
            	} else if ((bin2Num<10 | bin1Num<10) & (bin2Num+bin1Num)>=10) {
                    # Combine second and last bin to create binary variable
                    incrementCounter("cont.ordcattry.binsbinary")
                    cat("Combine last two bins and treat as binary || ",
                        file=varlogfile, append=TRUE)
                    phenoBinned[which(phenoBinned==2)] <- 1
                    # Binary, so logistic regression.
                    thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], phenoBinned)
                    data_to_add <- binaryLogisticRegression(varName, varType,
                        thisdatanew, varlogfile)
                    return(data_to_add)
                } else {
                    # Skip - not possible to create binary variable because combining
                    # bins would still be too small.
                    cat("SKIP 2 bins are too small(2) || ", file=varlogfile, append=TRUE)
                    incrementCounter("cont.ordcattry.smallbins2")
                    return(NULL)
                }
            }
        }
    } else {
        cat("IRNT || ", file=varlogfile, append=TRUE)
        incrementCounter("cont.main")

        # Check there are at least 500 examples
        numNotNA <- length(which(!is.na(phenoAvg)))
        if (numNotNA < opt$contnacutoff) {
            cat("CONTINUOUS-SKIP-", opt$contnacutoff, " (", numNotNA, ") || ",sep="",
                file=varlogfile, append=TRUE)
            incrementCounter(paste("cont.main.", opt$contnacutoff, sep=""))
            return(NULL)
        } else {
            # Inverse rank normal transformation
            phenoIRNT <- irnt(phenoAvg)

            # Do regression (use standardised geno values)
            geno <- scale(thisdata[,"geno"])
            confounders <- thisdata[,2:numPreceedingCols]

            # BEGIN TRYCATCH
            tryCatch(
            {
                incrementCounter("success.continuous")
                return(list(phenoIRNT, varName))

                # END TRYCATCH
            }, error = function(e) {
                cat("ERROR:", varName, gsub("[\r\n]", "", e))
                incrementCounter("continuous.error")
            })
        }
    }
}

irnt <- function(pheno) {
    set.seed(1234)
    numPhenos <- length(which(!is.na(pheno)))
    quantilePheno <- (rank(pheno, na.last = "keep", ties.method = "random") - 0.5) / numPhenos
    phenoIRNT <- qnorm(quantilePheno)	
    return(phenoIRNT)
}
