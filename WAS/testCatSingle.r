# Performs variable processing for categorical (single) fields, namely:
# 1) Reassigning values as specified in data coding information file
# 2) Reordering categories for ordered fields
# 3) Replacing missing codes - we assume values < 0 are missing for categorical (single) variables
# 4) Remove values with <10 cases
# 5) Determine correct test to perform, either binary, ordered or unordered.

testCategoricalSingle <- function(varName, varType, thisdata, varlogfile)
{
    cat("CAT-SINGLE || ", file=varlogfile, append=TRUE)
    pheno <- thisdata[,phenoStartIdx:ncol(thisdata)]

    # Assert variable has only one column
    if (!is.null(dim(pheno))) stop("More than one column for categorical single")

    pheno <- reassignValue(pheno, varName, varlogfile)

    # Get data code info - whether this data code is ordinal or not 
    # and any reordering
    dataPheno <- vl$phenoInfo[which(vl$phenoInfo$FieldID==varName),]
    dataCode <- dataPheno$DATA_CODING

    # Get data coding information.
   	dataCodeRow <- which(vl$dataCodeInfo$dataCode==dataCode)
    
    if (length(dataCodeRow)==0) {
        cat("ERROR: No row in data coding info file || ",
            file=varlogfile, append=TRUE)
        return(NULL)
    }

    dataDataCode <- vl$dataCodeInfo[dataCodeRow,]
    ordered <- dataDataCode$ordinal
    order <- as.character(dataDataCode$ordering)

    # Reorder variable values into increasing order
    # (we do this now as this may convert variable to binary rather than ordered)
    pheno <- reorderOrderedCategory(pheno, order, varlogfile)

    # If data code has a default_value then recode NA's to this value 
    # for participants with value in default_related_field
    # this is used where there is no zero option e.g. field 100200
    defaultValue <- dataDataCode$default_value
    defaultRelatedID <- dataDataCode$default_related_field
    pheno <- setDefaultValue(pheno, defaultValue, defaultRelatedID, varlogfile)

    # All categories coded as < 0 we assume are `missing' values
    pheno <- replaceMissingCodes(pheno)

    # Remove categories if < opt$mincategorysize examples
    pheno <- testNumExamples(pheno, varlogfile)

    uniqVar <- unique(na.omit(pheno))
    uniqVar <- sort(uniqVar)

    if (length(uniqVar) <= 1) {
    	cat("SKIP (only one value) || ", file=varlogfile, append=TRUE)
    	incrementCounter("catSin.onevalue")
        return(NULL)
    } else if (length(uniqVar) == 2) {		
        cat("CAT-SINGLE-BINARY || ", file=varlogfile, append=TRUE)
        incrementCounter("catSin.case3")

        # Binary so logistic regression
        phenoFactor <- factor(pheno)
        thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], phenoFactor)
        data_to_add <- binaryLogisticRegression(varName, varType, thisdatanew, varlogfile)
        
        return(data_to_add)
    } else {
        # > 2 categories
        if (is.na(ordered)) {
            cat(" ERROR: 'ordered' not found in data code info file")
            return(NULL)	
        } else {
            # Unordered
            if (ordered == 0) {
                cat("CAT-SINGLE-UNORDERED || ", file=varlogfile, append=TRUE)
                incrementCounter("catSin.case2")

                thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], pheno)
                data_to_add <- testCategoricalUnordered(varName, varType, thisdatanew,
                    varlogfile)
                return(data_to_add)
            } else if (ordered == 1) {
                # Ordered
                cat("ordered || ", file=varlogfile, append=TRUE)
                incrementCounter("catSin.case1")

                # Reorder variable values into increasing order
                thisdatanew <- cbind.data.frame(thisdata[,1:numPreceedingCols], pheno)
                data_to_add <- testCategoricalOrdered(varName, varType, thisdatanew,
                    varlogfile, order)
                return(data_to_add)
            } else if (ordered == -2) {
                cat(" EXCLUDED or BINARY variable: Should not get here in code. ")
                incrementCounter("catSin.binaryorexcluded")
                return(NULL)
            } else {
                print(paste("ERROR", varName, varType, dataCode))
                return(NULL)
            }
        }
    }
}

# Values are reordered and assigned values 1:N for N categories
reorderOrderedCategory <- function(pheno, order, varlogfile)
{
    # New pheno of NAs (all values not in order are assumed to be NA)
    if (!is.na(order) && nchar(order)>0) {
        # Make empty pheno		
        pheno2 <- rep(NA,length(pheno))
        # Get ordering
        orderParts <- unlist(strsplit(order,"\\|"))
        
        # Go through values in correct order and set value
        # from 1 to the number of values
        count <- 1
        for (i in orderParts) {
            idx <- which(pheno==i)
            pheno2[idx] <- count
            count <- count+1
        }

        cat("reorder ",order," || ", sep="", file=varlogfile, append=TRUE)
        return(pheno2)

    } else {
        return(pheno)
    }
}

# Sets default value for people with no value in pheno, but with a value in the
# field specified in the default_value_related_field column in the data coding info file.
# the default value is specified in the default_value column in the data coding info file.
setDefaultValue <- function(pheno, defaultValue, defaultRelatedID, varlogfile)
{
    if (!is.na(defaultValue) && nchar(defaultValue) > 0) {
        # Remove people who have no value for indicator variable
        indName <- paste("x", defaultRelatedID, "_0_0", sep="")
     	cat("Default related field:", indName, "|| ",
            file=varlogfile, append=TRUE)
            
        if(length(indName) == 0){
            cat("Related field ", indName, " not exist in the data frame. SKIPPING...|| ",
                file=varlogfile, append=TRUE)
        } else {
            indicatorVar <- data[,indName]
            # Remove participants with NA value in this related field
    	    indicatorVar <- replaceNaN(indicatorVar)
            # Check if there are already examples with default value and if so display warning
            numWithDefault <- length(which(pheno==defaultValue))

            if (numWithDefault > 0) 
                cat("(WARNING: already", numWithDefault, "values with default value) ")

            # Set default value in people who have no value in the pheno but do 
            # have a value in the default_value_related_field
    	    defaultIdxs <- which(!is.na(indicatorVar) & is.na(pheno))
            pheno[defaultIdxs] <- defaultValue
       	    cat("default value", defaultValue, "set, N=", length(defaultIdxs), "|| ",
                file=varlogfile, append=TRUE)

        }
    }
    return(pheno)
}
