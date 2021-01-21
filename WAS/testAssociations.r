# Tests the association of a field, determined by its field type
testAssociations <- function(currentVar, currentVarShort, thisdata, varlogfile, catmultcutoff)
{
    # Call file for variable type.
    tryCatch(
    {
        # Retrieve whether phenotype is excluded etc
        idx <- which(vl$phenoInfo$FieldID == currentVarShort)
        # Check if variable info is found for this field
        if (length(idx) == 0) {
            cat(currentVar, "|| Variable could not be found in pheno info file.\n",
                file=varlogfile, append=TRUE)			
            incrementCounter("notinphenofile")
            return(NULL)
        } else {
            # Get info from variable info file
            excluded <- vl$phenoInfo$EXCLUDED[idx]
            catSinToMult <- vl$phenoInfo$CAT_SINGLE_TO_CAT_MULT[idx]
            fieldType <- vl$phenoInfo$ValueType[idx]

            if (fieldType == "Integer") {
                # INTEGER
                cat(currentVar, "|| ", sep="", file=varlogfile, append=TRUE)
                if (excluded != "") {
                	cat("Excluded integer:", as.character(excluded), "||\n", file=varlogfile,
                        append=TRUE)
                	incrementCounter("excluded.int")
                    return(NULL)
                } else {
                    incrementCounter("start.int")
                    data_to_add <- testInteger(currentVarShort, "INTEGER",
                                               thisdata, varlogfile)
                    cat("\n", file=varlogfile, append=TRUE)
                    return(data_to_add)
                }
            } else if (fieldType == "Continuous") {
                # CONTINUOUS
                cat(currentVar, "|| ", sep="", file=varlogfile, append=TRUE)
                if (excluded != "") {
                    cat("Excluded continuous:", as.character(excluded), "||\n", file=varlogfile, append=TRUE)
                    incrementCounter("excluded.cont")
                    return(NULL)
                } else {
                    incrementCounter("start.cont")
                    data_to_add <- testContinuous(currentVarShort, "CONTINUOUS",
                                                  thisdata, varlogfile)
                    cat("\n", file=varlogfile, append=TRUE)
                    return(data_to_add)
                }   

            } else if (fieldType == "Categorical single" && catSinToMult == "") {
        	    # CAT SINGLE
                cat(currentVar, "|| ", sep="", file=varlogfile, append=TRUE)
                if (excluded != "") {
                    cat("Excluded cat-single:", as.character(excluded), "||\n", file=varlogfile,
                        append=TRUE)
                    incrementCounter("excluded.catSin")
                    return(NULL)
                } else {
                    incrementCounter("start.catSin")
                    # Changed input to below from currentVarShort to currentVar
                    # and evaluate visit number due to small set of edge cases that 
                    # error out due to starting at visit 2, not 0.
                    data_to_add <- testCategoricalSingle(currentVar, "CAT-SIN",
                                                         thisdata, varlogfile)
                    cat('\n', file=varlogfile, append=TRUE)
                    return(data_to_add)
                }

        	} else if (fieldType == "Categorical multiple" || catSinToMult != "") {
                # CAT MULTIPLE
                cat(currentVar, "|| ", sep="", file=varlogfile, append=TRUE)
                if (excluded != "") {
                    cat("Excluded cat-multiple:", as.character(excluded), "||\n",
                        file=varlogfile, append=TRUE)
                    incrementCounter("excluded.catMul")
                    return(NULL)
                } else {

                    if (catSinToMult != "") {
                        cat("cat-single to cat-multiple || ", sep="",
                            file=varlogfile, append=TRUE)
                        incrementCounter("catSinToCatMul")
                    }
                    
                    ### Start of Modification by Longmanz ###
                    # need to catch the error of sapply()! Modified on 2021 Jan 21
                    error_status <- try(testCategoricalMultiple(currentVarShort, "CAT-MUL", thisdata, varlogfile))
                    
                    if(class(error_status) == "try-error"){
                       cat("testAssociations.r: An \'sapply\' error occurs, please check!.")
                       return(NULL)
                     } else {
                       data_to_add <- testCategoricalMultiple(currentVarShort, "CAT-MUL", thisdata, varlogfile)
                     }
                    ### End of Modification ###
                    
                }
                cat("\n", file=varlogfile, append=TRUE)
                return(data_to_add)
            } else {
                cat("VAR MISSING - likely a 'text' field", currentVarShort, "\n",
                    file=varlogfile, append=TRUE)
                return(NULL)
            }
        }
    }, error = function(e) {
    	print(paste("ERROR:", currentVar, e, '\n'))
    })
    # Return the data that was used in the test.
    return(data_to_add)
}
