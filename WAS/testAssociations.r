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
                	cat("Excluded integer:", excluded, "||\n", file=varlogfile,
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
                    cat("Excluded continuous:", excluded, "||\n", file=varlogfile, append=TRUE)
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
                    cat("Excluded cat-single:", excluded, "||\n", file=varlogfile,
                        append=TRUE)
                    incrementCounter("excluded.catSin")
                    return(NULL)
                } else {
                    incrementCounter("start.catSin")
                    data_to_add <- testCategoricalSingle(currentVarShort, "CAT-SIN",
                                                         thisdata, varlogfile)
                    cat('\n', file=varlogfile, append=TRUE)
                    return(data_to_add)
                }

        	} else if (fieldType == "Categorical multiple" || catSinToMult != "") {
                # CAT MULTIPLE
                cat(currentVar, "|| ", sep="", file=varlogfile, append=TRUE)
                if (excluded != "") {
                    cat("Excluded cat-multiple:", excluded, "||\n",
                        file=varlogfile, append=TRUE)
                    incrementCounter("excluded.catMul")
                    return(NULL)
                } else {

                    if (catSinToMult != "") {
                        cat("cat-single to cat-multiple || ", sep="",
                            file=varlogfile, append=TRUE)
                        incrementCounter("catSinToCatMul")
                    }
                    
                    data_to_add <- testCategoricalMultiple(currentVarShort, "CAT-MUL",
                                                           thisdata, varlogfile)
                }
                cat("\n", file=varlogfile, append=TRUE)
                return(data_to_add)
            } else {
                cat("VAR MISSING", currentVarShort, "\n",
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
