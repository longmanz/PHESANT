# Reassigns values as specified in data coding info file.
reassignValue <- function(pheno, varName, varlogfile)
{
    # Get data code info - whether this data code is ordinal or not 
    # and any reordering and reassignments
    dataPheno <- vl$phenoInfo[which(vl$phenoInfo$FieldID==varName),]
    dataCode <- dataPheno$DATA_CODING

    # Not all variables will have a data code info row
    dataCodeRow <- which(vl$dataCodeInfo$dataCode==dataCode)

    if (length(dataCodeRow)==0) {
        return(pheno)
    } else if (length(dataCodeRow) > 1) {
        cat("WARNING: >1 ROWS IN DATA CODE INFO FILE || ",
            file=varlogfile, append=TRUE)
        return(pheno)
    }

    dataDataCode <- vl$dataCodeInfo[dataCodeRow,]
    reassignments <- as.character(dataDataCode$reassignments)
    return(reassignValue2(pheno, reassignments, varlogfile))
}

# Reassigns values in pheno, as specified in reassignments argument.
reassignValue2 <- function(pheno, reassignments, varlogfile)
{
    # Can be NA if row not included in data coding info file
    if (!is.na(reassignments) && nchar(reassignments)>0)
    {
        reassignParts <- unlist(strsplit(reassignments, "\\|"))
        cat("reassignments:", reassignments, "|| ", file=varlogfile, append=TRUE)

        # Do each reassignment
        reassignment_list <- vector("list", length(reassignParts))
        j <- 1
        for (i in reassignParts) {
            reassignParts <- unlist(strsplit(i,"="))
            # Matrix version
            reassignment_list[[j]]$idx <- which(pheno == reassignParts[1], arr.ind = TRUE)
            reassignment_list[[j]]$reassignment <- strtoi(reassignParts[2])
            j <- j+1
        }

        for (j in 1:length(reassignment_list)) {
            pheno[reassignment_list[[j]]$idx] <- reassignment_list[[j]]$reassignment
        }

        # Do each reassignment
        # for (i in reassignParts) {
        #     reassignParts <- unlist(strsplit(i,"="))
        #     # Matrix version
        #     idx <- which(pheno == reassignParts[1], arr.ind = TRUE)
        #     pheno[idx] <- strtoi(reassignParts[2])
        # }

        if(!is.null(dim(pheno))) {
            pNum <- as.data.frame(matrix(as.numeric(unlist(pheno)), ncol=ncol(pheno)))
            colnames(pNum) <- colnames(pheno)
        } else {
            pNum <- as.data.frame(as.numeric(unlist(pheno)))
        }

        # See if type has changed (this happens for field 216 (X changed to -1))
        # as.numeric will set non numeric to NA so we know if it's ok to do this by seeing if there are extra NA's after the conversion
        # pNum <- as.numeric(unlist(pheno))
        isNum <- length(which(is.na(pheno), arr.ind = TRUE)) == 
                 length(which(is.na(pNum), arr.ind = TRUE))
        print(isNum)
        if (isNum) pheno <- pNum
    }
    # print("hello")
	return(pheno)
}
