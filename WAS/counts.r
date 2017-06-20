# Adds given value to counter, that are used to count how many variables go down each route in the data flow
addToCounts <- function(countName, num)
{
    idx <- which(counters$name==countName)

    if (length(idx)==0) {
        # Counter does not exist so add with countValue 1
        counters <<- rbind(counters, data.frame(name=countName, countValue=num))
    } else {
        # Add to counter that already exists
        counters$countValue[idx] <<- counters$countValue[idx]+num
    }
}

# Saves the counters stored in count variables, to a file in results directory
saveCounts <- function() 
{
    countFile <- paste(opt$resDir, "variable-flow-counts-",
                       opt$varTypeArg,".txt", sep="")

    # Sort on counter name
    sortIdx <- order(as.character(counters[,"name"]))
    counters <<- counters[sortIdx,]

    write.table(counters, file=countFile, sep=",", quote=FALSE, row.names=FALSE)
}

# Increments counters used to count how many variables go down each route in the data flow
incrementCounter <- function(countName)
{
    idx <- which(counters$name == countName)

    if (length(idx) == 0) {
        # Counter does not exist so add with countValue 1
        counters <<- rbind(counters, data.frame(name=countName, countValue=1))
    } else {
        # Increment counter that already exists
        counters$countValue[idx] <<- counters$countValue[idx]+1
    }
}

