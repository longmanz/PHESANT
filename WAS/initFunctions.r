# Load the required R files
loadSource <- function() {
    source("loadData.r")
    source("reassignValue.r")
    source("validateInput.r")
    source("testNumExamples.r")
    source("binaryLogisticRegression.r")
    source("equalSizedBins.r")
    source("fixOddFieldsToCatMul.r")
    source("replaceMissingandNaN.r")
    source("testAssociations.r")
    source("testCatMultiple.r")
    source("testCatSingle.r")
    source("testContinuous.r")
    source("testInteger.r")
    source("testCategoricalOrdered.r")
    source("testCategoricalUnordered.r")
    source("counts.r")
}

# Init the counters used to determine how many variables took each path 
# in the variable processing flow.
initCounters <- function() {
	counters <- data.frame(name=character(),
						   countValue=integer(),
                           stringsAsFactors=FALSE)
	return(counters)
}

# Load the variable information and data code information files
initVariableLists <- function()
{
	phenoInfo <- read.table(opt$variablelistfile, sep="\t", header=1, comment.char="", quote="")
	dataCodeInfo <- read.table(opt$datacodingfile,sep=",", header=1)
	vars <- list(phenoInfo=phenoInfo, dataCodeInfo=dataCodeInfo)
	return(vars)
}
