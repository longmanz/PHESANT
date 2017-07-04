# Splits the pheno into 3 bins with the cut points between values rather at the exact value for the quantile
equalSizedBins <- function(phenoAvg, varlogfile)
{
    # Equal sized bins 
    q <- quantile(phenoAvg, probs=c(1/3, 2/3), na.rm=TRUE)
    minX <- min(phenoAvg, na.rm=TRUE)
    maxX <- max(phenoAvg, na.rm=TRUE)
    phenoBinned <- phenoAvg

    cat("Split into three bins: ", file=varlogfile, append=TRUE)

    if (q[1] == minX) {
        
        # Edge case - quantile value is lowest value
        # Assign min value as cat1
        idx1 <- which(phenoAvg == q[1])
        phenoBinned[idx1] <- 0
        cat("0: [", q[1], ",", q[1], "],", file=varlogfile, append=TRUE, sep="")

        # Divide remaining values into cat2 and cat3
        phenoAvgRemaining <- phenoAvg[which(phenoAvg!=q[1])]
        qx <- quantile(phenoAvgRemaining, probs=c(0.5), na.rm=TRUE)
        minXX <- min(phenoAvgRemaining, na.rm=TRUE)
        maxXX <- max(phenoAvgRemaining, na.rm=TRUE)

        if (qx[1] == minXX) {
            # Edge case again - quantile value is lowest value
            idx2 <- which(phenoAvg == qx[1])
            cat(" 1: [", qx[1], ",", qx[1], "],", file=varlogfile, append=TRUE, sep="")

            idx3 <- which(phenoAvg > qx[1])
            cat(" 2: >", qx[1], " ||", file=varlogfile, append=TRUE, sep="")

        } else if (qx[1] == maxXX) {
            # Edge case again - quantile value is max value
            idx2 <- which(phenoAvg < qx[1] & phenoAvg > q[1])
            cat(" 1: (", q[1], ",", qx[1], "),", file=varlogfile, append=TRUE, sep="")
            idx3 <- which(phenoAvg == qx[1])
            cat(" 2: [",qx[1],",", qx[1], "] ||", file=varlogfile, append=TRUE, sep="")
        } else {
            idx2 <- which(phenoAvg < qx[1] & phenoAvg > q[1])
            cat(" 1: (", q[1], ",", qx[1], "),", file=varlogfile, append=TRUE, sep="")
            idx3 <- which(phenoAvg >= qx[1])
            cat(" 2: >=", qx[1], " ||", file=varlogfile, append=TRUE, sep="")
        }

        phenoBinned[idx2] <- 1
        phenoBinned[idx3] <- 2

    } else if (q[2] == maxX) {
    	
        # Edge case - quantile value is highest value
    	# Assign max value as cat3
    	idx3 <- which(phenoAvg == q[2])
    	phenoBinned[idx3] <- 2

    	# Divide remaining values into cat1 and cat2
    	phenoAvgRemaining <- phenoAvg[which(phenoAvg!=q[2])]
        qx <- quantile(phenoAvgRemaining, probs=c(0.5), na.rm=TRUE)
    	minXX <- min(phenoAvgRemaining, na.rm=TRUE)
        maxXX <- max(phenoAvgRemaining, na.rm=TRUE)

    	if (qx[1] == minXX) {
            # edge case again - quantile value is lowest value
            idx1 <- which(phenoAvg == qx[1])
            cat("0: [", qx[1], ",", qx[1],"],", file=varlogfile, append=TRUE, sep="") 
            idx2 <- which(phenoAvg > qx[1] & phenoAvg < q[2])
            cat(" 1: (", qx[1], ",", q[2], "),", file=varlogfile, append=TRUE, sep="") 
        } else if	(qx[1] == maxXX) {
            # edge case again - quantile value is max value
            idx1 <- which(phenoAvg < qx[1])
            cat("0: <", qx[1], ",", file=varlogfile, append=TRUE, sep="") 
            idx2 <- which(phenoAvg == qx[1])
            cat(" 1: [", qx[1], ",", qx[1], "] ", file=varlogfile, append=TRUE, sep="") 
        } else {
            idx1 <- which(phenoAvg < qx[1])  
            cat("0: <", qx[1], ",", file=varlogfile, append=TRUE, sep="") 
            idx2 <- which(phenoAvg >= qx[1] & phenoAvg < q[2])
            cat(" 1: [", qx[1], ",", q[2], "),", file=varlogfile, append=TRUE, sep="") 
    	}

        phenoBinned[idx1] <- 0
        phenoBinned[idx2] <- 1

        cat(" 2: [", q[2], ",", q[2],"] ||", file=varlogfile, append=TRUE, sep="")        

    } else if (q[1] == q[2]) {
        
        # Both quantiles correspond to the same value so set 
        # cat1 as < this value, cat2 as exactly this value and
        # cat3 as > this value
        
        phenoBinned <- phenoAvg
        idx1 <- which(phenoAvg < q[1])
        idx2 <- which(phenoAvg == q[2])
        idx3 <- which(phenoAvg > q[2])
        cat("0: <", q[1], ", 1: [", q[2], ",", q[2], "], 2: >", q[2]," ||", file=varlogfile, append=TRUE, sep="")    
        phenoBinned[idx1] <- 0
        phenoBinned[idx2] <- 1
        phenoBinned[idx3] <- 2

    } else {
        
        # Standard case - split the data into three roughly equal parts where
        # cat1 < q1, cat2 between q1 and q2, and cat3 >= q2

        phenoBinned <- phenoAvg
        idx1 <- which(phenoAvg < q[1])
        idx2 <- which(phenoAvg >= q[1] & phenoAvg < q[2])
        idx3 <- which(phenoAvg >= q[2])
        cat("0: <", q[1],", 1: [", q[1], ",", q[2], "), 2: >=", q[2]," ||", file=varlogfile, append=TRUE, sep="")    

        phenoBinned[idx1] <- 0
        phenoBinned[idx2] <- 1
        phenoBinned[idx3] <- 2
    }

    cat("cat N: ", length(idx1), ", ", length(idx2), ", ",
        length(idx3), " || ", sep="", file=varlogfile, append=TRUE)

    return(phenoBinned)
}
