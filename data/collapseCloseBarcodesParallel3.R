#
# Need to develop a script to run on SciNet that will go through 
# a matrix of barcode data, and collapse barcode records that
# have a Hamming distance <= 2.
#
# June 5, 2019
# Kevin R Brown
######################################

# Parameters passed in:
# 1.  Input filename
# 2.  Output dir
args <- commandArgs(trailingOnly = T);

if (length(args) < 2) {
	print("Usage:  Rscript collapseCloseBarcodes.R <inFile> <outfile>");
	quit();
}

iFile <- args[1];
oFile <- args[2];

# Load data
data <- read.csv(iFile,header=T,sep="\t",stringsAsFactors=F);

# Define collapsing function.
mergeRows <- function(x) {

   m <- max(unlist(x[,-1]));
	
   # Only collapse rows if the MAX value is > some constant
   # MAX > 3
   if (( nrow(x) > 1) & (m > 1)) {
      # Find column that had max value
      # col <- which.max(apply(x[,-1],2,max))+1;
      col <- 2;

      # Reorder matrix by column with max value     ### MATRIX ALREADY ORDERED BEFORE CALL TO mergeRows() - DON'T REDO!
#		x <- x[ order(x[,col],decreasing=T), ];
      bc <- x[1,1];   # Retain most prominent barcode
		
      x2 <- sum(x[,-1]);
      x2 <-  data.frame(Barcode=bc,x2,stringsAsFactors=F)
		
   } else {
      x2 <- x;
   }

    colnames(x2) <- colnames(x);
    x2;
	
}

# Reorder data matrix
o <- order(data[,2],decreasing=T);
data <- data[ o,];

grps <- rep(0,nrow(data));
grp <- 1;

library(stringdist);
for (i in 1:nrow(data)) {

	if ((i %% 100) == 0) {
		print(i);
	}

    if (grps[i] != 0) {
    	next;
    }
    
    idxs <- which( grps == 0 );
    y <- stringdist(data[idxs[1],1],data[idxs,1],method="hamming",useBytes=T,nthread=16);
    ratios <- data[idxs,2] / data[idxs[1],2];
#    rows <- idxs[ y <= 2 ];
    rows <- idxs[ y == 0 | (y == 1 & ratios < 0.125) | (y == 2 & ratios < 0.025) ];
    grps[rows] <- grp;
    grp <- grp + 1;

#    print(paste("Remaining rows:",sum( grps == 0 )));    
}

#tmp <- do.call("rbind",tapply(1:nrow(data),grps,FUN=function(x,data) { mergeRows(data[x,],paste(oDir,"collapsing_test_log_version3.txt",sep="/")); },data=data ));

require(doMC);
registerDoMC(16);
tmp <- foreach(grp=unique(grps), .combine=rbind) %dopar% {
   mergeRows( data[grps == grp, ] );
}

write.table(tmp,oFile,sep="\t",row.names=F,quote=F);

