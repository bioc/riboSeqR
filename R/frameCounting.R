.oldframeCounting <-
function(riboDat, fastaCDS, lengths = 25:30)
  {
    message("Calling frames...", appendLF = FALSE)

    if(!("frame" %in% names(values(fastaCDS)))) fastaCDS$frame <- (start(fastaCDS) - 1) %% 3
    
    fr0GR <- fastaCDS[fastaCDS$frame == 0]
    fr1GR <- fastaCDS[fastaCDS$frame == 1]
    fr2GR <- fastaCDS[fastaCDS$frame == 2]

    frameCalls <- lapply(riboDat@riboGR, function(gral, lengths) {
#      gral$frame <- (start(gral) - 2) %% 3
      gral$frame <- (start(gral) - 1 - 0) %% 3
      if(length(fr0GR) > 0) fr0 <- .tableOverlaps(fr0GR, gral, lengthRange = lengths) else fr0 <- new("riboCoding", hits = array(dim=c(0,1,3,length(lengths))), unqHits = array(dim=c(0,1,3,length(lengths))))
      message(".", appendLF = FALSE)
      gral$frame <- (start(gral) - 1 - 1) %% 3
      if(length(fr1GR) > 0) fr1 <- .tableOverlaps(fr1GR, gral, lengthRange = lengths) else fr1 <- new("riboCoding", hits = array(dim=c(0,1,3,length(lengths))), unqHits = array(dim=c(0,1,3,length(lengths))))
      message(".", appendLF = FALSE)
      gral$frame <- (start(gral) - 1 - 2) %% 3
      if(length(fr2GR) > 0) fr2 <- .tableOverlaps(fr2GR, gral, lengthRange = lengths) else fr2 <- new("riboCoding", hits = array(dim=c(0,1,3,length(lengths))), unqHits = array(dim=c(0,1,3,length(lengths))))
      message(".", appendLF = FALSE)
      rc <- new("riboCoding",
          CDS = c(fr0@CDS, fr1@CDS, fr2@CDS),
          hits = do.call("abind", args = list(list(fr0@hits, fr1@hits, fr2@hits), along = 1)),
          unqHits = do.call("abind", args = list(list(fr0@unqHits, fr1@unqHits, fr2@unqHits), along = 1)))
      rc
    }, lengths = lengths)
    message(".done!", appendLF = TRUE)
    fCs <- new("riboCoding")
    fCs@hits <- do.call("abind", args = list(lapply(frameCalls, function(x) x@hits), along = 2))
    fCs@unqHits <- do.call("abind", args = list(lapply(frameCalls, function(x) x@unqHits), along = 2))
    fCs@CDS <- frameCalls[[1]]@CDS
    fCs@replicates <- riboDat@replicates
    fCs
  }


frameCounting <- function(riboDat, cds, lengths = 25:30, offset5p = 0, offset3p = 0)
  {
      message("Calling frames...", appendLF = FALSE)

      getHits <- function(rgr, offset5p, offset3p) {
          message(".", appendLF = FALSE)
          rgr <- rgr[width(rgr) %in% lengths]
          mcds <- cds; start(mcds) <- pmax(0, start(mcds) - offset5p); end(mcds) <- pmax(start(mcds), end(mcds) - offset3p)
          fo <- findOverlaps(mcds, rgr)
          dat <- cbind.data.frame(cds = factor(queryHits(fo), levels = 1:length(cds)),
                                  frame = factor((start(rgr)[subjectHits(fo)] - start(cds[queryHits(fo)])) %% 3, levels = 0:2),
                                  length = factor(width(rgr)[subjectHits(fo)], levels = lengths))
                           
          
          list(hits = table(dat), unqHits = table(dat[!duplicated(rgr)[subjectHits(fo)],]))
      }

      riboCounts <- lapply(riboDat@riboGR, getHits, offset5p = offset5p, offset3p = offset3p)
      #rnaCounts <- lapply(riboDat@rnaGR, getHits, offset5p = offset5p, offset3p = offset3p)
      message("done!")
      
      fCs <- new("riboCoding")
      fCs@hits <- aperm(do.call("abind", args = list(lapply(riboCounts, function(x) x$hits), along = 0)), c(2,1,3,4))
      fCs@unqHits <- aperm(do.call("abind", args = list(lapply(riboCounts, function(x) x$unqHits), along = 0)), c(2,1,3,4))
      if(length(dim(fCs@hits)) == 3) fCs@hits <- array(fCs@hits, c(dim(fCs@hits)[1], 1, dim(fCs@hits)[2:3]))
      if(length(dim(fCs@unqHits)) == 3) fCs@unqHits <- array(fCs@unqHits, c(dim(fCs@unqHits)[1], 1, dim(fCs@unqHits)[2:3]))

      dimnames(fCs@hits) <- list(NULL, names(riboDat@riboGR), 0:2, lengths)
      dimnames(fCs@unqHits) <- list(NULL, names(riboDat@riboGR), 0:2, lengths)
      
      fCs@CDS <- cds
      fCs@replicates <- riboDat@replicates      
      
      fCs
   
  }      

