libScales <- function(rC, riboDat, lengths, frames, method = "edgeR") {
    riboC <- sliceCounts(rC, lengths = lengths, frames = frames)
    if(length(riboDat@rnaGR) > 0) {
        rnaC <- rnaCounts(riboDat, rC@CDS)
        rnaLS = getLibsizes(data = rnaC, replicates = riboDat@replicates, estimationType = "edgeR")
    } else rnaLS <- NULL
    riboLS = getLibsizes(data = riboC, replicates = riboDat@replicates, estimationType = "edgeR")    
    list(riboLS = riboLS, rnaLS = rnaLS)
}
