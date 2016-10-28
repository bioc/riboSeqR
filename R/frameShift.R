frameShift <- function(...) {
  stop("This function is obsolete. Please use readingFrame instead.")
}

readingFrame <- function(rC, lengths = 26:30) {    

    frameCounts <- lapply(1:ncol(rC@hits), function(ii) {
        frameCounts <- sapply(lengths, function(length)
            sapply(0:2, function(frame) sum(rC@hits[,ii,as.character(frame), as.character(length)])))

        colnames(frameCounts) <- lengths
        rownames(frameCounts) <- 0:2
        
        frameCounts <- rbind(frameCounts, frame.ML = apply(frameCounts, 2, function(x) which.max(x)) - 1)
    })
                                        #frameShift <- which.max(frameCounts) - 1
                                        #c(frameShift, max(frameCounts) / sum(frameCounts))
    
  
  #colnames(frameShifts) <- as.character(lengths)
  #rownames(frameShifts) <- c("frameShift", "weighting")
                                        #frameShifts
    if(is.list(frameCounts) & length(frameCounts) == 1) frameCounts <- frameCounts[[1]]
    frameCounts
}


plotFS <- function(fS, lengths, legend.text = c("Frame 0", "Frame 1", "Frame 2"), ...)
    {

        if(!is.list(fS)) fS <- list(fS)
        if(is.null(names(fS))) names(fS) <- 1:length(fS)
        
        for(ff in 1:length(fS)) {
            if(!missing(lengths)) colsel <- which(colnames(fS[[ff]]) %in% as.character(lengths)) else colsel <- 1:ncol(fS[[ff]])   
            barplot(fS[[ff]][1:3, colsel], beside = TRUE, col = rainbow(3, s = 0.7), border = rainbow(3, s = 0.7), legend.text = legend.text, main = names(fS)[ff], ...)
        }
    
#    if(legend) legend(x = "topright", bty = "n", , fill = rainbow(3, s = 0.7), border = rainbow(3, s = 0.7), )
  }
  
