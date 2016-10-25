logoContext <- function(cds, ...) {
    baseTab <- apply(
        do.call("rbind", lapply(strsplit(cds$context, ""), function(x) c(rep(NA, 7 -length(x)), x)))
      , 2, function(x) table(factor(x, levels = c("A", "C", "G", "T"))))
    p <- makePWM(t(t(baseTab) / colSums(baseTab)))
    seqLogo(p, ...)
    invisible(baseTab)
}
