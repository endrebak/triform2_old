
test.genome <- function(min.z=MIN.Z,
												min.shift=MIN.SHIFT,
												min.width=MIN.WIDTH) {
  INFO <<- NULL
  for(chr in CHRS) {
    cat("\n\tProcessing",chr,"... "); flush.console()
    INFO <<- rbind(INFO, test.chr(
									 chr=chr,
									 min.z=min.z,
									 min.shift=min.shift,
									 min.width=min.width))
    cat("done\n\t\tfound",N.PEAKS,"peaks")
  }
	INFO <<- INFO[order(-INFO$NLP),]
	nlps <- unique(INFO$NLP)
	sizes <- sapply(nlps,function(nlp) sum(INFO$NLP == nlp))
	indices <- sapply(nlps,function(nlp) sum(INFO$NLP >= nlp))

	nlrs <- mapply(function(nlp, j) {
		m <- sum(INFO$MAX.NLP >= nlp)	# Tarone modification for discrete nlp
		b.y <- log10(sum((1:m)^-1))		# discrete Benjamini-Yekutieli offset
		nls <- nlp + log10(j/m)				# discrete Benjamini-Hochberg adjustment
		max(nls-b.y,0)								# discrete Benjamini-Yekutieli adjustment
	}, nlp=nlps, j=indices)

	M <- length(nlrs)
	nlqs <- numeric(M)
	for(i in 1:M) nlqs[i] <- max(nlrs[i:M])		# step-up procedure
	nlqss <- unlist(mapply(function(nlq,size) rep(nlq,size),
									nlq=nlqs, size=sizes))

	INFO <<- cbind(QVAL=10^-nlqss, NLQ=nlqss, INFO)

  cat("\n\nSaving results ... "); flush.console()
  write.table(INFO,file="Triform.srf.xls",col.names=NA,quote=FALSE,sep="\t")
	save(INFO,file="Triform.srf.info.RData")
  cat("Finished.\n\n")
}
