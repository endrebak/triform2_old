import pandas as pd
from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py


def compute_fdr(dfs):

    new_dfs = []
    for i, df in enumerate(dfs):
        if df is None:
            continue

        new_dfs.append(df)

    df = r["do.call"]("rbind", new_dfs)
    if r["is.null"](df)[0]:
        return pd.DataFrame()

    _compute_fdr = r("""function(INFO) {
    INFO <- INFO[order(-INFO$NLP),]
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

    cbind(QVAL=10^-nlqss, NLQ=nlqss, INFO)
    }""")

    df = ri2py(_compute_fdr(df))
    idx = df.index.get_level_values(0).to_series().astype(str).str.split(
        ":",
        expand=True).ix[:, 0].to_frame()

    df = idx.join(df).reset_index(drop=True)
    df.columns = ["CHROMOSOME"] + list(df.columns)[1:]
    df = df.set_index("CHROMOSOME START END".split())

    return df.sort_values(["QVAL", "FORM"])
