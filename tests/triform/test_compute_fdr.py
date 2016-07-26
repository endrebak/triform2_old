import pytest

import pandas as pd
import numpy as np

from io import StringIO


@pytest.fixture
def input_data():

    data = """    NLP  MAX.NLP      LOC WIDTH    START      END CVG SURL SURR FORM
chrY:57420858-57420942:1 3.923426 6.317243 57420916    85 57420858 57420942  10  0    2    1
chrY:10636518-10636576:2 2.869699 2.869699 10636545    59 10636518 10636576   9  0   16    2
chrY:10631247-10631396:3 4.039029 4.039029 10631334   150 10631247 10631396  14   24    0    3"""

    df = pd.read_table(data, index_col=0)

    return data


@pytest.fixture
def expected_result():
    data = """QVAL	NLQ	NLP	MAX.NLP	LOC	WIDTH	START	END	CVG	SURL	SURR	FORM
chrY:57420858-57420942:1	0.000178922590521532	3.74733482254352	3.9234260815992	6.31724273434332	57420916	85	57420858	57420942	10	0	2	1
chrY:10631247-10631396:3	0.000178922590521532	3.74733482254352	4.03902854326486	4.03902854326486	10631334	150	10631247	10631396	14	24	0	3
chrY:10636518-10636576:2	0.00247481305798851	2.60645760115479	2.86969903592937	2.86969903592937	10636545	59	10636518	10636576	9	0	16	2"""

    df = pd.read_table(data, index_col=0)

    return df


@pytest.mark.unit
def compute_fdr(df):
    """
    """

# INFO <<- INFO[order(-INFO$NLP),]
# nlps <- unique(INFO$NLP)
# sizes <- sapply(nlps,function(nlp) sum(INFO$NLP == nlp))
# indices <- sapply(nlps,function(nlp) sum(INFO$NLP >= nlp))

# nlrs <- mapply(function(nlp, j) {
# 	m <- sum(INFO$MAX.NLP >= nlp)	# Tarone modification for discrete nlp
# 	b.y <- log10(sum((1:m)^-1))		# discrete Benjamini-Yekutieli offset
# 	nls <- nlp + log10(j/m)				# discrete Benjamini-Hochberg adjustment
# 	max(nls-b.y,0)								# discrete Benjamini-Yekutieli adjustment
# }, nlp=nlps, j=indices)

# M <- length(nlrs)
# nlqs <- numeric(M)
# for(i in 1:M) nlqs[i] <- max(nlrs[i:M])		# step-up procedure
# nlqss <- unlist(mapply(function(nlq,size) rep(nlq,size),
# 								nlq=nlqs, size=sizes))

# INFO <<- cbind(QVAL=10^-nlqss, NLQ=nlqss, INFO)


def test_compute_fdr(input_data, expected_result):

    result = compute_fdr(input_data)
    assert result == expected_result
