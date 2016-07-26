from rpy2.robjects import r

from rpy2.robjects.packages import importr
importr("S4Vectors")
bioc = importr("BiocGenerics")

length_of_rle = r("function (x) sum(runLength(x))")


def make_treatment_control_same_length(treatment, control):

    treatment_result_dict = dict()
    control_result_dict = dict()
    for chromosome in treatment:

        chr_treatment, chr_control = treatment[chromosome], control[chromosome]

        treatment_maxlen = r["max"](r["sapply"](chr_treatment.values(),
                                                length_of_rle))
        control_maxlen = r["max"](r["sapply"](chr_control.values(),
                                              length_of_rle))

        maxlen = r["max"](treatment_maxlen, control_maxlen)

        lapply = r('function(cvg, maxlen) c(cvg,Rle(0,maxlen-length(cvg)))')

        treatment_result = {k: lapply(v, maxlen)
                            for (k, v) in chr_treatment.items()}

        control_result = {k: lapply(v, maxlen)
                          for (k, v) in chr_control.items()}

        treatment_result_dict[chromosome] = treatment_result
        control_result_dict[chromosome] = control_result

    return treatment_result_dict, control_result_dict
