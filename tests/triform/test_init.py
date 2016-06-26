## Use function below to create csv data of all the files! Will need it to crate
# test data for the rest of the sotware

# df <- data.frame(seqnames=seqnames(gr),
#   starts=start(gr)-1,
#   ends=end(gr),
#   names=c(rep(".", length(gr))),
#   scores=c(rep(".", length(gr))),
#   strands=strand(gr))

# write.table(df, file="foo.bed", quote=F, sep="\t", row.names=F, col.names=F)

# import pytest
# from rpy2.robjects import r, pandas2ri
# importr("GenomicRanges")

# from triform.init import _init

# @pytest.fixture
# def expected_result():

#     pass

# @pytest.fixture
# def input_data():

#     create_expected_result = r("""function() {
#   minus = "start end width
# 11743011 11743110 100
# 11936576 11936675 100
# 11764644 11764743 100"

#   plus = "start end width
# 10641343 10641442 100
# 12314811 12314910 100
# 57396665 57396764 100
# 11778648 11778747 100
# 10595033 10595132 100
# 16891529 16891628 100
# 11924940 11925039 100"

#     pcon <- textConnection(plus)
#     mcon <- textConnection(minus)

#     pdf = read.delim(pcon, sep=" ", header=1)
#     mdf = read.delim(mcon, sep=" ", header=1)

#     pi = IRanges(start=pdf$start, end=pdf$end, width=pdf$width)
#     mi = IRanges(start=mdf$start, end=mdf$end, width=mdf$width)

#     covers_expected_result = coverage(IRangesList("-"=mi, "+"=pi))
#     lsize = list("+"=7, "-"=3)
#     covers = list(SIZE=lsize, CVG=covers_expected_result)
# }
# """)
#     return create_expected_result()

# @pytest.mark.current
# def test_init(input_data, expected_result, args_tests):
#     result = _init(input_data, False, args_tests)

#     print("result")
#     print(result)
#     print("expected_result")
#     print(expected_result)

#     all_equals = r["isTRUE"](r["all.equal"](expected_result, result))
#     assert all_equals
