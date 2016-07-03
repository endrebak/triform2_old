## Use function below to create csv data of all the files! Will need it to crate
# test data for the rest of the sotware

# want to convert list of 6 chrcovers to csv so can store for tests

import pandas as pd
import pytest
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
importr("GenomicRanges")
importr("S4Vectors")

from triform.init import init
from triform.helper_functions import df_to_rle, rle_to_df


@pytest.fixture
def expected_result_chip(init_result_chip):
    return {k: pd.read_table(v, sep=" ") for k, v in init_result_chip.items()}


@pytest.fixture
def expected_result_input(init_result_input):
    return {k: pd.read_table(v, sep=" ") for k, v in init_result_input.items()}


@pytest.mark.current
def test_init(expected_result_chip, expected_result_input, input_data_control,
              input_data_treatment):
    print(input_data_treatment.items(), 'input_data_treatment')
    print(input_data_treatment, 'input_data_treatment')
    assert 0


@pytest.fixture
def input_data_control(chromosome_cover_results_input):
    control = {k: df_to_rle(r["read.table"](v,
                                            sep=" "))
               for (k, v) in chromosome_cover_results_input.items()}
    return control


@pytest.fixture
def input_data_treatment(chromosome_cover_results_chip):
    treatment = {k: df_to_rle(r["read.table"](v,
                                              sep=" "))
                 for (k, v) in chromosome_cover_results_chip.items()}
    return treatment


@pytest.fixture
def expected_result_short(run_length_encodings):
    return {k: pd.read_table(v,
                             sep=" ",
                             header=0)
            for (k, v) in run_length_encodings.items()}


@pytest.fixture
def input_data_short():

    create_input_data = r("""function() {
  minus = "start end width
11743011 11743110 100
11936576 11936675 100
11764644 11764743 100"

  plus = "start end width
10641343 10641442 100
12314811 12314910 100
57396665 57396764 100
11778648 11778747 100
10595033 10595132 100
16891529 16891628 100
11924940 11925039 100"

    pcon <- textConnection(plus)
    mcon <- textConnection(minus)

    pdf = read.delim(pcon, sep=" ", header=1)
    mdf = read.delim(mcon, sep=" ", header=1)

    pi = IRanges(start=pdf$start, end=pdf$end, width=pdf$width)
    mi = IRanges(start=mdf$start, end=mdf$end, width=mdf$width)

    covers_expected_result = coverage(IRangesList("-"=mi, "+"=pi))
    lsize = list("+"=7, "-"=3)
    covers = list(SIZES=lsize, CVG=covers_expected_result)
}
""")
    return create_input_data()


@pytest.mark.unit
def test_init_short(input_data, expected_result, args):
    result = init(input_data, False, args)

    # convert result to python for easier comparison
    names = list(r["names"](result))
    result_dict = {}
    for name in names:
        l = pd.Series(pandas2ri.ri2py_intvector(r["runLength"](result.rx2(
            name))))
        v = pd.Series(pandas2ri.ri2py_intvector(r["runValue"](result.rx2(
            name))))
        df = pd.concat([v, l], axis=1)
        df.columns = ["values", "lengths"]
        result_dict[name] = df.astype(int)

    # check that the results are as expected
    assert len(result_dict) == len(expected_result)
    assert set(result_dict.keys()) == set(expected_result.keys())
    for k in result_dict:
        print("result_dict[k]")
        print(result_dict[k])
        print("expected_result[k]")
        print(expected_result[k])
        assert result_dict[k].equals(expected_result[k])
