import pytest

from collections import OrderedDict

from triform.preprocess.bed_to_chromosome_dfs import bed_to_chromosome_dfs


@pytest.fixture
def expected_result():
    return OrderedDict([(
        'chr22',
        'chr22\t10592905\t10592930\t-\nchr22\t12968115\t12968140\t-\nchr22\t11926812\t11926837\t-\nchr22\t10633359\t10633384\t+\nchr22\t12128682\t12128707\t-\nchr22\t11942644\t11942669\t+\nchr22\t11938603\t11938628\t+\nchr22\t10646972\t10646997\t-\nchr22\t11251733\t11251758\t-\nchr22\t12968115\t12968140\t-\n'
    ), ('chrY',
        'chrY\t10641380\t10641405\t+\nchrY\t12314848\t12314873\t+\nchrY\t11743048\t11743073\t-\nchrY\t57396702\t57396727\t+\nchrY\t11778685\t11778710\t+\nchrY\t10595070\t10595095\t+\nchrY\t16891566\t16891591\t+\nchrY\t11936613\t11936638\t-\nchrY\t11924977\t11925002\t+\nchrY\t11764681\t11764706\t-\n'
        )])


@pytest.mark.unit
def test_bed_to_chromosome_dfs(args_tests, expected_result):
    result = bed_to_chromosome_dfs("examples/multiple_chromos.bed", args_tests)
    print("result")
    for k, v in result.items():
        print(k)
        print(v)

    print("expected_result")
    for k, v in expected_result.items():
        print(k)
        print(v)

    assert result == expected_result
