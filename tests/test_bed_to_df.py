import pytest

from subprocess import check_output

from joblib import Parallel, delayed

from triform.config.genomes import create_genome_size_dict


@pytest.fixture
def expected_result():
    # yapf: disable
    return {'chr22': 'chr22\t10592905\t10592930\t-\nchr22\t12968115\t12968140\t-\nchr22\t11926812\t11926837\t-\nchr22\t10633359\t10633384\t+\nchr22\t12128682\t12128707\t-\nchr22\t11942644\t11942669\t+\nchr22\t11938603\t11938628\t+\nchr22\t10646972\t10646997\t-\nchr22\t11251733\t11251758\t-\nchr22\t12968115\t12968140\t-\n', 'chrY': 'chrY\t10641380\t10641405\t+\nchrY\t12314848\t12314873\t+\nchrY\t11743048\t11743073\t-\nchrY\t57396702\t57396727\t+\nchrY\t11778685\t11778710\t+\nchrY\t10595070\t10595095\t+\nchrY\t16891566\t16891591\t+\nchrY\t11936613\t11936638\t-\nchrY\t11924977\t11925002\t+\nchrY\t11764681\t11764706\t-\n'}
# yapf: enable


def test_name(args_tests, expected_result):
    result = bed_to_chromosome_dfs("examples/multiple_chromos.bed", args_tests)
    result == expected_result


@pytest.mark.unit
def bed_to_chromosome_dfs(bed_file, args):

    chromosomes = create_genome_size_dict(args.genome)

    text_dfs = Parallel(n_jobs=args.number_cores)(
        delayed(_bed_to_chromosomes)(bed_file, chromosome)
        for chromosome in chromosomes)

    text_dfs = {c: df for (c, df) in zip(chromosomes, text_dfs) if df != ""}
    return text_dfs


def _bed_to_chromosomes(bed_file, chromosome):
    command = "grep -E '^{chromosome}\\b' {bed_file} | cut -f 1-3,6".format(
        chromosome=chromosome,
        bed_file=bed_file)
    output = check_output(command, shell=True)

    return output
