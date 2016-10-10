import logging
from subprocess import call
from os.path import dirname, join

import pandas as pd

from joblib import Parallel, delayed

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
pandas2ri.activate()

from triform.matrix.find_read_midpoints import find_read_midpoints
from triform.matrix.fdr_table_to_bins import fdr_table_to_bins


def _give_sequences_same_length(seqs):

    max_pos = 0
    for seq in seqs:
        max_pos = max(r["seqlengths"](seq)[0], max_pos)

    add_seqlength = r("function(x, l) { seqlengths(x) <- c(l); x}")

    new_maxlengths = []
    for seq in seqs:
        new_maxlengths.append(add_seqlength(seq, max_pos))

    return new_maxlengths


def create_matrix(treatment_iranges, control_iranges, fdr_table, args):

    enriched_bins = fdr_table_to_bins(fdr_table, args.matrix_bin_size)
    logging.info("Creating treatment matrix.")
    treatment_matrixes = find_read_midpoints(treatment_iranges, args)

    logging.info("Creating control matrix.")
    control_matrixes = find_read_midpoints(control_iranges, args)
    chromosomes = set(treatment_matrixes).intersection(control_matrixes)

    directory = dirname(args.matrix)
    if directory:
        call("mkdir -p {}".format(directory), shell=True)

    logging.info("Merging matrixes")
    matrix_outfiles = Parallel(n_jobs=args.number_cores)(
        delayed(_create_matrix)(
            treatment_matrixes[chromosome], control_matrixes[chromosome],
            enriched_bins[chromosome], directory, chromosome)
        for chromosome in chromosomes)

    matrix_outfiles = sorted(matrix_outfiles)

    matrix_outfiles_str = " ".join(matrix_outfiles)
    logging.info("Writing matrix to {}".format(args.matrix))
    call("cat {} > {}".format(matrix_outfiles_str, args.matrix), shell=True)
    call("rm {}".format(matrix_outfiles_str), shell=True)


def _create_matrix(treatment_matrixes, control_matrixes, enriched_bins, outdir,
                   chromosome):

    matrixes = treatment_matrixes + control_matrixes

    fnames = [f for f, _ in matrixes]
    seqs = [s for _, s in matrixes]

    if len(seqs) > 1:
        seqs = _give_sequences_same_length(seqs)

    u = r["unique"](r["c"](*[r["granges"](s) for s in seqs]))

    _add_element_metadata = r("""function(gr, u, name) {
        mcols(u)[match(gr, u), name] = 0
        mcols(u)[match(gr, u), name] = elementMetadata(gr)
        u
        }""")

    for name, seq in zip(fnames, seqs):
        u = _add_element_metadata(seq, u, name)

    coords = r(
        """function(u) data.frame(Chromosome=seqnames(u), Start=start(u), End=end(u))""")(
            u)
    counts = r["as.data.frame"](r["elementMetadata"](u))
    count_df = r["cbind"](coords, counts)
    cdf = ri2py(count_df).set_index(["Chromosome"])
    cdf = cdf.fillna(0).astype(int)
    cdf = cdf.set_index(["Start", "End"], append=True)

    if enriched_bins.empty:
        enriched_bins = pd.DataFrame(index=cdf.index)
        enriched_bins.insert(0, "Enriched", 0)

    cdf = cdf.join(enriched_bins).fillna(0)
    cdf = cdf.reset_index("End", drop=True).reset_index()

    new_columns = list(cdf.columns)
    new_columns[1] = "Bin"
    cdf.columns = new_columns
    cdf = cdf.set_index("Chromosome Bin Enriched".split())
    cdf = cdf[~cdf.index.duplicated()]

    outpath = join(outdir, chromosome + ".gz")
    cdf.to_csv(outpath, sep=" ", compression="gzip")

    return outpath

# ...........................................................................
# /local/home/endrebak/anaconda3/lib/python3.5/site-packages/pandas/core/frame.py in join(self=                              data.animal.Exp2_9...                     1

# [8738 rows x 4 columns], other=Empty DataFrame
# Columns: []
# Index: [], on=None, how='left', lsuffix='', rsuffix='', sort=False)
#    4380         -------
#    4381         joined : DataFrame
#    4382         """
#    4383         # For SparseDataFrame's benefit
#    4384         return self._join_compat(other, on=on, how=how, lsuffix=lsuffix,
# -> 4385                                  rsuffix=rsuffix, sort=sort)
#         rsuffix = ''
#         sort = False
#    4386
#    4387     def _join_compat(self, other, on=None, how='left', lsuffix='', rsuffix='',
#    4388                      sort=False):
#    4389         from pandas.tools.merge import merge, concat

# ...........................................................................
# /local/home/endrebak/anaconda3/lib/python3.5/site-packages/pandas/core/frame.py in _join_compat(self=                              data.animal.Exp2_9...                     1

# [8738 rows x 4 columns], other=Empty DataFrame
# Columns: []
# Index: [], on=None, how='left', lsuffix='', rsuffix='', sort=False)
#    4394             other = DataFrame({other.name: other})
#    4395
#    4396         if isinstance(other, DataFrame):
#    4397             return merge(self, other, left_on=on, how=how,
#    4398                          left_index=on is None, right_index=True,
# -> 4399                          suffixes=(lsuffix, rsuffix), sort=sort)
#         lsuffix = ''
#         rsuffix = ''
#         sort = False
#    4400         else:
#    4401             if on is not None:
#    4402                 raise ValueError('Joining multiple DataFrames only supported'
#    4403                                  ' for joining on index')

# ...........................................................................
# /local/home/endrebak/anaconda3/lib/python3.5/site-packages/pandas/tools/merge.py in merge(left=                              data.animal.Exp2_9...                     1

# [8738 rows x 4 columns], right=Empty DataFrame
# Columns: []
# Index: [], how='left', on=None, left_on=None, right_on=None, left_index=True, right_index=True, sort=False, suffixes=('', ''), copy=True, indicator=False)
#      34           suffixes=('_x', '_y'), copy=True, indicator=False):
#      35     op = _MergeOperation(left, right, how=how, on=on, left_on=left_on,
#      36                          right_on=right_on, left_index=left_index,
#      37                          right_index=right_index, sort=sort, suffixes=suffixes,
#      38                          copy=copy, indicator=indicator)
# ---> 39     return op.get_result()
#         op.get_result = <bound method _MergeOperation.get_result of <pandas.tools.merge._MergeOperation object>>
#      40 if __debug__:
#      41     merge.__doc__ = _merge_doc % '\nleft : DataFrame'
#      42
#      43

# ...........................................................................
# /local/home/endrebak/anaconda3/lib/python3.5/site-packages/pandas/tools/merge.py in get_result(self=<pandas.tools.merge._MergeOperation object>)
#     212     def get_result(self):
#     213         if self.indicator:
#     214             self.left, self.right = self._indicator_pre_merge(
#     215                 self.left, self.right)
#     216
# --> 217         join_index, left_indexer, right_indexer = self._get_join_info()
#         join_index = undefined
#         left_indexer = undefined
#         right_indexer = undefined
#         self._get_join_info = <bound method _MergeOperation._get_join_info of <pandas.tools.merge._MergeOperation object>>
#     218
#     219         ldata, rdata = self.left._data, self.right._data
#     220         lsuf, rsuf = self.suffixes
#     221

# ...........................................................................
# /local/home/endrebak/anaconda3/lib/python3.5/site-packages/pandas/tools/merge.py in _get_join_info(self=<pandas.tools.merge._MergeOperation object>)
#     334         left_ax = self.left._data.axes[self.axis]
#     335         right_ax = self.right._data.axes[self.axis]
#     336
#     337         if self.left_index and self.right_index:
#     338             join_index, left_indexer, right_indexer = \
# --> 339                 left_ax.join(right_ax, how=self.how, return_indexers=True)
#         left_ax.join = <bound method Index.join of MultiIndex(levels=[[...           names=['Chromosome', 'Start', 'End'])>
#         right_ax = Index([], dtype='object')
#         self.how = 'left'
#     340         elif self.right_index and self.how == 'left':
#     341             join_index, left_indexer, right_indexer = \
#     342                 _left_join_on_index(left_ax, right_ax, self.left_join_keys,
#     343                                     sort=self.sort)

# ...........................................................................
# /local/home/endrebak/anaconda3/lib/python3.5/site-packages/pandas/indexes/base.py in join(self=MultiIndex(levels=[['chrY'], [310831, 386851, 69...
#            names=['Chromosome', 'Start', 'End']), other=Index([], dtype='object'), how='left', level=None, return_indexers=True)
#    2440             # have the same levels/names so a simple join
#    2441             if self.names == other.names:
#    2442                 pass
#    2443             else:
#    2444                 return self._join_multi(other, how=how,
# -> 2445                                         return_indexers=return_indexers)
#         return_indexers = True
#    2446
#    2447         # join on the level
#    2448         if level is not None and (self_is_mi or other_is_mi):
#    2449             return self._join_level(other, level, how=how,

# ...........................................................................
# /local/home/endrebak/anaconda3/lib/python3.5/site-packages/pandas/indexes/base.py in _join_multi(self=MultiIndex(levels=[['chrY'], [310831, 386851, 69...
#            names=['Chromosome', 'Start', 'End']), other=Index([], dtype='object'), how='left', return_indexers=True)
#    2532         other_names = [n for n in other.names if n is not None]
#    2533         overlap = list(set(self_names) & set(other_names))
#    2534
#    2535         # need at least 1 in common, but not more than 1
#    2536         if not len(overlap):
# -> 2537             raise ValueError("cannot join with no level specified and no "
#    2538                              "overlapping names")
#    2539         if len(overlap) > 1:
#    2540             raise NotImplementedError("merging with more than one level "
#    2541                                       "overlap on a multi-index is not "

# ValueError: cannot join with no level specified and no overlapping names
# ___________________________________________________________________________
# Error in job run_triform while creating output files data/chip_caller/10/9h_PolII.regions, data/chip_caller/10/9h_PolII_matrix.csv.gz.
# RuleException:
# CalledProcessError in line 32 of /local/home/endrebak/code/programmable_epigenetics/rules/triform/triform.rules:
# Command 'triform2 -cpu 25 --matrix-bin-size 10 -m data/chip_caller/10/9h_PolII_matrix.csv.gz -b data/chip_caller/10_9h_PolII/ -t data/animal/Exp1_9h_PolII.bed data/animal/Exp2_9h_PolII.bed -c data/animal/Exp1_9h_Input.bed data/animal/Exp2_9h_Input.bed -o data/chip_caller/10/9h_PolII.regions' returned non-zero exit status 1
#   File "/local/home/endrebak/code/programmable_epigenetics/rules/triform/triform.rules", line 32, in __rule_run_triform
