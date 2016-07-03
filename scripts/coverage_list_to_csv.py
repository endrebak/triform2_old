from os.path import abspath
from rpy2 import robjects as ro
from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr

importr("S4Vectors")

from sys import argv

rdata = abspath(argv[1])
print(rdata)

cmd = "load('" + rdata + "')"
print(cmd)
r(cmd)
chromosome_covers_list = r["CVG"]

ranges_to_bed = r("""ranges_to_bed <- function(gr, outpath){
  values = runValue(gr)
  lengths = runLength(gr)

  df = cbind(values, lengths)
  names(df) = c("Values", "Lengths")

  write.table(df, file=outpath, quote=F, sep=" ", row.names=F, col.names=T)
}
""")

names = r["names"](chromosome_covers_list)
for n, l in zip(names, chromosome_covers_list):
    if "backgr" in n and "CENTER" not in n:
        continue

    print("name:", n)
    name = n.lower().replace(".", "_") + ".csv"
    ranges_to_bed(l, name)
