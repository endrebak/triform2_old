library(IRanges)

"Converts the output from test.init to csv data."

load("temp/cvg_no_maxlen.RData")

files = names(CVG)

outfolder = "tests/test_data/"

for (f in files){
  print(f)
  if (!grepl("CENTER", f) && grepl("backgr", f)) next()
  rle = CVG[[f]]
  print(rle)

  rv = runLength(rle)
  rl = runValue(rle)
  df = data.frame(lengths=rv, values=rl)

  outfile = paste0(outfolder, f, ".csv")
  print(outfile)
  write.table(df, outfile, sep=" ")
}
