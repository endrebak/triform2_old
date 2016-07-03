library(IRanges)

"Converts the output from makeChromosomeCoverFiles to csv data."

files = c("chrcovers/chrY_srf_huds_Gm12878_rep1.RData",
          "chrcovers/chrY_backgr_huds_Gm12878_rep1.RData",
          "chrcovers/chrY_srf_huds_Gm12878_rep2.RData",
          "chrcovers/chrY_backgr_huds_Gm12878_rep2.RData")

outfolder = "tests/test_data/"

for (f in files) {
  print(f)
  load(f)

  cvg.pos = covers$CVG$"+"
  rv.pos = runLength(cvg.pos)
  rl.pos = runValue(cvg.pos)
  pos = data.frame(lengths=rv.pos, values=rl.pos)
  outfile.pos = gsub(".RData", "_pos.csv", gsub("chrcovers/", outfolder, f))
  write.table(pos, outfile.pos, sep=" ")

  cvg.neg = covers$CVG$"-"
  rv.neg = runLength(cvg.neg)
  rl.neg = runValue(cvg.neg)
  neg = data.frame(lengths=rv.neg, values=rl.neg)
  outfile.neg = gsub(".RData", "_neg.csv", gsub("chrcovers/", outfolder, f))
  write.table(neg, outfile.neg, sep=" ")

}
