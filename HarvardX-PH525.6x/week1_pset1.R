# Genome-wide visualization of peak heights

library(ERBS)
library(ggbio)
library(GenomeInfoDb)
data("GM12878")
# the visualization of GM12878 ESRRA peaks
seqlevels(GM12878, force=T) = paste0("chr", 1:22)
autoplot(GM12878, layout="karyogram", aes(colour=log(peak)))

## An integrative visualization
# combines information on two cell lines and the measured peak values
data("HepG2")
HepG2$cell = "HepG2"
GM12878$cell = "Bcell"
tot = c(GM12878, HepG2)
tot$peak10 = tot$peak/10 # copes with automatic scale of y axis
seqlevels(tot, force=TRUE) = paste0("chr", 1:22)

library(scales)
p = autoplot(seqinfo(tot))
p = p + layout_karyogram(tot, aes(fill=cell, colour=cell), geom="rect") +
  scale_colour_manual(values = alpha(c("green", "red"), .1)) +
  scale_fill_manual(values = alpha(c("green", "red"), .1))
p + layout_karyogram(tot, aes(x=start, y=peak10), ylim=c(15,30),
                     geom="point", color="blue", size=.8)
stot = split(tot, as.character(seqnames(tot)))
w = sapply(stot, function(x) sum(width(x)))
sort(w/seqlengths(tot)[names(w)])


## multitrack visualization function
showz = function (sym = "ESRRA", radius = 1e+05) 
{
  require(GenomicRanges)
  require(ggbio)
  require(erma)
  require(ERBS)
  es = genemodel(sym)
  data(HepG2)
  data(GM12878)
  hsub = subsetByOverlaps(HepG2, es + radius)
  gsub = subsetByOverlaps(GM12878, es + radius)
  tracks(gene = es, hepNarr = autoplot(hsub), gmNarr = autoplot(gsub), 
         title = sym)
}
p = showz()
p

## incremental zooming into/out of display with ggbio
p+zoom(2)

## using the debugger with a troublesome gene
debug(showz)
showz("OCM", radius=2e6)




