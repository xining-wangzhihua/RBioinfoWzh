RPhylogenyZM

r funtions and scripts for phylogeny, written when the authour worked in nwipbcslzm.
the authour was a novice, so the robustness and effectiveness of the codes are not guaranteed.
there is no help file, nor documentation.

major functions are:
get.haplotype.and.genotype.from.dnasp.nex().r
generate.input.for.permut().r
generate.arp.input.for.arlequin().r

names of lower-level functions (usually) begin with "."

the package dependencies of these functions are:
ape, magrittr, seqinr, stringr, tibble
BiocManager, Biobase, Biostrings

codes to install these dependencies:
with(
  data = list(a = c("ape", "gdata", "magrittr", "readxl", "seqinr", "stringr", "tibble", "XML", "xml2"),
              b = c("Biobase", "Biostrings")),
  expr = {
    for(i in c(a,"BiocManager")) if(!require(i, character.only=TRUE)){install.packages(pkgs=i);}
    if(!require(BiocManager, character.only=FALSE)){stop("bug");}
    for(i in b) if(!require(i, character.only=TRUE)){BiocManager::install(pkgs=i);}
})
