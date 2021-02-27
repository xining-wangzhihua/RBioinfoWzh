RPhylogenyZM

This readme file is obsolete.

R funtions and scripts for phylogeny.
Since the authour is a novice, the robustness and effectiveness of the codes are not guaranteed.
There is no help files, nor documentation.

major functions are:
get.haplotype.and.genotype.from.dnasp.nex().r
generate.input.for.permut().r
generate.arp.input.for.arlequin().r

the package dependencies of these functions are:
ape, magrittr, seqinr, stringr, tibble,
BiocManager, Biobase, Biostrings,
etc

codes to install these dependencies:
with(
  data = list(a = c("ape", "gdata", "magrittr", "readxl", "seqinr", "stringr", "tibble", "XML", "xml2"),
              b = c("Biobase", "Biostrings")),
  expr = {
    for(i in c(a,"BiocManager")) if(!require(i, character.only=TRUE)){install.packages(pkgs=i);}
    if(!require(BiocManager, character.only=FALSE)){stop("bug");}
    for(i in b) if(!require(i, character.only=TRUE)){BiocManager::install(pkgs=i);}
})
