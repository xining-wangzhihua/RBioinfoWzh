RPhylogenyZM

R funtions and scripts for phylogeny.
The robustness and effectiveness of the codes are not guaranteed.

major functions are:
step*().r and {1~5}*().r

the package dependencies of these functions are (not completed):
magrittr, tibble, dplyr, stringr, purrr,
BiocManager, IRanges, Biostrings, ape, seqinr, pegas.

codes to install these dependencies (not completed):
with(
  data = list(a = c("magrittr", "tibble", "dplyr", "stringr", "purrr"),
              b = c("IRanges", "Biostrings")),
  expr = {
    for(i in c(a,"BiocManager")) if(!require(i, character.only=TRUE)){install.packages(pkgs=i);}
    if(!require(BiocManager, character.only=FALSE)){stop("bug");}
    for(i in b) if(!require(i, character.only=TRUE)){BiocManager::install(pkgs=i);}
})
