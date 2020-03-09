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
{
  install.packages(pkgs = c("ape", "magrittr", "seqinr", "stringr", "tibble"))
  install.packages(pkgs = "BiocManager")
  BiocManager::install(pkgs = c("Biobase", "Biostrings"))
}
