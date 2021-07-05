RBioinfoWzh

1.
function definations end with "().r", scripts end with "_not_parenthesis.r", expressions begin with "expr "
major functions begin with "step?."
functions that end with "2" are similar functions to original ones
functions that start with "widget." are used to achieve small functions. They are mainly used in scripts, and may be obsolete or replaced quickly.

2.
r codes to install these dependencies:
stop("need change"); with(
  data <- list(a=c("ape", "gdata", "magrittr", "readxl", "seqinr", "stringr", "tibble", "XML", "xml2"),
               b=c("Biobase", "Biostrings")),
  expr <- {
    #use lapply instead of for loops to install packages
    for(i in c(a,"BiocManager")) if(!require(i, character.only=TRUE)){install.packages(pkgs=i);}
    if(!require(BiocManager, character.only=FALSE)){stop("bug");}
    for(i in b) if(!require(i, character.only=TRUE)){BiocManager::install(pkgs=i);}
})

3.
to use RBioinfoWzh, put the folder in the directory where you install r, and source the script "script - RBioinfoWzh.r"
for example, if you install r at "D:/softwares/R-3.6.1", then move the directory there, and source "D:/softwares/R-3.6.1/RBioinfoWzh/script - RBioinfoWzh.r"

5.
for step2.GetCheckedSequences()
#ratio = 0 means there is no ambiguity
#ratio = 1 means which peak of the two bases are higher is unpredictable
#ratio = 2 means generally speaking the area size under one base peak is twice as much as the area size under another peak, i.e. the primary peak is obviously higher than the secondary one, and the secondary peak also exist obviously
#ratio = 3 means generally speaking, the primary peak is obviously higher than the secondary one, but the residue represented by the secondary peak may exist or not
#ratio = 4 usually is not used at first, and preserved for revision/collation when checking the rare haplotypes in the result of haplotype estimation. (sometimes ratio = 3 can also be used for revision/collation during haplotype estimation)
#ratio = 5, 6, 7, etc means the secondary peaks is getting more lower

6.
three types of functions: dependency, standalone functions (widget, etc), functions with dependencies.
