#a function to determine "substitution type"
#return a vector, where each item is one of c("ts", "tv", "indel")
#similar to strataG::subType(), but better than it. and the "strataG" package is obsolete and archieved

#get dependency
if(FALSE){
  require(magrittr)
  require(gtools)
  source(file=R.home(component="RBioinfoZm/AnyEmpty().r"))
}

warning("currently subType2() contains a bug for more than 2 residues. the algorithm ",
        "only permute the linear processes, which should be network processes. for ",
        "example, for residues (a c g t -), the algorithm only enumerate the processes ",
        "in a line: \"a\" --> \"c\" --> \"g\" --> \"t\" --> \"-\", etc. ",
        "in fact, the real processes may be:\n", "\"a\" --> \"c\" --> \"g\"\n",
        "  \\       \\\n", "   \\       --> \"-\"\n", "    --> \"t\"\n",
        "learn the permutation of nodes in networks and revise the codes")

subType2 <- function(x){
  #control the input of x-------------------------------------------------------
  if(!is.character(x)){stop("invalid x");}
  if(AnyEmpty(x)){stop("invalid x");}
  x %<>% base::tolower() %>% base::unique() %>% base::sort()
  if(!all(x %in% c("a", "c", "g", "t", "-"))){stop("invalid x");}
  #give a waninging if length(x) > 2--------------------------------------------
  if(length(x) > 2){
    warning("when length(unique(x)) > 2, it's usually impossible to determine ",
            "what processes have happened. for example, if (a, c, -) occur at ",
            "one site, the process may be: \"a\" --indel--> \"-\" --indel--> ",
            "\"c\", or it may be: \"a\" --transvertion--> \"c\" --indel--> ",
            "\"-\". subType2() will only return un-avoidable processes, e.g. ",
            "\"indel\" in the former example.")
  }
  #main manipulation------------------------------------------------------------
  f_ans <- function(datum){
    if(length(datum) != 2){stop("invalid input");}
    datum <- c(a="purine", c="pyrimidine", g="purine", t="pyrimidine", `-`="hyphen") %>%
      .[datum] %>% unname() %>% base::unique()
    if("hyphen" %in% datum){result <- "indel";}else{
      #see K80 model description in ?ape::dist.dna()
      if(length(datum) == 1){result <- "ts";}
      if(length(datum) == 2){result <- "tv";}
    }
    return(result)
  }
  if(length(x) > 1){
    x <- gtools::permutations(n=length(x), r=length(x), v=x, repeats.allowed=FALSE)
    result <- base::vector(mode="list", length=nrow(x))
    for(i in 1:length(result)){
      result[[i]] <- x[i,] %>% {base::rbind(utils::head(., n=-1), .[-1])} %>%
        base::apply(MARGIN=2, FUN=f_ans)
    }; remove(i);
    result %<>% lapply(FUN=base::unique)
    if(length(result) > 1) for(i in 2:length(result)){
      result[[1]] %<>% base::intersect(x=., y=result[[i]])
    }; i <- 0; remove(i);
    result %<>% .[[1]]
  }else{
    result <- character()
  }
  remove(f_ans, x)
  #anchor-----------------------------------------------------------------------
  return(result)
}

#here are the test codes
if(FALSE){
  x <- c("A", "C", "g", "t", "-") %>% sample(x=., size=5, replace=TRUE); x;
  x %>% subType2()
}
