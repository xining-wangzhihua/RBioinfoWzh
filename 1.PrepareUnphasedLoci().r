#! /usr/bin/Rscript
# license: BSD 3-Clause License
# author: https://github.com/ywd5
# version: 20210314

warning("ensure \"if there are only 2 sequences, a.primary and a.secondary\" won't generate bugs")

if("" == "get dependencies"){
  require(magrittr); require(tibble); require(dplyr); require(stringr);
  require(purrr)#; require(rlang);
  require(IRanges); require(Biostrings); require(ape); require(pegas);
}

if("" == "test codes"){
  readDNAStringSet(filepath="../test/its ps.fasta") %>% PrepareUnphasedLoci()
}

PrepareUnphasedLoci <- function(x, suffix_pattern="\\.primary$|\\.secondary$|\\.hap[[:digit:]]+$"){
  # suffix_pattern can to be "", which will leave sequence names unchanged
  # control the input of suffix_pattern-------------------------------------------------------------
  if(!is.character(suffix_pattern)){stop("suffix_pattern must be a string");}
  if(length(suffix_pattern) != 1){stop("length(suffix_pattern) must be 1");}
  if(is.na(suffix_pattern)){stop("suffix_pattern can't be NA");}
  # control the input of x--------------------------------------------------------------------------
  if(!is(object=x, class2="DNAStringSet")){stop("x must be a DNAStringSet object defined in Biostrings");}
  x %>% DNAMultipleAlignment() %>% {dim(.) == 0} %>% {if(any(.)){stop("x can't be empty");}}
  unname(x) %>% unique() %>% {length(.) == 1} %>%
    {if(.){stop("all sequences in x are same. no need to analyse");}}
  names(x) %>% {if(is.null(.)){stop("sequence names must be speicified");}}
  names(x) %>% {if(anyNA(.)){stop("sequence names must be speicified");}}
  names(x) %<>% stringr::str_squish(string=.)
  names(x) %>% {any(. == "")} %>% {if(.){stop("sequence names must be specified");}}
  # sequence names can contain duplicated items, i.e. looks like the suffixes are already removed
  # control the relationship between x and suffix_pattern-------------------------------------------
  if(suffix_pattern != ""){
    base::grepl(pattern=suffix_pattern, x=names(x)) %>%
      {if(!all(.)){stop("some sequence names don't match suffix_pattern");}}
    names(x) %<>% base::sub(pattern=suffix_pattern, replacement="", x=.)
  }
  names(x) %>% table() %>% {any(. != 2)} %>% {if(.){stop("each sequence must have 2 copies");}}
  remove(suffix_pattern)
  x %<>% {.[stringr::str_order(x=names(.), numeric=TRUE)];}
  # get info----------------------------------------------------------------------------------------
  warning("write test codes here to check the behavior of ape::sge.sites()")
  # here ambiguities aren't split, i.e. M(AC)/M aren't treated as segregating sites, (M/A do).
  info <- x %>% ape::as.DNAbin() %>% ape::seg.sites(x=., strict=TRUE) %>%
    {list(name=unique(names(x)), locus_position=., not_varied_loci=.,
          nrow=as.integer(length(x) / 2), ncol=length(.))}
  info$not_varied_loci %<>% IRanges::IRanges(start=., end=.) %>%
    {list(IRanges::gaps(x=IRanges::reduce(x=.), start=1, end=IRanges::nchar(x[[1]])), .);} %>%
    mapply(FUN=function(xx, yy){setNames(object=xx, nm=rep(yy, times=length(xx)));},
           xx=., yy=c("not_varied", "varied"), SIMPLIFY=FALSE, USE.NAMES=FALSE) %>%
    {c(.[[1]], .[[2]])} %>% BiocGenerics::sort() %>% Biostrings::extractAt(x=x[[1]], at=.) %>%
    as.character(use.names=TRUE) %>%
    {base::replace(x=., list=(names(.) == "varied"), values="varied")} %>% unname()
  # split x into x1 and x2--------------------------------------------------------------------------
  ans <- base::split(x=1:length(x), f=names(x)) %>% base::simplify2array(x=.) %>% base::t(x=.) %>%
    {tibble(id=rownames(.), primary=.[,1], secondary=.[,2])} %>% {.[stringr::str_order(x=.$id, numeric=TRUE),]}
  x1 <- x[ans$primary]; x2 <- x[ans$secondary];
  remove(x, ans)
  any(names(x1) != info$name) %>% {if(.){stop("names' order bug");}}
  # leave snp loci only in x1 and x2----------------------------------------------------------------
  # info$locus_position %>% R.utils::seqToIntervals() %>% {IRanges::IRanges(start=.[,"from"], end=.[,"to"])}
  ans <- info$locus_position %>% {IRanges::IRanges(start=., end=.)} %>% IRanges::reduce()
  # x1 %>% as.DNAbin() %>% as.matrix() %>% .[,info$locus_position] %>% as.character() %>% base::toupper()
  x1 %<>% extractAt(x=., at=ans) %>% lapply(FUN=unlist) %>% DNAStringSet() %>% as.matrix()
  x2 %<>% extractAt(x=., at=ans) %>% lapply(FUN=unlist) %>% DNAStringSet() %>% as.matrix()
  remove(ans)
  # prepare a loci object (defined in pegas)--------------------------------------------------------
  if("" == "runs slower"){
    unphased_loci <- matrix(data="", nrow=info$nrow, ncol=info$ncol * 2,
                            dimnames=list(info$name, paste0("col", rep(info$locus_position, each=2))))
    ans <- (1:info$ncol) * 2; unphased_loci[,ans - 1] <- x1; unphased_loci[,ans] <- x2; remove(ans);
    unphased_loci %<>% pegas::alleles2loci(x=., ploidy=2)
  }
  if("" != "runs faster"){
    unphased_loci <- mapply(FUN=base::paste, purrr::array_tree(array=x1, margin=2),
                            purrr::array_tree(array=x2, margin=2), sep="/", SIMPLIFY=TRUE,
                            USE.NAMES=FALSE) %>% base::`rownames<-`(x=., value=info$name) %>%
      base::`colnames<-`(x=., value=paste0("col", info$locus_position)) %>%
      as.data.frame(stringsAsFactors=TRUE) %>% pegas::as.loci(x=., allele.sep="/")
  }
  remove(x1, x2)
  identical(x=rownames(unphased_loci), y=info$name) %>% {if(!.){stop("rownames are dropped");}}
  grepl(pattern="^col[[:digit:]]+$", x=names(unphased_loci)) %>% {if(!all(.)){stop("colnames are dropped");}}
  # return unphased_loci----------------------------------------------------------------------------
  unphased_loci <- list(unphased_loci=unphased_loci, locus_position=info$locus_position,
                        not_varied_loci=info$not_varied_loci)
  remove(info)
  return(unphased_loci)
}
