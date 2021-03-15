#! /usr/bin/Rscript
# license: BSD 3-Clause License
# author: https://github.com/ywd5
# version: 20210314

warning("ensure \"if there are only 2 sequences, a.primary and a.secondary\" won't generate bugs")

# loci object can store DNA sequences, microsatellite, AFLP, RFLP, etc data
# so it's prefered to DNAStringSet object
# snpStats::loci2SnpMatrix(x)
# require(PopGenReport)
# require(genepop)

if("" == "get dependencies"){
  require(magrittr); require(tibble); require(dplyr); require(stringr); require(purrr);
  require(pegas);
}

if("" == "test codes"){
  WritePhaseInput(x=unphased_loci$unphased_loci,
                  locus_position=unphased_loci$locus_position)
}

WritePhaseInput <- function(x, locus_position=integer(), known=character(),
                            is_microsatellite=base::rep(FALSE, times=ncol(x))){
  # if a locus is ssr, but you want to by pass the "stepwise mutation model", specify it as not ssr
  # the arguement "known" represents "segments whose phases are known"
  # define internal functions-----------------------------------------------------------------------
  f_sub <- function(xx, pending=TRUE){
    # PHASE support a locu being "?" (for non-ssr) or "-1" (for ssr),
    # where PHASE will guess its phase according to data elsewhere
    # test data: xx <- c("A", "C", "G")
    # test data: xx <- c("A", "C")
    # test data: xx <- c("99")
    # test data: xx <- c("0", "99") # this also need to be modified
    # test data: xx <- c("0", "99", "9") # this also need to be modified
    # test data: xx <- c("?", "A", "C")
    # test data: xx <- c("?", "A", "C", "G")
    # test data: xx <- c("-1", "0", "99")
    # test data: xx <- c("-1", "0", "99", "9")
    # test data: xx <- c("0", "99", "A")
    # test data: xx <- c("?", "-1", "0")
    # test data: xx <- c("?", "-1", "0", "99")
    # test data: xx <- c("?", "-1", "A")
    # test data: xx <- c("?", "-1", "0", "A")
    # test data: xx <- c("?", "0", "A")
    # test data: xx <- c("?", "0", "99", "A")
    if((!is.character(xx)) | (length(xx) == 0) | (anyDuplicated(xx) != 0)){stop("invalid input");}
    if(("-1" %in% xx) & ("?" %in% xx)){
      warning("a locus contains \"?\" and \"-1\", which is frustrating for PHASE")
    }
    if(length(pending) != 1){stop("invalid input");}
    #
    all_integer <- base::suppressWarnings(expr=as.integer(xx)) %>% base::anyNA(x=.) %>% {!.}
    result <- list(original=character(), integer=character())
    if(all_integer) if(length(xx) == 1L){pending <- FALSE;}
    if(!all_integer) if(length(xx) < 3L){pending <- FALSE;}
    if(!all_integer) if((length(xx) == 3L) & ("?" %in% xx)){pending <- FALSE;}
    if(pending) if(all_integer){
      if("-1" %in% xx){result$original %<>% c(., "-1"); result$integer %<>% c(., "-1"); xx %<>% setdiff(x=., y="-1");}
      result$original %<>% c(., xx); result$integer %<>% {c(., as.character(0:(length(xx) - 1)))};
    }
    if(pending) if(!all_integer){
      if("?" %in% xx){result$original %<>% c(., "?"); result$integer %<>% c(., "-1"); xx %<>% setdiff(x=., y="?");}
      result$original %<>% c(., xx); result$integer %<>% {c(., as.character(0:(length(xx) - 1)))};
    }
    remove(pending, all_integer, xx)
    #
    if(length(result$original) > 79L){
      warning("usually PHASE (version 2.1.1) doesn't support more than 79 different ",
              "types of alleles at a locus. the user may get a weird error message from ",
              "PHASE. \"something like K=..., KMAX=50, please recompile PHASE\"")
    }
    return(result)
  }
  # control the input of x--------------------------------------------------------------------------
  # check structure
  if(!is(object=x, class2="loci")){stop("x must be a loci object defined in package pegas");}
  if(any(dim(x) == 0)){stop("x is empty");}
  attr(x=x, which="locicol", exact=TRUE) %>% {length(.) == 0} %>%
    {if(.){stop("x doesn't contain any loci data");}}
  if(nrow(x) == 1){stop("only 1 segment. no need to analyse");}
  unique(x) %>% {nrow(.) == 1} %>% {if(.){stop("all segments are same. no need to analyse");}}
  # check names
  rownames(x) %>% is.null() %>% {if(.){stop("rownames(x) must be specified");}}
  rownames(x) %>% anyNA() %>% {if(.){stop("rownames(x) must be specified");}}
  rownames(x) %<>% stringr::str_squish(string=.)
  rownames(x) %>% {anyDuplicated(.) != 0} %>% {if(.){stop("rownames(x) can't contain duplicated items");}}
  rownames(x) %>% {any(. == "")} %>% {if(.){stop("rownames(x) must be specified");}}
  # check each item (except meta data) is factor, and then ploidies are 2
  attr(x=x, which="locicol", exact=TRUE) %>% x[.] %>% lapply(FUN=is.factor) %>% unname(obj=.) %>%
    unlist() %>% {if(!all(.)){stop("generally, all variables in x should be factors");}}
  getPloidy(x=x) %>% {any(. != 2L)} %>% {if(.){stop("isn't diploid, not supported");}}
  # remove irrelevant columns in x------------------------------------------------------------------
  x %<>% {.[attr(x=., which="locicol", exact=TRUE)]}
  # control the input of locus_position-------------------------------------------------------------
  if(!is.numeric(locus_position)){stop("locus_position must be integers");}
  locus_position %<>% as.integer()
  if(length(locus_position) != 0){
    if(length(locus_position) != ncol(x)){stop("length(locus_position) must be 0 or equal ncol(x)");}
    if(anyNA(locus_position)){stop("locus_position must be integers");}
    # if(any(locus_position < 1)){stop("locus_position can't be smaller than 1");}
  }
  # control the input of known----------------------------------------------------------------------
  known %>% is.character() %>% {if(!.){stop("known must be a character vector");}}
  if(length(known) == 1) if(file.exists(known)){
    message("Assuming \"known\" specifies a file name, where there is information.")
    warning("consider the confict to specifing \"known\" file input directly")
    known %<>% readLines(con=., encoding="UTF-8") %>% stringr::str_trim(string=., side="both") %>%
      grep(pattern="^#", x=., value=TRUE, invert=TRUE) %>% {.[. != ""]}
  }
  if(length(known) != 0){
    known %>% anyNA() %>% {if(.){stop("known can't contain NA");}}
    known %<>% stringr::str_trim(string=., side="both") %>% {.[. != ""]}
    known %>% {anyDuplicated(.) != 0} %>% {if(.){stop("known can't contain duplicated items");}}
    known %>% {any(. == "")} %>% {if(.){stop("known can't contain empty strings");}}
  }
  # control the relationship between known and x----------------------------------------------------
  if(length(known) > 0) if(!all( known %in% rownames(x) )){
    stop("known don't match with rownames(x)")
  }
  # control the input of is_microsatellite----------------------------------------------------------
  if(!is.logical(is_microsatellite)){stop("is_microsatellite must be a logical vector");}
  if(length(is_microsatellite) != ncol(x)){stop("length(is_microsatellite) must equal ncol(x)");}
  if(anyNA(is_microsatellite)){stop("is_microsatellite can't contain NA");}
  # get info----------------------------------------------------------------------------------------
  info <- rep(FALSE, times=ncol(x)) %>%
    list(name=rownames(x), known=known, locus_name=colnames(x), locus_position=locus_position,
         nrow=nrow(x), ncol=ncol(x), true_ssr=is_microsatellite, pseudo_ssr=., sub="", output="")
  remove(locus_position, is_microsatellite, known)
  # create a directory to store the input of phase--------------------------------------------------
  info$output <- strftime(x=Sys.time(),format="./r_phase, %Y%m%d-%H%M%S")
  if(dir.exists(info$output)){
    stop("WritePhaseInput() wants to create the directory \"", info$output,
         "\", but there is already a directory with that name")
  }
  dir.create(info$output)
  if(!dir.exists(info$output)){stop("directory creation error");}
  # split x into 2 matrices-------------------------------------------------------------------------
  ans <- lapply(X=x, FUN=as.character) %>% lapply(FUN=base::strsplit, split="/", fixed=TRUE) %>%
    lapply(FUN=simplify2array) %>% simplify2array()
  x1 <- ans[1,,]
  x2 <- ans[2,,]
  remove(x, ans)
  # get info$sub and info$pseudo_ssr----------------------------------------------------------------
  # substitute "0/1/A/C/G/T/hap1/hap2" with "1/2/3" at multi-allele loci is necessary
  # because phase only accept 2-allele snp loci and microsatellite loci, the user must make
  # their multi-allele loci looks like microsatellite loci, i.e. represent them with by integers
  #
  # if don't consider "stepwise mutation model", substitute "100/200/999/etc" ssr data with "1/2/3" is also necessary
  # because a bug in PHASE, it won't work if "max(ssr) - min(ssr) > 79"
  info$sub <- mapply(FUN=base::c, purrr::array_tree(array=x1, margin=2), purrr::array_tree(array=x2, margin=2),
                     SIMPLIFY=FALSE, USE.NAMES=FALSE) %>% lapply(FUN=unique) %>% lapply(FUN=sort) %>%
    mapply(FUN=f_sub, xx=., pending=!info$true_ssr, SIMPLIFY=TRUE, USE.NAMES=FALSE) %>% base::t(x=.) %>%
    as_tibble() %>% dplyr::mutate(col_index=1:info$ncol, col_name=info$locus_name, .before=1L)
  remove(f_sub)
  info$pseudo_ssr <- (base::lengths(info$sub$original) != 0)
  if(any(info$true_ssr & info$pseudo_ssr)){stop("bug: calculation of info$sub");}
  # substitute info$sub$original with info$sub$integer----------------------------------------------
  i_set <- which(info$pseudo_ssr); if(length(i_set) > 0L) for(i in i_set){
    ans <- setNames(object=info$sub$integer[[i]], nm=info$sub$original[[i]])
    x1[,i] %<>% ans[.] %>% unname()
    x2[,i] %<>% ans[.] %>% unname()
    remove(ans)
  }; i <- 0L; remove(i_set, i);
  # create substitute_with_integers.tsv-------------------------------------------------------------
  info$sub %>% dplyr::mutate(original=lapply(X=original, FUN=paste0, collapse=","),
                             integer=lapply(X=integer, FUN=paste0, collapse=",")) %>%
    dplyr::mutate(original=unlist(original), integer=unlist(integer)) %>%
    utils::write.table(x=., file=paste0(info$output, "/substitute_with_integers.tsv"), append=FALSE,
                       quote=TRUE, sep="\t", row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")
  # create input_d.txt------------------------------------------------------------------------------
  # ifelse(test=info$true_ssr, yes="0", no="1")
  base::rep("1", times=info$ncol) %>% base::replace(x=., list=info$true_ssr, values="0") %>%
    base::paste(., collapse=" ") %>% base::paste0(., "\n") %>%
    base::cat(., file=paste0(info$output, "/input_d.txt"), append=FALSE, sep="")
  # create input_k.txt------------------------------------------------------------------------------
  base::strrep("*", times=info$ncol) %>% base::rep(., times=info$nrow) %>%
    base::replace(x=., list=(info$name %in% info$known), values=base::strrep("0", times=info$ncol)) %>%
    cat(., file=paste0(info$output, "/input_k.txt"), append=FALSE, sep="\n")
  # create input.txt--------------------------------------------------------------------------------
  tulip <- list(individual_number=as.character(info$nrow), locus_number=as.character(info$ncol),
                position="", snp_or_microsatellite="")
  if(length(info$locus_position) != 0){
    tulip$position <- paste0(info$locus_position, collapse=" ") %>% paste0("P ", .)
  }
  tulip$snp_or_microsatellite <- rep("S", times=info$ncol) %>%
    base::replace(x=., list=(info$true_ssr | info$pseudo_ssr), values="M") %>% paste0(collapse=" ")
  cat("\n", tulip$individual_number, "\n", tulip$locus_number, "\n\n",
      tulip$position, "\n", tulip$snp_or_microsatellite, "\n",
      file=paste0(info$output, "/input.txt"), append=FALSE, sep="")
  remove(tulip)
  mapply(FUN=base::list, info$name, purrr::array_tree(array=x1, margin=1),
         purrr::array_tree(array=x2, margin=1), SIMPLIFY=FALSE, USE.NAMES=FALSE) %>%
    lapply(FUN=function(xx){lapply(X=xx, FUN=paste0, collapse=" ");}) %>%
    lapply(FUN=function(xx){paste0("\n", xx[1], "\n", xx[2], "\n", xx[3], "\n")}) %>%
    unlist() %>% cat(., file=paste0(info$output, "/input.txt"), append=TRUE, sep="")
  # return result-----------------------------------------------------------------------------------
  remove(x1, x2)
  message("\ninput files for PHASE are written to: \"", info$output, "\"\n")
  return(info$output)
}
