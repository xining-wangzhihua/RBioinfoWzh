#the result is always aligned

#get dependency
if(FALSE){
  require(magrittr); require(tibble); require(dplyr); require(stringr);
  require(Biostrings); require(ape);
  source(file=R.home(component="RBioinfoZm/as.DNAStringSet.DNAbin().r"))
  source(file=R.home(component="RBioinfoZm/widget.write.seg.sites().r")); require(R.utils);
  source(file=R.home(component="RBioinfoZm/widget.Align().r")); require(msa); require(ips);
}

#test codes
if(FALSE){
  x <- "./its/"
  suffix_pattern <- ", pairwise aligned\\.txt$"; ignore_ratio_larger_than=2L;
  step3.GetCheckedSequences(x=x, suffix_pattern=suffix_pattern,
                               ignore_ratio_larger_than=ignore_ratio_larger_than)
  remove(x, suffix_pattern, ignore_ratio_larger_than)
}

step3.GetCheckedSequences <- function(
  x="a_directory", suffix_pattern=", pairwise aligned\\.txt$", ignore_ratio_larger_than=2L
){
  irlt <- ignore_ratio_larger_than; remove(ignore_ratio_larger_than);
  #control the input of x-------------------------------------------------------
  if(!is.character(x)){stop("x must be a string, specifing a directory");}#this line isn't necessary
  if(length(x) != 1){stop("length(x) must be 1");}
  if(!dir.exists(paths=x)){stop(x, " doesn't exist");}
  if(!file.info(x)$isdir){stop(x, " isn't a directory");}#this line isn't necessary
  if(grepl(pattern="/$", x=x)){x <- base::sub(pattern="/$", replacement="", x=x);}
  #control the input of suffix_pattern------------------------------------------
  if(!is.character(suffix_pattern)){stop("suffix_pattern must be a character string");}
  if(length(suffix_pattern) != 1){stop("length(suffix_pattern) must be 1");}
  if(is.na(suffix_pattern)){stop("suffix_pattern can't be NA");}
  #control the input of ignore_ratio_larger_than--------------------------------
  if(!is.numeric(irlt)){stop("ignore_ratio_larger_than must be a number");}
  if(length(irlt) != 1){stop("length(ignore_ratio_larger_than) must be 1");}
  irlt <- as.integer(irlt)
  if(is.na(irlt)){stop("ignore_ratio_larger_than must be an integer");}
  if(!(irlt %in% 0L:9L)){stop("ignore_ratio_larger_than must be one of 0 ~ 9");}
  #get file list and id, and read in files--------------------------------------
  file_name <- list.files(path=x, pattern=suffix_pattern, all.files=TRUE, full.names=TRUE, recursive=TRUE)
  if(length(file_name) == 0){stop("suffix_pattern error, can't find files");}
  file_name <- tibble(id=base::sub(pattern=suffix_pattern, replacement="", x=file_name), file=file_name) %>%
    dplyr::mutate(id=base::basename(id)) %>% {.[str_order(x=.$id, numeric=TRUE),]}
  remove(suffix_pattern)
  if(anyDuplicated(file_name$id) != 0){stop("there can't be duplicated sequence IDs, based on the first part of file names");}
  x <- lapply(X=file_name$file, FUN=base::readLines) %>% lapply(FUN=function(va){va[va != ""]})
  #if file contents are valid (each file contain 8 lines)----------------
  if(any(lengths(x) != 8)){
    file_name$file[lengths(x) != 8] %>% paste0(collapse="\n") %>%
      stop("the following files don't contain 8 lines:\n", .)
  }
  #if file contents are valid (line width are same, row names are valid)--------
  x_nchar <- lapply(X=x, FUN=base::nchar) %>% lapply(FUN=base::unique)
  if(any(lengths(x_nchar) != 1)){
    file_name$file[lengths(x_nchar) != 1] %>% paste0(collapse="\n") %>%
      stop("line widths within the following files aren't all equal:\n", .)
  }
  x_nchar %<>% base::unlist()
  ans <- c("forward_primary_sequence", "forward_secondary_sequence", "forward_ratio",
           "consistency", "backward_primary_sequence", "backward_secondary_sequence",
           "backward_ratio", "f_or_b_is_used") %>%
    {stringr::str_pad(string=., width=max(nchar(.)) + 2, side="right", pad=" ")}
  ans_nchar <- nchar(ans[1])
  if(any(x_nchar <= ans_nchar)){
    file_name$file[x_nchar <= ans_nchar] %>% paste0(collapse="\n") %>%
      stop("lines in the following files are too short:\n", .)
  }
  tempo <- lapply(X=x, FUN=base::substr, start=1, stop=ans_nchar) %>%
    lapply(FUN=base::identical, y=ans) %>% base::unlist() %>% {!.}
  if(any(tempo)){
    file_name$file[tempo] %>% paste0(collapse="\n") %>%
      stop("the leading parts of each line are invalid in the following files:\n", .)
  }
  remove(tempo)
  x %<>% lapply(FUN=function(va, vb){base::substr(x=va, start=vb + 1, stop=nchar(va))}, vb=ans_nchar)
  remove(x_nchar, ans, ans_nchar)
  #convert x[[i]] from character strings to tibble------------------------------
  x %<>% lapply(FUN=base::strsplit, split="") %>% lapply(FUN=as_tibble, .name_repair="minimal") %>%
    lapply(FUN=stats::setNames, nm=c("fp", "fs", "fr", "consistency", "bp", "bs", "br", "f_or_b"))
  #if file contents are valid (no weird characters)-----------------------------
  ans <- list(fp=setdiff(x=DNA_ALPHABET, y=c("+", ".")), fs=setdiff(x=DNA_ALPHABET, y="+"),
              fr=as.character(0:9), consistency=c("_", "*"), bp="", bs="", br="",
              f_or_b=c("f", "b")); ans$bp <- ans$fp; ans$bs <- ans$fs; ans$br <- ans$fr;
  ans <- lapply(X=x, FUN=function(datum, la){
    result <- rep(TRUE, times=length(datum))
    for(i in 1:length(datum)){result[i] <- base::all(datum[[i]] %in% la[[i]])}; remove(i);
    return(!result)
  }, la=ans) %>% lapply(FUN=base::any) %>% base::unlist()
  if(any(ans)){
    file_name$file[ans] %>% paste0(collapse="\n") %>%
      stop("the following files contain weird characters:\n", .)
  }
  remove(ans)
  #if file contents are valid (secondary_sequence is consistent with ratio)-----
  ans <- lapply(X=x, FUN=function(datum){
    datum <- list(s=c(datum$fs, datum$bs), r=c(datum$fr, datum$br))
    datum <- base::any((datum$s == ".") != (datum$r == "0"))
    return(datum)
  }) %>% base::unlist()
  if(any(ans)){
    file_name$file[ans] %>% paste0(collapse="\n") %>%
      stop("f/r_secondary_sequence aren't consistent with f/r_ratio in the following files:\n", .)
  }
  remove(ans)
  #if file contents are valid (primary_sequence != secondary_sequence)----------
  ans <- lapply(X=x, FUN=function(dfa){any(dfa$fp == dfa$fs) | any(dfa$bp == dfa$bs)}) %>% base::unlist()
  if(any(ans)){
    file_name$file[ans] %>% paste0(collapse="\n") %>%
      stop("some secondary_sequence are same with primary_sequence in the following files:\n", .)
  }
  remove(ans)
  #if file contents are valid (not checked)-------------------------------------
  if(FALSE){"not checked: primary_sequence is ambiguity, and secondary_sequence is \".\"";}
  if(FALSE){
    "not checked: when primary_sequence and secondary_sequence are both ambiguity, ";
    "they can't overlap each other, i.e. M(AC)/R(AG) sholdn't allowed";
  }
  #get checked sequences (remove gaps -/./0) according to f_or_b_is_used--------
  fa <- function(ma){
    apply(X=ma, MARGIN=1, FUN=function(va){if(va[8] == "f"){va[1:3]}else{va[5:7]};})
  }
  fb <- function(dfa, max_ratio){
    should_ignore <- (dfa$r > max_ratio); dfa$s[should_ignore] <- "."; dfa$r[should_ignore] <- 0L; remove(should_ignore);
    is_gap <- (dfa$p == "-") & (dfa$s == "."); dfa <- dfa[!is_gap,]; remove(is_gap);#nrow(dfa) can't be 0 after removing empty nucleotides
    return(dfa)
  }
  x <- lapply(X=x, FUN=base::as.matrix) %>% lapply(FUN=fa) %>% lapply(FUN=base::t) %>%
    lapply(FUN=as_tibble, .name_repair="minimal") %>%
    lapply(FUN=setNames, nm=c("p", "s", "r")) %>%
    lapply(FUN=dplyr::mutate, r=as.integer(r)) %>%
    lapply(FUN=fb, max_ratio=irlt)
  remove(irlt, fa, fb)
  ans <- lapply(X=x, FUN=base::nrow) %>% base::unlist() %>% {. == 0}; if(any(ans)){
    file_name$file[ans] %>% paste0(collapse="\n") %>%
      stop("after removing gap residues, the following files contain no residues:\n", .)
  }; remove(ans);
  #get primary, secondary and merged sequences----------------------------------
  f_m <- function(dfa){
    #remove -/[ACGT]/[1-9], -/./0 is already removed
    dfa %<>% {.[.$p != "-",];}; if(nrow(dfa) == 0){stop("after removing gaps, a merged sequence is empty");}
    dfa$s[dfa$s == "-"] <- ""
    dfa$s[dfa$s == "."] <- ""
    dfa <- paste0(dfa$p, dfa$s) %>% Biostrings::mergeIUPACLetters() %>% paste0(collapse="")
    return(dfa)
  }
  f_ps <- function(dfa){
    is_same <- (dfa$s == "."); dfa$s[is_same] <- dfa$p[is_same]; remove(is_same);
    dfa <- c(primary=paste0(dfa$p, collapse=""), secondary=paste0(dfa$s, collapse=""))
    return(dfa)
  }
  x_m <- lapply(X=x, FUN=f_m) %>% setNames(nm=file_name$id) %>% base::unlist() %>% DNAStringSet()
  x_ps <- lapply(X=x, FUN=f_ps) %>% setNames(nm=file_name$id) %>% base::unlist() %>% DNAStringSet()
  remove(x, f_m, f_ps)
  if(anyDuplicated(names(x_ps)) != 0){warnings("there are duplicated IDs after appending .primary/.secondary");}
  #align sequences--------------------------------------------------------------
  x_m %<>% widget.Align()
  x_ps %<>% widget.Align()
  #output files-----------------------------------------------------------------
  result <- base::strftime(x=Sys.time(), format="./step3.GetCheckedSequences() %Y%m%d-%H%M%S")
  result <- c("merged sequences.fasta", "primary and secondary sequences.fasta",
              "merged sequences, segregating sites.txt",
              "primary and secondary sequences, segregating sites.txt",
              "processed files.log") %>% paste0(result, ", ", .)
  ape::as.DNAbin(x_m) %>% ape::write.FASTA(x=., file=result[1])
  ape::as.DNAbin(x_ps) %>% ape::write.FASTA(x=., file=result[2])
  ans <- widget.write.seg.sites(xss=x_m); if(!file.rename(from=ans, to=result[3])){
    warning("file rename error"); result[3] <- ans;
  }; remove(ans);
  ans <- widget.write.seg.sites(xss=x_ps); if(!file.rename(from=ans, to=result[4])){
    warning("file rename error"); result[4] <- ans;
  }; remove(ans);
  cat(file_name$file, sep="\n", file=result[5])
  #clean environment and return result------------------------------------------
  message("\"A-\" \"CG-\" etc, (if there is any), can't be represented as ",
          "ambiguity codes in merged sequences, which may generate bugs")
  remove(x_m, x_ps, file_name)
  return(result)
}
