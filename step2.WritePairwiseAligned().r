#this function is used to pairwise-align bi-directionally sequenced segments
#the input should be a directory (folder), and output is a kind of file (8 rows)
#the author (not strongly) recommend you to backup data before running (it's always a good practice to backup your data)

# similar to ?Biostrings::writePairwiseAlignments
warning("change the name WritePairwiseAligned() to writePairwiseAlignments2()")

#get dependency
if(FALSE){
  require(magrittr); require(tibble); require(dplyr); require(stringr);
  require(Biostrings); require(DECIPHER); require(seqinr);
}

#test codes
if(FALSE){
  x <- "./its/"
  forward_suffix_pattern <- ", its1a\\.seq$"
  backward_suffix_pattern <- ", its4\\.seq$"
  step2.WritePairwiseAligned(x=x, forward_suffix_pattern=forward_suffix_pattern,
                             backward_suffix_pattern=backward_suffix_pattern)
  remove(x, forward_suffix_pattern, backward_suffix_pattern)
}

step2.WritePairwiseAligned <- function(
  x="a_directory", forward_suffix_pattern=", forward_primer_name\\.seq$",
  backward_suffix_pattern=", backward_primer_name\\.seq$"
){
  fsp <- forward_suffix_pattern; bsp <- backward_suffix_pattern;
  remove(forward_suffix_pattern, backward_suffix_pattern)
  #control the input of x-------------------------------------------------------
  if(!is.character(x)){stop("x must be a string, specifing a directory");}#this line isn't necessary
  if(length(x) != 1){stop("length(x) must be 1");}
  if(!dir.exists(paths=x)){stop(x, " doesn't exist");}
  if(!file.info(x)$isdir){stop(x, " isn't a directory");}#this line isn't necessary
  if(grepl(pattern="/$", x=x)){x <- base::sub(pattern="/$", replacement="", x=x);}
  #control the input of forward_suffix_pattern and backward_suffix_pattern------
  if(!is.character(fsp)){stop("forward_suffix_pattern must be a character string");}
  if(!is.character(bsp)){stop("backward_suffix_pattern must be a character string");}
  if(length(fsp) != 1){stop("length(forward_suffix_pattern) must be 1");}
  if(length(bsp) != 1){stop("length(backward_suffix_pattern) must be 1");}
  if(is.na(fsp)){stop("forward_suffix_pattern can't be NA");}
  if(is.na(bsp)){stop("backward_suffix_pattern can't be NA");}
  #get .seq files in the directory----------------------------------------------
  x <- list(f=list.files(path=x, pattern=fsp, all.files=TRUE, full.names=TRUE, recursive=TRUE),
            b=list.files(path=x, pattern=bsp, all.files=TRUE, full.names=TRUE, recursive=TRUE))
  if(length(x$f) == 0){stop("forward_suffix_pattern error, can't find files");}
  if(length(x$b) == 0){stop("forward_suffix_pattern error, can't find files");}
  x$f <- tibble(id=base::sub(pattern=fsp, replacement="", x=x$f), f=x$f)
  x$b <- tibble(id=base::sub(pattern=bsp, replacement="", x=x$b), b=x$b)
  remove(fsp, bsp)
  if(!base::setequal(x=x$f$id, y=x$b$id)){
    c(base::setdiff(x=x$f$id, y=x$b$id), base::setdiff(x=x$b$id, y=x$f$id)) %>%
      paste0(collapse="\n") %>% stop("the following files don't have paired sequences:\n")
  }
  x <- dplyr::inner_join(x=x$f, y=x$b, by="id") %>% dplyr::arrange(id) %>%
    dplyr::mutate(aligned=paste0(id, ", pairwise aligned.txt"))
  #read in .seq files-----------------------------------------------------------
  y <- DNAStringSet(base::rep("", times=nrow(x))) %>% list(f="", b="", af=., ab=.)
  y$f <- x$f %>% lapply(FUN=base::readLines, warn=FALSE) %>%
    lapply(FUN=base::paste0, collapse="") %>% base::unlist() %>% DNAStringSet()
  if(any(BiocGenerics::width(x=y$f) == 0)){
    x$f[BiocGenerics::width(x=y$f) == 0] %>% paste0(collapse="\n") %>% stop("the following .seq files are empty:\n", .)
  }
  y$b <- x$b %>% lapply(FUN=base::readLines, warn=FALSE) %>%
    lapply(FUN=base::paste0, collapse="") %>% base::unlist() %>% DNAStringSet()
  if(any(BiocGenerics::width(x=y$b) == 0)){
    x$b[BiocGenerics::width(x=y$b) == 0] %>% paste0(collapse="\n") %>% stop("the following .seq files are empty:\n", .)
  }
  #reverse and complement the backward sequences--------------------------------
  y$b %<>% Biostrings::reverseComplement()
  #align forward and backward sequences-----------------------------------------
  subs_matr <- Biostrings::nucleotideSubstitutionMatrix()
  #subs_matr <- Biostrings::DNA_ALPHABET %>% base::setdiff(x=., y=c("-", "+", ".")) %>%
  #  {matrix(data=0, nrow=length(.), ncol=length(.), dimnames=list(., .))}
  #diag(subs_matr) <- 5
  for(i in 1:nrow(x)){
    #base::max(substr_matr) is 1
    ans <- pairwiseAlignment(pattern=y$f[[i]], subject=y$b[[i]], substitutionMatrix=subs_matr,
                             gapOpening=base::max(subs_matr) * 2, gapExtension=base::max(subs_matr) / 2)
    y$af[[i]] <- Biostrings::alignedPattern(ans)[[1]]
    y$ab[[i]] <- Biostrings::alignedSubject(ans)[[1]]
    remove(ans)
  }; remove(i);
  rm(subs_matr)
  if(!all(DECIPHER::RemoveGaps(y$af) == y$f)){stop("the pairwise alignments aren't global");}
  if(!all(DECIPHER::RemoveGaps(y$ab) == y$b)){stop("the pairwise alignments aren't global");}
  #get the consistency information of aligned sequences-------------------------
  f_consistency <- function(datum){
    #?Biostrings::compareStrings; ?stringi::stri_compare; ??SNP;
    datum <- seqinr::s2c(datum[["af"]]) != seqinr::s2c(datum[["ab"]])
    datum <- base::replace(x=base::rep("_", times=length(datum)), list=datum, values="*")
    datum <- paste0(datum, collapse="")
    return(datum)
  }
  x %<>% dplyr::select(aligned) %>% dplyr::rename(file=aligned) %>%
    dplyr::mutate(af=as.character(y$af), ab=as.character(y$ab), consistency="")
  if(is.list(x$af)){stop("behaviour of as.character() for XStringSet has changed");}
  x$consistency <- apply(X=x, MARGIN=1, FUN=f_consistency)
  remove(y, f_consistency)
  #output-----------------------------------------------------------------------
  tulip <- c("forward_primary_sequence", "forward_secondary_sequence",
             "forward_ratio", "consistency", "backward_primary_sequence",
             "backward_secondary_sequence", "backward_ratio", "f_or_b_is_used") %>%
    stringr::str_pad(string=., width=base::max(base::nchar(.)) + 2, side="right", pad=" ")
  for(i in 1:nrow(x)){
    nchar(x$af[i]) %>% {c(base::strrep(x=".", times=.), base::strrep(x="0", times=.))} %>%
      {c(x$af[i], ., x$consistency[i], x$ab[i], ., "")} %>%
      paste0(tulip, .) %>% cat(., file=x$file[i], append=FALSE, sep="\n")
  }; remove(i);
  remove(tulip)
  #return file names------------------------------------------------------------
  return(x$file)
}
