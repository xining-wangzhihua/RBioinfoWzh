#get dependency
if(FALSE){
  require(magrittr); require(tibble); require(dplyr);
  require(Biostrings)
}

#test codes
if(FALSE){
  xss <- "./its/result of zm.step3.GetCheckedSequences()/primary and secondary sequences.fasta" %>%
    readDNAMultipleAlignment(filepath=., format="fasta")
  primers <- "../step 1, look for samples, extract DNA, PCR/primer sequences.fasta" %>%
    readDNAStringSet(filepath=.); names(primers) %<>% base::sub(pattern="^ ", replacement="", x=.);
  widget.RepairPrimerPair(xss=xss, forward_primer=primers$ITS1a, backward_primer=primers$ITS4)
  remove(xss, primers)
}

widget.RepairPrimerPair <- function(xss, forward_primer, backward_primer, mismatch_frequency=0.3){
  fp <- forward_primer; bp <- backward_primer; mf <- mismatch_frequency;
  remove(forward_primer, backward_primer, mismatch_frequency)
  #control the input of xss-----------------------------------------------------
  if(methods::is(object=xss, class2="DNAStringSet")){xss %<>% DNAMultipleAlignment();}
  if(!methods::is(object=xss, class2="DNAMultipleAlignment")){
    stop("x must be a DNAStringSet or DNAMultipleAlignment object")
  }
  if(any(dim(xss) == 0)){stop("xss can't be empty");}
  #control the input of fp, bp--------------------------------------------------
  if(methods::is(object=fp, class2="DNAString")){fp %<>% as.character();}
  if(methods::is(object=bp, class2="DNAString")){bp %<>% as.character();}
  if(!is.character(fp)){stop("fp must be a character string or DNAString");}
  if(!is.character(bp)){stop("bp must be a character string or DNAString");}
  if(length(fp) != 1){stop("length(fp) must be 1");}
  if(length(bp) != 1){stop("length(bp) must be 1");}
  if(nchar(fp) == 0){stop("nchar(fp) can't be 0");}
  if(nchar(bp) == 0){stop("nchar(bp) can't be 0");}
  #control the input of mismatch_frequency--------------------------------------
  if(!is.numeric(mf)){stop("mismatch_frequency must be a number");}
  if(mf <= 0){stop("mismatch_frequency must be larger than 0");}
  #control the relationship among xss, forward_primer, backward_primer----------
  if(ncol(xss) <= nchar(fp) + nchar(bp)){
    stop("sequence length must be longer than the sum of primer pair")
  }
  #reverse and complement backward_primer---------------------------------------
  bp %<>% DNAString() %>% reverseComplement() %>% as.character()
  #get consensus string of xss--------------------------------------------------
  iupac <- IUPAC_CODE_MAP %>% enframe() %>% dplyr::mutate(value=paste0(value, "-")) %>%
    tibble::add_row(name="-", value="-") %>% deframe() %>% c(IUPAC_CODE_MAP, .)
  xss_consensus <- consensusString(x=xss, ambiguityMap=iupac, threshold=0.2) %>% DNAString()
  remove(iupac)
  #get forward/backward range---------------------------------------------------
  tulip <- list(fr="", br="")
  #
  tulip$fr <- round(nchar(fp) * mf) %>% base::max(3L, .) %>% base::min(nchar(fp), .) %>%
    matchPattern(pattern=fp, subject=xss_consensus, max.mismatch=., with.indels=TRUE) %>%
    IRanges::ranges() %>% {.[order(IRanges::start(.), decreasing=FALSE)]}
  if(length(tulip$fr) == 0){stop("can't find forward pattern, try a bigger mismatch_frequency");}
  if(length(tulip$fr) > 1) if(tulip$fr[1] %over% tulip$fr[-1]){
    #overlapsAny(query=ans[1], subject=ans[-1])
    stop("too many, and overlaped, ranges found, try a smaller mismatch_frequency")
  }
  tulip$fr %<>% .[1]; IRanges::start(tulip$fr) <- 1;
  tulip$fr %>% as.character() %>% message("the forward range is: ", .)
  #
  tulip$br <- round(nchar(bp) * mf) %>% base::max(3L, .) %>% base::min(nchar(bp), .) %>%
    matchPattern(pattern=bp, subject=xss_consensus, max.mismatch=., with.indels=TRUE) %>%
    IRanges::ranges() %>% {.[order(IRanges::end(.), decreasing=TRUE)]}
  if(length(tulip$br) == 0){stop("can't find backward pattern, try a bigger mismatch_frequency");}
  if(length(tulip$br) > 1) if(tulip$br[1] %over% tulip$br[-1]){
    stop("too many, and overlaped, ranges found, try a smaller mismatch_frequency")
  }
  tulip$br %<>% .[1]; IRanges::end(tulip$br) <- nchar(xss_consensus);
  tulip$br %>% as.character() %>% message("the backward range is: ", .)
  #
  remove(mf, xss_consensus)
  xss %<>% Biostrings::unmasked()
  #ensure range flanks are consistent-------------------------------------------
  tulip$fr %>% IRanges::flank(width=5, start=FALSE) %>% extractAt(x=xss, at=.) %>%
    unlist() %>% unname() %>% as.character() %>% unique() %>%
    {if(length(.) != 1){stop("sequences at the start flanks aren't consistent");}}
  tulip$br %>% IRanges::flank(width=5, start=TRUE) %>% extractAt(x=xss, at=.) %>%
    unlist() %>% unname() %>% as.character() %>% unique() %>%
    {if(length(.) != 1){stop("sequences at the end flanks aren't consistent");}}
  #replace forward/backward part with fr/br-------------------------------------
  #attention: replace bp has to be ahead of fp, for replacing fp may change sequence length
  xss %<>% replaceAt(at=tulip$br, value=bp) %>% replaceAt(at=tulip$fr, value=fp)
  remove(fp, bp, tulip)
  #return result----------------------------------------------------------------
  return(xss)
}
