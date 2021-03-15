#a function to get segregating sites and their flanks
#to facilitate (assist) checking an alignment result

if("" != "get dependency"){
  require(magrittr); require(tibble); require(stringr);
  require(Biostrings); require(ape);
  require(R.utils); require(IRanges);
}

widget.write.seg.sites <- function(xss){
  #control the input of xss-----------------------------------------------------
  if(!methods::is(object=xss, class2="XStringSet")){stop("xss must be an XStringSet");}
  xss %<>% DNAMultipleAlignment(x=.)
  if(any(dim(xss) == 0)){stop("xss can't be empty");}
  if(is.null(rownames(xss))){rownames(xss) <- paste0("NO.", 1:nrow(xss));}
  #get locations of snp and their flanks----------------------------------------
  tulip <- xss %>% as.DNAbin() %>% ape::seg.sites(strict=TRUE) %>% R.utils::seqToIntervals() %>%
    {IRanges::IRanges(start=.[,"from"], end=.[,"to"])} %>% list(before=0L, snp=., after=0L)
  tulip$before <- flank(x=tulip$snp, width=1L, start=TRUE, both=FALSE)
  tulip$after <- flank(x=tulip$snp, width=1L, start=FALSE, both=FALSE)
  tulip %<>% {c(.$before, .$snp, .$after)} %>%
    IRanges::reduce(drop.empty.ranges=TRUE, min.gapwidth=2L) %>%
    IRanges::restrict(start=1, end=ncol(xss))
  #get a matrix-----------------------------------------------------------------
  xss %<>% methods::as(Class="DNAStringSet") %>% Biostrings::extractAt(at=tulip) %>%
    lapply(FUN=as.character) %>% as_tibble(.name_repair="minimal") %>% as.matrix() %>%
    base::t() %>% {cbind(rownames(.), .)} %>% rbind(c("#id", as.character(tulip)), .)
  rownames(xss) <- NULL
  remove(tulip)
  xss %<>% apply(MARGIN=2, FUN=function(datum){
    str_pad(string=datum, width=max(nchar(datum)) + 2, side="both", pad=" ")
  })
  #write result-----------------------------------------------------------------
  result <- base::strftime(x=Sys.time(), format="./show.seg.sites() %Y%m%d-%H%M%S.txt")
  write.table(x=xss, file=result, append=FALSE, quote=FALSE, sep="",
              row.names=FALSE, col.names=FALSE)
  return(result)#remove(xss, result)
}
