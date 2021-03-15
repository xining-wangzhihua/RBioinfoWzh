#a function similar to pegas::nuc.div(),
#pegas::nuc.div() will count the nucleotide difference between "a" & "-" as 0.
#nuc.div2() will treat the difference between "c" & "-" as normal, i.e. 1.

#get dependency
if(FALSE){
  require(magrittr)
  require(ape)
  require(pegas)
  source(file=R.home(component="RBioinfoZm/AnyEmpty().r"))
}

nuc.div2 <- function(x){
  #control the input of x-------------------------------------------------------
  if(is(object=x, class2="DNAbin")){
    if(AnyEmpty(x)){stop("DNAbin structure error");}
    if(is.list(x)) if(length(unique(lengths(x))) != 1){stop("DNAbin structure error");}
    x %<>% as.matrix() %>% pegas::haplotype()
  }
  if(!is(object=x, class2="haplotype")){
    stop("x must be a DNAbin object defined in ape, or haplotype object defined in pegas")
  }
  if(!is.matrix(x)){
    warning("please revise codes. when converting DNAbin to haplotype, it should be ",
            "automatically converted to matrix format.")
    x %<>% as.matrix()
  }
  wt <- attr(x=x, which="index") %>% lengths()
  x %<>% as.character()
  if(length(setdiff(x=unique(c(x)), y=c("a", "c", "g", "t", "-"))) != 0){
    stop("x contain invalid characters other than a, c, g, t, -.")
  }
  #main manipulation------------------------------------------------------------
  if(nrow(x) > 1){
    tulip <- utils::combn(x=1:nrow(x), m=2)
    ans_n_mismatch <- tulip %>% base::apply(MARGIN=2, FUN=function(datum, b){
      length(which(b[datum[1],] != b[datum[2],]))
    }, b=x)
    ans_wt <- tulip %>% base::apply(MARGIN=2, FUN=function(datum, b){
      b[datum[1]] * b[datum[2]]
    }, b=wt)
    tulip <- sum(ans_n_mismatch * ans_wt) / ncol(x) * 2 / sum(wt) / (sum(wt) - 1)#attention here!!!
    remove(ans_n_mismatch, ans_wt)
  }else{
    tulip <- 0
  }
  #clean environment and return results-----------------------------------------
  remove(x, wt)
  return(tulip)#remove(tulip)
}

#here are the test codes
if(FALSE){
  c("a", "c", "g", "t", "-") %>% rbind(.) %>% as.DNAbin() %>% nuc.div2()
  c("a", "c", "g", "t", "-") %>% rbind(., ., .) %>% as.DNAbin() %>% nuc.div2()
}
