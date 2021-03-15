# this function is just a wrapper of msa::msa() and ips::mafft()
# ips::mafft() don't implement the algorithm in r, it just call "mafft" command from shell

if("" == "get dependency"){
  require(magrittr)
  require(Biostrings); require(msa); require(ips); require(ape);
  source(file=R.home(component="RBioinfoZm/as.DNAStringSet.DNAbin().r"))
}

if("" == "test codes"){
  c(a="ACCGGT", b="ACCGT", a="ACGGT") %>% DNAStringSet() %>% widget.Align(xss=.)
}

widget.Align <- function(xss){
  #control the input of xss-----------------------------------------------------
  if(!methods::is(object=xss, class2="XStringSet")){stop("xss must be an XStringSet object");}
  if(length(xss) == 0){stop("length(xss) can't be 0");}
  if(any(nchar(xss) == 0)){stop("some sequence lengths in xss is 0");}
  #backup names(xss)------------------------------------------------------------
  original_name <- names(xss)
  hash_name <- (1:length(xss)) %>% as.hexmode() %>% as.character()
  names(xss) <- hash_name
  #align via ips::mafft or msa::msa()-------------------------------------------
  mafft_path <- Sys.which("mafft") %>% .[1] %>% normalizePath(path=., winslash="/", mustWork=FALSE)
  v_message <- paste0(" aligning via ", if(mafft_path == ""){"msa::msa()"}else{"ips::mafft()"}, "\n")
  cat("begin", v_message, sep="")
  if(mafft_path == ""){
    #msa::msa() is too slow, and is an expedient (ad hoc) method
    xss %<>% msa::msa(substitutionMatrix=nucleotideSubstitutionMatrix()) %>%
      methods::as(object=., Class="DNAStringSet")
  }else{
    #ips::mafft(options="--fmodel --inputorder --nuc") won't work
    xss %<>% as.DNAbin() %>% ips::mafft(x=., exec=paste0(mafft_path, " --fmodel --inputorder"),
                                        method="localpair", maxiterate=1000L, op=1.53, ep=0.123) %>% as.DNAStringSet.DNAbin()
  }
  cat("end", v_message, sep="")
  remove(mafft_path, v_message)
  #set names(xss) back----------------------------------------------------------
  xss %<>% .[hash_name] %>% setNames(nm=original_name)
  remove(original_name, hash_name)
  return(xss)
}
