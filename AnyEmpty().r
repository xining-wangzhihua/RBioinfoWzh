#get dependency
if(FALSE){
  require(magrittr); require(stringr); require(tibble);
  require(Biostrings); require(ape);
}

AnyEmpty <- function(x){
  ae <- function(datum){base::UseMethod(generic="ae");}
  ae.NULL <- function(datum){return(TRUE);}
  ae.logical <- function(datum){
    result <- FALSE
    if(length(datum) == 0){result <- TRUE;}
    if(!result) if(anyNA(datum)){result <- TRUE;}
    return(result)
  }
  ae.integer <- ae.logical
  ae.numeric <- ae.logical
  ae.complex <- ae.logical
  ae.character <- function(datum){
    result <- FALSE
    if(length(datum) == 0){result <- TRUE;}
    if(!result) if(anyNA(datum)){result <- TRUE;}
    if(!result){datum %<>% stringr::str_squish();}
    if(!result) if(any(datum == "")){result <- TRUE;}
    return(result)
  }
  ae.list <- function(datum){
    return(length(datum) == 0)
  }
  ae.array <- function(datum){
    result <- FALSE
    if(any(dim(datum) == 0)){result <- TRUE;}
    if(!result) if(anyNA(datum)){result <- TRUE;}
    return(result)
  }
  ae.matrix <- ae.array
  ae.data.frame <- ae.array
  ae.XString <- function(datum){
    return(nchar(datum) == 0)
  }
  ae.XStringSet <- function(datum){
    result <- FALSE
    if(length(datum) == 0){result <- TRUE;}
    if(!result) if(any(nchar(datum) == 0)){result <- TRUE;}
    return(result)
  }
  ae.DNAbin <- function(datum){
    result <- FALSE
    if(is.matrix(datum)){
      #why not as.list(datum) and merge this condition with if(is.list(datum)){}:
      #matrix(data="a", nrow=0, ncol=0) %>% as.DNAbin() %>% as.list() is an error
      if(any(dim(datum) == 0)){result <- TRUE;}
      if(!result) if(anyNA(as.character(datum))){result <- TRUE;}
    }else{
      datum %<>% as.list()
      if(length(datum) == 0){stop("bug. this condition shouldn't happen");}
      if(any(lengths(datum) == 0)){result <- TRUE;}
      if(!result) if(anyNA(unlist(as.character(datum)))){result <- TRUE;}
    }
    return(result)
  }
  ae.AAbin <- ae.DNAbin
  x <- ae(datum=x)
  return(x)
}

#here are the test codes
if(FALSE){
  NULL %>% AnyEmpty()
  vector() %>% AnyEmpty()
  logical() %>% AnyEmpty()
  TRUE %>% AnyEmpty()
  integer() %>% AnyEmpty()
  NA_integer_ %>% AnyEmpty()
  0L %>% AnyEmpty()
  numeric() %>% AnyEmpty()
  NA_real_ %>% AnyEmpty()
  1.1 %>% AnyEmpty()
  complex() %>% AnyEmpty()
  NA_complex_ %>% AnyEmpty()
  (1 + 2i) %>% AnyEmpty()
  character() %>% AnyEmpty()
  NA_character_ %>% AnyEmpty()
  c("a", " ") %>% AnyEmpty()
  c("a", NA_character_) %>% AnyEmpty()
  c("a", "b") %>% AnyEmpty()
  matrix(data="", nrow=0, ncol=0) %>% AnyEmpty()
  matrix(data="", nrow=2, ncol=0) %>% AnyEmpty()
  matrix(data="", nrow=0, ncol=2) %>% AnyEmpty()
  matrix(data=c("", NA_character_), nrow=1, ncol=2) %>% AnyEmpty()
  matrix(data="", nrow=2, ncol=3) %>% AnyEmpty()
  array(data=logical(), dim=c(0, 0, 0)) %>% AnyEmpty()
  array() %>% AnyEmpty()
  array(data=1:24, dim=c(4, 3, 2)) %>% AnyEmpty()
  data.frame() %>% AnyEmpty()
  matrix(data="", nrow=2, ncol=0) %>% as.data.frame() %>% AnyEmpty()
  matrix(data="", nrow=0, ncol=2) %>% as.data.frame() %>% AnyEmpty()
  matrix(data=c("", NA_character_), nrow=1, ncol=2) %>% as.data.frame() %>% AnyEmpty()
  data.frame(x=1:6) %>% AnyEmpty()
  tibble() %>% AnyEmpty()
  matrix(data="", nrow=2, ncol=0) %>% as_tibble(.name_repair="minimal") %>% AnyEmpty()
  matrix(data="", nrow=0, ncol=2) %>% as_tibble(.name_repair="minimal") %>% AnyEmpty()
  matrix(data=c("", NA_character_), nrow=1, ncol=2) %>% as_tibble(.name_repair="minimal") %>% AnyEmpty()
  tibble(x=1:6) %>% AnyEmpty()
  DNAString() %>% AnyEmpty()
  DNAString("acgt") %>% AnyEmpty()
  RNAStringSet() %>% AnyEmpty()
  RNAStringSet(c("acgu", "")) %>% AnyEmpty()
  RNAStringSet(c("acgu", "cgu")) %>% AnyEmpty()
  matrix(data="", nrow=0, ncol=0) %>% as.DNAbin() %>% AnyEmpty()
  matrix(c("a", NA_character_), nrow=1) %>% as.DNAbin() %>% AnyEmpty()
  matrix(c("a", "c"), nrow=1) %>% as.DNAbin() %>% AnyEmpty()
  list("a", character()) %>% as.DNAbin() %>% AnyEmpty()
  list("a", c("c", NA_character_)) %>% as.DNAbin() %>% AnyEmpty()
  list("a", c("c", "g")) %>% as.DNAbin() %>% AnyEmpty()
  c("A", NA_character_, "g") %>% as.DNAbin() %>% AnyEmpty()
}
