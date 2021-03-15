#a simple function similar to ape::read.nexus(), but not only work for phylogenetic datum

#get dependency
if(FALSE){
  require(magrittr); require(stringr);
  require(IRanges)
}

message("\nthere are packages in language cpp or others, which can read files in nexus ",
        "format. please learn cpp and the package Rcpp, and then revise the codes\n")
message("\ncomments (which begin with #) are not removed, for \"#\" may locate in ",
        "quotations, e.g. \"\'a#b\'\", \"(a#b)\", etc, which makes code too complex\n")

read.nexus2 <- function(file_name){
  #control the input of file_name-----------------------------------------------
  if(!is.character(file_name)){stop("file_name must be a string");}
  if(length(file_name) != 1){stop("length(file_name) must be 1");}
  if(!file.exists(file_name)){stop("the file specified by file_name doesn't exist");}
  x <- readLines(con=file_name)
  if(length(x) > 0){
    x %<>% stringr::str_trim(side="both") %>% {.[. != ""];}
  }
  p_begin <- "^BEGIN +([[:alnum:]_]+)( *;|)$"
  p_end <- "^END( +[[:alnum:]_]+ *| *)(;|)$"
  #split the function in 3 conditions-------------------------------------------
  if(length(x) == 0){status = 0L;}
  if(length(x) > 0){
    if(any(grepl(pattern=p_begin, x=x))){status <- 2;}else{status <- 1;}
  }
  if(status == 0){result <- list();}
  if(status == 1){result <- list(x);}
  if(status == 2){
    #get document structure
    msu <- grep(pattern=p_begin, x=x)
    spbu <- grep(pattern=p_end, x=x)
    if(length(msu) != length(spbu)){
      stop("number of BEGIN patterns must equal with numbers of END patterns")
    }
    if(!all(msu < spbu)){stop("all BEGIN patterns must locate before END patterns");}
    if( !all(msu[-1] > head(spbu, n=-1)) ){
      stop("some (BEGIN text END) sections overlap with each other")
    }
    msu <- IRanges::IRanges(start=msu, end=spbu)
    spbu <- IRanges::IRanges(start=1, end=length(x)) %>% IRanges::setdiff(x=., y=msu)
    names(msu) <- rep("hit", times=length(msu))
    names(spbu) <- rep("gap", times=length(spbu))
    itmo <- c(msu, spbu) %>% sort()#c(msu, spbu) %>% {.[order(.)];}
    if(base::is.unsorted( start(itmo) )){
      stop("now sort() doesn't work for IRange objects. please revise the code")
    }
    itmo <- data.frame(hit=(names(itmo) == "hit"), start=start(itmo), end=end(itmo))
    remove(msu, spbu)
    #split the document according to the structure
    result <- list(character()) %>% rep(times=nrow(itmo)) %>% setNames(nm=rep("", times=nrow(itmo)))
    for(i in 1:nrow(itmo)) if(itmo$hit[i]){
      ans <- x[itmo$start[i]:itmo$end[i]]
      msu <- sub(pattern=p_begin, replacement="\\1", x=ans[1])
      spbu <- sub(pattern=p_end, replacement="\\1", x=tail(ans, n=1)) %>%
        stringr::str_trim(side="both")
      if(spbu != "") if(msu != spbu){
        stop("END pattern names, if there are, must be same as BEGIN pattern names")
      }
      names(result)[i] <- msu
      result[[i]] <- ans[-c(1, length(ans))]
      remove(ans, msu, spbu)
    }; remove(i);
    for(i in 1:nrow(itmo)) if(!itmo$hit[i]){
      result[[i]] <- x[ itmo$start[i]:itmo$end[i] ]
    }; remove(i);
    remove(itmo)
  }
  #return the result------------------------------------------------------------
  remove(file_name, x, p_begin, p_end, status)
  return(result)
}

#here are the test codes
if(FALSE){
  file_name <- "demo phase output.txt"
  file_name <- "demo dnasp output.nex"
}
