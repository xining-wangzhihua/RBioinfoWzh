#
#last edited at 20200712

#a function similar to numbers::mGCD(),
#but calculates approximate GCDs, return a numeric vector,
#and may runs for a long time for its bad algorithm

#get dependency
if(FALSE){
  require(magrittr); require(tibble); require(dplyr);
  require(parallel)
}

message("\nmGCD_approximate() is very slow. change the algorithm and use rust\n")
warning("\nmGCD_approximate isn't stable, and isn't fully tested\n")

mGCD_approximate <- function(x){
  #control the input of x-------------------------------------------------------
  if(!is.integer(x)){
    stop("class(x) must be integer, not numeric. please use as.integer() first.")
  }
  if(length(x) == 0){stop("invalid x");}
  x %<>% stats::na.omit() %>% {.[. != 0];} %>% base::abs()
  if(length(x) == 0){stop("invalid x");}
  #convert x from a vector to tibble--------------------------------------------
  x %<>% base::table() %>% {tibble(x=as.integer(names(.)), weight=unname(.))}
  #get info array---------------------------------------------------------------
  if(nrow(x) != 1){
    divisor <- seq(from=1.1, to=ceiling(max(x$x) / 2), by=0.1)
    f_get_parameter <- function(dvs, dvdd){
      if(length(dvs) != 1){stop("invalid input");}
      qtt <- dvdd %/% dvs
      rmd <- dvdd %% dvs
      ans <- rmd > (dvs / 2)
      qtt[ans] <- qtt[ans] + 1
      rmd[ans] <- dvs - rmd[ans]
      remove(ans)
      offset_dvdd <- base::signif(x=rmd / dvdd, digits=3)
      offset_dvs <- base::signif(x=rmd / dvs, digits=3)
      return(cbind(qtt, offset_dvdd, offset_dvs))
    }
    if(length(divisor) > 5 * 10^5){
      a_cluster <- base::ceiling(length(divisor) / 10^5) %>%
        base::min(., detectCores(logical=TRUE) * 2) %>%
        parallel::makeCluster(spec=., outfile="")
      chunk_size <- ceiling(length(divisor) / 10^5 / length(a_cluster)) %>%
        {length(divisor) / length(a_cluster) / .} %>% ceiling()
      system.time(expr={
        info <- parSapply(cl=a_cluster, X=divisor, FUN=f_get_parameter, dvdd=x$x,
                          simplify="array", chunk.size=chunk_size)
      })
      stopCluster(cl=a_cluster); remove(a_cluster, chunk_size);
    }else{
      system.time(expr={
        info <- sapply(X=divisor, FUN=f_get_parameter, dvdd=x$x, simplify="array")
      })
    }
    dimnames(info) <- list(x=as.character(x$x),
                           parameter=c("quotient", "offset_dvdd", "offset_dvs"),
                           divisor=as.character(divisor))
    remove(divisor, f_get_parameter)
  }
  if(nrow(x) == 1) if(x$x != 1){
    info <- array(data=0, dim=c(1, 3, 1))
    info[,1,] <- 1
    dimnames(info) <- list(x=as.character(x$x),
                           parameter=c("quotient", "offset_dvdd", "offset_dvs"),
                           divisor=as.character(x$x))
  }
  if(nrow(x) == 1) if(x$x == 1){
    info <- array(data=0, dim=c(1, 3, 0))
    dimnames(info) <- list(x=as.character(x$x),
                           parameter=c("quotient", "offset_dvdd", "offset_dvs"),
                           divisor=character())
  }
  #analyse info array-----------------------------------------------------------
  if(dim(info)[3] > 0){
    info %<>% .[,"offset_dvs",,drop=FALSE] %>%
      apply(MARGIN=3, FUN=function(a){all(a < 0.4)}) %>% info[,,.,drop=FALSE]
  }
  if(dim(info)[3] > 0){
    info %<>% .[,"offset_dvdd",,drop=FALSE] %>%
      apply(MARGIN=3, FUN=function(a){all(a[a!=1]<0.1)}) %>% info[,,.,drop=FALSE]
  }
  if(dim(info)[3] > 0){
    f_ans <- function(datum, b){
      datum <- base::unname(datum[,"quotient"]) / b
      datum <- base::sort(datum[datum != 0])
      if(length(datum) == 0){stop("this may not be a bug, but a rare and intreasting situation");}
      if(length(datum) > 1){
        datum <- c(continuous=base::max((datum[-1] - head(datum, n=-1)) / head(datum, n=-1)),
                   range=(max(datum) - min(datum)) / min(datum))
      }else{datum <- c(continuous=0, range=0);}
      return(datum)
    }
    tulip <- info[,"quotient",,drop=FALSE] %>% apply(MARGIN=3, FUN=f_ans, b=x$x)
    tulip <- ((tulip["continuous",] < 0.05) & (tulip["range",] < 0.1))
    info %<>% .[,,tulip,drop=FALSE]
    remove(tulip, f_ans)
  }
  if(dim(info)[3] > 1){
    tulip <- dimnames(info)[[3]] %>% as.numeric()
    f_tulip <- function(i, datum, l){
      if(i %% 1000 == 0){cat(i, ", ");}
      msu <- datum[i]; spbu <- datum[(i + 1):l]; spbu <- spbu[(spbu %/% msu) > 1];
      return(all((spbu %% msu) > (msu * 0.1)))
    }
    i <- which(tulip > (max(tulip) / 2 * 1.11)) %>% .[1] %>% {.:1}
    dahlia <- rep(TRUE, times=length(tulip))
    cat("\nthe slowest part begins\n")
    system.time(expr={
      dahlia[i] <- base::sapply(X=i, FUN=f_tulip, datum=tulip, l=length(tulip), simplify="array")
    })
    cat("\nthe slowest part ends\n")
    info %<>% .[,,dahlia,drop=FALSE]
    remove(tulip, f_tulip, i, dahlia)
  }
  #split info according to quotient---------------------------------------------
  if(dim(info)[3] > 0){
    tulip <- info[,"quotient",,drop=FALSE] %>% base::unique(MARGIN=3)
    dahlia <- rep(x=list(NULL), times=dim(tulip)[3]); for(i in 1:length(dahlia)){
      dahlia[[i]] <- tulip[,"quotient",i,drop=TRUE] %>% base::unname()
    }
    remove(tulip)
    tulip <- base::apply(X=info, MARGIN=3,
                         FUN=function(a, b){match(list(unname(a[,"quotient",drop=TRUE])), b)},
                         b=dahlia) %>% base::unname()
    for(i in 1:length(dahlia)){
      dahlia[[i]] <- info[,,tulip == i,drop=FALSE]
    }
    info <- dahlia
    remove(tulip, dahlia, i)
  }else{info <- list();}
  #get GCD----------------------------------------------------------------------
  dahlia <- tibble(GCD=rep(0, times=length(info)), is_skew="", number_of_alternatives=0L)
  if(length(info) > 0) for(i in 1:length(info)){
    tulip <- info[[i]] %>% dimnames() %>% .[[3]] %>% as.numeric()
    tulip_mean <- info[[i]][,"quotient",1,drop=TRUE] %>% base::unname() %>% {. / x$x} %>%
      {.[. != 0]} %>% base::mean() %>% {1 / .}
    dahlia$GCD[i] <- base::abs(tulip - tulip_mean) %>% base::round(digits=3) %>%
      {which(. == min(.))} %>% tulip[.] %>% base::max()
    dahlia$number_of_alternatives[i] <- length(tulip)
    remove(tulip, tulip_mean)
  }; i <- 0; remove(i);
  if(length(info) == 0){
    warning("\ncan't find approximate GCD. an emptye tibble is returned\n")
  }
  info <- dahlia; remove(dahlia);
  #get info$is_skew-------------------------------------------------------------
  f_ans <- function(datum){
    msu <- datum$expected; spbu <- base::range(datum$observed);
    datum <- (msu < spbu[1]) | (msu > spbu[2])
    return(datum)
  }
  if(nrow(info) > 0) for(i in 1:nrow(info)){
    info$is_skew[i] <- (x$x / info$GCD[i]) %>% base::split(x=., f=base::round(.)) %>%
      {.[names(.) != 0]} %>% {tibble(expected=as.numeric(names(.)), observed=unname(.))} %>%
      apply(MARGIN=1, FUN=f_ans) %>% base::table(dnn=NULL) %>% unclass() %>% c() %>%
      c(., `TRUE`=0, `FALSE`=0) %>% {split(x=unname(.), f=names(.))} %>% lapply(FUN=sum) %>%
      {c(.[["TRUE"]], .[["TRUE"]] + .[["FALSE"]])} %>% {paste0(.[1], "/", .[2])}
  }; i <- 0; remove(i);
  remove(f_ans)
  #return results---------------------------------------------------------------
  info %<>% dplyr::arrange(dplyr::desc(GCD))
  return(info)
}

#here are the test codes
if(FALSE){
  #rnorm(n=10^8, mean=5, sd=0.05) %>% range()
  x <- c(rnorm(n=10, mean=0, sd=0.05), rnorm(n=10, mean=10, sd=0.05),
         rnorm(n=10, mean=15, sd=0.05), NA, NaN) %>%
    {. * 10^4} %>% base::round() %>% as.integer()
  x <- c(rnorm(n=10, mean=0, sd=0.05), rnorm(n=10, mean=10, sd=0.05),
         rnorm(n=10, mean=15, sd=0.05), NA, NaN) %>%
    {. * 10^3} %>% base::round() %>% as.integer()
  x <- c(rnorm(n=10, mean=0, sd=0.05), rnorm(n=10, mean=12, sd=0.05),
         rnorm(n=10, mean=16, sd=0.05), NA, NaN) %>%
    {. * 10^2} %>% base::round() %>% as.integer()
  x <- c(0, 1, NA, NaN) %>% sample(size=20, replace=TRUE) %>% as.integer()
  x <- c(0, 3, NA, NaN) %>% sample(size=20, replace=TRUE) %>% as.integer()
  x <- c(1, 18) %>% sample(size=20, replace=TRUE) %>% as.integer()
  x <- c(16, 18) %>% sample(size=20, replace=TRUE) %>% as.integer()
  x <- c(17, 18) %>% sample(size=20, replace=TRUE) %>% as.integer()
  x <- 17:31 %>% sample(size=50, replace=TRUE) %>% sort() %>% as.integer()
  #
  hist(x, breaks=100)
  system.time(expr={
    ans <- mGCD_approximate(x=x)
  })
  remove(x, mGCD_approximate, ans)
}
