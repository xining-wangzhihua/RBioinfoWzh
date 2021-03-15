#a function similar to pegas::MMD(), and more stable, more pretty,
#but accept x as a numeric vector. which is c(dist.dna(x))

#get dependency
if(FALSE){
  require(magrittr); require(tibble); require(dplyr);
  require(ggplot2)
}

MMD2 <- function(x, output_directory=tempfile(fileext="_MMD2")){
  #control the input of x-------------------------------------------------------
  if(!is.numeric(x)){stop("invalid x");}
  if(length(x) == 0){stop("invalid x");}
  if(anyNA(x)){stop("invalid x");}
  if(any(x < 0)){stop("invalid x");}
  #control the input of output--------------------------------------------------
  if(!is.character(output_directory)){stop("invalid output");}
  if(length(output_directory) != 1){stop("invalid output");}
  if(anyNA(output_directory)){stop("invalid output");}
  output_directory %<>% base::normalizePath(winslash="/", mustWork=FALSE)
  if(!dir.exists(output_directory)){dir.create(output_directory);}
  #plot-------------------------------------------------------------------------
  x_backup <- x
  f_ans <- function(datum, n){
    if(length(datum) != 1){stop("invalid input");}
    if(length(n) == 0){stop("invalid input");}
    #result <- (datum^i) / (datum + 1)^(i + 1)#function used in pegas::MMD()
    result <- (datum / (datum + 1))^n / (datum + 1)
    return(result)
  }
  for(i in c("smooth", "moderate", "steep")) for(j in 10:30){
    #change the scale of x
    x <- x_backup * j / diff(range(x_backup))
    #prepare data
    msu <- switch(i, smooth=2, moderate=4.5, steep=7) %>% {sqrt(var(x)) / .} %>%
      stats::density(x=x, bw=., kernel="gaussian", from=0, to=max(x) * 1.2) %>%
      {tibble(x=.$x, y=.$y, group="msu", alpha="placeholder")}
    spbu <- base::seq(from=0, to=max(x) * 1.2, length.out=2^9) %>%
      {tibble(x=., y=f_ans(datum=mean(x), n=.), group="spbu", alpha="placeholder")}
    itmo <- dplyr::bind_rows(msu, spbu)
    remove(msu, spbu)
    #plot
    result <- ggplot() +
      ggplot2::geom_histogram(mapping=aes(x=distance, y=after_stat(density)),
                              size=0, colour=NA_character_,#colour="#FFFFFF00" is also okay
                              fill=grDevices::gray(level=0.5, alpha=0.6),
                              data=tibble(distance=x), bins=j, show.legend=FALSE) +
      ggplot2::geom_line(mapping=aes(x=x, y=y, group=group, colour=group,
                                     linetype=group, alpha=alpha),
                         show.legend=c(colour=TRUE, linetype=TRUE, alpha=FALSE),
                         data=itmo, size=1) +
      ggplot2::scale_colour_manual(values=c(msu="blue", spbu="red"), name="",
                                   labels=c("observed", "expected")) +#name should be NULL
      ggplot2::scale_linetype_manual(values=c(msu="solid", spbu="dashed"), name="",
                                     labels=c("observed", "expected")) +#name should be NULL
      ggplot2::scale_alpha_manual(values=c(placeholder=0.7)) +#should contain parameter "name=NULL"
      ggplot2::theme(text=element_text(size=30),
                     legend.position=c(0.72, 0.85),
                     #legend.position="right",
                     legend.title=element_blank(),#maybe omitted when setting "name=NULL" instead of "name=\"\""
                     #legend.spacing=grid::unit(x=0, units="mm"),#inheritation error
                     legend.spacing.y=grid::unit(x=0, units="mm"),
                     legend.margin=ggplot2::margin(t=3, r=3, b=3, l=3, unit="mm"),
                     legend.key=element_rect(fill=gray(level=0.9)),
                     legend.key.height=grid::unit(x=10, units="mm"),
                     legend.key.width=grid::unit(x=30, units="mm"),
                     plot.margin=grid::unit(x=c(5, 5, 5, 5), units="mm"))
    ggsave(filename=paste0(output_directory, "/", i, ", ", j, " bins.png"),
           plot=result, device="png", width=20, height=15, units="cm", dpi=300)
    remove(itmo, result)
  }; remove(i, j);
  remove(f_ans)
  #clean environment and invisible(NULL)----------------------------------------
  remove(x, x_backup, output_directory)
  invisible(NULL)
}

#here are the test codes
if(FALSE){
  x <- runif(n=20, min=10, max=60) %>% base::rep(., times=sample(x=2:6, size=20, replace=TRUE)) %>%
    base::round(digits=2)
  MMD2(x=x)
}
