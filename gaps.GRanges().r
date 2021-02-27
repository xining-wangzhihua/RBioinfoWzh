# license: BSD 3-Clause License
# author: https://github.com/ywd5
# version: 20200227

# This is a function similar to IRanges::gaps().
# Although IRanges::gaps() already has a method for GRanges object, its behaviour is weird.
# The author don't know how to use S4 object-orientation in r, so the function's name looks like S3 method.

if("" == "get dependency"){
  require(magrittr); require(tibble);
  require(GenomicRanges); require(GenomeInfoDb);
  require(BiocGenerics); require(S4Vectors); require(IRanges);
}

gaps.GRanges <- function(x, start=1L, end=GenomeInfoDb::seqlengths(x), is_circular=GenomeInfoDb::isCircular(x)){
  # get dependency----------------------------------------------------------------------------------
  f1 <- function(v, g){
    #v <- 1:2; g <- GRanges(seqnames="a", ranges=IRanges(start=c(1, 1, 2), end=4:6), type=1:3, gene=4:6);
    lapply(X=v, FUN=IRanges::findOverlaps, subject=IRanges::ranges(g)) %>%
      lapply(FUN=S4Vectors::to) %>% lapply(FUN=function(i, xx){xx[i];}, xx=g) %>%
      lapply(FUN=S4Vectors::mcols) %>% lapply(FUN=function(xx, i){xx[i];}, i=c("type", "gene")) %>%
      lapply(FUN=function(xx){c(type=paste0(xx$type, collapse=","), name=paste0(xx$gene, collapse=","));}) %>%
      simplify2array() %>% base::t() %>% as_tibble()
  }
  # control the input of x--------------------------------------------------------------------------
  if(!is(object=x, class2="GRanges")){stop("x isn't GRanges");}
  S4Vectors::mcols(x) %<>% .[c("type", "gene")]
  GenomeInfoDb::seqnames(x) %>% levels.Rle() %>% {if(length(.) != 1){stop("all sequence names must be same");}}
  ans <- x[mcols(x)$type == "gene"]; tempo <- ans %in% IRanges::reduce(ans); if(!all(tempo)){
    # IRanges::setdiff(ans, IRanges::reduce(ans))
    # !(ans %in% IRanges::reduce(ans))
    ans[!tempo] %>% {.[order(IRanges::ranges(.))];} %>%
      {cbind(BiocGenerics::start(IRanges::ranges(.)), BiocGenerics::end(IRanges::ranges(.)),
             BiocGenerics::strand(.), S4Vectors::mcols(.))} %>%
      lapply(FUN=as.character) %>% simplify2array() %>% apply(MARGIN=1, FUN=paste0, collapse="\t") %>%
      c("the following gene regions overlap with each other or are within other gene regions:",
        "start\tend\tstrand\ttype\tgene", .) %>% paste0(., "\n") %>% warning()
  }; remove(ans, tempo);
  # control the input of start, end, is_circular----------------------------------------------------
  info <- list(start=as.integer(start), end=as.integer(end), is_circular=as.logical(is_circular)) %>%
    lapply(FUN=base::unname)
  remove(start, end, is_circular)
  anyNA(info) %>% {if(.){stop("start, end or is_circular can't be NA");}}
  lengths(info) %>% {if(any(. != 1)){stop("currently more than 1 genomes in x isn't supported");}}
  # preparation-------------------------------------------------------------------------------------
  y_backup <- IRanges::gaps(x=IRanges::ranges(x), start=info$start, end=info$end)
  l <- length(y_backup)
  y <- tibble(name="", up_pending=TRUE, up_boundary=S4Vectors::start(y_backup) - 1L, up_type="", up_name="",
              down_pending=TRUE, down_boundary=S4Vectors::end(y_backup) + 1L, down_type="", down_name="")
  # manipulate gaps on head and tail----------------------------------------------------------------
  gap_position <- c(head=(y$up_boundary[1] == info$start - 1L), tail=(y$down_boundary[l] == info$end + 1L))
  if(info$is_circular){
    if(all(gap_position)){y$up_boundary[1] <- y$up_boundary[l]; y$down_boundary[l] <- y$down_boundary[1];}
    if(gap_position[1] & (!gap_position[2])){y$up_boundary[1] <- info$end;}
    if((!gap_position[1]) & gap_position[2]){y$down_boundary[l] <- info$start;}
  }else{
    if(gap_position[1]){y$up_pending[1] <- FALSE; y$up_boundary[1] <- NA_integer_;}
    if(gap_position[2]){y$down_pending[l] <- FALSE; y$down_pending[l] <- NA_integer_;}
  }
  remove(gap_position)
  # manipulate gaps in the body---------------------------------------------------------------------
  y[y$up_pending, c("up_type", "up_name")] <- f1(v=y$up_boundary[y$up_pending], g=x)
  y[y$down_pending, c("down_type", "down_name")] <- f1(v=y$down_boundary[y$down_pending], g=x)
  y$name <- paste0("inter (", y$up_name, ") (", y$down_name, ")")
  # return the result-------------------------------------------------------------------------------
  y <- GRanges(seqnames=levels.Rle(GenomeInfoDb::seqnames(x)), ranges=y_backup, strand="*",
               type="gap", name=y$name, seqinfo=GenomeInfoDb::seqinfo(x),
               upstream_boundary=y$up_boundary, upstream_type=y$up_type, upstream_name=y$up_name,
               downstream_boundary=y$down_boundary, downstream_type=y$down_type, downstream_name=y$down_name)
  remove(f1, info, l, x, y_backup)
  return(y)
}
