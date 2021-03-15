#! /usr/bin/Rscript
# license: BSD 3-Clause License
# author: https://github.com/ywd5
# version: 20210314

if("" == "get dependency"){
  require(magrittr); require(tibble); require(dplyr); require(stringr); require(purrr);
  require(pegas)
  require(questionr)
  require(digest)#; require(openssl);
}

if("" == "test codes"){
  ans <- list.files(pattern="^output\\.xml$", full.names=TRUE, recursive=TRUE) %>%
    grep(pattern="first run", x=., value=TRUE)
  tempo <- tempfile(); phased_loci <- vector(mode="list", length=length(ans)); for(i in 1:length(ans)){
    ans[i] %>% read_xml() %>% xml_child(search="phased_loci") %>% xml_contents() %>%
      as.character() %>% cat(., file=tempo, append=FALSE, sep="")
    phased_loci[[i]] <- pegas::read.loci(file=tempo, header=TRUE, row.names=1L)
  }; file.remove(tempo); remove(tempo, i);
  #
  load(file="../can delete, unphased_loci.rdata")
  locus_position <- unphased_loci$locus_position
  not_varied_loci <- unphased_loci$not_varied_loci
  remove(unphased_loci)
  # pegas::is.snp(phased_loci) # this code will treate triple-phase locus as not snp
  # phased_loci %>% pegas::rarefactionplot()
  # ?pegas::haplotype()
  # ?pegas::haplotype.loci()
  # haplo %>% ape::as.DNAbin()
  ExamineRarePhasedHaplotypes(phased_loci=phased_loci, locus_position=locus_position,
                              not_varied_loci=not_varied_loci)
}

message("ExamineRarePhasedHaplotypes() may be able to tackle (deal with, cope with, handle) ",
        "polyploid, but currently only support diploid")
# this function can work with snp and microsatellite loci

ExamineRarePhasedHaplotypes <- function(
  phased_loci, locus_position,
  not_varied_loci=base::rep("varied", times=length(locus_position))
){
  # control the input-------------------------------------------------------------------------------
  if(!is(object=phased_loci, class2="loci")){stop("phased_loci must be a loci object");}
  if(any(dim(phased_loci) == 0)){stop("phased_loci is empty");}
  pegas::is.phased(phased_loci) %>% {if(!all(.)){warning("seperators of haplotypes should be |, not /");}}
  if("" != "ensure data are from diploid"){
    getPloidy(phased_loci) %>% {any(. != 2)} %>% {if(.){stop("only diploid data are supported");}}
  }
  #
  if(!is.numeric(locus_position)){stop("locus_position must be integers");}
  locus_position %<>% as.integer() %>% base::sort()
  if(anyNA(locus_position)){stop("locus_position must be integers");}
  if(any(locus_position < 1L)){stop("locus_position can't be smaller than 1");}
  if(anyDuplicated(locus_position) != 0){stop("locus_position can't contain duplicated items");}
  #
  if(!is.character(not_varied_loci)){stop("not_varied_loci must be a character vector");}
  if(length(not_varied_loci) == 0){stop("not_varied_loci is empty");}
  if(anyNA(not_varied_loci)){stop("not_varied_loci can't contain NA");}
  # control the relationship among inputs-----------------------------------------------------------
  (ncol(phased_loci) != length(locus_position)) %>% {if(.){stop("locus number are consistent");}}
  not_varied_loci %>% {which(. == "varied")} %>% {length(.) != length(locus_position)} %>%
    {if(.){stop("locus_position doesn't match with not_varied_loci");}}
  # get info, remove locus_position-----------------------------------------------------------------
  info <- phased_loci %>% {list(name=rownames(.), locus_position=locus_position, haplo_freq=integer(),
                                n_indiv=nrow(.), n_locus=ncol(.), n_haplo=0L)}
  remove(locus_position)
  # convert phased_loci from loci into 2 matrices--------------------------------------------------------------
  if("" == "method 1, use pegas::haplotype()"){
    # phased_loci %>% {pegas::haplotype(x=., locus=1:ncol(.), compress=TRUE)}
    ans <- {1:nrow(phased_loci)} %>% {. * 2} %>% {list(. - 1, .)}
    phased_loci %<>% {pegas::haplotype(x=., locus=1:ncol(.), compress=FALSE)} %>% base::t(x=.) %>%
      {mapply(FUN=function(xx, ii){xx[ii,]}, xx=list(.), ii=ans, SIMPLIFY=FALSE, USE.NAMES=FALSE)} %>%
      lapply(FUN=base::`rownames<-`, value=rownames(phased_loci))
    remove(ans)
  }
  if("" != "method 2, use pegas::loci2alleles()"){
    ans <- {1:ncol(phased_loci)} %>% {. * 2} %>% {list(. - 1, .)}
    phased_loci %<>% loci2alleles() %>%
      {base::`colnames<-`(x=., value=base::sub(pattern="\\.(1|2)$", replacement="", x=colnames(.)))} %>%
      {mapply(FUN=function(xx, ii){xx[,ii]}, xx=list(.), ii=ans, SIMPLIFY=FALSE, USE.NAMES=FALSE)}
    remove(ans)
  }
  # split phased_loci into haplo and info_indiv, get info$haplo_freq, info$n_haplo------------------------
  phased_loci %<>% lapply(FUN=base::`colnames<-`, value=NULL) %>% lapply(FUN=purrr::array_tree, margin=1)
  #
  haplo <- phased_loci %>% lapply(FUN=base::unname) %>% {c(.[[1]], .[[2]])} %>%
    {list(all=., unique=base::unique(.))} # create haplo in list format
  haplo <- base::match(x=haplo$all, table=haplo$unique) %>% base::table() %>%
    as.character() %>% setNames(object=haplo$unique, nm=.) %>% simplify2array() %>%
    base::t(x=.) # get haplo occurrence times, change haplo from list into matrix
  str_rank <- function(xx){
    base::unique(xx) %>% stringr::str_sort(x=., numeric=TRUE) %>% base::match(x=xx, table=.)
  }
  haplo %<>% apply(X=., MARGIN=2, FUN=str_rank) %>% purrr::array_tree(array=., margin=2) %>%
    c(list(-as.integer(rownames(haplo))), .) %>% base::do.call(what=base::order, args=.) %>%
    haplo[.,] # sort haplo according to frequency and alphanumeric order
  remove(str_rank)
  info$haplo_freq <- rownames(haplo) %>% as.integer()
  info$n_haplo <- length(info$haplo_freq)
  #
  info_indiv <- tibble(name=info$name, genotype="",
                       haplotype_1=base::match(x=phased_loci[[1]], table=purrr::array_tree(array=haplo, margin=1)),
                       haplotype_2=base::match(x=phased_loci[[2]], table=array_tree(array=haplo, margin=1)))
  info_indiv[info_indiv$haplotype_1 > info_indiv$haplotype_2,3:4] %<>% base::rev(x=.)
  #
  haplo %<>% {base::`rownames<-`(x=., value=paste0("hap", 1:nrow(.)));}
  info_indiv$haplotype_1 %<>% {rownames(haplo)[.]}
  info_indiv$haplotype_2 %<>% {rownames(haplo)[.]}
  info_indiv %<>% dplyr::mutate(genotype=paste0(haplotype_1, "/", haplotype_2))
  #
  remove(phased_loci)
  # create output_directory-------------------------------------------------------------------------
  output_directory <- base::tempfile() %>% base::normalizePath(winslash="/", mustWork=FALSE)
  base::dir.create(path=output_directory, recursive=FALSE)
  dir.exists(paths=output_directory) %>% {if(!.){stop("can't create directory: ", output_directory);}}
  if(!grepl(pattern="/$", x=output_directory)){output_directory <- paste0(output_directory, "/");}
  # prepare an array for .xml-----------------------------------------------------------------------
  rare_haplo <- c("position", "residue", "residue_frequency", "other_residue_frequencies",
    "excised_frequency") %>% list(five_items=., locus=NULL, haplotype=rownames(haplo)) %>%
    array(data="", dim=c(5L, info$n_locus, info$n_haplo), dimnames=.)
  rare_haplo[1,,] <- info$locus_position %>% as.character() %>% list() %>% rep(times=info$n_haplo) %>% simplify2array()
  rare_haplo[2,,] <- base::t(haplo)
  #
  f_frequency <- function(xx, weight){
    # xx <- sample(letters[1:6], size=20, replace=TRUE)
    # weight <- sample(1:6, size=20, replace=TRUE)
    if(!is.character(xx) | length(xx) == 0){stop("invalid input");}
    if(!is.integer(weight) | length(weight) != length(xx)){stop("invalid input");}
    if(any(weight < 1)){stop("invalid input");}
    xx <- questionr::wtd.table(x=xx, weights=weight, digits=10) %>% # wtd.table() will sort in lexical order
      {setNames(object=as.integer(.), nm=names(.));}
    ans <- rep("0", times=length(xx))
    if(length(xx) > 1) for(i in 1:length(xx)){
      ans[i] <- xx[-i] %>% {paste0(unname(.), names(.), collapse=",");}
    }else{
      # Do nothing. Only happen when there's only 1 locus.
    }; i <- 0; remove(i);
    xx <- list(xx, setNames(object=ans, nm=names(xx)))
    return(xx) # remove(xx, weight, ans)
  }
  haplo_minus <- haplo %>% lapply(X=1:ncol(haplo), FUN=function(ii, xx){
    # apply(X=xx[,-ii], MARGIN=1, FUN=paste0, collapse="")
    xx <- purrr::array_tree(array=xx[,-ii], margin=1)
    xx <- base::match(x=xx, table=base::unique(xx))
    xx <- paste0("_", xx)
    return(xx)
  }, xx=.) %>% simplify2array()
  rare_haplo[3,,] <- haplo %>% apply(X=., MARGIN=2, FUN=f_frequency, weight=info$haplo_freq) %>%
    lapply(FUN=function(xx){xx[[1]]}) %>%
    mapply(FUN=function(xx, ii){unname(xx[ii]);}, xx=., ii=purrr::array_tree(array=haplo, margin=2),
           SIMPLIFY=TRUE, USE.NAMES=FALSE) %>% base::t(x=.)
  rare_haplo[4,,] <- haplo %>% apply(X=., MARGIN=2, FUN=f_frequency, weight=info$haplo_freq) %>%
    lapply(FUN=function(xx){xx[[2]]}) %>%
    mapply(FUN=function(xx, ii){unname(xx[ii]);}, xx=., ii=purrr::array_tree(array=haplo, margin=2),
           SIMPLIFY=TRUE, USE.NAMES=FALSE) %>% base::t(x=.)
  rare_haplo[5,,] <- haplo_minus %>% apply(X=., MARGIN=2, FUN=f_frequency, weight=info$haplo_freq) %>%
    lapply(FUN=function(xx){xx[[1]]}) %>%
    mapply(FUN=function(xx, ii){unname(xx[ii]);}, xx=., ii=purrr::array_tree(array=haplo_minus, margin=2),
           SIMPLIFY=TRUE, USE.NAMES=FALSE) %>% base::t(x=.)
  remove(f_frequency, haplo_minus)
  # get info_geno, modify info_indiv$genotype, get info_haplo------------------------------------------
  info_geno <- info_indiv %>% {base::split(x=.$name, f=.$genotype)} %>%
    {tibble(name=names(.), haplotype_1="", haplotype_2="", freq=lengths(.),
            occurrence_in_individual=lapply(X=unname(.), FUN=stringr::str_sort, numeric=TRUE))}
  ans <- info_geno$name %>% base::strsplit(split="/") %>% simplify2array()
  info_geno$haplotype_1 <- ans[1,]
  info_geno$haplotype_2 <- ans[2,]
  remove(ans)
  info_geno %<>% .[c("haplotype_1", "haplotype_2")] %>%
    lapply(FUN=function(xx){base::match(x=xx, table=stringr::str_sort(x=unique(xx), numeric=TRUE));}) %>%
    {base::order(.[[1]], .[[2]])} %>% info_geno[.,] %>% dplyr::mutate(name=paste0("gen", 1:length(name)))
  #
  # if the haplotype names contain "/", there may be bugs.
  info_indiv$genotype <- info_geno %>% {paste0(.$haplotype_1, "/", .$haplotype_2)} %>%
    setNames(object=info_geno$name, nm=.) %>% .[info_indiv$genotype] %>% unname()
  #
  min_r <- rare_haplo[3,,] %>% apply(MARGIN=2, FUN=as.integer) %>% apply(MARGIN=2, base::min) %>% unname()
  max_e <- rare_haplo[5,,] %>% apply(MARGIN=2, FUN=as.integer) %>% apply(MARGIN=2, base::max) %>% unname()
  occurrence_in_i <- info_indiv %>% {base::split(x=c(.$name, .$name), f=c(.$haplotype_1, .$haplotype_2))} %>%
    {.[stringr::str_order(x=names(.), numeric=TRUE)]} %>% lapply(FUN=base::unique) %>%
    lapply(FUN=stringr::str_sort, numeric=TRUE)
  occurrence_in_g <- info_geno %>% {base::split(x=c(.$name, .$name), f=c(.$haplotype_1, .$haplotype_2))} %>%
    {.[stringr::str_order(x=names(.), numeric=TRUE)]} %>% lapply(FUN=base::unique) %>%
    lapply(FUN=stringr::str_sort, numeric=TRUE)
  info_haplo <- occurrence_in_i %>% {base::replace(x=., list=lengths(.) > 6L, values=list("..."))} %>%
    lapply(FUN=paste0, collapse=", ") %>% unlist() %>% unname() %>%
    tibble(name=rownames(haplo), freq=info$haplo_freq, occurrence_in_individual=occurrence_in_i,
           occurrence_in_genotype=occurrence_in_g, min_residue_freq=min_r, max_excised_freq=max_e, file=.) %>%
    dplyr::mutate(file=paste0(freq, " ", min_residue_freq, " ", max_excised_freq, "; ", name, "; ", file, ";.xml")) %>%
    dplyr::mutate(file=paste0(output_directory, file))
  remove(min_r, max_e, occurrence_in_i, occurrence_in_g)
  # fill not_varied_loci in the array---------------------------------------------------------------
  dahlia <- array(data="", dim=c(5, length(not_varied_loci), info$n_haplo), dimnames=dimnames(rare_haplo))
  dahlia[2,,] <- not_varied_loci %>%
    purrr::map_if(.x=., .p=(. != "varied"), .f=stringr::str_trunc, width=15, side="center") %>%
    unlist() %>% list() %>% rep(times=info$n_haplo) %>% simplify2array()
  dahlia[,not_varied_loci == "varied",] <- rare_haplo
  rare_haplo <- dahlia
  remove(dahlia)
  # prepare the array for pretty print--------------------------------------------------------------
  # rare_haplo %>% apply(MARGIN=c(2, 3), FUN=function(xx){
  #   stringr::str_pad(string=xx, width=max(nchar(xx)) + 2, side="both")
  # }) %>% {base::`dimnames<-`(x=., value=base::replace(x=dimnames(rare_haplo), list=1L, values=dimnames(rare_haplo)[1]))}
  for(i in 1:length(not_varied_loci)) for(j in 1:info$n_haplo){
    rare_haplo[,i,j] %<>% {stringr::str_pad(string=., width=max(nchar(.)) + 2, side="both")}
  }; remove(i, j);
  dimnames(rare_haplo)[[1]] %<>% stringr::str_pad(string=., width=max(nchar(.)) + 2, side="right") %>% paste0(., ": ")
  # export and remove the array---------------------------------------------------------------------
  for(i in 1:nrow(info_haplo)){
    cat("<root>\n",
        "<!--file name's meaning: haplotype_frequency minimum_residue_frequency maximum_excised_frequency; ",
        "haplotype_name; occurrence_in_individuals_or_..._if_too_many;.xml-->\n",
        "<haplotype_name>", info_haplo$name[i], "</haplotype_name>\n",
        "<haplotype_frequency>", info_haplo$freq[i], "</haplotype_frequency>\n",
        "<occurrence_in_individual>", paste0(info_haplo$occurrence_in_individual[[i]], collapse=", "), "</occurrence_in_individual>\n",
        "<minimum_residue_frequency>", info_haplo$min_residue_freq[i], "</minimum_residue_frequency>\n",
        "<maximum_excised_frequency>", info_haplo$max_excised_freq[i], "</maximum_excised_frequency>\n",
        "<array>\n", file=info_haplo$file[i], append=FALSE, sep="")
    write.table(x=rare_haplo[,,i], file=info_haplo$file[i], append=TRUE, quote=FALSE,
                sep="", row.names=TRUE, col.names=FALSE)
    cat("</array>\n", "</root>\n", file=info_haplo$file[i], append=TRUE, sep="")
  }; remove(i);
  remove(rare_haplo)
  # export haplo, haplo_full, info_haplo, info_geno, info_indiv, etc to info.xml---------------------------
  haplo_full <- not_varied_loci %>% list() %>% rep(times=info$n_haplo) %>% simplify2array() %>%
    base::t(x=.) %>% base::`rownames<-`(x=., value=rownames(haplo))
  haplo_full[,not_varied_loci == "varied"] <- haplo
  info_haplo %<>% dplyr::select(name, freq, occurrence_in_individual, occurrence_in_genotype) %>%
    dplyr::mutate(occurrence_in_individual=lapply(X=occurrence_in_individual, FUN=paste0, collapse=", "),
                  occurrence_in_genotype=lapply(X=occurrence_in_genotype, FUN=paste0, collapse=", ")) %>%
    dplyr::mutate(occurrence_in_individual=unname(unlist(occurrence_in_individual)),
                  occurrence_in_genotype=unname(unlist(occurrence_in_genotype)))
  info_geno %<>% dplyr::mutate(occurrence_in_individual=lapply(X=occurrence_in_individual, FUN=paste0, collapse=", ")) %>%
    dplyr::mutate(occurrence_in_individual=unname(unlist(occurrence_in_individual)))
  file_ans <- paste0(output_directory, "info.xml")
  cat("<root>\n",
      "<n_individual>", info$n_indiv, "</n_individual>\n",
      "<n_genotype>", nrow(info_geno), "</n_genotype>\n",
      "<n_haplotype>", info$n_haplo, "</n_haplotype>\n",
      "<locus_position>", paste0(info$locus_position, collapse=" "), "</locus_position>\n",
      "<haplotype>\n", file=file_ans, append=FALSE, sep="")
  write.table(x=haplo, file=file_ans, append=TRUE, sep="\t", row.names=TRUE, col.names=FALSE, fileEncoding="UTF-8")
  cat("</haplotype>\n", "<haplotype_full>\n", file=file_ans, append=TRUE, sep="")
  write.table(x=haplo_full, file=file_ans, append=TRUE, sep="\t", row.names=TRUE, col.names=FALSE, fileEncoding="UTF-8")
  cat("</haplotype_full>\n", "<info_individual>\n", file=file_ans, append=TRUE, sep="")
  suppressWarnings(expr={
    write.table(x=info_indiv, file=file_ans, append=TRUE, sep="\t", row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")
  })
  cat("</info_individual>\n", "<info_genotype>\n", file=file_ans, append=TRUE, sep="")
  suppressWarnings(expr={
    write.table(x=info_geno, file=file_ans, append=TRUE, sep="\t", row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")
  })
  cat("</info_genotype>\n", "<info_haplotype>\n", file=file_ans, append=TRUE, sep="")
  suppressWarnings(expr={
    write.table(x=info_haplo, file=file_ans, append=TRUE, sep="\t", row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")
  })
  cat("</info_haplotype>\n", "</root>\n", file=file_ans, append=TRUE, sep="")
  remove(file_ans)
  # calculate hash value of output_directory--------------------------------------------------------
  # file_name <- list.files(path=output_directory, full.names=TRUE, no..=TRUE)
  # file_sha1 <- vector(mode="character", length=length(file_name))
  # for(i in 1:length(file_name)){
  #   file_con <- base::file(description=file_name[i], open="rb")
  #   file_sha1[i] <- as.character(openssl::sha1(x=file_con))
  #   base::close(con=file_con)
  #   remove(file_con)
  # }; remove(i);
  old_files <- output_directory %>% list.files(path=., full.names=TRUE, no..=TRUE)
  output_directory <- old_files %>% stringr::str_sort(numeric=TRUE) %>%
    lapply(FUN=digest::digest, algo="crc32", serialize=FALSE, file=TRUE) %>%
    unlist() %>% paste0(collapse="") %>% digest::digest(algo="crc32", serialize=TRUE) %>%
    paste0("./", info$n_haplo, "_", ., strftime(x=Sys.time(), format="_%H%M%S/")) %>%
    base::normalizePath(winslash="/", mustWork=FALSE)
  dir.create(path=output_directory, recursive=FALSE)
  if(!dir.exists(output_directory)){stop("can't create directory: ", output_directory);}
  file.copy(from=old_files, to=output_directory, overwrite=TRUE) %>%
    {if(!all(.)){stop("error when coping files to new directory");}} # try file.rename() in revision
  remove(old_files)
  # clean environment and return output_directory---------------------------------------------------
  remove(not_varied_loci, haplo, haplo_full, info, info_indiv, info_geno, info_haplo)
  message("\nresults are written to this directory:\n", output_directory, "\n")
  return(output_directory)
}
