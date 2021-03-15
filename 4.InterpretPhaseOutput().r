#! /usr/bin/Rscript
# license: BSD 3-Clause License
# author: https://github.com/ywd5
# version: 20210314

if("" == "get dependency"){
  require(magrittr); require(stringr); require(tibble); require(dplyr); require(purrr);
  require(pegas)
  source(file=R.home(component="RBioinfoZm/read.nexus2().r")); require(IRanges);
  # source(file=R.home(component="RBioinfoZm/as.DNAStringSet.DNAbin().r"))
}

if("" == "test codes"){
  x <- list.dirs(path="C:/Users/wangz/OneDrive/github/RBioinfoZm/new phase", full.names=TRUE,
            recursive=FALSE) %>% {.[grep(pattern="^first run -MR.+[[:digit:]]$", x=basename(.))]} %>%
    paste0(., "/output.nexus") %>% .[1]
  x %<>% InterpretPhaseOutput()
  x %<>% InterpretPhaseOutput(substitute_with_integers=NA_character_)
}

if("" == "test codes for reading loci inside .xml"){
  ans <- tempfile()
  xml2::read_xml(x="../first run -MR0, 1/output.xml") %>%
    xml2::xml_child(x=., search="phased_loci") %>% xml_contents() %>%
    as.character() %>% cat(., file=ans, sep="", append=FALSE)
  pegas::read.loci(file=ans, header=TRUE, row.names=1L)
  identical(.Last.value, phased_loci$phased_loci)
  file.remove(ans); remove(ans);
}

message("add a part to determine which result is the best one")
warning("using to much %>%, apply(), lapply(), simplify2array(), etc may have bugs, when ",
        "data is small (of length 0, 1, 2 etc). the robustness isn't tested.")

InterpretPhaseOutput <- function(
  x, substitute_with_integers=c("./substitute_with_integers.tsv", NA_character_)[1],
  output_file=base::sub(pattern="\\.(nexus|txt)$", replacement=".xml", x=x)
){
  # control the input-------------------------------------------------------------------------------
  if(!is.character(x)){stop("x should be a file name");}
  if(length(x) != 1L){stop("x should be a file name");}
  if(is.na(x)){stop("x can't be NA");}
  if(!file.exists(x)){stop("file specified by x doesn't exist");}
  if(length(substitute_with_integers) != 1){stop("length(substitute_with_integers) must be 1");}
  substitute_with_integers %<>% as.character()
  if(!is.na(substitute_with_integers)) if(!file.exists(substitute_with_integers)){
    warning("The file specified by substitute_with_integers doesn't exist. It's not used.")
    substitute_with_integers <- NA_character_
  }
  if(!is.character(output_file)){stop("output_file should be a file name");}
  if(length(output_file) != 1L){stop("output_file should be a file name");}
  if(is.na(output_file)){stop("output_file can't be NA");}
  if(file.exists(output_file)){warning("file specified by output_file already exists. it will be overwritten.");}
  # read in x---------------------------------------------------------------------------------------
  x %<>% read.nexus2(file_name=.)
  c("INPUT_SUMMARY", "LIST_SUMMARY", "BESTPAIRS_SUMMARY", "BESTPAIRS1") %>% {. %in% names(x)} %>%
    {!all(.)} %>% {if(.){stop("Structure of PHASE output has changed. Please revise codes.");}}
  c("", "COMMAND_LINE", "OUTFILE_LIST", "INPUT_SUMMARY", "", "LIST_SUMMARY", "", "BESTPAIRS_SUMMARY",
    "", "BESTPAIRS1", "", "BESTPAIRS2", "", "PHASEPROBS") %>% identical(., names(x)) %>%
    {if(!.){warning("Structure of PHASE output has some trivial (subtle) changes. Please revise codes.");}}
  x %<>% .[c("INPUT_SUMMARY", "LIST_SUMMARY", "BESTPAIRS_SUMMARY", "BESTPAIRS1")]
  if(length(x$BESTPAIRS_SUMMARY) * 3L != length(x$BESTPAIRS1)){stop("PHASE output is disrupted (chaotic)");}
  # in old codes, where x_old[[i]] is x, you can run "x_old %<>% unique()"
  # x$INPUT_SUMMARY # c("Number of Individuals", "Number of Loci", "Positions of Loci")
  # utils::head(x$LIST_SUMMARY, n=4) # x$LIST_SUMMARY structure is "hap_index hap_seq ?probability?"
  # utils::head(x$BESTPAIRS_SUMMARY, n=4) # x$BESTPAIRS_SUMMARY structure is "individual: (hap_index, hap_index)"
  # utils::head(x$BESTPAIRS1, n=4) # x$BESTPAIRS1 structure is "0 individual", "hap_seq", "hap_seq"
  x %<>% setNames(nm=c("position", "hap", "index", "seq"))
  # get x$position----------------------------------------------------------------------------------
  # "Positions of loci" will always appear in the output.
  # When positions aren't designated, it will be 0, 1000, 2000, etc.
  if(length(x$position) != 3){stop("Structure of PHASE output has changed. Please revise codes.");}
  x$position %<>% base::strsplit(split=": ", fixed=TRUE) %>% simplify2array()
  (dim(x$position) != c(2L, 3L)) %>% {if(any(.)){stop("Structure of PHASE output has changed. Please revise codes.");}}
  (x$position[1,] != c("Number of Individuals", "Number of Loci", "Positions of loci")) %>%
    {if(any(.)){stop("Structure of PHASE output has changed. Please revise codes.");}}
  x$position %<>% .[2,3] %>% stringr::str_trim(string=., side="both") %>%
    base::strsplit(split=" ") %>% .[[1]] %>% as.integer()
  anyNA(x$position) %>% {if(.){stop("PHASE output is disrupted (chaotic)");}}
  # change character vector into matrix-------------------------------------------------------------
  p_hap <- "^([[:digit:]]+) +([^ ][ [:alnum:]?.+-]+[^ ]|[^ ]{1,2}) +[[:digit:].]+$"
  if(!all(grepl(pattern=p_hap, x=x$hap))){stop("regular expression don't match output format");}
  x$hap %<>% stringr::str_match(string=., pattern=p_hap) %>% .[,-1]
  remove(p_hap)
  if(any(as.character(1:nrow(x$hap)) != x$hap[,1])){stop("PHASE output is disrupted (chaotic)");}
  if(anyDuplicated(x$hap[,2]) != 0){stop("PHASE output is disrupted (chaotic)");}
  ans <- x$hap[,2] %>% stringr::str_squish(string=.) %>% strsplit(split=" ")
  ans %>% lengths() %>% unique() %>% {if(length(.) != 1L){stop("PHASE output is disrupted (chaotic)");}}
  x$hap[,2] <- ans %>% simplify2array() %>% purrr::array_tree(array=., margin=1) %>%
    purrr::map_if(.x=., .p=function(xx){anyNA(suppressWarnings(as.integer(xx)));},
                  .f=function(xx){base::strsplit(x=xx, split="");}) %>%
    lapply(FUN=function(xx){lapply(X=xx, FUN=paste0, collapse=" ");}) %>%
    simplify2array() %>% apply(MARGIN=1, FUN=paste0, collapse=" ")
  remove(ans)
  #
  p_index <- "^([^:]+): \\(([[:digit:]]+),([[:digit:]]+)\\)$"
  if(!all(grepl(pattern=p_index, x=x$index))){stop("regular expression don't match output format");}
  x$index %<>% stringr::str_match(string=., pattern=p_index) %>% .[,-1]
  remove(p_index)
  #
  x$seq <- ((1:nrow(x$index)) * 3) %>% {cbind(x$seq[. - 2], x$seq[. - 1], x$seq[.])}
  p_seq <- "^0 ([^ ]+)$"
  if(!all(grepl(pattern=p_seq, x=x$seq[,1]))){stop("regular expression don't match output format");}
  x$seq[,1] %<>% base::sub(pattern="^0 ([^ ]+)$", replacement="\\1", x=.)
  remove(p_seq)
  p_seq <- c("\\(([^()]+)\\)", "\\[([^]\\[]+)\\]")
  x$seq[,2] %<>% base::gsub(pattern=p_seq[1],  replacement="\\1", x=.) %>%
    base::gsub(pattern=p_seq[2], replacement="\\1", x=.) %>% stringr::str_squish(string=.)
  x$seq[,3] %<>% base::gsub(pattern=p_seq[1], replacement="\\1", x=.) %>%
    base::gsub(pattern=p_seq[2], replacement="\\1", x=.) %>% stringr::str_squish(string=.)
  remove(p_seq)
  # in old codes, where x_old[[i]] is x, you can run "x_old %<>% unique()"
  # check the relationship among x[[i]]-------------------------------------------------------------
  (c(x$seq[,-1]) %in% x$hap[,2]) %>% {if(!all(.)){stop("PHASE output is disrupted (chaotic)");}}
  (x$index[,1] != x$seq[,1]) %>% {if(any(.)){stop("PHASE output is disrupted (chaotic)");}}
  base::match(x=x$index[,2], table=x$hap[,1]) %>% x$hap[.,2] %>% {x$seq[,2] != .} %>%
    {if(any(.)){stop("PHASE output is disrupted (chaotic)");}}
  base::match(x=x$index[,3], table=x$hap[,1]) %>% x$hap[.,2] %>% {x$seq[,3] != .} %>%
    {if(any(.)){stop("PHASE output is disrupted (chaotic)");}}
  if("" == "check haplotype order"){
    ans <- as.integer(x$index[,2]) > as.integer(x$index[,3]); if(any(ans)){
      warning("in the haplotype composition part in phase results, haplotypes with small ",
              "index should be put before haplotypes represented by larger index")
      x$index[ans,2:3] <- x$index[ans,3:2]
    }; remove(ans);
  }
  x$hap[1,2] %>% base::strsplit(split=" ") %>% .[[1]] %>% {length(.) != length(x$position)} %>%
    {if(.){stop("PHASE output is disrupted (chaotic)");}}
  # remove x$seq, modify x$hap, change x$index from matrix into tibble------------------------------
  x %<>% .[c("position", "hap", "index")]
  if(any(as.character(1:nrow(x$hap)) != x$hap[,1])){stop("PHASE output is disrupted (chaotic)");}
  x$hap %<>% .[,2] %>% base::strsplit(split=" ") %>% simplify2array() %>% base::t(x=.)
  x$index %<>% as_tibble(.name_repair="minimal") %>%
    setNames(nm=c("individual", "primary", "secondary")) %>%
    dplyr::mutate(primary=as.integer(primary), secondary=as.integer(secondary))
  # if(!is.na(substitute_with_integers)), substitute integers with original letters-----------------
  if(!is.na(substitute_with_integers)){
    substitute_with_integers %<>% read.table(file=., header=TRUE, sep="\t", stringsAsFactors=FALSE,
                                             fileEncoding="UTF-8") %>% as_tibble()
    if(any(dim(substitute_with_integers) == 0)){stop("the file substitute_with_integers is disrupted (chaotic)");}
    (names(substitute_with_integers) != c("col_index", "col_name", "original", "integer")) %>%
      {if(any(.)){stop("the file substitute_with_integers is disrupted (chaotic)");}}
    #
    if(ncol(x$hap) != nrow(substitute_with_integers)){stop("x don't match with substitute_with_integers");}
    ans <- which(substitute_with_integers$original != ""); if(length(ans) > 0) for(i in ans){
      x$hap[,i] <- substitute_with_integers[i, 3:4] %>% unlist() %>% strsplit(split=",") %>%
        {.$original[base::match(x=x$hap[,i], table=.$integer)]}
    }; i <- 0; remove(ans, i);
    #
    x %<>% {c(.["position"], list(locus_name=substitute_with_integers$col_name), .["hap"], .["index"])}
  }
  remove(substitute_with_integers)
  # prepare a loci object---------------------------------------------------------------------------
  phased_loci <- mapply(FUN=function(xx, p, s){paste0(xx[p,], "|", xx[s,])}, p=x$index$primary,
                        s=x$index$secondary, xx=list(x$hap), SIMPLIFY=TRUE, USE.NAMES=FALSE) %>%
    base::t(x=.) %>% base::`colnames<-`(x=., value=x$locus_name) %>% # set rownames for a matrix would make each column a named vector
    as.data.frame(stringsAsFactors=TRUE) %>%
    base::`rownames<-`(x=., value=x$index$individual) %>% pegas::as.loci()
  # export x to .xml, return phased_loci------------------------------------------------------------
  cat("<root>\n", "<locus_position>", paste0(x$position, collapse=" "), "</locus_position>\n",
      "<locus_name>", paste0(x$locus_name, collapse=", "), "</locus_name>\n",
      "<haplotype>\n", file=output_file, sep="", append=FALSE)
  write.table(x=x$hap, file=output_file, append=TRUE, quote=TRUE, sep="\t",
              row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8")
  cat("</haplotype>\n", "<haplotypes_of_individual>\n", file=output_file, sep="", append=TRUE)
  suppressWarnings(expr={
    write.table(x=x$index, file=output_file, append=TRUE, quote=TRUE, sep="\t",
                row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")
  })
  cat("</haplotypes_of_individual>\n", "<phased_loci>\n", file=output_file, sep="", append=TRUE)
  suppressWarnings(expr={
    pegas::write.loci(x=phased_loci, file=output_file, append=TRUE, col.names=NA)
  })
  cat("</phased_loci>\n", "</root>\n", file=output_file, sep="", append=TRUE)
  remove(output_file)
  phased_loci %<>% list(phased_loci=., locus_position=x$position)
  remove(x)
  return(phased_loci)
}
