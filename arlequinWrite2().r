#this is a stable version
#last edited at 20200617

#similar to strataG::arlequinWrite()

#get dependency
if(FALSE){
  require(magrittr); require(stringr);
}

arlequinWrite2 <- function(haplotype_sequences, info, file,
                           sorted_population_names=NA, sorted_haplotype_names=NA,
                           ignore_individual=FALSE, with_structure=TRUE){
  #control the input of haplotype_sequences-------------------------------------
  if(!is(object=haplotype_sequences, class2="DNAStringSet")){
    stop("haplotype_sequences must be a DNAStringSet")
  }
  if(length(haplotype_sequences) == 0){stop("length(haplotype_sequences) can't be 0");}
  if(length(unique(nchar(haplotype_sequences))) != 1){stop("all sequences length must equal");}
  if(nchar(haplotype_sequences[[1]]) == 0){stop("sequences must contain at least 1 residue");}
  #
  if(is.null(names(haplotype_sequences))){stop("names(haplotype_sequences) must be specified");}
  if(anyNA(names(haplotype_sequences))){stop("names(haplotype_sequences) must be specified");}
  names(haplotype_sequences) %<>% stringr::str_squish()
  if(any(names(haplotype_sequences) == "")){stop("names(haplotype_sequences) must be specified");}
  if(anyDuplicated(names(haplotype_sequences))){stop("names(haplotype_sequences) can't contain duplicated names")}
  #control the input of info----------------------------------------------------
  if(!is.list(info)){stop("info must be a list");}
  if(!identical(names(info), c("main", "genotype_component"))){
    stop("the 3 elements of info must be main, genotype_component respectively")
  }
  if(!identical(names(info$main), c("population", "individual", "genotype"))){
    stop("names(info$main) must be population, individual, genotype respectively")
  }
  #
  if(length(info$genotype_component) < 2L){stop("genotype must contain at least 1 allele");}
  if(names(info$genotype_component)[1] != "genotype"){
    stop("names(info$genotype_component)[1] must be genotype")
  }
  if(!all( grepl(pattern="^allele", x=names(info$genotype_component)[-1]) )){
    stop("names(info$genotype_component)[-1] must all begin with allele")
  }
  if(anyDuplicated(info$genotype_component$genotype)){
    stop("info$genotype_component$genotype can't contain duplicated items")
  }
  message("\nask the author to resive codes. other control of parameters are needed\n")
  #control the input of file----------------------------------------------------
  if(!is.character(file)){stop("file must be a character string");}
  if(length(file) != 1){stop("file must be a character string");}
  file %<>% base::normalizePath(winslash="/", mustWork=FALSE)
  #control the input of with_structure------------------------------------------
  if((!identical(with_structure, TRUE)) & (!identical(with_structure, FALSE))){
    stop("with_structure must be TRUE or FALSE")
  }
  #control the input of ignore_individual------------------------------------------
  if((!identical(ignore_individual, TRUE)) & (!identical(ignore_individual, FALSE))){
    stop("ignore_individual must be TRUE or FALSE")
  }
  #control sorted_haplotype_names----------------------------------------------
  if(length(sorted_haplotype_names) == 1) if(is.na(sorted_haplotype_names)){
    sorted_haplotype_names <- stringr::str_sort(x=names(haplotype_sequences), numeric=TRUE)
  }
  if(length(haplotype_sequences) != length(sorted_haplotype_names)){
    stop()
  }
  if(!setequal(x=names(haplotype_sequences), y=sorted_haplotype_names)){
    stop()
  }
  #control sorted_population_names----------------------------------------------
  if(length(sorted_population_names) == 1) if(is.na(sorted_population_names)){
    sorted_population_names <- info$main$population %>% base::unique() %>% str_sort(numeric=TRUE)
  }
  if(length(unique(info$main$population)) != length(sorted_population_names)){
    stop()
  }
  if(!setequal(x=info$main$population, y=sorted_population_names)){
    stop()
  }
  #control the relationship among parameters------------------------------------
  if(!setequal(x=info$main$genotype, y=info$genotype_component$genotype)){
    stop()
  }
  if(!setequal(x=names(haplotype_sequences), y=unname(unlist(info$genotype_component[-1])))){
    stop()
  }
  #rename parameters------------------------------------------------------------
  hs <- haplotype_sequences; remove(haplotype_sequences);
  file_name <- file; remove(file);
  ohn <- sorted_haplotype_names; remove(sorted_haplotype_names);
  opn <- sorted_population_names; remove(sorted_population_names);
  #modify parameters------------------------------------------------------------
  info <- base::merge(x=info$main, y=info$genotype_component, by="genotype") %>%
    .[c("population", "individual", "genotype", names(info$genotype_component)[-1])] %>%
    {.[str_order(.$individual, numeric=TRUE),];} %>% {split(x=., f=.$population);} %>% .[opn]
  for(i in 1:length(info)){rownames(info[[i]]) <- NULL;}
  hs <- hs[ohn]
  #get general information------------------------------------------------------
  l <- list(ploidy=length(info[[1]]) - 3, npop=length(opn))
  connection_name <- file(description=file_name, open="wt", blocking=FALSE)
  #export to file---------------------------------------------------------------
  cat("[Profile]\n",
      "  Title=\"created by r at ", strftime(x=Sys.time(), format="%Y%m%d"), "\"\n",
      "  NbSamples=", l$npop, "\n",
      "  DataType=DNA\n",
      "  GenotypicData=", if(ignore_individual){"0"}else{"1"}, "\n",
      if(ignore_individual){""}else{"  GameticPhase=1\n"},
      "  LocusSeparator=WHITESPACE\n",
      "  MissingData=\'?\'\n",
      "  Frequency=ABS\n",
      "\n",
      "[Data]\n",
      "\n",
      file=connection_name, append=FALSE, sep="")
  #
  if(!ignore_individual){
    cat("  [[Samples]]\n", file=connection_name, append=TRUE, sep="")
    individual_nchar <- lapply(X=info, FUN=function(datum){nchar(datum$individual);}) %>%
      unlist() %>% max() %>% {. + 2;}
    for(i in 1:l$npop){
      info[[i]]$individual %<>% {str_pad(string=., width=individual_nchar, side="right");}
    }
  }
  if(!ignore_individual) for(i in 1:l$npop){
    cat(strrep(" ", times=4), "SampleName=\"", names(info[i]), "\"\n",
        strrep(" ", times=4), "SampleSize=", nrow(info[[i]]), "\n",
        strrep(" ", times=4), "SampleData={\n",
        file=connection_name, append=TRUE, sep="")
    for(j in 1:nrow(info[[i]])){
      cat(strrep(" ", times=6), info[[i]]$individual[j], "1  ", as.character(hs[ info[[i]][j,3+1] ]), "\n",
          file=connection_name, append=TRUE, sep="")
      if(l$ploidy > 1) for(k in 2:l$ploidy){
        cat(strrep(" ", times=6 + individual_nchar + 3), as.character(hs[ info[[i]][j,3+k] ]), "\n",
            file=connection_name, append=TRUE, sep="")
      }
    }
    cat(strrep(" ", times=4), "}\n", file=connection_name, append=TRUE, sep="")
  }
  if(!ignore_individual){k <- 0L; remove(individual_nchar, j, k);}
  #
  if(ignore_individual){
    cat("  [[HaplotypeDefinition]]\n",
        strrep(" ", times=4), "HaplListName=\"haplotype sequences\"\n",
        strrep(" ", times=4), "HaplList={\n",
        file=connection_name, append=TRUE, sep="")
    ans <- as.character(hs)
    names(ans) %<>% {str_pad(string=., width=max(nchar(.)) + 2, side="right");}
    for(i in 1:length(ans)){
      cat(strrep(" ", times=6), names(ans[i]), unname(ans[i]), "\n", file=connection_name, append=TRUE, sep="")
    }; remove(i);
    remove(ans)
    cat(strrep(" ", times=4), "}\n",
        "\n",
        "  [[Samples]]\n",
        file=connection_name, append=TRUE, sep="")
    haplotype_nchar <- lapply(X=info, FUN=function(datum){nchar(unlist(datum[4:length(datum)]));}) %>%
      unlist() %>% max() %>% {. + 2L;}
    for(i in 1:l$npop) for(j in 1:l$ploidy){
      info[[i]] [[j + 3]] %<>% {str_pad(string=., width=haplotype_nchar, side="right");}
    }; remove(j);
    remove(haplotype_nchar)
  }
  if(ignore_individual) for(i in 1:l$npop){
    info[[i]] %<>% .[4:(3+l$ploidy)] %>% unlist() %>% unname() %>% table() %>% c() %>%
      {.[str_order(names(.), numeric=TRUE)];}
    cat(strrep(" ", times=4), "SampleName=\"", names(info[i]), "\"\n",
        strrep(" ", times=4), "SampleSize=", sum(info[[i]]), "\n",
        strrep(" ", times=4), "SampleData={\n",
        paste0(paste0(strrep(" ", times=6), names(info[[i]]), unname(info[[i]]), "\n"), collapse=""),
        strrep(" ", times=4), "}\n",
        file=connection_name, append=TRUE, sep="")
  }
  #
  if(with_structure){
    cat("\n","  [[Structure]]\n",
        strrep(" ", times=4), "StructureName=\"a group of ", l$npop, " populations\"\n",
        strrep(" ", times=4), "NbGroups=1\n",
        strrep(" ", times=4), "Group={\n",
        paste0(strrep(" ", times=6), "\"", opn, "\"\n"),
        strrep(" ", times=4), "}\n",
        file=connection_name, append=TRUE, sep="")
  }
  base::close(con=connection_name)
  remove(hs, info, file_name, with_structure, ignore_individual, ohn, opn, connection_name, l, i);
  invisible(NULL)
}
