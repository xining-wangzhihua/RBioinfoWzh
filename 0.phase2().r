#! /usr/bin/Rscript
# license: BSD 3-Clause License
# author: https://github.com/ywd5
# version: 20210314

# this function is a wrapper of the following components:
# PrepareUnphasedLoci(), WritePhaseInput(), RunPhase(), remove extra phase output,
# ReadPhaseOutput(), ExamineRarePhasedHaplotypes(),
# convert haplotype sequences from .tsv to .fasta

# PHASE only accept diploid data, so do PrepareUnphasedLoci(), WritePhaseInput(),
# InterpretPhaseOutput() and ExamineRarePhasedHaplotypes()

# this function only accept DNAStringSet as input (to use loci, omit PrepareUnphasedLoci() manually)

message("When estimating haplotypes, remember to check if the DNA segement is repeated in genome. ",
        "If it's repeated, it's better to distinguish among paralogues (antonym of orthologue).")

warning("this functions isn't tested as a whole from head to tail")

if("" == "get dependency"){
  require(magrittr); require(tibble); require(dplyr); require(stringr); require(purrr);
  require(xml2)
  require(Biostrings); require(pegas);
  #
  stop("get dependency functions manually")
}

if("" == "test codes"){
  x <- readDNAStringSet(filepath="../test/its ps.fasta")
  smaller_for_faster <- c(0.005, 1)
}

phase2 <- function(x, suffix_pattern="\\.primary$|\\.secondary$|\\.hap[[:digit:]]+$",
                   known=character(), smaller_for_faster=c(5, 200)){
  # the arguement "known" represents "segments whose phases are known"
  # control the input-------------------------------------------------------------------------------
  # inputs are pooly controled here, and are strictly controled in dispatched functions
  if(!is(object=x, class2="DNAStringSet")){stop("x must be a DNAStringSet object");}
  if(!is.character(suffix_pattern)){stop("suffix_pattern must be a string");}
  if(!is.character(known)){stop("known must be a character vector");}
  if(!is.numeric(smaller_for_faster)){stop("smaller_for_faster must be 2 numbers");}
  # call PrepareUnphasedLoci(), remove(x, suffix_pattern)-------------------------------------------
  unphased_loci <- PrepareUnphasedLoci(x=x, suffix_pattern=suffix_pattern)
  remove(x, suffix_pattern)
  # call WritePhaseInput(), remove(known)-----------------------------------------------------------
  output_directory <- WritePhaseInput(x=unphased_loci$unphased_loci,
                                      locus_position=unphased_loci$locus_position, known=knwon)
  setwd(dir=output_directory)
  cat("<root>\n", "<unphased_loci>\n", file="./unphased_loci.xml", append=FALSE, sep="")
  suppressWarnings(expr={
    pegas::write.loci(x=unphased_loci$unphased_loci, file="./unphased_loci.xml", append=TRUE, col.names=NA)
  })
  cat("</unphased_loci>\n",
      "<locus_position>", paste0(unphased_loci$locus_position, collapse=" "), "</locus_position>\n",
      "<not_varied_loci>\n", paste0(unphased_loci$not_varied_loci, collapse="\n"), "\n</not_varied_loci>\n",
      "</root>\n", file="./unphased_loci.xml", append=TRUE, sep="")
  remove(known)
  # call RunPhase(), remove(smaller_for_faster)-----------------------------------------------------
  phased_loci <- RunPhase(smaller_for_faster=smaller_for_faster)
  remove(smaller_for_faster)
  list.files(path=".", pattern="\\.nexus_[[:alpha:]]+$", full.names=TRUE, recursive=TRUE) %>%
    {if(length(.) > 0) if(!all(file.remove(.))){warning("error when removing extra PHASE outputs");}}
  # call InterpretPhaseOutput()---------------------------------------------------------------------
  phased_loci %<>% lapply(FUN=InterpretPhaseOutput) %>%
    lapply(FUN=function(xx){xx$phased_loci;}) %>%
    tibble(phased_loci=., old_directory=base::dirname(phased_loci)) %>%
    {.[match(x=base::unique(.$phased_loci), table=.$phased_loci),]}
  # (currently omitted) ensure original and phased loci are consistent------------------------------
  # original phases may contain "?" or "-1", which represent missing values for PHASE to guess,
  # so setdiff(origianl_loci, phased_loci) at each locus each individual may contain "?" or "-1"
  # message("\nPHASE has guessed some missing values\n")
  # call ExamineRarePhasedHaplotypes()--------------------------------------------------------------
  phased_loci %<>% .$phased_loci %>%
    mapply(FUN=ExamineRarePhasedHaplotypes, phased_loci=.,
           locus_position=list(unphased_loci$locus_position),
           not_varied_loci=list(unphased_loci$not_varied_loci)) %>%
    dplyr::mutate(.data=phased_loci, new_directory=.)
  mapply(FUN=file.copy, from=lapply(X=phased_loci$old_directory, FUN=list.files, full.names=TRUE, no..=TRUE),
         to=phased_loci$new_directory, SIMPLIFY=TRUE, USE.NAMES=FALSE) %>%
    {if(!all(.)){stop("can't copy output.nexus and output.xml");}}
  phased_loci %<>% .$new_directory
  remove(unphased_loci)
  ans <- tempfile(); for(i in 1:length(phased_loci)){
    # get haplotype sequences from info.xml
    paste0(phased_loci[i], "/info.xml") %>% xml2::read_xml() %>%
      xml2::xml_child(x=., search="haplotype_full") %>% xml2::xml_contents(x=.) %>%
      as.character() %>% cat(., file=ans, append=FALSE, sep="")
    read.table(file=ans, header=FALSE, sep="\t", row.names=1L, colClasses="character") %>%
      as.matrix() %>% apply(MARGIN=1, FUN=paste0, collapse="") %>% DNAStringSet() %>%
      writeXStringSet(x=., filepath=paste0(phased_loci[i], "/haplotype_full.fasta"), width=100)
  }; file.remove(ans); remove(ans, i);
  phased_loci %>% paste0("/haplotype_full.fasta") %>% lapply(FUN=readDNAStringSet) %>%
    lapply(FUN=uniqueLetters) %>% unlist() %>% base::unique() %>%
    {. %in% c("A", "C", "G", "T", "-")} %>% {!all(.)} %>%
    {if(.){warning("Some phased sequences have ambiguous genetic codes. This may be ",
                   "attributed (ascribed, imputed) to ape::seg.sites(), which treat only ",
                   "the 2nd but not 4th site in AACST / ACCST as segregating sites.");}}
  remove(phased_loci)
  # return result-----------------------------------------------------------------------------------
  setwd(dir="../") # set back after some operations
  message("\noutput files are written to: ", getwd(), "\n")
  return(output_directory)
}
