mut_counter <- function(full_orf_aln, 
                        genetic_code_path = './data/genetic_code.txt', 
                        sn_table_path = './data/SN_codonTable.txt') {
  ### This function take blah blah blah
  
  # get number of lines of the input file
  con <- full_orf_aln
  n_lines <- length(readLines(con))
  
  # Loading constant data tables
  SN_table <- 
    read.table(file = sn_table_path, header = TRUE, sep = '\t')
  
  genetic_code <- 
    read.table(file = genetic_code_path, header = TRUE, sep = '\t')

  # creating empty output table
  mut_table <- genetic_code
  mut_table$s <- 0
  mut_table$n <- 0
  
  # reading first sequence as a reference
  message("Analyzing the first sequence as a reference.")
  to_skip <- 1 # in order to read lines one by one
  ref <- 
    read_codons(file = full_orf_aln, n = 1, skip = to_skip)
  # counting codons in REF
  ref_count <- codon_counter(ref, sn_table = SN_table)
  
  # reading the following sequences
  message("Analyzing the following sequences")
  while(to_skip < (n_lines - 1)) {
  to_skip <- to_skip + 2
  cat('- sequence number:', to_skip %/% 2, '\n')
  curr_seq <- 
    read_codons(file = full_orf_aln, n = 1, skip = to_skip)
  # finding mismatches positions
  curr_mismatches <- 
    find_mismatches(seq_1 = ref, seq_2 = curr_seq)
  
  if (length(curr_mismatches) == 0) {
    next()
  }
  for (i in 1:length(curr_mismatches)) {
    # print(ref[curr_mismatches[i]])
    # print(curr_seq[curr_mismatches[i]])
    curr_mut <- syn_or_nonsyn(codon_ref = ref[curr_mismatches[i]], 
                  codon_mut = curr_seq[curr_mismatches[i]], 
                  gen_code = genetic_code)
    mut_table$s[mut_table$codon == curr_mut$ref[1]] <- 
      mut_table$s[mut_table$codon == curr_mut$ref[1]] + curr_mut$s[1]   
    mut_table$n[mut_table$codon == curr_mut$ref[1]] <- 
     mut_table$n[mut_table$codon == curr_mut$ref[1]] + curr_mut$n[1]              
  }
}
  mut_table <- merge(ref_count, mut_table, by = 'codon')
  full_output <- compute_kaKs(mutation_table = mut_table)
  order_wanted <- c('codon', 'AA','count', 'S', 'N', 's', 'n', 'sS', 'nN')
  full_output[[1]] <- full_output[[1]][ , order_wanted]
  cat('Ka/Ks is', round(full_output[[2]], 3),'\n')
  return(full_output)
}


n.readLines <- function(fn,n,comment="#", skip=0, header=FALSE)
{
  ### This is part of the reader package from Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
  # read at least 'n' lines of a file, skipping lines and ignoring any starting with comment
  if(!file.exists(fn)) { warning("file doesn't exist"); return(NULL) }
  if(!is.character(comment)) { warning("illegal comment char, reverting to #"); comment <- "#" }
  rl <- 0; cc <- 0 + {if(is.numeric(skip)) skip else 0 }
  while(rl<n) { 
    test.bit <- readLines(fn,n+cc)
    if(skip>0 & length(test.bit>1)) { test.bit <- test.bit[-(1:(min((length(test.bit)-1),skip)))] }
    cmnt <- which(substr(test.bit,1,1)==comment)
    rl <- n+cc-length(cmnt)
    cc <- cc + length(cmnt)
  }
  if(length(cmnt)>0) { test.bit <- test.bit[-cmnt] } 
  if(length(test.bit)>1 & header) { test.bit <- test.bit[-1] }
  return(test.bit)
}

read_codons <- function(file, n = 1,skip) {
  codons <- 
    n.readLines(fn = file, n = n, skip = skip, header = FALSE)
  codons <- toupper(codons)
  # if (grepl(pattern = '[^ATGC]', x = codons)) {
  #   stop('ARGH ! INVALID NUCLEOTIDES, WE ARE DOOOOOMED')
  # }
  codons <- 
    gsub(pattern = '(.{3})', replacement = '\\1;', x = codons)
  codons <-
    strsplit(x = codons, split = ';', fixed = TRUE)[[1]]
  return(codons)
}

codon_counter <- function(ref, sn_table = SN_table) {
  # count codon occurence
  codon_count_table <- data.frame(codon=unique(ref), count=0)
  for (i in 1:length(ref)) {
    codon.i <- ref[i]
    codon_count_table$count[match(codon.i, codon_count_table$codon)] <- 
      codon_count_table$count[match(codon.i, codon_count_table$codon)] + 1
  }
  codon_count_table$S <- 
    sn_table$S[match(codon_count_table$codon, sn_table$codon)]
  codon_count_table$N <- 
    sn_table$N[match(codon_count_table$codon, sn_table$codon)]
  return (codon_count_table)
}

find_mismatches <- function(seq_1, seq_2) {
  if (length(seq_1) != length(seq_2)) {
    stop('OHHHH NOOOOOO,
ONE OF THEM IS TOO LONG AND THE OTHER TOO SHORT, 
WHY G*D WHY HAVE YOU FORSAKEN US !')
  }
  mis <- which(seq_1 != seq_2)
  # print(mis)
  return(mis)
}

syn_or_nonsyn <- 
  function(codon_ref, codon_mut, gen_code = genetic_code) {
    in_codons <- c(codon_ref, codon_mut)
    result <- list(ref = codon_ref, mut = 'xxx', s = 0, n = 0)
    
    codon_1_invalid <- grepl(pattern = '[^ATGC]', x = in_codons[1])
    codon_2_invalid <- grepl(pattern = '[^ATGC]', x = in_codons[2])
    if (codon_1_invalid | codon_2_invalid) {
      return(result)
    }
    
    if (length(in_codons) != 2) {
      stop('ARGHHH it is not only two codons, we are all going to die !')
    }
    
    if (codon_ref == codon_mut) {
      stop('Biology 101 much ?')
    }
    out_aa <- 
      gen_code$AA[sapply(X = in_codons, FUN = function(x) grep(pattern = x, x = gen_code$codon))]
    
    
    if (out_aa[1] == out_aa[2]) {
      result <- list(ref = codon_ref, 
                     mut = codon_mut, 
                     aa_ref = out_aa[1], 
                     aa_mut = out_aa[2], 
                     s = 1, n = 0)
    } else if (out_aa[1] != out_aa[2] ) {
      result <- list(ref = codon_ref, 
                     mut = codon_mut, 
                     aa_ref = out_aa[1], 
                     aa_mut = out_aa[2], 
                     s = 0, n = 1)
    }
  return(result)
  }

compute_kaKs <- function(mutation_table) {
  # ATG/M and TTG/W are unique, hence the divide by 0
  mutation_table$sS <- mutation_table$s / (mutation_table$S * mutation_table$count)
  mutation_table$sS[is.na(mutation_table$sS)] <- 0
  mutation_table$nN <- mutation_table$n / mutation_table$N
  kaks <- sum(mutation_table$nN) / sum(mutation_table$sS * mutation_table$count)
  return(list(mutation_table, ka_ks = kaks))
}

#### test ----


system.time({
core_coli <-
  mut_counter("./data/core_gene_alignment_test.aln")
})

system.time({
  core_coli_consensus <-
    mut_counter("./data/core_gene_alignment_test_consensus.aln")
})


# 
# rm(glou_2)
# glou_2 <- mut_counter('data/rpoB_align_1l.aln')
# View(glou_2[[1]])
# 
# glou_3 <- mut_counter(full_orf_aln = 'data/tetC.align')
