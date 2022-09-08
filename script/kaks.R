kaks <- function(full_orf_aln,
                        genetic_code_path = './data/genetic_code.txt',
                        sn_table_path = './data/SN_codonTable.txt',
                        n_core = 1) 
  {
  ### This function take blah blah blah

  # Loading constant data tables
  library(foreach)
  library(doParallel)
  SN_table <- 
    read.table(file = sn_table_path, header = TRUE, sep = '\t')
  
  genetic_code <- 
    read.table(file = genetic_code_path, header = TRUE, sep = '\t')
  # getting sequence names and indexes
  all_comb <- comb_tab(full_orf_aln)
  # Creating the output list
  # output_list <-
  # vector(mode = 'list', length = nrow(all_comb))
  # my_clust <- makeCluster(spec = n_core, type = 'FORK')
  # registerDoParallel(cl = my_clust, cores = n_core)
  registerDoParallel(cores = n_core)
  # for (i in 1:nrow(all_comb)) {
  output_list <- 
  foreach(i = 1:nrow(all_comb)) %dopar% {
    seq_names <- 
      c(all_comb$seq_1[i], all_comb$seq_2[i])
    seq_index <- 
      c(all_comb$index_1[i], all_comb$index_2[i])
    # Checking current names and indexes    
  if (seq_names[1] == seq_names[2]) {
    stop('Identical sequences names !')
  }
  if (seq_index[1] == seq_index[2]) {
    stop('Identical sequences names !')
  }
  
  if (length(seq_names) < 2 || length(seq_index) < 2) {
    stop('Not enough sequence indexes or name')
  }
  
  if (length(seq_names) > 2 || length(seq_index) > 2) {
    message('Too many sequence indexes or name !')
    message('Only the first two will be used')
  }

  # creating empty output table
  mut_table <- genetic_code
  mut_table$seq_1 <- seq_names[1] # as.factor(seq_names[1])
  mut_table$seq_2 <- seq_names[2] # as.factor(seq_names[2])
  mut_table$s <- 0
  mut_table$n <- 0
  
  cat('Step', i, '/', nrow(all_comb),'\n')
  cat('Computing sequence', seq_names[1],' vs ', seq_names[2], '\n')
  
  # reading first sequence as a reference
  
  ref <- 
    read_codons(file = full_orf_aln, n = 1, to_skip = seq_index[1])
  # counting codons in REF
  ref_count <- codon_counter(ref, sn_table = SN_table)
  # reading the second sequences
  
  curr_seq <- 
    read_codons(file = full_orf_aln, n = 1, to_skip = seq_index[2])
  # finding mismatches positions
  curr_mismatches <- 
    find_mismatches(seq_1 = ref, seq_2 = curr_seq)
  
  if (length(curr_mismatches) == 0) {
    cat('The sequences', seq_names[1],' and ', seq_names[2], ' are identical !\n')
    return(NA)
  }
  
  for (j in 1:length(curr_mismatches)) {
    # print(ref[curr_mismatches[i]])
    # print(curr_seq[curr_mismatches[i]])
    curr_mut <- syn_or_nonsyn(codon_ref = ref[curr_mismatches[j]], 
                  codon_mut = curr_seq[curr_mismatches[j]], 
                  gen_code = genetic_code)
    mut_table$s[mut_table$codon == curr_mut$ref[1]] <- 
      mut_table$s[mut_table$codon == curr_mut$ref[1]] + curr_mut$s[1]   
    mut_table$n[mut_table$codon == curr_mut$ref[1]] <- 
     mut_table$n[mut_table$codon == curr_mut$ref[1]] + curr_mut$n[1]              
  }
  mut_table <- merge(ref_count, mut_table, by = 'codon')
  curr_output <- compute_kaKs(mutation_table = mut_table)
  order_wanted <- c('seq_1', 'seq_2','codon', 'AA','count', 'S', 'N', 's', 'n', 'Ks', 'Ka')
  curr_output <- curr_output[ , order_wanted]
  # cat('Ka/Ks is', round(curr_output[[2]], 3),'\n')
  cat('\n...\n\n')
  return(curr_output) # to get the result into output_list
  }
  
  output_df <- do.call(what = rbind, output_list)
  output_df <- output_df[complete.cases(output_df), ]
  kaks_df <- final_KaKs(df = output_df)
  
  return(list(output_df, kaks_df))
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

comb_tab <- function(full_orf_aln) {
  n_lines <- length(readLines(full_orf_aln))
  seq.names <- integer(0)
  for (i in seq(0,(n_lines-1),2)) {
    seq.names <- c(seq.names, n.readLines(full_orf_aln, n=1, skip=i))
  }
  seq.names <- gsub(">", "", seq.names) 		
  
  index.combn <- combn(seq(1,n_lines,2),2)
  seq.combn <- combn(seq.names,2)
  
  result <- data.frame(cbind(t(seq.combn), t(index.combn)))
  
  colnames(result) <- c("seq_1", "seq_2", "index_1", "index_2")
  result$index_1 <- as.integer(result$index_1)
  result$index_2 <- as.integer(result$index_2)
  return(result)
}

read_codons <- function(file, n = 1, to_skip) {
  codons <- 
    n.readLines(fn = file, n = n, skip = to_skip, header = FALSE)
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
  mutation_table$Ks <- mutation_table$s / (mutation_table$S * mutation_table$count)
  mutation_table$Ks[is.na(mutation_table$Ks)] <- 0
  mutation_table$Ka <- mutation_table$n / (mutation_table$N * mutation_table$count)
  if (sum(mutation_table$Ks) == 0) {
    cat('No synonymous mutation detected !\n')
  }
  # mutation_table$KaKs <- 
  #   sum(mutation_table$Ka) / sum(mutation_table$Ks)
  # return(list(mutation_table, ka_ks = kaks))
  return(mutation_table)
}

final_KaKs <- function(df) {
  df$id <- 
      paste0(df$seq_1, df$seq_2)
  uni_id <- 
    unique(df$id)
  n_uni <- length(uni_id)
  out_df <- data.frame(seq_1 = rep(NA, n_uni), 
                       seq_2 = rep(NA, n_uni), 
                       KaKs = rep(NA, n_uni))
  for (i in 1:n_uni) {
    curr_df <- 
      df[df$id == uni_id[i], ]
    
    out_df$seq_1[i] <- curr_df$seq_1[1]
    out_df$seq_2[i] <- curr_df$seq_2[1]
    out_df$KaKs[i] <-  sum(curr_df$Ka) / sum(curr_df$Ks)
  }
  return(out_df)
}
