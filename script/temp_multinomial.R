#### test ----

source('./script/kaks.R')

system.time({
  p91 <-
    kaks('./data/i_P91.aln', n_core = 16)
})

system.time({
  p92 <-
    kaks('./data/i_P92.aln', n_core = 16)
})

system.time({
  p22 <-
    kaks('./data/i_P22.aln', n_core = 16)
})

system.time({
  p2 <-
    kaks('./data/i_P2.aln', n_core = 33)
})

summary(p91[[2]]$KaKs)
ecdf91 <- ecdf(p91[[2]]$KaKs)
ecdf91(1)

summary(p92[[2]]$KaKs)
ecdf92 <- ecdf(p92[[2]]$KaKs)
ecdf92(1)

summary(p22[[2]]$KaKs)
ecdf22 <- ecdf(p22[[2]]$KaKs)
ecdf22(1)

summary(p2[[2]]$KaKs)
ecdf2 <- ecdf(p2[[2]]$KaKs)
ecdf2(1)


multi_simul <- function(df, n_draw = 1e1) {
  


codons <- 
  unique(df$codon)

df_out <- data.frame(codon = codons)
for (i in 1:length(codons)) {
  curr_codon <- codons[i]
  df_out$count[df_out$codon == curr_codon] <- sum(df$count[df$codon == curr_codon])
  df_out$S[df_out$codon == curr_codon] <- df$S[df$codon == curr_codon][1]
  df_out$N[df_out$codon == curr_codon] <- df$N[df$codon == curr_codon][1]
  df_out$s[df_out$codon == curr_codon] <- sum(df$s[df$codon == curr_codon])
  df_out$n[df_out$codon == curr_codon] <- sum(df$n[df$codon == curr_codon])
  df_out$Ks[df_out$codon == curr_codon] <- sum(df$Ks[df$codon == curr_codon])
  df_out$Ka[df_out$codon == curr_codon] <- sum(df$Ka[df$codon == curr_codon])
}
df_out$prop <- df_out$count / sum(df_out$count)

# Generate mutations randomly weighted by codon proportions
tot_s <- sum(df_out$s)
tot_n <- sum(df_out$n)
prop_S <- (df_out$count * df_out$S) / sum(df_out$count * df_out$S)
prop_N <- (df_out$count * df_out$N) / sum(df_out$count * df_out$N)

sim_kaks <- rep(NA, n_draw)

multi_s <- 
  rmultinom(n = n_draw, size = tot_s, prob = prop_S)

multi_n <- 
  rmultinom(n = n_draw, size = tot_n, prob = prop_N)


multi_Ks <- 
  apply(X = multi_s, MARGIN = 2, FUN = function(x) {
  x / (df_out$S)
})  

multi_Ka <- 
  apply(X = multi_n, MARGIN = 2, FUN = function(x) {
    x / (df_out$N)
  })

multi_KaKs <- 
  sapply(X = 1:n_draw, FUN = function(x) sum(multi_Ka[,x], na.rm = TRUE) / sum(multi_Ks[,x], na.rm = TRUE))
  
  return(multi_KaKs)
}

multi_simul(df = p92[[1]])
