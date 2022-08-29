codon_counter <- function (full_ORF.aln) {
	
	# codon list
	g <- read.table("genetic_code.txt", header=TRUE)$codon
	
	# first sequence as the reference
	s <- readLines("core_gene_alignment_test.aln", 2) [2]
	s <- strsplit(s, "")
	
	# count codon occurence
	codon_count_table <- data.frame(codon=g, count=0)
	for (i in seq(1, length(s[[1]]), 3)) {
		codon.i <- toupper(paste(s[[1]][i:(i+2)], collapse=""))
		codon_count_table[match(codon.i, codon_count_table$codon), 2] <- codon_count_table[match(codon.i, codon_count_table$codon), 2] + 1
		print(i)
		}
	
	return (codon_count_table)
	}
	

# readchar