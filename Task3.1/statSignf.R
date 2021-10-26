# Create an R function to show the statistical significance of a pairwise alignment.

statSignf <- function(seq1, seq2, kind.seq, kind.align, subst.matrix, gap.scores, num.shuffles, shuffle.seq, plot){
  # Input parameters: 
  # seq1: Sequence 1 in Fasta format
  # seq2: Sequence 2 in Fasta format
  # kind.seq: Kind of sequence: Protein or DNA
  # kind.align: Kind of pairwise alignment: “local” or “global”
  # subst.matrix: Substitution matrix: PAMn, BLOSUMn, …
  # gap.scores: Gap score in vector format c(Open penalty, extended penalty).
  # num.shuffles: Number of suffles: N
  # shuffle.seq: Sequence using for shuffling: 1 or 2
  # plot: boolean value. plot = TRUE plots histogram. plot = FALSE does not plot. 
  
  # Output: 
  # Score of the original sequences. 
  # Histogram plot of scores obtained for the num.shuffles amount of permutations. 
  # p-value of the 
  
  
  # Should probably have some form of asserts on the input parameters. 
  
  # Load the dependencies inside function for now, not sure if this should be done when making the package. 
  library(seqinr)
  library(Biostrings)
  data(subst.matrix)
  
  # Should first calculate the score of the original sequences. 
  OG.score <- pairwiseAlignment(AAString(seq1), AAString(seq1), substitutionMatrix = subst.matrix, 
                                gapOpening = gap.scores[1], gapExtension = gap.scores[2], 
                                scoreOnly = TRUE)
  
  
  # This calculates scores with random permutations of one of the sequences, in order to analyse statistical significance. 
  scores <- rep(0, length.out = num.shuffles)
  for (i in 1:num.shuffles){ # change the loop to a replicate() later (or an apply() or something fitting).
    if (shuffle.seq = 1){
      x <- sample(seq1, length(seq1), replace = F)
    } else { # shuffle.seq must be 2 - should or could be asserted earlier. 
      x <- sample(seq2, length(seq2), replace = F)
    }
    scores[i] <- pairwiseAlignment(AAString(x), AAString(seq2),
                        substitutionMatrix = subst.matrix, 
                        gapOpening = gap.scores[1], gapExtension = gap.scores[2], 
                        scoreOnly = TRUE)
  }
  
  if (plot){
    hist(scores)
  }
  
  # WHAT IS THE BEST MANNER IN TESTING THIS FUNCTION? HOW CAN I TEST IT (E.G WITH WHAT TYPE OF DATA)
  # MAYBE I SHOULD USE THE EXAMPLE WE SAW IN CLASS (USE THE SAME PROTEINS), IN ORDER TO SEE THAT I GET
  # THE SAME OUTPUT!
  
  # Next: Estimate a Gumbel distribution to the histogram and calculate p-value! 
  
}
