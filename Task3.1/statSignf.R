# Create an R function to show the statistical significance of a pairwise alignment.

statSignf <- function(seq1, seq2, kind.seq, kind.align, subst.matrix, gap.scores, num.shuffles, shuffle.seq, plot){
  # Input parameters: 
  # seq1: Sequence 1 in Fasta format.
  # seq2: Sequence 2 in Fasta format.
  # kind.seq: Kind of sequence: Protein or DNA.
  # kind.align: Kind of pairwise alignment: “local” or “global”.
  # subst.matrix: Substitution matrix: PAMn, BLOSUMn, etc.
  # gap.scores: Gap score in vector format c(Open penalty, extended penalty).
  # num.shuffles: Number of suffles.
  # shuffle.seq: Sequence using for shuffling: 1 or 2.
  # plot: boolean value. plot = TRUE plots histogram and Gumbel. plot = FALSE does not plot. 
  
  # Output: 
  # Score of the original sequences. 
  # Histogram plot of scores obtained for the num.shuffles amount of permutations. 
  # p-value of the original score. 
  
  
  # Should probably have some form of asserts on the input parameters. 
  
  # Load the dependencies inside function for now, not sure if this should be done when making the package. 
  library(seqinr)
  library(Biostrings)
  library(evir) # Could use this to estimate a Gumbel, but will use the formulas given in the pdf in bioinfo-course. 
  data(subst.matrix)
  
  # I am assuming that the sequences already are in fasta format, such that I can use functions like "length" and "names".
  # If this is not the case, I can use a "readDNAStringSet" here, to fix this, such that the rest of the program should work (hopefully).
  
  m <- length(seq1)
  n <- length(seq2)
  
  # Should first calculate the score of the original sequences. 
  OG.score <- pairwiseAlignment(AAString(seq1), AAString(seq1), substitutionMatrix = subst.matrix, 
                                gapOpening = gap.scores[1], gapExtension = gap.scores[2], 
                                scoreOnly = TRUE)
  
  
  # This calculates scores with random permutations of one of the sequences, in order to analyse statistical significance. 
  scores <- rep(0, length.out = num.shuffles)
  for (i in 1:num.shuffles){ # change the loop to a replicate() later (or an apply(), sapply() etc).
    if (shuffle.seq = 1){
      x <- sample(seq1, m, replace = F)
    } else { # shuffle.seq must be 2 - should or could be asserted earlier. 
      x <- sample(seq2, n, replace = F)
    }
    
    # Need to add an option for local or global alignment here as well!
    scores[i] <- pairwiseAlignment(AAString(x), AAString(seq2),
                        substitutionMatrix = subst.matrix, 
                        gapOpening = gap.scores[1], gapExtension = gap.scores[2], 
                        scoreOnly = TRUE)
  }
  
  # Calculated mean and standard error of the scores. 
  score.mean <- mean(scores)
  score.se <- sd(scores)
  
  # Estimate the parameters of the Gumbel distribution. 
  lambda.hat <- 1.2825/score.se
  u.hat <- score.mean - 0.45*score.se
  
  # Estimate K.
  K.hat <- exp(lambda.hat * u.hat)/(m*n)
  
  # Standardize the score. 
  standard.score <- lambda.hat*OG.score - ln(K.hat*m*n)
  
  # p-value estimated from Gumbel.
  p.gumbel <- 1 - exp(-exp(-standard.score))
  
  # p-value counted directly from scores.
  p.counted <- sum(scores>OG.score)/num.shuffles
  
  gumbel.dist <- function(x){
    lambda.hat*exp(-lambda.hat*(x-u.hat)-exp(-lambda.hat*(x-u.hat)))
  }
  
  if (plot){
    x <- seq(from = min(scores)-1, to = max(scores) + 1, along.with = scores) # make x-axis for Gumbel. 
    hist(scores) # plot the histogram. 
    lines(x, gumbel.dist(x), col = 2) # plot the estimated Gumbel alongside it. 
    abline(x = OG.score, col = 3, lty = 2) # Add line for the original score on the plot. 
  }
  
  # Print the summary of the function (could be formatted more nicely later probably).
  cat("The names of the sequences are ", names(seq1), " and ", names(seq2), ". \n")
  cat("These sequences are of type ", kind.seq, ". \n")
  cat("The type of alignment that has been done here is ", kind.align, ". \n")
  cat("The substitution matrix used is ", subst.matrix, ". \n")
  cat("The gap scores are ", gap.scores[1], " for open penalty and ", gap.scores[2], " for extended penalty. \n")
  cat("The number of shuffles done are ", num.shuffles, ". \n")
  if (shuffle.seq = 1){
    cat("The shuffling was done on the first sequence mentioned above. \n")
  } else { # Just as earlier, the other option for shuffle.seq has to be 2!
    cat("The shuffling was done on the second sequence mentioned above. \n")
  }
  
  cat("The origian score is ", OG.score, ". \n")
  cat("The parameters of the Gumbel distribution are ", expression(hat(lambda)), "=", lambda.hat, 
            " and ", expression(hat(u)), "=", u.hat, ". \n")
  cat("The p-value estimated using the estimated Gumbel distribution is ", p.gumbel, ". \n")
  cat("The p-value estimated empirically by counting in the histogram is ", p.counted, ". \n")
  
  # Make a summary() function of the scores obtained by shuffling. 
  
  cat("The estimated K is ", expression(hat(K)), "=", K.hat, ". \n")
  cat("The standardized score is S' = ", standard.score, ". \n")
  
  # WHAT IS THE BEST MANNER IN TESTING THIS FUNCTION? HOW CAN I TEST IT (E.G WITH WHAT TYPE OF DATA)
  # MAYBE I SHOULD USE THE EXAMPLE WE SAW IN CLASS (USE THE SAME PROTEINS), IN ORDER TO SEE THAT I GET
  # THE SAME OUTPUT!
  
}
