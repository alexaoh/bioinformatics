# Create an R function to show the statistical significance of a pairwise alignment.

statSignf <- function(seq1, seq2, kind.seq = "Protein", kind.align = "global", subst.matrix = "BLOSUM62",
                      gap.scores = c(3, 1), num.shuffles = 1000, shuffle.seq = 1, plot = TRUE){

  # Should probably have some form of asserts on the input parameters.

  # Load the dependencies inside function for now.
  library(seqinr)
  library(Biostrings)
  data(list = subst.matrix)

  # Split the fastas into one name-part and one sequence-part
  string1.split <- strsplit(toString(AAStringSet(seq1)), "\n")[[1]]
  string2.split <- strsplit(toString(AAStringSet(seq2)), "\n")[[1]]

  # Get the name-part of the Fasta-string-inputs.
  name1 <- string1.split[1]
  name2 <- string2.split[1]

  # Get the sequence part of the Fasta-string-inputs.
  s1 <- string1.split[2]
  s2 <- string2.split[2]

  # Gets the length of the strings.
  m <- nchar(s1)
  n <- nchar(s2)

  if (kind.align == "global"){
    OG.score <- pairwiseAlignment(s1, s2, substitutionMatrix = subst.matrix,
                                  gapOpening = gap.scores[1], gapExtension = gap.scores[2],
                                  scoreOnly = TRUE, type = "global")
  } else {
    OG.score <- pairwiseAlignment(s1, s2, substitutionMatrix = subst.matrix,
                                  gapOpening = gap.scores[1], gapExtension = gap.scores[2],
                                  scoreOnly = TRUE, type = "local")
  }

  # This calculates scores with random permutations of one of the sequences, in order to analyse statistical significance.
  scores <- rep(0, length.out = num.shuffles)
  s1.splitted <- strsplit(s1,split="")[[1]] # Makes sampling possible.
  s2.splitted <- strsplit(s2,split="")[[1]] # Makes sampling possible.

  for (i in 1:num.shuffles){ # should make this more efficient eventually, since it is very slow! Not now though.
    if (shuffle.seq == 1){
      x <- paste(sample(s1.splitted, m, replace = F), collapse = "")
      x2 <- s2
    } else { # shuffle.seq must be 2 - should or could be asserted earlier.
      x <- paste(sample(s2.splitted, n, replace = F), collapse = "")
      x2 <- s1
    }

    if (kind.align == "global"){
      scores[i] <- pairwiseAlignment(x, x2,
                                     substitutionMatrix = subst.matrix,
                                     gapOpening = gap.scores[1], gapExtension = gap.scores[2],
                                     scoreOnly = TRUE, type = "global")
    } else {
      scores[i] <- pairwiseAlignment(x, x2,
                                     substitutionMatrix = subst.matrix,
                                     gapOpening = gap.scores[1], gapExtension = gap.scores[2],
                                     scoreOnly = TRUE, type = "local")
    }

  }

  # Calculated mean and standard error of the scores.
  score.mean <- mean(scores)
  score.se <- sd(scores)

  # More summary statistics from 'scores'.
  score.min <- min(scores)
  score.max <- max(scores)
  score.median <- median(scores)
  score.1st.quant <- quantile(scores, 0.25)
  score.3rd.quant <- quantile(scores, 0.75)

  # Estimate the parameters of the Gumbel distribution.
  lambda.hat <- 1.2825/score.se
  u.hat <- score.mean - 0.45*score.se

  # Estimate K.
  K.hat <- exp(lambda.hat * u.hat)/(m*n)

  # Standardize the score.
  standard.score <- lambda.hat*OG.score - log(K.hat*m*n)

  # p-value estimated from Gumbel.
  p.gumbel <- 1 - exp(-exp(-standard.score))

  # p-value counted directly from scores.
  p.counted <- sum(scores>OG.score)/num.shuffles

  gumbel.dist <- function(x){
    lambda.hat*exp(-lambda.hat*(x-u.hat)-exp(-lambda.hat*(x-u.hat)))
  }

  if (plot){
    par(mar = c(4, 4, 6, 4))
    x <- seq(from = min(scores)-1, to = max(scores) + 1, along.with = scores) # make x-axis for Gumbel.
    gumbel <- gumbel.dist(x) # Calculated Gumbel curve.
    hist(scores, freq = F, ylim = c(0, max(gumbel))) # plot the histogram.
    lines(x, gumbel, col = "blue") # plot the estimated Gumbel alongside it.
    abline(v = OG.score, col = 2, lty = 2) # Add line for the original score on the plot.
    legend("topleft", legend = c("Gumbel", "Original Score"), col = c("blue", 2), lty = 1:2)
    # mtext(paste0("Names: ", substr(name1, start = 2, stop = nchar(name1)-1)
    #              , " and ", substr(name2, start = 2, stop = nchar(name2)-1)), side = 3)
    mtext(paste0("Gap Scores: -", gap.scores[1]," for open penalty and -", gap.scores[2], " for extended penalty. \n"))
    #mtext(paste0(num.shuffles, " shuffles. \n"), side = 4, adj = 1, padj = 0)
    #mtext(paste0("Original Score: ", OG.score, "\n"))
    mtext(paste0("Parametric p-value ", round(p.gumbel, 3), " and counted p-value: ", round(p.counted, 3)))
  }

  # Print the summary of the function (could be formatted more nicely later probably).
  cat("The names of the sequences are '", substr(name1, start = 2, stop = nchar(name1)-1), "' and '",
      substr(name2, start = 2, stop = nchar(name2)-1), "'. \n", sep = "")
  cat("These sequences are of type '", kind.seq, "'. \n", sep = "")
  cat("The type of alignment that has been done here is '", kind.align, "'. \n", sep = "")
  cat("The substitution matrix used is '", subst.matrix, "'. \n", sep = "")
  cat("The gap scores are -", gap.scores[1], " for open penalty and -", gap.scores[2], " for extended penalty. \n", sep = "")
  cat("The number of shuffles done are ", num.shuffles, ". \n", sep = "")
  if (shuffle.seq == 1){
    cat("The shuffling was done on the first sequence mentioned above. \n")
  } else { # Just as earlier, the other option for shuffle.seq has to be 2!
    cat("The shuffling was done on the second sequence mentioned above. \n")
  }

  cat("The original score is ", OG.score, ". \n", sep = "")
  cat("The estimations of the parameters of the Gumbel distribution are: Scale parameter lambda = ", lambda.hat,
      " and mode u = ", u.hat, ". \n", sep = "")
  cat("The p-value estimated using the estimated Gumbel distribution is ", p.gumbel, ". \n", sep = "")
  cat("The p-value estimated empirically by counting in the histogram is ", p.counted, ". \n", sep = "")

  cat("The estimated K is ", K.hat, ". \n", sep = "")
  cat("The standardized score is S' = ", standard.score, ". \n", sep = "")

  # Summary.
  cat("Summary of the", num.shuffles, "scores calculated after shuffling: \n")
  cat("Min.\t1st Qu.\tMedian\tMean\t3rd Qu.\tMax.\n")
  cat(paste(score.min, score.1st.quant, score.median, score.mean, score.3rd.quant, score.max,sep="\t"))
  cat("\n")

}
