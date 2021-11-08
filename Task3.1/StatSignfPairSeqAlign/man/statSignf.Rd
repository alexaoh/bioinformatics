\name{statSignf}
\alias{statSignf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistical significance of a pairwise alignment of sequences.
}
\description{
Function to check statistical significance of a pairwise alignment of sequences. Does so by permutation of one of the sequences (without replacement), before re-alignment with the other sequence. Estimates p-values non-parametrically (by counting) and by fitting a Gumbel distribution to the obtained scores after permutation.
}
\usage{
statSignf(seq1, seq2, kind.seq = "Protein", kind.align = "global", subst.matrix = "BLOSUM62",
              gap.scores = c(3, 1), num.shuffles = 1000, shuffle.seq = 1, plot = TRUE)
}

\arguments{
  \item{seq1}{
    Sequence 1 in Fasta format (as a string). An example is shown in "Details".
  }
  \item{seq2}{
    Sequence 2 in Fasta format (as a string). An example is shown in "Details".
  }
  \item{kind.seq}{
    Kind of sequence: "Protein" or "DNA". Default is "Protein".
  }
  \item{kind.align}{
    Kind of pairwise alignment: “local” or “global”. Default is "global".
  }
  \item{subst.matrix}{
    Substitution matrix: PAMn, BLOSUMn, etc. Default is "BLOSUM62".
  }
  \item{gap.scores}{
    Gap scores in vector format c(Open penalty, extended penalty). Default is c(3,1).
  }
  \item{num.shuffles}{
    Number of shuffles for permutation. Default is 1000.
  }
  \item{shuffle.seq}{
    Sequence using for shuffling: 1 or 2. Default is 1.
  }
  \item{plot}{
    Boolean value. plot = TRUE plots histogram and Gumbel. plot = FALSE does not plot. Default is TRUE.
  }

}
\details{
  It is assumed that the function takes the sequences in the type of format shown in "Examples" (important with a newline). The name of the   sequence should be specified on the first line, with the sequence added in the next line.
}
\value{
The function returns a vast summary of the input variables used when executing the function. It returns the names of the sequences aligned, the kind of sequences, the type of alignment done, the substitution matrix used, the gap scores, the number of permutations, the score of the original alignment of the two sequences and which of the two sequences that has been permuted. Moreover, it returns the estimates of the parameters for the Gumbel distribution, the p-value from the sequence alignment parametrically and non-parametrically, the estimated K, the standardized score and a summary table of the permuted scores (minimum, maximum, mean, median, 0.25 quantile and 0.75 quantile). Lastly, if plot = TRUE, a histogram of the scores is plotted together with the fitted Gumbel distribution function and the original score. Some of the most important information is repeated as text in the plot.
}
\references{
Course of bioinformatics at UPC.
}
\author{
Alex-Alex-Helena
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Two random sequences in the specified format.
seq1 <- ">random sequence 1 consisting of 20 residues.
KMMIDIHWGMWWYEYMMCLD"
seq2 <- ">random sequence 1 consisting of 20 residues.
DVYRVCQNVFRYHHFCKRTI"

# Simple alignment.
statSignf(seq1, seq2)

# More specific alignment.
statSignf(seq1, seq2, num.shuffles = 5000,
      subst.matrix = "PAM30", kind.align = "local")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
