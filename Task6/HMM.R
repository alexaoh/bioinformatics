# Tasca 6 en bioinfo. 

generate.sequences <- function(Tr, Em, n, p){
    # Tr = transition matrix.
    # Em = emission matrix.
    # n = length of generated sequence.
    # p = initial state probabilities. 
    
  symbols <- c("A", "C", "G", "T")
  states <- c(1, 2) # 1 = GC and 2 = AT
  
  # x are hidden states and y are observed symbols. 
  x <- y <- rep(NA, length.out = n)
  
  x[1] <- sample(states, size = 1, prob = p)
  y[1] <- sample(symbols,size = 1,replace=T,prob = Em[,x[1]])
  
  for(i in 1:(n-1)){
    x[i+1] <- sample(states, size = 1, prob = Tr[x[i],])
    y[i+1] <- sample(symbols, size = 1, prob = Em[, x[i]])
  }
  return(cbind("Hidden" = x, "Obs" = y))
  
}

# Testing when developing.
Tr <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)
rownames(Tr) <- c("AT", "CG")
Em <- matrix(c(0.1, 0.4, 0.4, 0.1, 0.4, 0.1, 0.1, 0.4), nrow = 4)
rownames(Em) <- c("A", "C", "G", "T")
colnames(Em) <- c("AT", "CG")
p <- c(0.5,0.5)
n <- 10000

df <- generate.sequences(Tr, Em, n, p)

# Show sequence of states and nucleotides together:
df[1:50, ]

# Make graph. THIS IS CODE FROM LAST WEEK, I THINK WE CAN USE THE SAME CODE TO MAKE THE GRAPHS THIS WEEK!
make.for.graphing <- function(sequence.string){
  df.seq1 <- rep(NA, length = 15*25-8)
  seq1 <- concatenated
  seq1.vec <- strsplit(seq1, "")[[1]]
  for (i in 1:((15*25)-7)){
    prods <- rep(NA, length = 8)
    for(j in 1:8){
      prods[j] <- rel.freq[seq1.vec[j+i-1], j] 
    }
    df.seq1[i] <- prod(prods)
  }
  return(log(df.seq1))
}
