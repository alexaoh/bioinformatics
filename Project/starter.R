my.knit = knitr::knit("main.Rnw")
## document.tex is the latex file that will be compiled by the two following command:

system(paste0("pdflatex ", "main.tex")) 
system(paste0("pdflatex ", "main.tex")) 
