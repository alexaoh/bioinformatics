render.report = function(path, fasta_name){
  rmarkdown::render(
    "/home/ajo/gitRepos/bioinformatics/Task1/SeqAna.Rmd", params = list(
      path = path, 
      name_fasta = fasta_name
    )
  )
}
# Should put a output filename and path in the Rmd also, such that the html does not get overwritten each time.
# Should also make a parameter for html or pdf. 
render.report("/home/ajo/gitRepos/bioinformatics/Task1", "sequence3.fa")
render.report("/home/ajo/gitRepos/bioinformatics/Task1", "sequence4.fa")
