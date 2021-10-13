render.report = function(path, fasta_name, output_filename, output_fileformat){
  rmarkdown::render(
    "/home/ajo/gitRepos/bioinformatics/Task1/SeqAna.Rmd", params = list(
      path = path, 
      name_fasta = fasta_name
    ),
    output_file = output_filename,
    output_format = if (output_fileformat == "all") "all" else paste(output_fileformat, "document", sep = "_") 
    # Can specify only one output format or "all" for both html and pdf. Does not give the correct names to both pdf and html now, not sure why. 
  )
}

render.report("/home/ajo/gitRepos/bioinformatics/Task1", "sequence3.fa", "SeqAna3", "pdf")
render.report("/home/ajo/gitRepos/bioinformatics/Task1", "sequence3.fa", "SeqAna3", "html")

render.report("/home/ajo/gitRepos/bioinformatics/Task1", "sequence4.fa", "SeqAna4", "pdf")
render.report("/home/ajo/gitRepos/bioinformatics/Task1", "sequence4.fa", "SeqAna4", "html")

render.report("/home/ajo/gitRepos/bioinformatics/Task1", "gi32141095_N_1.fa", "SeqAna1", "pdf")
render.report("/home/ajo/gitRepos/bioinformatics/Task1", "gi32141095_N_1.fa", "SeqAna1", "html")

render.report("/home/ajo/gitRepos/bioinformatics/Task1", "gi32141095_N_0.fa", "SeqAna0", "pdf")
render.report("/home/ajo/gitRepos/bioinformatics/Task1", "gi32141095_N_0.fa", "SeqAna0", "html")