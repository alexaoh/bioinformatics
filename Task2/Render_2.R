render.report = function(path, fasta_name, output_filename, output_fileformat){
  rmarkdown::render(
    "/home/ajo/gitRepos/bioinformatics/Task2/SeqAna_2.Rmd", params = list(
      path = path, 
      name_fasta = fasta_name
    ),
    output_file = output_filename,
    output_format = if (output_fileformat == "all") "all" else paste(output_fileformat, "document", sep = "_") 
    # Can specify only one output format or "all" for both html and pdf. 
  )
}

render.report("/home/ajo/gitRepos/bioinformatics/Task2", "ENSG00000114374.fa", "SeqAna_2_ENSG00000114374", "pdf")
render.report("/home/ajo/gitRepos/bioinformatics/Task2", "ENSG00000114374.fa", "SeqAna_2_ENSG00000114374", "html")
