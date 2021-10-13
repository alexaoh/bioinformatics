render.report = function(path, fasta_name, output_filename, output_fileformat){
  rmarkdown::render(
    "/home/ajo/gitRepos/bioinformatics/Task2/SeqAna_2.Rmd", params = list(
      path = path, 
      name_fasta = fasta_name
    ),
    output_file = paste(output_filename,format(Sys.Date(),'%Y-%m-%d'), sep = "_"),
    output_format = if (output_fileformat == "all") "all" else paste(output_fileformat, "document", sep = "_") 
    # "all" does not give the correct names to both pdf and html here either. 
  )
}

render.report("/home/ajo/gitRepos/bioinformatics/Task2/fastas", "ENSG00000114374.fa", "ENSG00000114374", "html")
render.report("/home/ajo/gitRepos/bioinformatics/Task2/fastas", "ENSG00000114374.fa", "ENSG00000114374", "pdf")

render.report("/home/ajo/gitRepos/bioinformatics/Task2/fastas", "ENSG00000169800.fa", "ENSG00000169800", "html")
render.report("/home/ajo/gitRepos/bioinformatics/Task2/fastas", "ENSG00000169800.fa", "ENSG00000169800", "pdf")

render.report("/home/ajo/gitRepos/bioinformatics/Task2/fastas", "ENSG00000169953.fa", "ENSG00000169953", "html")
render.report("/home/ajo/gitRepos/bioinformatics/Task2/fastas", "ENSG00000169953.fa", "ENSG00000169953", "pdf")
