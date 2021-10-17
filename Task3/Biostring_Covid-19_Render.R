setwd("/home/ajo/gitRepos/bioinformatics/Task3")
mypath_in <- "./Fastas"
mypath_out <-"./Reports"

# Create a folder if is not exists
if(!dir.exists(mypath_out)) dir.create(mypath_out)


# Get the param values from the CSV-file.
parameters <- read.csv2("Listado_sec.csv", header = T)

for (i in 1:nrow(parameters)){
  # assemble the output file name
  row <- parameters[i,]
  outfile.name <- paste0("Biostring_Covid-19_",row[5])
  fn <-file.path(mypath_out, outfile.name)
  
  #Process                                    
  rmarkdown::render("/home/ajo/gitRepos/bioinformatics/Task3/Biostring_Covid-19_Tasca.Rmd", params = list(
                                                                 gene_seq=row[1]$sec_genica, 
                                                                 start_vec=eval(parse(text = row[2]$pos_start)), 
                                                                 end_vec=eval(parse(text = row[3]$pos_end)),
                                                                 protein_seq=row[4]$sec_prot,
                                                                 name_protein=row[5]$nombre.CDS,
                                                                 path_in=mypath_in,
                                                                 path_out=mypath_out
                                                                 ),
                    output_format =  c("html_document","pdf_document"),
                    output_file = c(paste0(fn, output_format ='.html'),
                                    paste0(fn, output_format ='.pdf')),
                    encoding = "UTF8")
  
}
