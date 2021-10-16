setwd("/home/ajo/gitRepos/bioinformatics/Task3")
mypath_in <- "./Fastas"
mypath_out <-"./Reports"

# Create a folder if is not exists
if(!dir.exists(mypath_out)) dir.create(mypath_out)


# Get the param values from the CSV-file.
parameters <- read.csv2("Listado_sec.csv", header = T)

for (i in nrow(parameters)){
  # assemble the output file name
  row <- parameters[i,]
  outfile.name <- paste0("Biostring_Covid-19_",row[5])
  fn <-file.path(mypath_out, outfile.name)
  
  #Process                                    
  rmarkdown::render("/home/ajo/gitRepos/bioinformatics/Task3/Biostring_Covid-19_Tasca.Rmd", params = list(
                                                                 gene_seq=row[1], 
                                                                 start_vec=row[2], 
                                                                 end_vec=row[3],
                                                                 protein_seq=row[4],
                                                                 name_protein=row[5],
                                                                 path_in=mypath_in,
                                                                 path_out=mypath_out
                                                                 ),
                    output_format =  c("html_document","pdf_document"),
                    output_file = c(paste0(fn, output_format ='.html'),
                                    paste0(fn, output_format ='.pdf')),
                    encoding = "UTF8")
  
}

# assemble the output file name
#row <- strsplit(parameters[2,], split="\t")[[1]] # In case of header = F when getting values from CSV-file.
row <- parameters[2,]
outfile.name <- paste0("Biostring_Covid-19_",row[5])
fn <-file.path(mypath_out, outfile.name)

#Process                                    
rmarkdown::render("/home/ajo/gitRepos/bioinformatics/Task3/Biostring_Covid-19_Tasca.Rmd", params = list(
  gene_seq=row[1], 
  start_vec=row[2]$pos_start, 
  end_vec=row[3]$pos_end,
  protein_seq=row[4],
  name_protein=row[5],
  path_in=mypath_in,
  path_out=mypath_out
),
output_format =  c("html_document","pdf_document"),
output_file = c(paste0(fn, output_format ='.html'),
                paste0(fn, output_format ='.pdf')),
encoding = "UTF8")
