library(data.table)
library(openxlsx)

wb <- createWorkbook()
encode_dir<-"fastq/kundaje_encode/manorm"

for (subd in list.dirs(encode_dir, full.names=FALSE, recursive=FALSE)) {

  fpaths <- Sys.glob(file.path(encode_dir,subd,"*HCT116*all_MAvalues.xls"))
  message(sprintf('checking %s',subd))
  
  if (length(fpaths)==1){
    message(sprintf('adding %s',fpaths[[1]]))
    dt <- fread(fpaths[[1]],header=T)
    sheet_name <- substring(subd,first=3)
    addWorksheet(wb,sheet_name)
    writeData(wb, sheet=sheet_name,x=dt)
    message(sprintf('[%s] updated',sheet_name))
  }
}

saveWorkbook(wb, file.path(encode_dir,"Supplementary_Table_S5.xlsx"), overwrite = TRUE)

