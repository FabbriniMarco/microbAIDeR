addSheet <- function(path,sheet.name, addition, col.save=TRUE, row.save=TRUE, overwrite = TRUE){
  check_and_load_package("openxlsx")
  wb <- openxlsx::loadWorkbook(path)
  if ( sheet.name %in% openxlsx::getSheetNames(path) & overwrite ){ openxlsx::removeWorksheet(wb, sheet.name) }
  openxlsx::addWorksheet(wb,sheet.name)
  openxlsx::writeData(wb,sheet.name, as.data.frame(addition), colNames = col.save, rowNames = row.save )
  openxlsx::saveWorkbook(wb,path,overwrite = TRUE)
}