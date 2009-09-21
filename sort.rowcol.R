"sort.rowcol" <- function(M) {
  roworder <- order(rowSums(M))
  colorder <- order(colSums(M))
  M <- M[roworder,colorder] }