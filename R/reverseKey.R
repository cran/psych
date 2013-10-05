#written by David Stanley
reverseKey=function(itemNames,itemNames2ReverseKey) {
  lengthItems=length(itemNames)
  lengthItemsUnique=length(unique(itemNames))
  
  if (lengthItems!=lengthItemsUnique) {
    print("Item names are not unique: reverseKey will not work")
    return(NA)
  }
  
  keyValues=rep(1,lengthItems)
  matchKey=!is.na(match(itemNames,itemNames2ReverseKey))
  keyValues[matchKey]=-1
  return(keyValues)
}