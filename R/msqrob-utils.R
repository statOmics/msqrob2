#' Smallest unique protein groups
#'
#' @description For a given vector of protein group names, outputs the names of those protein groups for which none of its member proteins is present in a smaller protein group.
#' @param proteins A vector of characters or factors containing single proteins and/or protein groups (i.e. proteins separated by a separator symbol).
#' @param split The character string that is used to separate the indivudual protein names in each protein group.
#' @return A character vector containing the names of the protein groups for which none of its proteins is present in a smaller protein group.
#' @examples #TODO
#' @export
smallestUniqueGroups <- function(proteins, split=";"){

  b <- strsplit(x=as.character(unique(proteins)),split=split,fixed=TRUE)

  erbij <- vector()

  j <- 1
  while(length(b)!=0)
  {

    erbij <- c(erbij,sapply(b[sapply(b, length)==j], function(x) paste(x, collapse=split)))

    a <- unlist(b[sapply(b, length)==j])
    b <- b[sapply(b, length)>j]

    if(length(b)!=0){
      welke <- vector()
      for(i in 1:length(b))
      {
        welke[i] <- !any(b[[i]] %in% a)
      }

      b <- b[welke]
      j <- j+1
    }
  }

  erbij <- unlist(erbij)
  return(erbij)
}
