
extract.centroids <- function(lgl){
    # grabs centroids from a plot so it's not repeated for each input file
    library(data.table)
    lgm <- lgl[[1]][[1]][[2]][,c(1,6,7)]
    lgm <- as.data.table(lgm)
    setnames(lgm, c('cell','x','y'))
    # default center point; overwritten if it is different
    lgm[,dist := sqrt((x-50.25)^2 + (y-50.25)^2)]

    return(lgm)
}
