# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function to obtain block diagonal structure through thresholding
#'
#' Given a correlation matrix \eqn{R}, this function discovers the block-diagonal structure by thresholding the absolute values of the entries of the correlation matrix at \eqn{c}. We create an adjacency matrix with the elements being 1 if and only if the corresponding member of the correlation matrix has an absolute value \eqn{\ge c}. This is equivalent to performing a single linkage hierarchical clustering on the variables, with the distance matrix given by \eqn{1 - |R|} and cutting the tree at height \eqn{1-c}.
#' 
#' @param R a \eqn{p \times p} correlation matrix
#' @param c a threshold
#' @param fig a logical input; if `TRUE`, plots `R`, the adjacency matrix obtained by thresholding `R`, and the undirected graph corresponding to the adjacency matrix; default value is `FALSE`
#' @return A \eqn{p} length integer vector whose \eqn{i^{th}} element denotes the group \eqn{i^{th}} variable belongs to. If `fig = TRUE`, returns the plots for `R`, the adjacency matrix obtained by thresholding `R`, and the undirected graph corresponding to the adjacency matrix.
#' @examples
#' # Simulates a 10 x 5 X from N(0, I)
#' set.seed(1)
#' X <- matrix(rnorm(50), 10, 5)
#' # Compute the correlation matrix of X.
#' corX <- cor(X)
#' # Compute the block diagonal structure at c=0.5
#' block_diag(R=corX, c=0.5)
#' @export 
block_diag <- function(R, c, fig = FALSE){
  dis_R <- 1 - abs(R)
  test <- stats::as.dist(dis_R, diag = TRUE)
  clust_result <- stats::hclust(test, method = "single")
  if(!fig){
    return(stats::cutree(clust_result, h = (1 - c)))
  }
  if(fig){
    adj <- 2*(abs(R) > 0.5)/2
    g1 <- igraph::graph_from_adjacency_matrix(adj , diag = F, mode = "undirected")  
    graphics::par(mfrow=c(1,3))
    graphics::image(t(abs(R))[,ncol(abs(R)):1], main="Absolute sample correlation", col=rev(grDevices::heat.colors(10)), xaxt='n', yaxt='n')
    graphics::image(t(adj)[,ncol(adj):1], main="Adjacency Matrix", col=rev(grDevices::heat.colors(10)), xaxt='n', yaxt='n')
    igraph::plot.igraph(g1, vertex.size = 30, main= "Graph")
    return(stats::cutree(clust_result, h = (1 - c)))
  }
}
