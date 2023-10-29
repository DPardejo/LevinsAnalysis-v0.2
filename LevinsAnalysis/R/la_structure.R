#' Community Network Structure
#'
#' Computes and enumerates the value of multiple structural properties from the
#' Community Network matrix \code{J}.
#'
#' @param J a square matrix describing the qualitative direct casual effect
#'   between the variables of a Community Network. See \code{Note} for more
#'   details.
#' @param fname the name for the exporting file. By default is \code{NA} and no
#'   file is exported. See \code{Details} for more information on use.
#' @export
#' @importFrom  igraph graph_from_adjacency_matrix
#' @importFrom igraph all_simple_paths
#' @details \code{fname} should be used giving the name intended for the
#'   exporting files noting that a suffix \code{"_cm_strt.txt"} will be added to
#'   the file name. See \code{Examples} for pratical use.
#'
#' @note \itemize{ \item The \code{J} matrix should be formatted as follows:
#'   \itemize{ \item should be written by indicating the qualitative direct
#'   casual effect from j (rows) variables on i (columns) variables; \item the
#'   signs of qualitative direct casual effect allowed are: \code{1} for
#'   positive, \code{-1} for negative and \code{0} for null.} \item One
#'   requisite for Community Network matrix stability is the sign of the the
#'   \code{determinant}. For matrices wih even number of variables, the
#'   \code{determinant} must be positive, while for matrices with odd number of
#'   variables the \code{determinant} must be negative.}
#'
#' @return  It returns a print of the Community Network matrix \code{J} and of
#'   the value for the following strutural properties : No. of Nodes of the
#'   Community Matrix, Total no. of Nodal Interactions, No. of Nodal
#'   Interactions minus self-loops, Total no. of Positive Interactions, No. of
#'   Positive Interactions minus self-loops, Total no. of Negative Interactions,
#'   No. of Negative Interactions minus self-loops, No. of Self-loop
#'   interactions, Number of possible paths of the Community Network matrix
#'   \code{J}, Number of paths crossing each variable of the Community Network
#'   matrix \code{J}, Trophic Level, Determinant of the Community
#'   Network matrix.
#'
#'   If \code{fname} is indicated, the results will be exported to a
#'   \code{fname.txt} file.
#'
#'   An \code{invisible} list with the results is generated. It can be assigned
#'   to an object.
#' @keywords array, attribute, math, print
#' @concept interaction, levins, loop, network, path, stability
#' @seealso \code{\link{cm_interaction}}, \code{\link{cm_paths}},
#'   \code{\link{cm_stability}}, \code{\link{J}}, \code{\link{igraph}},
#'   \code{\link[igraph]{all_simple_paths}},
#'   \code{\link[igraph]{graph_from_adjacency_matrix}}
#'
#' @examples
#' ## To calculate the structural indices of M without printing the results to a .txt file
#' cm_structure(M, fname = NA)
#' ## or (to assign the invisible results list)
#' results<-cm_structure(M)
#'
#' ## To calculate the structural indices of M, printing the results to a "fname.txt" file
#' cm_structure(M, fname = "M")
#'
#' @author Daniel Pereira \email{dsldlf@@unife.it}, Antonio Bodini, Marco Scotti
#'   \email{marcoscot@@gmail.com}
#' @references \itemize{ \item Our manual; \item Puccia, C. J. and Levins, R.
#'   (1986) Qualitative Modeling of Complex Systems: An Introduction to Loop
#'   Analysis and Time Averaging. Cambridge: Harvard University Press; \item
#'   Levins, R. and Schultz, B. B. (1996) Effects of density dependence,
#'   feedback and environmental sensitivity on correlations among predators,
#'   prey and plant resources: Models and practical implications. Journal of
#'   Animal Ecology, 65(6),802-812}.
#'
#'
la_structure <- function(J, paths = "n", fname = NA){
	##
	## error message if the input is not a square numeric matrix
	if(is.matrix(J) == FALSE | is.numeric(J) == FALSE){
		cat("\nerror: the input must be a square numeric matrix")
		cat("\n\n")
		invisible(NULL)
	}
	else{
		## number of rows and columns
		n_rows <- nrow(J)
		n_cols <- ncol(J)
		##
		## square matrix check
		{
		if(n_rows == n_cols){
			##
			## number of nodes
			n_nodes <- sqrt(length(J))
			##
			## maximum number of interactions possible
			n_nodes2 <- n_nodes^2
			##
			## transpose of the matrix
			t_J <- t(J)
			##
			## setting all non-zero links as equal to 1
			Mp <- ifelse(J < 0, 1, J)
			##
			## community matrix
			cat("\ncommunity matrix:\n")
			print(J)
			cat("\n")
			##
			## number of nodes
			cat(paste("number of nodes: ", n_nodes, sep = ""))
			cat("\n\n")
			##
			## creation of the matrix j0d that corresponds to J without self-loops
			j0d <- J
			diag(j0d) <- 0
			##
			## number of interactions
			t_inter <- length(which(J != 0))
			cat(paste("number of interactions: ", t_inter, sep = ""))
			cat("\n\n")
			##
			## number of interactions without self-loops
			pt_inter <- length(which(j0d != 0))
			cat(paste("number of interactions without self-loops: ", pt_inter, sep = ""))
			cat("\n\n")
			##
			## number of positive interactions
			tp_inter <- length(which(J > 0))
			cat(paste("number of positive interactions: ", tp_inter, sep = ""))
			cat("\n\n")
			##
			## number of positive interactions without self-loops
			pp_inter <- length(which(j0d > 0))
			cat(paste("number of positive interactions without self-loops: ", pp_inter, sep = ""))
			cat("\n\n")
			##
			## number of negative interactions
			tn_inter <- length(which(J < 0))
			cat(paste("number of negative interactions: ", tn_inter, sep = ""))
			cat("\n\n")
			##
			## number of negative interactions without self-loops
			pn_inter <- length(which(j0d < 0))
			cat(paste("number of negative interactions without self-loops: ", pn_inter, sep = ""))
			cat("\n\n")
			##
			## number of self-loops (irrespective of the sign)
			ts_inter <- length(which(diag(J) != 0))
			cat(paste("number of self-loops: ", ts_inter, sep = ""))
			cat("\n\n")
			##
			{
			if(paths == "y"){
				## list of all simple paths starting from any node
				## (i.e. any variable) in the community matrix
				list_paths_Mp <- all_paths_raw(Mp)
				##
				## number of paths
				n_paths <- length(list_paths_Mp)
				cat(paste("number of paths: ", n_paths, sep = ""))
				cat("\n\n")
				##
				## number of paths crossing each node
				unl_list_paths_mp <- unlist(list_paths_Mp)
				n_paths_cross <- unlist(lapply(colnames(J),
				function(x)length(which(unl_list_paths_mp == x))))
				names(n_paths_cross) <- colnames(J)
				cat("paths crossing:\n")
				print(n_paths_cross)
				cat("\n")
				}
			else{
				n_paths <- noquote("value not calculated")
				cat(paste("number of paths: ", n_paths, sep = ""))
				cat("\n\n")
				n_paths_cross <- noquote("vector not generated")
				cat("paths crossing: ", n_paths_cross, sep = "")
				cat("\n\n")
				}
			}
			##
			## calculate determinant
			det_J <- det(t_J)
			##
			## partial feeding matrix
			{
			if(det_J != 0){
				T_j <- ifelse(J < 0, 0, J)
				IN_j <- apply(T_j, 2, sum)
				G_j <- matrix(rep(0, n_nodes^2), n_nodes)
				for(i in 1:n_nodes){
					for(j in 1:n_nodes){
						if(T_j[i,j] != 0)G_j[i,j] <- T_j[i,j]/IN_j[j]
					}
				}
				ID_j <- diag(1, n_nodes)
				t_level <- apply(solve(ID_j - G_j), 2, sum)
				names(t_level) <- colnames(J)
				}
			else{
				t_level <- noquote("vector not calculated (singular matrix)")
				}
			}
			##
			## trophic levels
			{
			if(det_J != 0){
				cat("trophic levels of nodes:\n")
				print(t_level)
				cat("\n")
				}
			else{
				cat("trophic levels of nodes: ", t_level, sep = "")
				cat("\n\n")
				}
			}
			##
			## print determinant
			cat(paste("determinant: ", det_J, sep = ""))
			cat("\n\n")
			##
			## list of results
			list_results <- as.list(rep(NA, 13))
			names(list_results) <- c("community matrix",
			"number of nodes",
			"number of interactions", "number of interactions without self-loops",
			"number of positive interactions", "number of positive interactions without self-loops",
			"number of negative interactions", "number of negative interactions without self-loops",
			"number of self-loops (irrespective of the sign)",
			"number of paths", "paths crossing each node",
			"trophic levels of nodes",
			"determinant")
			##
			list_results[[1]] <- J
			list_results[[2]] <- n_nodes
			list_results[[3]] <- t_inter
			list_results[[4]] <- pt_inter
			list_results[[5]] <- tp_inter
			list_results[[6]] <- pp_inter
			list_results[[7]] <- tn_inter
			list_results[[8]] <- pn_inter
			list_results[[9]] <- ts_inter
			list_results[[10]] <- n_paths
			list_results[[11]] <- n_paths_cross
			list_results[[12]] <- t_level
			list_results[[13]] <- det_J
			##
			## the list of results is saved in a text file
			{
			if(is.character(fname) == TRUE){
				file_name_txt <- paste(fname, "_la_structure.txt", sep = "")
				sink(file_name_txt)
				print(list_results)
				sink()
				cat(paste("results are in the file \"", file_name_txt, "\"", sep = ""))
				cat("\n\n")
				##
				if(paths != "y" & paths != "n"){
					cat("error: paths must indicate the proper keyword (i.e. either \"y\" or \"n\")")
					cat("\n\n")
					}
				}
			##
			else{
				if(length(fname) != 1){
					cat("error: the filename must be of type \"character\"")
					cat("\n\n")
					}
				else if(is.na(fname) == FALSE){
					cat("error: the filename must be of type \"character\"")
					cat("\n\n")
				}
				##
				if(paths != "y" & paths != "n"){
					cat("error: paths must indicate the proper keyword (i.e. either \"y\" or \"n\")")
					cat("\n\n")
					}
				}
			##
			}
			## invisible object for the list of results
			invisible(list_results)
		}
		else{
			cat("\nerror: the input must be a square numeric matrix")
			cat("\n\n")
			invisible(NULL)
			}
		}
	}
}
