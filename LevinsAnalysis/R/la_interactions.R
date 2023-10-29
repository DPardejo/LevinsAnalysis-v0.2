#' Summary of simple paths
#'
#' \code{cm_interactions} returns the summary of the simple paths
#' between all pairs of nodes in the graph corresponding to the
#' community matrix \code{J}.
#'
#' @param J the square matrix describing the qualitative direct
#'	causal effects between the variables of the community matrix.
#'	See \code{Note} for more details.
#' @param NtN the vector that indicates the nodes from which
#'	(first element) and/or to which (second element) simple
#'	paths should be considered. See \code{Details} for additional
#'	information.
#' @param PL selects the path lengths for which results have to
#'	be returned: \code{"all"} (default), \code{"sum"}, or vector
#'	with integer, strictly positive numbers shorter than the
#'	maximum path length in the graph. See \code{Details} for
#'	further information.
#' @param fname the name of the file where the results are stored.
#'	No file is exported when NA (default) is selected.
#'
#' @export
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph all_simple_paths
#' @importFrom igraph V
#' @importFrom igraph as_ids
#'
#' @details \itemize{ \item The argument \code{NtN} allows the
#'	user to choose whether: \itemize{ \item shortest paths
#'	between all pairs of nodes are returned (in the outcome the
#'	simple paths are classified based on their length and also
#'	the sum of all simple paths is returned by default); for
#'	generating this output the vector used for the argument
#'	\code{NtN} should be composed of two \code{NA} elements:
#'  \code{NtN} = \code{c(NA,NA)} (default option); \item
#'	all shortest paths that start from one specific node are
#'	returned; if this option is selected, the first element of
#'	the vector used for the argument \code{NtN} must indicate
#'	either the name of the node or the corresponding row index,
#'	while the second element must be equal to \code{NA}:
#'	\code{NtN} = \code{c(m,NA)}, where \code{m} is an integer,
#'	strictly positive number;
#'
#'
#'
#'
#'   the user to choose the level of information to be printed, by use of a
#'   vector, as follows: \itemize{ \item for printing the interaction
#'   information between all nodes (regardless of path lenght) plus a sum
#'   matrix, the vector \code{NtN} should be empty (\code{c()}; default); \item
#'   for priting the interation information for all paths begining in node
#'   \code{i} (regardless of path lenght) plus a sum matrix, the vector
#'   \code{NtN} should be of the form \code{c(i)}; \item for printing the
#'   interaction information between two specific \code{i} and \code{j} nodes
#'   (regardless of path lenght) plus a sum matrix, the vector \code{NtN} should
#'   be of the form \code{c(i,j)}. \item see \code{Note} for more information
#'   and \code{Examples} for practical use.} \item Through the argument
#'   \code{PL} the function allows the user to choose the path length from which
#'   he wants the information to be printed, by choosing the appropriate option
#'   as follows: \itemize{ \item \code{"all"} for printing the information for
#'   all path lengths plus a sum matrix (default); \item \code{t} for printing
#'   the information for paths of length \code{t}, where \code{t} is the path
#'   length desired to be printed (note that the length of a \code{simple path}
#'   is equal to the number of variables of the path minus 1); \item
#'   \code{"sum"} to print only the sum matrix. \item see \code{Note} for more
#'   information and \code{Examples} for practical use.} \item \code{fname}
#'   should be used giving the name intended for the exporting files noting that
#'   a suffix \code{"_cm_inter.txt"} will be added to the file name. See
#'   \code{Examples} for practical use.}
#'
#' @note \itemize{ \item The \code{J} matrix should be formatted as follows:
#'   \itemize{ \item should be written by indicating the qualitative direct
#'   casual effect from j (rows) variables on i (columns) variables; \item the
#'   signs of qualitative direct casual effect allowed are: \code{1} for
#'   positive, \code{-1} for negative and \code{0} for null.} \item Arguments
#'   \code{NtN} and \code{PL} can't be used together, i.e, if one or two
#'   variables are passed on the argument \code{NtN}, then all the information
#'   will be retrieved for those nodes, passing the argument \code{PL} as
#'   \code{"all"} even if this is changed by the user in the function call; \item See
#'   \code{Examples} for practical use.}
#'
#' @return List of matrices that summarize the number of pathways from \code{j}
#'   (row) nodes to \code{i} (column) nodes; each matrix corresponds to pathways
#'   of different lengths; plus a sum matrix of all interactions; note that it
#'   counts only \code{simple paths}, that is, paths that pass only once on each
#'   node.
#'
#'   If \code{fname} is indicated, the results will be exported to a
#'   \code{fname.txt} file.
#'
#'   An \code{invisible} list with the results is generated. It can be assigned
#'   to an object.
#'
#' @keywords arith, array, attribute, list, math, print
#' @concept interaction, levins, loop, network, path, summary
#'
#' @seealso \code{\link{cm_paths}}, \code{\link{cm_structure}}, \code{\link{J}},
#'   \code{\link{igraph}}, \code{\link[igraph]{all_simple_paths}},
#'   \code{\link[igraph]{graph_from_adjacency_matrix}}
#'
#' @examples
#' ## To calculate the interaction level of all nodes of M without printing the results to a .txt file
#' cm_interaction(M, NtN = c(), PL = "all", fname = NA)
#' ## or (to assign the invisible results list)
#' results<-cm_interaction(M, NtN = c(), PL = "all", fname = NA)
#'
#' ## To calculate the interaction level of all nodes of M, printing the results to a "fname.txt" file
#' cm_interaction(M, NtN = c(), PL = "all", fname = "M")
#'
#' ## To calculate the interaction level of all paths begining in the node "A"
#' cm_interaction(M, NtN = c("A"), PL = "all", fname = NA)
#'
#' ## To calculate the interaction level of all paths begining in the node "A" and ending in node "E"
#' cm_interaction(M, NtN = c("A", "E"), PL = "all", fname = NA)
#'
#' ## To calculate the interaction level of all paths of length 4
#' cm_interaction(M, NtN= c(), PL = 4, fname = NA)
#'
#' ## To calculate the sum matrix of the interaction level of all nodes of M
#' cm_interaction(M, NtN= c(), PL = "sum", fname = NA)
#'
#' @author Daniel Pereira \email{dsldlf@@unife.it}, Antonio Bodini, Marco Scotti
#'   \email{marcoscot@@gmail.com}
#' @references \itemize{ \item Our manual; \item Puccia, C. J. and Levins, R.
#'   (1986) Qualitative Modeling of Complex Systems: An Introduction to Loop
#'   Analysis and Time Averaging. Cambridge: Harvard University Press; \item
#'   Levins, R. and Schultz, B. B. (1996) Effects of density dependence,
#'   feedback and environmental sensitivity on correlations among predators,
#'   prey and plant resources: Models and practical implications. Journal of
#'   Animal Ecology, 65(6),802-812.}
#'
#'
la_interactions <- function(J, PL = "all", NtN = c(NA,NA), fname = NA){
	##
	counter <- 0
	##
	## error message if the input is not a square numeric matrix
	{
	if(is.matrix(J) == FALSE | is.numeric(J) == FALSE){
		counter <- 1
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
			## number and name of nodes
			nodi <- nrow(J)
			n_nodi <- c(1:nodi)
			nomi <- colnames(J)
			##
			## list of all simple paths starting from any node
			## (i.e. any variable) in the community matrix
			Mp <- ifelse(J < 0, 1, J)
			list_paths_Mp <- all_paths_raw(Mp)
			##
			## number of nodes in each path
			path_var_count <- unlist(lapply(list_paths_Mp,length))
			##
			## maximum length of the simple paths
			uni_len <- unique(path_var_count) - 1
			mxpl <- max(uni_len)
			##
			## creation of empty matrices to store numbers of
			## simple paths between all pairs of nodes
			M_list_all_paths <- as.list(rep(NA,mxpl))
			for(i in 1:mxpl){
				Mpres_n <- matrix(rep(0,length(Mp)), sqrt(length(Mp)))
				colnames(Mpres_n) <- colnames(Mp)
				rownames(Mpres_n) <- rownames(Mp)
				M_list_all_paths[[i]] <- Mpres_n
			}
			##
			## count of simple paths with different lengths
			## between all pairs of nodes
			start_node <- unlist(lapply(list_paths_Mp,function(x)x[1]))
			lalp <- length(path_var_count)
			end_node <- unlist(lapply(list_paths_Mp,function(x)x[length(x)]))
			##
			for(i in 1:lalp){
				lotp <- path_var_count[i] - 1
				snn <- start_node[i]
				enn <- end_node[i]
				M_list_all_paths[[lotp]][snn,enn] <- M_list_all_paths[[lotp]][snn,enn] + 1
			}
			total_int <- Reduce("+", M_list_all_paths)
			##
			## [1] print matrices that summarize the number of pathways
			## from row nodes to column nodes
			## [2] each matrix corresponds to pathways of different lengths
			## (i.e. from n = 1 step to n = max number of steps)
			## [3] also a matrix summarizing the total number paths between
			## all pairs of nodes is generated
			## [4] the user can query either information on paths of a specific
			## length or visualize the summary matrix for all paths
			{
			if(length(which(is.na(NtN) == TRUE)) == 2 & length(NtN) == 2 & is.vector(NtN) == TRUE){
				if(PL[1] != "all" & PL[1] != "sum" & is.numeric(PL) == FALSE){
					counter <- 1
					cat("\nerror: PL must indicate the proper keyword (i.e. either \"all\" or \"sum\") or include strictly positive integers")
					cat("\n\n")
					list_results <- NULL
				}
				##
				if(PL[1] == "all"){
					cat(paste("\npaths by length"))
					cat("\n\n")
					for(i in 1:mxpl){
						cat(paste("paths of length ", i, sep = ""))
						cat("\n\n")
						print.noquote(M_list_all_paths[[i]])
						cat("\n\n")
					}
					cat(paste("total paths"))
					cat("\n\n")
					print.noquote(total_int)
					cat("\n")
					nameslist <- c(paste("paths of length ", c(1:mxpl), sep = ""),"total paths")
					list_results <- M_list_all_paths
					list_results[[(mxpl+1)]] <- total_int
					names(list_results) <- nameslist
				}
				##
				if(is.numeric(PL) == TRUE){
					if(is.matrix(PL) == TRUE | is.list(PL) == TRUE |
					check_positive_integers(PL) == FALSE){
						counter <- 1
						cat("\nerror: PL must indicate the proper keyword (i.e. either \"all\" or \"sum\") or include strictly positive integers")
						cat("\n\n")
						list_results <- NULL
					}
					else{
						over <- which(PL > mxpl)
						if(length(over) > 0)PL <- PL[-over]
						ll <- length(PL)
						{
						if(ll > 0){
							nameslist <- paste("paths of length ", PL, sep = "")
							list_results <- as.list(rep(NA,ll))
							for(i in 1:ll){
								if(i == 1){
									cat(paste("\npaths of length ", PL[i], sep = ""))
									cat("\n\n")
									print.noquote(M_list_all_paths[[PL[i]]])
									cat("\n")
									list_results[[i]] <- M_list_all_paths[[PL[i]]]
									}
								else{
									if(i == ll){
										cat(paste("\npaths of length ", PL[i], sep = ""))
										cat("\n\n")
										print.noquote(M_list_all_paths[[PL[i]]])
										cat("\n")
										list_results[[i]] <- M_list_all_paths[[PL[i]]]
										}
									else{
										cat(paste("\npaths of length ", PL[i], sep = ""))
										cat("\n\n")
										print.noquote(M_list_all_paths[[PL[i]]])
										cat("\n\n")
										list_results[[i]] <- M_list_all_paths[[PL[i]]]
										}
								}
							}
							names(list_results) <- nameslist
							}
						else{
							list_results <- NULL
							PL <- "max_length"
							}
						}
						if(length(over) > 0){
							if(ll > 0){
								cat("error: path length(s) required longer than the maximum path length in the network (i.e. ",
								mxpl, ")", sep = "")
								cat("\n\n")
							}
							else{
								counter <- 1
								cat("\nerror: path length(s) required longer than the maximum path length in the network (i.e. ",
								mxpl, ")", sep = "")
								cat("\n\n")
							}
						}
					}
				}
				##
				if(PL[1] == "sum"){
					cat(paste("\ntotal paths"))
					cat("\n\n")
					print.noquote(total_int)
					cat("\n")
					list_results <- list(total_int)
					names(list_results) <- "total paths"
				}
			}
			##
			else{
				if(length(NtN) == 2 & is.vector(NtN) == TRUE){
					if(length(which(is.na(NtN) == FALSE)) == 2){
						if((length(which(NtN[1] == nomi)) == 1 | length(which(NtN[1] == n_nodi)) == 1) &
						(length(which(NtN[2] == nomi)) == 1 | length(which(NtN[2] == n_nodi)) == 1)){
							k <- NtN[1]
							j <- NtN[2]
							cat(paste("\npaths between nodes", k, "and", j, sep = " "))
							cat("\n\n")
							for(i in 1:mxpl){
								cat(paste("paths of length ", i, sep = ""))
								cat("\n\n")
								print.noquote(M_list_all_paths[[i]][k,j])
								cat("\n\n")
							}
							cat("total paths")
							cat("\n\n")
							print.noquote(total_int[k,j])
							cat("\n")
							nameslist <- paste("paths of length ", c(1:mxpl), sep = "")
							M_list_paths_i_j <- as.list(rep(0, mxpl))
							names(M_list_paths_i_j) <- nameslist
							for(i in 1:mxpl){
								if(M_list_all_paths[[i]][k,j] > 0){
									M_list_paths_i_j[[i]] <- M_list_all_paths[[i]][k,j]
								}
							}
							nameslist2 <- c(nameslist,"total paths")
							list_results <- M_list_paths_i_j
							list_results[[(mxpl+1)]] <- total_int[k,j]
							names(list_results) <- nameslist2
						}
						else{
							counter <- 1
							cat("\nerror: NtN must indicate the names (or integer codes) of two network nodes")
							cat("\n\n")
							list_results <- NULL		
						}
					}
					##
					else{
						if(is.na(NtN[1]) == FALSE & is.na(NtN[2]) == TRUE){
							if(length(which(NtN[1] == nomi)) == 1 | length(which(NtN[1] == n_nodi)) == 1){
								cat(paste("\npaths from node", NtN[1], sep = " "))
								cat("\n\n")
								k <- NtN[1]
								for(i in 1:mxpl){
									cat(paste("paths of length ", i, sep = ""))
									cat("\n\n")
									print.noquote(M_list_all_paths[[i]][k,])
									cat("\n\n")
								}
								cat("total paths")
								cat("\n\n")
								print.noquote(total_int[k,])
								cat("\n")
								nameslist <- paste("paths of length ", c(1:mxpl), sep = "")
								M_list_paths_i <- as.list(rep(NA, mxpl))
								for(i in 1:mxpl){
									M_list_paths_i[[i]] <- M_list_all_paths[[i]][k,]
								}
								nameslist2 <- c(nameslist,"total paths")
								list_results <- M_list_paths_i
								list_results[[(mxpl+1)]] <- total_int[k,]
								names(list_results) <- nameslist2
							}
							else{
								counter <- 1
								cat("\nerror: NtN must indicate the names (or integer codes) of two network nodes")
								cat("\n\n")
								list_results <- NULL
							}
						}
						else{
							if(is.na(NtN[1]) == TRUE & is.na(NtN[2]) == FALSE){
								if(length(which(NtN[2] == nomi)) == 1 | length(which(NtN[2] == n_nodi)) == 1){
									cat(paste("\npaths to node", NtN[2], sep = " "))
									cat("\n\n")
									k <- NtN[2]
									for(i in 1:mxpl){
										cat(paste("paths of length ", i, sep = ""))
										cat("\n\n")
										print.noquote(M_list_all_paths[[i]][,k])
										cat("\n\n")
									}
									cat("total paths")
									cat("\n\n")
									print.noquote(total_int[,k])
									cat("\n")
									nameslist <- paste("paths of length ", c(1:mxpl), sep = "")
									M_list_paths_i <- as.list(rep(NA, mxpl))
									for(i in 1:mxpl){
										M_list_paths_i[[i]] <- M_list_all_paths[[i]][,k]
									}
									nameslist2 <- c(nameslist,"total paths")
									list_results <- M_list_paths_i
									list_results[[(mxpl+1)]] <- total_int[,k]
									names(list_results) <- nameslist2
								}
								else{
									counter <- 1
									cat("\nerror: NtN must indicate the names (or integer codes) of two network nodes")
									cat("\n\n")
									list_results <- NULL
								}
							}
						}
					}
				}
				else{
					counter <- 1
					cat("\nerror: NtN must indicate the names (or integer codes) of two network nodes")
					cat("\n\n")
					list_results <- NULL
				}
			}
			##
			}
			##
			## the list of results is saved in a text file
			{
			if(is.character(fname) == TRUE){
				if(counter == 0){
					file_name_txt <- paste(fname, "_la_interactions.txt", sep = "")
					{
					if(length(which(is.na(NtN) == TRUE)) == 2){
						sink(file_name_txt)
						cat("paths by length")
						cat("\n\n")
						print(list_results)
						sink()
						}
					else{
						if(length(which(is.na(NtN) == FALSE)) == 2){
							sink(file_name_txt)
							cat("paths between nodes ", NtN[1], " and ", NtN[2], sep = "")
							cat("\n\n")
							print(list_results)
							sink()
							}
						else{
							if(is.na(NtN[2]) == TRUE){
								sink(file_name_txt)
								cat("paths from node ", NtN[1], sep = "")
								cat("\n\n")
								print(list_results)
								sink()
								}
							else{
								sink(file_name_txt)
								cat("paths to node ", NtN[2], sep = "")
								cat("\n\n")
								print(list_results)
								sink()
								}
							}
						}
					}
					cat(paste("results are in the file \"", file_name_txt, "\"", sep = ""))
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
			}
			##
			}
			##
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
}
