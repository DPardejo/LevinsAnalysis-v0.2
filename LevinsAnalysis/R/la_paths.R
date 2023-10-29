#' Community Network Paths Information
#'
#' Generates a list describing each \code{simple path} from the Community
#' Network matrix \code{J} with their length, sign and strength information.
#'
#' @param J a square matrix describing the qualitative direct casual effect
#'   between the variables of a Community Network. See \code{Note} for more
#'   details.
#' @param SM a square matrix with the interaction strength between the variable
#'   of the Community Network matrix \code{J}. See \code{Note} for more
#'   details.
#' @param type_info selects which path info (if any) on path information should
#'   be in place for the printing of the paths paths information: \code{"all"}
#'   (default), \code{"P"}, \code{"N"}, \code{"PN"}, \code{"VarC"},
#'   \code{"StrgP"}, \code{"NodetN"}. See \code{Details} for more information on
#'   use.
#' @param fname the name for the exporting file. Indication of parameter is
#'   obligatory. See \code{Details} for more information on use.
#'
#' @export
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph all_simple_paths
#' @importFrom igraph V
#' @importFrom igraph as_ids
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom utils write.table
#'
#'
#' @details  \itemize{ \item Through the argument \code{type_info} the function
#'   allows the user to choose the path information to be returned
#'   for the printing, as follows: \itemize{ \item \code{"all"} indicates
#'   that there is no type_info in place and the info from all paths should be
#'   printed (default); \item \code{"P"} indicates that the info to be printed
#'   should be only from the paths with positive sign; \item \code{"N"}
#'   indicates that the info to be printed should be only from the paths with
#'   negative sign; \item \code{"PN"} indicates that the both the info from
#'   paths with positive and negative sign should be printed (in to separated
#'   files); \item \code{"VarC"} indicates that the info to be printed should be
#'   only from paths of higher, equal or lower number of variable as indicated
#'   in the  \code{VarC} vector. The vector should be constructed in the form
#'   \code{c("Up", x)}, as follows: \itemize{ \item \code{"Up", "Eq" or "Dn"}
#'   first element of the vector indicating that the print should be for paths
#'   with variable number higher, equal or below to the value indicated in the
#'   second element; \item \code{x} second element indicating the number of
#'   variables of the path.} \item \code{"StrgP"} indicated that the info to be
#'   printed should be only from paths of higher, equal or lower path strength
#'   as indicated in the \code{StrgP} vector. The vector should be constructed
#'   in the form \code{c("Up", x)}, as follows: \itemize{ \item \code{"Up" or
#'   "Dn"} first element of the vector indicating that the print should be for
#'   paths with strength higher or below to the value indicated in the second
#'   element; \item \code{x} second element indicating the strength of the
#'   path.} \item \code{"NodetN"} indicates that the info to be printed should
#'   be only from paths begin and/or ending in the nodes indicated in the
#'   \code{NodetN} vector. The vector should be constructed in the form
#'   \code{c(x, y)}, as follows: \itemize{ \item \code{x} first element of the
#'   vector indicating that the print should be for paths begining in the
#'   variable indicated; \item \code{x} second element of the vector indicating
#'   that the print should be for paths ending in the variable indicated.} \item
#'   see \code{Examples} for practical use.} \item \code{fname} should be used
#'   giving the name intended for the exporting files noting that the following
#'   suffixes will be added to the file name: \itemize{ \item
#'   \code{"_path_info.txt"} for \code{type_info
#'   = "all"}; \item \code{"_path_info_sign_p.txt"} for \code{type_info = "P"};  \item
#'   \code{"_path_info_sign_n.txt"} for
#'   \code{type_info = "N"};  \item \code{"_path_info_sign_p.txt"},
#'   \code{"_path_info_sign_n.txt"} for
#'   \code{type_info = "PN"};  \item \code{"_path_info_above_varc.txt"} for \code{type_info = "VarC"} where
#'   the first element of the vector \code{VarC} is \code{Up}; \item
#'   \code{"_path_info_equal_varc.txt"} for \code{type_info = "VarC"} where
#'   the first element of the vector \code{VarC} is \code{Eq}; \item
#'   \code{"_path_info_below_varc.txt"} for \code{type_info = "VarC"} where
#'   the first element of the vector \code{VarC} is \code{Dn}; \item
#'   \code{"_path_info_above_strength.txt"} for \code{type_info = "StrgP"}
#'   where the first element of the vector \code{StrgP} is \code{Up}; \item
#'   \code{"_path_info_below_strength.txt"} for \code{type_info = "StrgP"}
#'   where the first element of the vector \code{StrgP} is \code{Dn};\item
#'   \code{"_path_info_begin_nodetn.txt"} for \code{type_info = "NodetN"}
#'   where the first element of the vector \code{NodetN} is a variable of the
#'   Community Network matrix \code{J} and the second element is 0;\item
#'   \code{"_path_info_end_nodetn.txt"} for \code{type_info = "NodetN"} where
#'   the first element of the vector \code{NodetN} is 0 and the second element
#'   is a variable of the Community Network matrix \code{J}; \item
#'   \code{"_path_info_begin_end_nodetn.txt"} for \code{type_info = "NodetN"}
#'   where the both elements of the vector \code{NodetN} are variables of the
#'   Community Network matrix \code{J}; \item see \code{Note} for more
#'   information and \code{Examples} for practical use.}}
#'
#' @note \itemize{ \item The \code{J} matrix should be formatted as follows:
#'   \itemize{ \item should be written by indicating the qualitative direct
#'   casual effect from j (rows) variables on i (columns) variables; \item The
#'   signs of qualitative direct casual effect allowed are: \code{1} for
#'   positive, \code{-1} for negative and \code{0} for null.} \item The
#'   \code{SM} matrix should be formatted as follows: \itemize{ \item Should be
#'   written by indicating the qualitative direct casual effect from j (rows)
#'   variables on i (columns) variables; \item The strength values between each
#'   variable should be within the interval \code{[0,1]}, where \code{0}
#'   indicates no interaction, and \code{1} indicates maximum interaction
#'   strength.} \item Regarding the argument \code{fname}, if not specified the
#'   function will be stopped and an \code{error message} will be shown
#'   indicating to specify \code{fname} in the function call.}
#'
#' @return It prints a list of each path within the selected type_info option
#'   indictating the variables of the path, its strength, the number of nodes
#'   and its sign to \code{fname.txt} file.
#'
#'   An \code{invisible} list with the results is generated. It can be assigned
#'   to an object.
#'
#' @keywords arith, array, attribute, list, print
#' @concept length, levins, network, paths, sign, strength, summary,
#'
#' @seealso \code{\link{cm_interaction}}, \code{\link{cm_strength}},
#'   \code{\link{cm_structure}}, \code{\link{M}}, \code{\link{SM}},
#'   \code{\link{igraph}}, \code{\link[igraph]{all_simple_paths}},
#'   \code{\link[igraph]{graph_from_adjacency_matrix}}
#'
#' @examples ## To print the information from all paths of M, printing the results to "fname.txt" and "fname.xlsx" files
#' cm_paths(M, SM, type_info = "all", fname = "M")
#' ## or (to assign the invisible results list)
#' results <- cm_paths(M, SM, type_info = "all", fname = "M")
#'
#' ## To print the information from paths with positive sign
#' cm_paths(M, SM, type_info = "P", fname = "M")
#'
#' ## To print the information from paths with more than 5 variables
#' VarC <- c("Up", 3)
#' cm_paths(M, SM, type_info = "VarC", fname = "M")
#'
#' ## To print the information from paths with strength lower than 0.5
#' StrgP <- c("Dn", 0.5)
#' cm_paths(M, SM, type_info = "StrgP", fname = "M")
#'
#' ## To print the information from paths begining in node "A" and ending in node "F"
#' NodetN <- c("A", "E")
#' cm_paths(M, SM, type_info = "NodetN", fname = "M")
#'
#'
#'
#' @author Daniel Pereira \email{dsldlf@@unife.it}, Antonio Bodini, Marco Scotti
#'   \email{marcoscot@@gmail.com}
#' @references \itemize{ \item Our manual; \item and something from Loop
#'   analysis}
#'
#'
#'
la_paths <- function(J, SM, type_info = "all", path_info = NA, fname = "community"){
	##
	## error message if the input is not a square numeric matrix
	if(is.matrix(J) == FALSE | is.numeric(J) == FALSE){
		cat("\nerror: the input must be a square numeric matrix")
		cat("\n\n")
		return(invisible(NULL))
	}
	else{
		## number of rows and columns
		n_rows <- nrow(J)
		n_cols <- ncol(J)
		##
		## square matrix check
		{
		if(n_rows == n_cols){
			nomilista <- c("total number of matrices generated",
			"total number of stable matrices",
			"probability distribution",
			"average strength matrix")
			{
			if(length(names(SM)) == 0) l_ml <- 1
			else l_ml <- length(which(is.na(match(names(SM),nomilista)) == TRUE))
			}
			##
			{
			if(is.list(SM) == TRUE){
				if(l_ml != 0){
					cat("\nerror: SM must be a list generated by the function \"la_strength\"")
					cat("\n\n")
					return(invisible(NULL))
					}
				else{
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
					##
					{
					if(is.character(fname) == FALSE){
						cat("\nwarning: the name of the .txt file must be of type \"character\"")
						cat("\n")
						cat("the name was reset to default option (i.e. \"community\")")
						cat("\n")
						fname <- "community"
						}
					else if(is.character(fname) == TRUE & length(fname) > 1){
						cat("\nwarning: the name of the .txt file must be a single element of type \"character\"")
						cat("\n")
						cat("the name was reset to default option (i.e. \"community\")")
						cat("\n")
						fname <- "community"
						}
					}
					##
					##
					all_thresholds <- c("all", "P", "N", "PN", "VarC", "StrgP", "NodetN")
					{
					if(length(type_info) > 1 | is.vector(type_info) == FALSE){
						cat("\nwarning: type_info must be one of the following: \"all\", \"P\", \"N\", \"PN\", \"VarC\", \"StrgP\", \"NodetN\"")
						cat("\n")
						cat("type_info was reset to default option (i.e. \"all\")")
						cat("\n")
						type_info <- "all"
						}
					else{
						s_type_info <- which(all_thresholds == type_info)
						if(length(s_type_info) != 1){
							cat("\nwarning: type_info must be one of the following: \"all\", \"P\", \"N\", \"PN\", \"VarC\", \"StrgP\", \"NodetN\"")
							cat("\n")
							cat("type_info was reset to default option (i.e. \"all\")")
							cat("\n")
							type_info <- "all"
							}
						}
					}
					##
					##
					all_extra <- c("VarC", "StrgP", "NodetN")
					all_VarC <- c("Up","Eq","Dn")
					all_StrgP <- c("Up","Dn")
					all_NodetN <- rownames(J)
					##
					s_type_info_extra <- which(all_extra == type_info)
					if(length(s_type_info_extra) == 1){
						if(type_info == "VarC"){
							if(is.data.frame(path_info) == TRUE){
								if(ncol(path_info) == 2 & nrow(path_info) == 1 &
								is.character(path_info[1,1]) == TRUE &
								check_positive_integers_char(path_info[1,2]) == TRUE){
									s_VarC <- which(all_VarC == path_info[1,1])
									{
									if(length(s_VarC) == 1){
										VarC1 <- path_info[1,1]
										}
									else{
										VarC1 <- "wrong"
										cat("\nerror: path_info must be a two column, single row data frame where:\n")
										cat("(1) the option associated to type_info \"VarC\" must be one of the following: \"Up\", \"Eq\" or \"Dn\"\n")
										cat("(2) the number of path variables must be \"numeric\" (and integer) in the interval [2, (max_path_length + 1)]")
										cat("\n\n")
										return(invisible(NULL))
										}
									}
									##
									{
									if(path_info[1,2] > 1 & path_info[1,2] <= (mxpl+1)){
										VarC2 <- path_info[1,2]
										}
									else{
										VarC2 <- -1
										cat("\nerror: path_info must be a two column, single row data frame where:\n")
										cat("(1) the option associated to type_info \"VarC\" must be one of the following: \"Up\", \"Eq\" or \"Dn\"\n")
										cat("(2) the number of path variables must be \"numeric\" (and integer) in the interval [2, (max_path_length + 1)]")
										cat("\n\n")
										return(invisible(NULL))
										}
									}
								}
								else{
									VarC1 <- "wrong"
									VarC2 <- -1
									cat("\nerror: path_info must be a two column, single row data frame where:\n")
									cat("(1) the option associated to type_info \"VarC\" must be one of the following: \"Up\", \"Eq\" or \"Dn\"\n")
									cat("(2) the number of path variables must be \"numeric\" (and integer) in the interval [2, (max_path_length + 1)]")
									cat("\n\n")
									return(invisible(NULL))
									}
								}
							else{
								VarC1 <- "wrong"
								VarC2 <- -1
								cat("\nerror: path_info must be a two column, single row data frame where:\n")
								cat("(1) the option associated to type_info \"VarC\" must be one of the following: \"Up\", \"Eq\" or \"Dn\"\n")
								cat("(2) the number of path variables must be \"numeric\" (and integer) in the interval [2, (max_path_length + 1)]")
								cat("\n\n")
								return(invisible(NULL))
							}
						}
						##
						if(type_info == "StrgP"){
							if(is.data.frame(path_info) == TRUE){
								if(ncol(path_info) == 2 & nrow(path_info) == 1 &
								is.character(path_info[1,1]) == TRUE &
								is.numeric(path_info[1,2]) == TRUE){
									s_StrgP <- which(all_StrgP == path_info[1,1])
									{
									if(length(s_StrgP) == 1){
										StrgP1 <- path_info[1,1]
										}
									else{
										StrgP1 <- "wrong"
										cat("\nerror: path_info must be a two column, single row data frame where:\n")
										cat("(1) the option associated to type_info \"StrgP\" must be one of the following: \"Up\" or \"Dn\"\n")
										cat("(2) the strength of the path must be of type \"numeric\" in the interval [1e-6, 1]")
										cat("\n\n")
										return(invisible(NULL))
										}
									}
									##
									{
									if(path_info[1,2] > 0 & path_info[1,2] <= 1){
										StrgP2 <- path_info[1,2]
										}
									else{
										StrgP2 <- -1
										cat("\nerror: path_info must be a two column, single row data frame where:\n")
										cat("(1) the option associated to type_info \"StrgP\" must be one of the following: \"Up\" or \"Dn\"\n")
										cat("(2) the strength of the path must be of type \"numeric\" in the interval [1e-6, 1]")
										cat("\n\n")
										return(invisible(NULL))
										}
									}
								}
								else{
									StrgP1 <- "wrong"
									StrgP2 <- -1
									cat("\nerror: path_info must be a two column, single row data frame where:\n")
									cat("(1) the option associated to type_info \"StrgP\" must be one of the following: \"Up\" or \"Dn\"\n")
									cat("(2) the strength of the path must be of type \"numeric\" in the interval [1e-6, 1]")
									cat("\n\n")
									return(invisible(NULL))
									}
								}
							else{
								StrgP1 <- "wrong"
								StrgP2 <- -1
								cat("\nerror: path_info must be a two column, single row data frame where:\n")
								cat("(1) the option associated to type_info \"StrgP\" must be one of the following: \"Up\" or \"Dn\"\n")
								cat("(2) the strength of the path must be of type \"numeric\" in the interval [1e-6, 1]")
								cat("\n\n")
								return(invisible(NULL))
							}
						}
						##
						if(type_info == "NodetN"){
							if(is.data.frame(path_info) == TRUE){
								if(ncol(path_info) == 2 & nrow(path_info) == 1 &
								is.character(path_info[1,1]) == TRUE &
								is.character(path_info[1,2]) == TRUE){
									s_NodetN1 <- which(all_NodetN == path_info[1,1])
									{
									if(length(s_NodetN1) == 1){
										NodetN1 <- path_info[1,1]
										}
									else{
										NodetN1 <- "wrong1"
										cat("\nerror: path_info must be a two column, single row data frame where:\n")
										cat("(1) the first node name must match one of the rownames of the community matrix\n")
										cat("(2) the second node name must match one of the rownames of the community matrix")
										cat("\n\n")
										return(invisible(NULL))
										}
									}
									##
									s_NodetN2 <- which(all_NodetN == path_info[1,2])
									{
									if(length(s_NodetN2) == 1){
										NodetN2 <- path_info[1,2]
										}
									else{
										NodetN2 <- "wrong2"
										cat("\nerror: path_info must be a two column, single row data frame where:\n")
										cat("(1) the first node name must match one of the rownames of the community matrix\n")
										cat("(2) the second node name must match one of the rownames of the community matrix")
										cat("\n\n")
										return(invisible(NULL))
										}
									}	
								}
								else{
									NodetN1 <- "wrong1"
									NodetN2 <- "wrong2"
									cat("\nerror: path_info must be a two column, single row data frame where:\n")
									cat("(1) the first node name must match one of the rownames of the community matrix\n")
									cat("(2) the second node name must match one of the rownames of the community matrix")
									cat("\n\n")
									return(invisible(NULL))
									}
								}
							else{
								NodetN1 <- "wrong1"
								NodetN2 <- "wrong2"
								cat("\nerror: path_info must be a two column, single row data frame where:\n")
								cat("(1) the first node name must match one of the rownames of the community matrix\n")
								cat("(2) the second node name must match one of the rownames of the community matrix")
								cat("\n\n")
								return(invisible(NULL))
							}
						}
					}
					##
					##
					list_paths_Mp_vector <- unlist(lapply(list_paths_Mp,
					function(x)paste0(x, collapse = ", ")))
					##
					## number of variables in each path
					Path_var_count <- unlist(lapply(list_paths_Mp,function(x)length(x)))
					##
					## weight of each path
					Path_weight_J <- unlist(lapply(list_paths_Mp,
					function(v)prod(unlist(Map(f = function(x,y)SM[[4]][x,y],
					head(v,-1),tail(v,-1))))))
					##
					Path_weight_J <- ifelse(Path_weight_J < 0, -Path_weight_J, Path_weight_J)
					##
					## sign of each path
					Path_sign_J <- unlist(lapply(list_paths_Mp,
					function(v)prod(unlist(Map(f = function(x,y)J[x,y],
					head(v,-1), tail(v,-1))))))
					##
					## details concerning each path are stored in a data frame and printed in a text file
					Path_info <- data.frame(cbind(list_paths_Mp_vector, Path_weight_J,
					Path_var_count, Path_sign_J))
					colnames(Path_info) <- c("nodes", "path_strength", "#_nodes", "path_sign")
					##
					##
					## creation of various types of data frames (depending on "type_info" option)
					if(type_info == "all"){
						file_name_txt <- paste(fname, "_path_info_all.txt", sep = "")
						## file_name_xlsx <- paste(fname, "_path_info_all.xlsx", sep = "")
						{
						if(length(path_info) != 1){
							cat("\nwarning: when type_info = \"all\" the argument path_info is ignored\n")
							cat("type_info was reset to default option (i.e. NA)")
							cat("\n")
							}
						else{
							if(is.na(path_info) == FALSE){
								cat("\nwarning: when type_info = \"all\" the argument path_info is ignored\n")
								cat("type_info was reset to default option (i.e. NA)")
								cat("\n")
								}
							}
						}
						rownames(Path_info) <- NULL
						write.table(file = file_name_txt, Path_info, row.names = FALSE, quote = FALSE, sep = "\t")
						cat("\nproperties of all paths")
						cat("\n\n")
						print(Path_info)
						cat("\n")
						return(invisible(Path_info))
						## write.xlsx(Path_info, file = file_name_xlsx, sheetName = "all_paths", row.names = FALSE)
						cat("\ninfo on paths are in the file \"", file_name_txt, "\"", sep = "")
						cat("\n\n")
						}
					##
					##
					if(type_info == "P"){
						file_name_txt <- paste(fname, "_path_info_sign_P.txt", sep = "")
						## file_name_xlsx <- paste(fname, "_path_info_sign_P.xlsx", sep = "")
						{
						if(length(path_info) != 1){
							cat("\nwarning: when type_info = \"P\" the argument path_info is ignored\n")
							cat("type_info was reset to default option (i.e. NA)")
							cat("\n")
							}
						else{
							if(is.na(path_info) == FALSE){
								cat("\nwarning: when type_info = \"P\" the argument path_info is ignored\n")
								cat("type_info was reset to default option (i.e. NA)")
								cat("\n")
								}
							}
						}
						Paths_positive <- subset(Path_info, Path_sign_J == 1)
						rownames(Paths_positive) <- NULL
						write.table(file = file_name_txt, Paths_positive, row.names = FALSE, quote = FALSE, sep = "\t")
						cat("\nproperties of all positive paths")
						cat("\n\n")
						print(Paths_positive)
						cat("\n")
						return(invisible(Paths_positive))
						## write.xlsx(Paths_positive, file = file_name_xlsx, sheetName = "positive_paths", row.names = FALSE)
						cat("\ninfo on paths are in the file \"", file_name_txt, "\"", sep = "")
						cat("\n\n")
						}
					##
					##
					if(type_info == "N"){
						file_name_txt <- paste(fname, "_path_info_sign_N.txt", sep = "")
						## file_name_xlsx <- paste(fname, "_path_info_sign_N.xlsx", sep = "")
						{
						if(length(path_info) != 1){
							cat("\nwarning: when type_info = \"N\" the argument path_info is ignored\n")
							cat("type_info was reset to default option (i.e. NA)")
							cat("\n")
							}
						else{
							if(is.na(path_info) == FALSE){
								cat("\nwarning: when type_info = \"N\" the argument path_info is ignored\n")
								cat("type_info was reset to default option (i.e. NA)")
								cat("\n")
								}
							}
						}
						Paths_negative <- subset(Path_info, Path_sign_J == -1)
						rownames(Paths_negative) <- NULL
						write.table(file = file_name_txt, Paths_negative, row.names = FALSE, quote = FALSE, sep = "\t")
						cat("\nproperties of all negative paths")
						cat("\n\n")
						print(Paths_negative)
						cat("\n")
						return(invisible(Paths_negative))
						## write.xlsx(Paths_negative, file = file_name_xlsx, sheetName = "negative_paths", row.names = FALSE)
						cat("\ninfo on paths are in the file \"", file_name_txt, "\"", sep = "")
						cat("\n\n")
						}
					##
					##
					if(type_info == "PN"){
						file_name_p_txt <- paste(fname, "_path_info_sign_P.txt", sep = "")
						file_name_n_txt <- paste(fname, "_path_info_sign_N.txt", sep = "")
						{
						if(length(path_info) != 1){
							cat("\nwarning: when type_info = \"PN\" the argument path_info is ignored\n")
							cat("type_info was reset to default option (i.e. NA)")
							cat("\n")
							}
						else{
							if(is.na(path_info) == FALSE){
								cat("\nwarning: when type_info = \"PN\" the argument path_info is ignored\n")
								cat("type_info was reset to default option (i.e. NA)")
								cat("\n")
								}
							}
						}
						## file_name_xlsx <- paste(fname, "_path_info_sign_PN.xlsx", sep = "")
						Paths_positive <- subset(Path_info, Path_sign_J == 1)
						rownames(Paths_positive) <- NULL
						write.table(file = file_name_p_txt, Paths_positive, row.names = FALSE, quote = FALSE, sep = "\t")
						## write.xlsx(Paths_positive, file = file_name_xlsx, sheetName = "positive_paths", row.names = FALSE)
						Paths_negative <- subset(Path_info, Path_sign_J == -1)
						rownames(Paths_negative) <- NULL
						write.table(file = file_name_n_txt, Paths_negative, row.names = FALSE, quote = FALSE, sep = "\t")
						cat("\nproperties of all positive paths")
						cat("\n\n")
						print(Paths_positive)
						cat("\n\nproperties of all negative paths")
						cat("\n\n")
						print(Paths_negative)
						cat("\n")
						PN_properties <- as.list(rep(NA,2))
						PN_properties[[1]] <- Paths_positive
						PN_properties[[2]] <- Paths_negative
						names(PN_properties) <- c("positive_paths", "negative_paths")
						return(invisible(PN_properties))
						## write.xlsx(Paths_negative, file = file_name_xlsx, sheetName = "negative_paths", row.names = FALSE, append = TRUE)
						cat("\ninfo on paths are in the files \"", file_name_p_txt, "\" and \"", file_name_n_txt, "\"", sep = "")
						cat("\n\n")
						}
					##
					##
					if(type_info == "VarC"){
						if(VarC1 == "Up"){
							file_name_txt <- paste(fname, "_path_info_VarC_Up", VarC2, ".txt", sep = "")
							## file_name_xlsx <- paste(fname, "_path_info_VarC_Up", VarC2, ".xlsx", sep = "")
							Paths_above_varc <- subset(Path_info, Path_var_count >= VarC2)
							rownames(Paths_above_varc) <- NULL
							write.table(file = file_name_txt, Paths_above_varc, row.names = FALSE, quote = FALSE, sep = "\t")
							cat("\nproperties of all paths with ", VarC2, " or more nodes", sep = "")
							cat("\n\n")
							print(Paths_above_varc)
							cat("\n")
							return(invisible(Paths_above_varc))
							## write.xlsx(Paths_above_varc, file = file_name_xlsx, sheetName = "higher_than_VarC", row.names = FALSE)
							cat("\ninfo on paths are in the file \"", file_name_txt, "\"", sep = "")
							cat("\n\n")
							}
						##
						if(VarC1 == "Dn"){
							file_name_txt <- paste(fname, "_path_info_VarC_Dn", VarC2, ".txt", sep = "")
							## file_name_xlsx <- paste(fname, "_path_info_VarC_Dn", VarC2, ".xlsx", sep = "")
							Paths_below_varc <- subset(Path_info, Path_var_count <= VarC2)
							rownames(Paths_below_varc) <- NULL
							write.table(file = file_name_txt, Paths_below_varc, row.names = FALSE, quote = FALSE, sep = "\t")
							cat("\nproperties of all paths with no more than ", VarC2, " nodes", sep = "")
							cat("\n\n")
							print(Paths_below_varc)
							cat("\n")
							return(invisible(Paths_below_varc))
							## write.xlsx(Paths_below_varc, file = file_name_xlsx, sheetName = "lower_than_VarC", row.names = FALSE)
							cat("\ninfo on paths are in the file \"", file_name_txt, "\"", sep = "")
							cat("\n\n")
							}
						##
						if(VarC1 == "Eq"){
							file_name_txt <- paste(fname, "_path_info_VarC_Eq", VarC2, ".txt", sep = "")
							## file_name_xlsx <- paste(fname, "_path_info_VarC_Eq", VarC2, ".xlsx", sep = "")
							Paths_equal_varc <-	subset(Path_info, Path_var_count == VarC2)
							rownames(Paths_equal_varc) <- NULL
							write.table(file = file_name_txt, Paths_equal_varc, row.names = FALSE, quote = FALSE, sep = "\t")
							cat("\nproperties of all paths with exactly ", VarC2, " nodes", sep = "")
							cat("\n\n")
							print(Paths_equal_varc)
							cat("\n")
							return(invisible(Paths_equal_varc))
							## write.xlsx(Paths_equal_varc, file = file_name_xlsx, sheetName = "equal_to_VarC", row.names = FALSE)
							cat("\ninfo on paths are in the file \"", file_name_txt, "\"", sep = "")
							cat("\n\n")
							}
						}
					##
					##
					if(type_info == "StrgP"){
						if(StrgP1 == "Up"){
							file_name_txt <- paste(fname, "_path_info_StrgP_Up", StrgP2, ".txt", sep = "")
							## file_name_xlsx <- paste(fname, "_path_info_StrgP_Up", StrgP2, ".xlsx", sep = "")
							Paths_above_strength <- subset(Path_info, Path_weight_J >= StrgP2)
							rownames(Paths_above_strength) <- NULL
							write.table(file = file_name_txt, Paths_above_strength, row.names = FALSE, quote = FALSE, sep = "\t")
							cat("\nproperties of all paths with strength equal to or greater than ", StrgP2, sep = "")
							cat("\n\n")
							print(Paths_above_strength)
							cat("\n")
							return(invisible(Paths_above_strength))
							## write.xlsx(Paths_above_strength, file = file_name_xlsx, sheetName = "higher_than_StrgP", row.names = FALSE)
							cat("\ninfo on paths are in the file \"", file_name_txt, "\"", sep = "")
							cat("\n\n")
						}
						##
						if(StrgP1 == "Dn"){
							file_name_txt <- paste(fname, "_path_info_StrgP_Dn", StrgP2, ".txt", sep = "")
							## file_name_xlsx <- paste(fname, "_path_info_StrgP_Dn", StrgP2, ".xlsx", sep = "")
							Paths_below_strength <- subset(Path_info, Path_weight_J <= StrgP2)
							rownames(Paths_below_strength) <- NULL
							write.table(file = file_name_txt, Paths_below_strength, row.names = FALSE, quote = FALSE, sep = "\t")
							cat("\nproperties of all paths with strength equal to or less than ", StrgP2, sep = "")
							cat("\n\n")
							print(Paths_below_strength)
							cat("\n")
							return(invisible(Paths_below_strength))
							## write.xlsx(Paths_below_strength, file = file_name_xlsx, sheetName = "lower_than_StrgP", row.names = FALSE)
							cat("\ninfo on paths are in the file \"", file_name_txt, "\"", sep = "")
							cat("\n\n")
							}
						}
					##
					##
					if(type_info == "NodetN"){
						indica <- 0
						for(i in 1:length(list_paths_Mp)){
							if(list_paths_Mp[[i]][1] == NodetN1 & list_paths_Mp[[i]][length(list_paths_Mp[[i]])] == NodetN2){
								if(indica == 0){
									vector_paths_begin_end_NodetN <- i
									indica <- 1
									}
								else vector_paths_begin_end_NodetN <- c(vector_paths_begin_end_NodetN, i)
								}
							}
						l_lista_NN <- length(vector_paths_begin_end_NodetN)
						file_name_txt <- paste(fname, "_path_info_Node", NodetN1 ,"_to_Node", NodetN2,".txt", sep = "")
						## file_name_xlsx <- paste(fname, "_path_info_Node", NodetN1 ,"to_Node", NodetN2,".xlsx", sep = "")
						paths_begin_end_NodetN_vector <- rep(NA,l_lista_NN)
						giri <- 1
						for(i in 1:l_lista_NN){
							paths_begin_end_NodetN_vector[giri] <- paste0(list_paths_Mp[[vector_paths_begin_end_NodetN[i]]],
							collapse = ", ")
							giri <- giri + 1
							}
						corrisp_NN <- which(is.na(match(list_paths_Mp_vector,paths_begin_end_NodetN_vector)) == FALSE)
						Paths_begin_end_NodetN <- Path_info[corrisp_NN,]
						rownames(Paths_begin_end_NodetN) <- NULL
						write.table(file = file_name_txt, Paths_begin_end_NodetN, row.names = FALSE, quote = FALSE, sep = "\t")
						cat("\nproperties of all paths from node ", NodetN1, " to node ", NodetN2, sep = "")
						cat("\n\n")
						print(Paths_begin_end_NodetN)
						cat("\n")
						return(invisible(Paths_begin_end_NodetN))
						## write.xlsx(Paths_begin_end_NodetN, file = file_name_xlsx, sheetName = "NodetN", row.names = FALSE)
						cat("\ninfo on paths are in the file \"", file_name_txt, "\"", sep = "")
						cat("\n\n")
						}
					}
				}
			else{
				cat("\nerror: SM must be a list generated by the function \"la_strength\"")
				cat("\n\n")
				return(invisible(NULL))
				}
			}
		}
		else{
			cat("\nerror: the input must be a square numeric matrix")
			cat("\n\n")
			return(invisible(NULL))
			}
		}
	}
}
