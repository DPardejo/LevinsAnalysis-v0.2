#' Signed digraph of the community matrix
#'
#' \code{cm_graph} visualizes the community matrix \code{J} as a
#'	signed digraph
#'
#' @param J the square matrix describing the qualitative direct
#'	causal effects between the variables of the community matrix.
#'	See \code{Note} for more details.
#' @param fname the name of the \code{.dot} file generated.
#'	Default name is \code{community}. See \code{Details} for
#'	more information.
#' @param color selects the color option to visualize the signed
#'	digraph: \code{"bw"} (default), \code{"color"}, or
#'	\code{"greyscale"}. See \code{Details} for more information.
#' @param save_file defines the extension for saving the image file:
#'	\code{"pdf"} (default) or \code{"svg"}
#'
#' @export
#'
#' @importFrom rsvg rsvg_svg
#' @importFrom rsvg rsvg_pdf
#' @importFrom DiagrammeRsvg export_svg
#'
#' @details \itemize{ \item \code{fname} should be used to indicate
#'	the name of the \code{.dot} file created. It should be noted that
#'	the suffix \code{"_cm_graph.dot"} is automatically added to the
#'	file name. See \code{Note} for information and \code{Examples}
#'	for practical use; \code{color} allows the user to define the
#'	color option for graphical visualization. The \code{color}
#'	argument defines which options are available. The default option
#'	(\code{"bw"}) renders a black and white signed digraph. Other
#'	options are the use of random colors (\code{color}) and greyscale
#'	(\code{"greyscale"}). See \code{Examples} for practical use.}
#'
#' @note \itemize{ \item The matrix \code{J} should be formatted
#'	as follows: \itemize{ \item qualitative direct causal links
#'	must be from \code{i} (rows) variables to \code{j} (columns)
#'	variables; \item direct causal links takes the following
#'	coefficients: 1 for positive interactions, -1 for negative
#'	interactions, and 0 for no interaction; \item the positive
#'	effect of the variable \code{i} on the variable \code{j} is
#'	represented by an arrow-headed link from \code{i} to \code{j};
#'	\item the negative effect of the variable \code{i} on the
#'	variable \code{j} is represented by a circle-headed link
#'	from \code{i} to \code{j}.} \item If the argument \code{fname}
#'	is not specified the default name \code{"community"} is assigned.
#'	\item This function is a modification of the function
#'	\code{\link{graph.cm}} from the R package \code{\link{LoopAnalyst}}.
#'	This function was created on 22-11-2017.}
#'
#' @return The \code{.dot} file of the signed digraph corresponding to
#'	the community matrix \code{J} is created and saved. Such file
#'	can be opened with the command \code{DiagrammeR::grViz}. Also
#'	the image of the signed digraph in either \code{.pdf} (default)
#'	or \code{.svg} format is generated and saved using the same name
#'	of the \code{.dot} file.
#'
#' @keywords array, color, dplot, graphs, hplot, interface, print
#'
#' @concept dot, levins, loop, network, path, pdf, svg
#'
#' @seealso \code{\link{pred_graph}}, \code{\link{M}},
#'	\code{\link{DiagrammeRsvg}}, \code{\link{LoopAnalyst}},
#'	\code{\link[LoopAnalyst]{graph.cm}}, \code{\link{rsvg}}
#'
#' @examples ## To save "community.pdf" image of the signed
#'	## digraph for the community matrix M using the default color
#'	## scheme ("bw"); the signed digraph is also saved as
#'	## "community_cm_graph.dot"
#'	cm_graph(M)
#'	DiagrammeR::grViz("community_cm_graph.dot")
#'
#'	## To visualize the signed digraph of M with the greyscale
#'	## color scheme and save both ".dot" and ".pdf" files with
#'	## the name "graphM"
#'	cm_graph(M, fname = "graphM", color = "greyscale")
#'	DiagrammeR::grViz("community_cm_graph.dot")
#'
#'	## To save the colored version of the signed digraph of M
#'	## as ".svg" file, using default name (i.e. "community")
#'	cm_graph(M, color = "color", save_file = "svg")
#'	DiagrammeR::grViz("community_cm_graph.dot")
#'
#' @author Daniel Pereira \email{dsldlf@@unife.it}, Alexis Dinno
#'	\email{alexis.dinno@@pdx.edu}, Antonio Bodini, Marco Scotti
#'	\email{marcoscot@@gmail.com}
#'
#' @references \itemize{ \item Our manual \item
#'   \href{https://CRAN.R-project.org/package=LoopAnalyst}{\code{LoopAnalyst -
#'   CRAN}} \item \href{https://alexisdinno.com/LoopAnalyst/}{LoopAnalyst -
#'   Software}}
#'
#'
la_graph <- function(J, fname = "community", layout_g = "neato", coord = NA, save_file = "pdf",
	colors_nodes = NULL, colors_edges = NULL, col_option = NA){
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
			contatore <- 0
			{
			if(is.character(fname) == FALSE){
				cat("\nwarning: the name of the .dot and image files must be of type \"character\"\n")
				cat("the name was reset to default option (i.e. \"community\")\n\n")
				fname <- "community"
				contatore <- 1
				}
			else if(is.character(fname) == TRUE & length(fname) > 1){
				cat("\nwarning: the name of the .dot and image files must be a single element of type \"character\"\n")
				cat("the name was reset to default option (i.e. \"community\")\n\n")
				fname <- "community"
				contatore <- 1
				}
			}
			file_name_dot <- paste(fname, "_la_graph.dot", sep = "")
			##
			layout_g_options <- c("neato", "user")
			{
			if(length(layout_g) > 1){
				{
				if(contatore == 1)cat("warning: only one of the following graph layout options must be provided: \"neato\" or \"user\"\n")
				else cat("\nwarning: only one of the following graph layout options must be provided: \"neato\" or \"user\"\n")
				}
				cat("the graph layout was reset to default option (i.e. \"neato\")\n\n")
				layout_g <- "neato"
				contatore <- 1
				}
			else{
				layout_g_match <- which(is.na(match(layout_g_options, layout_g)) == FALSE)
				layout_g_match_l <- length(layout_g_match)
				if(layout_g_match_l == 0){
					{
					if(contatore == 1)cat("warning: only one of the following graph layout options must be provided: \"neato\" or \"user\"\n")
					else cat("\nwarning: only one of the following graph layout options must be provided: \"neato\" or \"user\"\n")
					}
					cat("the graph layout was reset to default option (i.e. \"neato\")\n\n")
					layout_g <- "neato"
					contatore <- 1
					}
				}
			}
			##
			if(layout_g == "user"){
				if(is.matrix(coord) == FALSE){
					{
					if(contatore == 1)cat("warning: user-defined coordinates must be provided through a rectangular numeric matrix with:\n")
					else cat("\nwarning: user-defined coordinates must be provided through a rectangular numeric matrix with:\n")
					}
					cat("(1) number of rows = number of community matrix rows\n")
					cat("(2) number of columns = 2\n")
					cat("(i.e. first column for x-values and second column for y-values)\n")
					cat("the graph layout was reset to default option (i.e. \"neato\")\n\n")
					layout_g <- "neato"
					coord <- NA
					contatore <- 1
					}
				else{
					if(nrow(coord) != nrow(J) | is.numeric(coord) == FALSE){
					{
					if(contatore == 1)cat("warning: user-defined coordinates must be provided through a rectangular numeric matrix with:\n")
					else cat("\nwarning: user-defined coordinates must be provided through a rectangular numeric matrix with:\n")
					}
					cat("(1) number of rows = number of community matrix rows\n")
					cat("(2) number of columns = 2\n")
					cat("(i.e. first column for x-values and second column for y-values)\n")
					cat("the graph layout was reset to default option (i.e. \"neato\")\n\n")
					layout_g <- "neato"
					coord <- NA
					contatore <- 1
					}
				}
			}
			##
			count_colors <- 0
			##
			if(length(colors_nodes) != 0){
				if(is.vector(colors_nodes) == FALSE | length(colors_nodes) != nrow(J) | is.numeric(colors_nodes) == FALSE){
					{
					if(contatore == 1)cat("warning: colors_nodes must be a numeric vector with length as the number of community matrix rows\n")
					else cat("\nwarning: colors_nodes must be a numeric vector with length as the number of community matrix rows\n")
					}
					cat("all node and label colors were reset to default value (i.e. \"black\")\n\n")
					colors_nodes <- NULL
					contatore <- 1
				}
				else count_colors <- 1
			}
			##
			count_colors_edges <- 0
			if(length(colors_edges) != 0){
				if(is.matrix(colors_edges) == TRUE & is.numeric(colors_edges) == TRUE){
					if(nrow(colors_edges) == nrow(J) & ncol(colors_edges) == ncol(J)){
						count_colors <- 1
						count_colors_edges <- 1
					}
					else{
						{
						if(contatore == 1)cat("warning: colors_edges must be a numeric square matrix with same size as the community matrix\n")
						else cat("\nwarning: colors_edges must be a numeric square matrix with same sizes as the community matrix\n")
						}
						cat("all edge colors were reset to default value (i.e. \"black\")\n\n")
						colors_edges <- NULL
						contatore <- 1
					}
				}
				else{
					{
					if(contatore == 1)cat("warning: colors_edges must be a numeric square matrix with same size as the community matrix\n")
					else cat("\nwarning: colors_edges must be a numeric square matrix with same sizes as the community matrix\n")
					}
					cat("all edge colors were reset to default value (i.e. \"black\")\n\n")
					colors_edges <- NULL
					contatore <- 1
				}
			}
			##
			## grey and green as possible color options
			col_option_all <- c("grey", "green")
			if(length(col_option) == 1 & is.na(col_option) == FALSE){
				if(length(which(col_option_all == col_option)) == 1){
					{
					if(col_option == "green")col_option_g <- "#009E73"
					else if(col_option == "grey")col_option_g <- "#C0C0C0"
					}
					count_colors <- count_colors + 1
				}
				else{
					{
					if(contatore == 1)cat("warning: the color used for nodes and edges must be either \"grey\" or \"green\"\n")
					else cat("\nwarning: the color used for nodes and edges must be either \"grey\" or \"green\"\n")
					}
					cat("all colors were reset to default value (i.e. \"black\")\n\n")
					col_option <- NA
					contatore <- 1
				}
			}
			##
			## construction of the graph for visualization
			if(length(colnames(J)) == 0 & length(rownames(J)) == 0){
				colnames(J) <- rownames(J) <- c(1:n_rows)
			}
			match_names <- length(which(is.na(match(colnames(J), rownames(J)))==TRUE))
			if(match_names != 0){
				{
				if(contatore == 1)cat("warning: row names and column names of the community matrix must be the same\n")
				else cat("\nwarning: row names and column names of the community matrix must be the same\n")
				}
				cat("the names were set to \"1\" ... \"", n_rows, "\"\n\n", sep = "")
				colnames(J) <- rownames(J) <- c(1:n_rows)
				contatore <- 1
			}
			tot_links <- length(which(J != 0))
			col1 <- rep(NA,tot_links)
			col2 <- rep(NA,tot_links)
			signs <- rep(NA,tot_links)
			conta <- 1
			{
			if(count_colors == 2){
				colors_nodes_v <- rep("black", n_rows)
				if(length(which(colors_nodes != 0))!=0)colors_nodes_v[which(colors_nodes != 0)] <- col_option_g
				##
				colors_edges_v <- rep("black", tot_links)
				for(i in 1:n_rows){
					for(j in 1:n_cols){
						if(J[i,j] > 0){
							col1[conta] <- i
							col2[conta] <- j
							signs[conta] <- "normal"
							if(count_colors_edges == 1){
								if(colors_edges[i,j] != 0)colors_edges_v[conta] <- col_option_g
							}
							conta <- conta + 1
						}
						else if(J[i,j] < 0){
							col1[conta] <- i
							col2[conta] <- j
							signs[conta] <- "odot"
							if(count_colors_edges == 1){
								if(colors_edges[i,j] != 0)colors_edges_v[conta] <- col_option_g
							}
							conta <- conta + 1
							}
						}
					}	
				}
			else{
				for(i in 1:n_rows){
					for(j in 1:n_cols){
						if(J[i,j] > 0){
							col1[conta] <- i
							col2[conta] <- j
							signs[conta] <- "normal"
							conta <- conta + 1
							}
						else if(J[i,j] < 0){
							col1[conta] <- i
							col2[conta] <- j
							signs[conta] <- "odot"
							conta <- conta + 1
							}
						}
					}
				}
			}
			##
			## edgelist and unique set of nodes
			df_g <- data.frame(from = col1, to = col2, stringsAsFactors = FALSE)
			nodes_u <- colnames(J)
			nodes_u_n <- unique(c(df_g$from, df_g$to))
			##
			{
			if(count_colors == 2){
				{
				if(layout_g == "user"){
					node_df <- create_node_df(n = n_rows, label = nodes_u, fontcolor = colors_nodes_v,
					style = "filled", fillcolor = "white", color = colors_nodes_v, shape = "circle",
					x = coord[,1], y = coord[,2])
					}
				else{
					node_df <- create_node_df(n = n_rows, label = nodes_u, fontcolor = colors_nodes_v,
					style = "filled", fillcolor = "white", color = colors_nodes_v, shape = "circle")
					}
				}
				##
				## edgelist dataframe
				edge_df <- create_edge_df(from = df_g$from, to = df_g$to,
				color = colors_edges_v, arrowhead = signs)
				}
			else{
				##
				## node dataframe
				{
				if(layout_g == "user"){
					node_df <- create_node_df(n = n_rows, label = nodes_u, fontcolor = "black",
					style = "filled", fillcolor = "white", color = "black", shape = "circle",
					x = coord[,1], y = coord[,2])
					}
				else{
					node_df <- create_node_df(n = n_rows, label = nodes_u, fontcolor = "black",
					style = "filled", fillcolor = "white", color = "black", shape = "circle")
					}
				}
				##
				## edgelist dataframe
				edge_df <- create_edge_df(from = df_g$from, to = df_g$to,
				color = "black", arrowhead = signs)
				}
			}
			##
			## creation of the graph object
			g <- create_graph(nodes_df = node_df, edges_df = edge_df)
			##
			## .dot file
			sink(file = file_name_dot)
			cat(generate_dot(g))
			sink()
			##
			## extension for saving the image file
			all_ext <- c("pdf","svg")
			{
			if(length(save_file) > 1 | is.vector(save_file) == FALSE){
				{
				if(contatore == 0)cat("\nwarning: the format of the image file must be one of the following: \"pdf\" or \"svg\"")
				else cat("warning: the format of the image file must be one of the following: \"pdf\" or \"svg\"")
				}
				cat("\n")
				cat("the format was reset to default (i.e. \"pdf\")")
				cat("\n\n")
				save_file <- "pdf"
				}
			else{
				s_ext <- which(all_ext == save_file)
				if(length(s_ext) != 1){
					{
					if(contatore == 0)cat("\nwarning: the format of the image file must be one of the following: \"pdf\" or \"svg\"")
					else cat("warning: the format of the image file must be one of the following: \"pdf\" or \"svg\"")
					}
					cat("\n")
					cat("the format was reset to default (i.e. \"pdf\")")
					cat("\n\n")
					save_file <- "pdf"
					}
				}
			}
			##
			## produce an svg file in the working directory
			{
			if(save_file == "svg"){
				file_name1 <- paste(fname, "_la_graph.", save_file, sep = "")
				## g %>% export_graph(file_name = file_name1, file_type = save_file)
				rsvg::rsvg_svg(charToRaw(DiagrammeRsvg::export_svg(grViz(file_name_dot))), file = file_name1)
				}
			else{
				##
				## produce a pdf file in the working directory
				if(save_file == "pdf"){
					file_name2 <- paste(fname, "_la_graph.", save_file, sep = "")
					## g %>% export_graph(file_name = file_name2, file_type = save_file)
					rsvg::rsvg_pdf(charToRaw(DiagrammeRsvg::export_svg(grViz(file_name_dot))), file = file_name2)
					}
				}
			}
			##
			## graph visualization in RStudio
			return(render_graph(g))
			}
		else{
			cat("\nerror: the input must be a square numeric matrix")
			cat("\n\n")
			invisible(NULL)
			}
		}
	}
}
