all_paths_raw <- function(G){
	GB <- ifelse(G != 0, 1, G)
	graphG <- graph_from_adjacency_matrix(GB, mode = c("directed"),
	weighted = NULL, diag = FALSE, add.colnames = NULL, add.rownames = NA)
	##
	pathsG <- unlist(lapply(V(graphG),
	function(x)all_simple_paths(graphG, from = x)), recursive = FALSE)
	##
	list_paths_G <- lapply(1:length(pathsG),
	function(x)as_ids(pathsG[[x]]))
	##
	invisible(list_paths_G)
}
