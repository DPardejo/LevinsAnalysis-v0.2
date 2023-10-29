random_nd_graph_notext <- function(J, fname = NA){
	##
	## number of rows and columns
	n_rows <- nrow(J)
	n_cols <- ncol(J)
	##
	## total number of links
	tot_links <- length(which(J != 0))
	##
	## total number of positive links
	tot_pos <- length(which(J > 0))
	##
	## total number of negative links
	tot_neg <- length(which(J < 0))
	##
	rotazioni <- 1
	##
	repeat{
		##
		## random graph creation
		g_gnm <- erdos.renyi.game(n_rows, tot_links, type = c("gnm"),
		directed = TRUE, loops = TRUE)
		##
		## edgelist of the random graph
		edgelist_g <- as_edgelist(g_gnm)
		##
		## test to verify the random graph is connected
		{
		if(is.connected(g_gnm) == TRUE & rotazioni <= 1000){
			##
			## links in the community matrix that represents the random graph
			M_g <- matrix(rep(0,n_rows^2), nrow = n_rows)
			##
			## preparing the links for the community matrix
			{
			if(tot_pos > 0)tot_pos_ele <- rep(1, tot_pos)
			else tot_pos_ele <- NULL
			}
			##
			{
			if(tot_neg > 0)tot_neg_ele <- rep(-1, tot_neg)
			else tot_neg_ele <- NULL
			}
			##
			tot_ele <- c(tot_pos_ele, tot_neg_ele)
			l_tot_tot <- length(tot_ele)
			tot_elements <- sample(x = tot_ele, size = l_tot_tot)
			##
			for(i in 1:nrow(edgelist_g))M_g[edgelist_g[i,1],edgelist_g[i,2]] <- tot_elements[i]
			##
			{
			if(length(rownames(J)) != n_rows)rownames(M_g) <- c(1:n_rows)
			else rownames(M_g) <- rownames(J)
			}
			{
			if(length(colnames(J)) != n_cols)colnames(M_g) <- c(1:n_cols)
			else colnames(M_g) <- colnames(J)
			}
			##
			return(M_g)
			##
			break
		}
		else{
			rotazioni <- rotazioni + 1
			if(rotazioni > 1000){
				return(NULL)
				break
				}
			}
		}
	}
}
