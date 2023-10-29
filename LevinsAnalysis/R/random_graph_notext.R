random_graph_notext <- function(J, fname = NA){
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
	## total number of links along the diagonal
	dia_tot <- length(which(diag(J) != 0))
	##
	## total number of positive links along the diagonal
	dia_pos <- length(which(diag(J) > 0))
	##
	## total number of negative links along the diagonal
	dia_neg <- length(which(diag(J) < 0))
	##
	## total number of links excluding the diagonal
	nondia_tot <- tot_links - dia_tot
	##
	## total number of positive links excluding the diagonal
	nondia_pos <- tot_pos - dia_pos
	##
	## total number of negative links excluding the diagonal
	nondia_neg <- tot_neg - dia_neg
	##
	rotazioni <- 1
	##
	repeat{
		##
		## random graph creation
		g_gnm <- erdos.renyi.game(n_rows, nondia_tot, type = c("gnm"),
		directed = TRUE, loops = FALSE)
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
			## preparing the links in the community matrix, except the diagonal
			{
			if(nondia_pos > 0)nondia_pos_ele <- rep(1, nondia_pos)
			else nondia_pos_ele <- NULL
			}
			##
			{
			if(nondia_neg > 0)nondia_neg_ele <- rep(-1, nondia_neg)
			else nondia_neg_ele <- NULL
			}
			##
			nondia_tot_ele <- c(nondia_pos_ele, nondia_neg_ele)
			l_nondia_tot <- length(nondia_tot_ele)
			nondia_elements <- sample(x = nondia_tot_ele, size = l_nondia_tot)
			##
			for(i in 1:nrow(edgelist_g))M_g[edgelist_g[i,1],edgelist_g[i,2]] <- nondia_elements[i]
			##
			## preparing the diagonal of the community matrix
			dia_zer <- n_rows - (dia_pos + dia_neg)
			##
			{
			if(dia_zer > 0)dia_zer_ele <- rep(0, dia_zer)
			else dia_zer_ele <- NULL
			}
			##
			{
			if(dia_pos > 0)dia_pos_ele <- rep(1, dia_pos)
			else dia_pos_ele <- NULL
			}
			##
			{
			if(dia_neg > 0)dia_neg_ele <- rep(-1, dia_neg)
			else dia_neg_ele <- NULL
			}
			##
			dia_tot_ele <- c(dia_zer_ele, dia_pos_ele, dia_neg_ele)
			l_dia_tot <- length(dia_tot_ele)
			dia_elements <- sample(x = dia_tot_ele, size = l_dia_tot)
			diag(M_g) <- dia_elements
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
