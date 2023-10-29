random_nd_graph <- function(J, fname = NA){
	##
	## error message if the input is not a square numeric matrix
	if(is.matrix(J) == FALSE | is.numeric(J) == FALSE){
		cat("\nerror: the input must be a square numeric matrix")
		cat("\n\n")
	}
	else{
		## number of rows and columns
		n_rows <- nrow(J)
		n_cols <- ncol(J)
		##
		## square matrix check
		{
		if(n_rows == n_cols){
			{
			if(length(fname) > 1){
				cat("\nwarning: the name of the text file must be a single element of type \"character\" or NA\n")
				cat("the name was reset to default option (i.e. NA)\n\n")
				fname <- NA
				}
			else if(is.character(fname) == FALSE & is.na(fname) == FALSE){
				cat("\nwarning: the name of the text file must be of type \"character\" or NA\n")
				cat("the name was reset to default option (i.e. NA)\n\n")
				fname <- NA
				}
			}
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
					if(is.na(fname) == FALSE){
						file_name_txt <- paste(fname, "_random_nd_graph.txt", sep = "")
						write.table(file = file_name_txt, M_g, row.names = FALSE, 
						col.names = FALSE, quote = FALSE, sep = "\t")
						cat("\nrandom graph saved in \"", file_name_txt, "\"", sep = "")
						cat("\n\n")
					}
					##
					return(M_g)
					##
					break
				}
				else{
					rotazioni <- rotazioni + 1
					if(rotazioni > 1000){
						cat("\nerror: the generation of a connected random graph failed after 1000 attempts\n")
						cat("please, verify the numeric community matrix used as input is connected\n\n")
						invisible(NULL)
						break
						}
					}
				}
			}
		}
		else{
			cat("\nerror: the input must be a square numeric matrix")
			cat("\n\n")
			invisible(NULL)
			}	
		}
	}
}
