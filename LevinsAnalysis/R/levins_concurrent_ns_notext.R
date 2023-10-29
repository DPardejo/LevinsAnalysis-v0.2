levins_concurrent_ns_notext <- function(OUT_LP, DF = NA){
	##
	## names of the nodes
	names_n <- colnames(OUT_LP[[8]])
	##
	## total number of nodes
	n_nodes <- length(names_n)
	##
	## total number of possible links (complete matrix with self-loops)
	n_nodes2 <- n_nodes^2
	##
	## percentages of positive responses
	per_p <- OUT_LP[[5]]
	##
	## percentages of negative responses
	per_m <- OUT_LP[[6]]
	##
	## percentages of (real) null responses (i.e. 0)
	per_z <- OUT_LP[[7]]
	##
	l_sgn <- 0
	##
	if(length(DF) == 2 & is.vector(DF) == TRUE)DF <- matrix(DF, nrow = 1)
	##
	## remove DF rows if names different from node names (or integer
	## numbers larger than network size) are in the second column
	remo_n1 <- which(is.na(suppressWarnings(as.numeric(DF[,2]))) == FALSE)
	if(length(remo_n1) != 0){
		num_to_nod <- which(is.na(match(c(1:n_nodes),remo_n1)) == FALSE)
		if(length(num_to_nod) != 0){
			num_to_nod_sele <- as.numeric(DF[num_to_nod,2])
			nod_from_num <- names_n[num_to_nod_sele]
			remo_n2 <- which(is.na(match(c(DF[-num_to_nod,2],nod_from_num),names_n)) == TRUE)
			remo <- remo_n2
			DF[num_to_nod,2] <- nod_from_num
		}
	}
	##
	in_sgn <- rep(1,nrow(DF))
	l_sgn <- length(in_sgn)
	nega <- which(DF[,1] == "-")
	if(length(nega) != 0)in_sgn[which(DF[,1] == "-")] <- -1
	in_var <- DF[,2]
	##
	levins_conc <- matrix(rep(0, n_nodes), (l_sgn + 2), n_nodes)
	colnames(levins_conc) <- names_n
	rownames(levins_conc) <- c(in_var, "mean", "signs")
	##
	levins_zeros <- matrix(rep(0, n_nodes), (l_sgn + 1), n_nodes)
	colnames(levins_zeros) <- names_n
	rownames(levins_zeros) <- c(in_var, "mean")
	##
	for(i in 1:l_sgn){
		{
		if(in_sgn[i] > 0)levins_conc[i,] <- round(per_p[in_var[i],] - per_m[in_var[i],], 3)
		else levins_conc[i,] <- round(per_m[in_var[i],] - per_p[in_var[i],], 3)
		}
		##
		levins_zeros[i,] <- round(per_z[in_var[i],], 3)
	}
	##
	{
	if(l_sgn == 1){
		levins_conc["mean",] <- levins_conc[in_var,]
		rownames(levins_conc)[1] <- paste(DF[1], DF[2], sep = "")
		##
		levins_zeros["mean",] <- levins_zeros[in_var,]
		}
	else{
		levins_conc["mean",] <- round(apply(levins_conc[in_var,], 2, mean),3)
		rownames(levins_conc)[1:l_sgn] <- apply(DF,1,function(x)paste(x[1], x[2], sep = ""))
		##
		levins_zeros["mean",] <- round(apply(levins_zeros[in_var,], 2, mean),3)
		}
	}
	##
	avrg_n <- as.numeric(levins_conc["mean",])
	sg_n <- rep(NA, n_nodes)
	avrg_z <- as.numeric(levins_zeros["mean",])
	for(v in 1:n_nodes){
		if(avrg_z[v] == 100)(sg_n[v] <- "0")
			else if(avrg_n[v] >= 50)(sg_n[v] <- "+")
				else if(avrg_n[v] < 50 & avrg_n[v] > 20)(sg_n[v] <- "?+")
					else if(avrg_n[v] <= 20 & avrg_n[v] >= -20)(sg_n[v] <- "0*")
						else if(avrg_n[v] < -20 & avrg_n[v] > -50)(sg_n[v] <- "?-")
							else if(avrg_n[v] <= -50)(sg_n[v] <- "-")
	}
	levins_conc["signs",] <- sg_n		
	levins_conc_list <- as.matrix(noquote(levins_conc))
	##
	## invisible object for the results
	invisible(levins_conc_list)
}
