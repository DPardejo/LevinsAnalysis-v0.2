levins_predictions_notext <- function(J, INT_MIN = NA, INT_MAX = NA, NT = 1000){
	##
	## number of rows and columns
	n_rows <- nrow(J)
	n_cols <- ncol(J)
	##
	## names of the nodes
	names_n <- colnames(J)
	##
	## defaults distribution type for random sampling
	RM = "uniform"
	##
	## transpose of J
	t_J <- t(J)
	##
	## total number of nodes
	n_nodes <- sqrt(length(J))
	##
	## total number of possible links (complete matrix with self-loops)
	n_nodes2 <- n_nodes^2
	##
	ntent <- (length(t_J[1,]) * NT)
	##
	## empty vector to store the determinant of all random matrices
	det_list <- rep(NA,ntent)
	##
	## empty list to store all random matrices that are stable
	## (in case of unstable matrices the corresponding element of the list is 0)
	weighted_list <- as.list(rep(0,ntent))
	##	
	## empty vector to record which random matrices are stable
	stable_list <- rep(0,ntent)
	##
	## empty vectors to count positive, negative or null responses
	n_plus <- as.vector(rep(0,n_nodes2), mode = "integer")
	n_min <- as.vector(rep(0,n_nodes2), mode = "integer")
	n_o <- as.vector(rep(0,n_nodes2), mode = "integer")
	##
	k <- 1
	##
	if(is.na(INT_MIN[1]) == TRUE)t_user_d_int_min <- matrix(rep(0,n_nodes2), nrow = n_nodes)
	else t_user_d_int_min <- INT_MIN
	##
	if(is.na(INT_MAX[1]) == TRUE)t_user_d_int_max <- matrix(rep(0,n_nodes2), nrow = n_nodes)
	else t_user_d_int_max <- INT_MAX
	##
	t_user_d_int_min_f <- t_user_d_int_min
	t_user_d_int_max_f <- t_user_d_int_max
	##
	for(i in 1:n_nodes){
		for(j in 1:n_nodes){
			if(t_user_d_int_min[i,j] == 0)t_user_d_int_min_f[i,j] <- 1e-6
			##
			if(t_user_d_int_max[i,j] == 0)t_user_d_int_max_f[i,j] <- 1
		}
	}
	
	
	##
	for(k in 1:ntent){
		##
		casuale <- matrix(rep(0,n_nodes2), nrow = n_nodes)
		for(i in 1:n_nodes){
			for(j in 1:n_nodes){
				casuale[i,j] <- runif(n = 1, min = t_user_d_int_min_f[i,j],
				max = t_user_d_int_max_f[i,j])
			}
		}
		num <- 0
		ww <- matrix(rep(0,n_nodes2), nrow = n_nodes)
		for(i in 1:n_nodes){
			for(j in 1:n_nodes){
				ww[i,j] <- t_J[i,j] * casuale[i,j]
			}
		}
		##
		det_ww <- round(det(ww), digits = 6)
		##
		## the criterion for stability is that the real part
		## of every eigenvalue of J must be negative
		eig_vre <- round(Re(eigen(ww)$values), digits = 6)
		##
		## the imaginary part of eigenvalues of J is not considered
		## eig_vim <- round(Im(eigen(ww)$values), digits = 6)
		##
		num <- length(which(eig_vre < 0))
		##
		if((num - n_nodes) == 0){
			##
			## if all the eigenvalues of [ww] have real part negative then the system is stable
			## here we store the determinant of stable matrices
			det_list[k] <- det(ww)
			##
			## the stable matrices are stored
			## in case of unstable matrices, the element of the list has value = 0
			weighted_list[[k]] <- ww
			##
			## vector with elements = 1 in case of stable matrices
			stable_list[k] <- 1
			##
			## counts of positive, negative or null responses
			inv_ww_vector <- as.vector(-round(ginv(ww), digits = 6))
			for(i in 1:n_nodes2){
				if(inv_ww_vector[i] > 0)(n_plus[i] <- n_plus[i] + 1)
				else if(inv_ww_vector[i] < 0)(n_min[i] <- n_min[i] + 1)
					else if(inv_ww_vector[i] == 0)(n_o[i] <- n_o[i] + 1)
			}
		}
	}
	##
	## number of stable matrices
	n_stable <- sum(stable_list)
	{
	if(n_stable > 0){
		##
		## percentages of positive, negative and null responses
		per_p <- round(matrix(c((n_plus * 100)/n_stable), nrow = n_nodes, byrow = TRUE), digits = 5)
		per_m <- round(matrix(c((n_min * 100)/n_stable), nrow = n_nodes, byrow = TRUE), digits = 5)
		per_o <- round(matrix(c((n_o * 100)/n_stable),	nrow = n_nodes, byrow = TRUE), digits = 5)
		##
		colnames(per_p) <- names_n
		rownames(per_p) <- names_n
		colnames(per_m) <- names_n
		rownames(per_m) <- names_n
		colnames(per_o) <- names_n
		rownames(per_o) <- names_n
		##
		v_p <- as.vector(per_p)
		v_m <- as.vector(per_m)
		v_o <- as.vector(per_o)
		v_d <- v_p - v_m
		##
		## table of predictions
		## the signs in the table of predictions depend on the % of
		## positive results found with the simulations
		tab <- rep(NA,n_nodes2)
		for(i in 1:n_nodes2){
			if(v_o[i] == 100)(tab[i] <- "0")
				else if(v_d[i] >= 50)(tab[i] <- "+")
					else if(v_d[i] < 50 & v_d[i] > 20)(tab[i] <- "?+")
						else if(v_d[i] <= 20 & v_d[i] >= -20)(tab[i] <- "0*")
							else if(v_d[i] < -20 & v_d[i] > -50)(tab[i] <- "?-")
								else if(v_d[i] <= -50)(tab[i] <- "-")
		}
		##
		tab_m <- matrix(c(tab), nrow = n_nodes, byrow = TRUE)
		table_of_predictions <- t(tab_m)
		colnames(table_of_predictions) <- names_n
		rownames(table_of_predictions) <- names_n	
		levins_pred <- as.matrix(noquote(table_of_predictions))
		##
		## list where results are saved
		list_results <- as.list(rep(NA,8))
		list_results[[1]] <- ntent
		list_results[[2]] <- n_stable
		list_results[[3]] <- RM
		list_results[[4]] <- J
		list_results[[5]] <- per_p
		list_results[[6]] <- per_m
		list_results[[7]] <- per_o
		list_results[[8]] <- levins_pred
		nameslist <- c("total number of matrices generated", "total number of stable matrices",
		"probability distribution", "community matrix", "(%) +", "(%) -", "(%) 0", "table of predictions")
		names(list_results) <- nameslist
		##
		}
	else{
		list_results <- NULL
		}
	}
	##
	invisible(return(list_results))
}
