levins_concurrent_ns <- function(OUT_LP, DF = NA, fname = NA){
	##
	## error message if the input is not a list object with the following eight elements:
	## [1] total number of matrices generated
	## [2] total number of stable matrices
	## [3] probability distribution
	## [4] community matrix
	## [5] (%) +
	## [6] (%) -
	## [7] (%) 0
	## [8] table of predictions
	##
	nomilista <- c("total number of matrices generated", "total number of stable matrices",
	"probability distribution", "community matrix", "(%) +", "(%) -", "(%) 0", "table of predictions")
	{
	if(length(names(OUT_LP)) == 0) l_ml <- 1
	else l_ml <- length(which(is.na(match(names(OUT_LP),nomilista))==TRUE))
	}
	##
	{
	if(is.list(OUT_LP) == TRUE){
		if(l_ml != 0){
			cat("\nerror: OUT_LP must be a list generated by the function \"levins_predictions\"")
			cat("\n\n")
			invisible(NULL)
		}
		##
		## concurrent predictions
		else{
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
			if(is.matrix(DF) == FALSE & is.vector(DF) == FALSE){
			DF <- "wrong"
			}
			##
			if(is.matrix(DF) == FALSE & is.vector(DF) == TRUE){
				if(length(DF) != 2)DF <- "wrong"
			}
			##
			if(is.matrix(DF) == TRUE & is.vector(DF) == FALSE){
				if(ncol(DF) != 2)DF <- "wrong"
			}
			##
			l_sgn <- 0
			##
			if(length(DF) == 2 & is.vector(DF) == TRUE & is.character(DF) == TRUE){
				DF <- matrix(DF, nrow = 1)
			}
			##
			DFin <- DF
			##
			## remove DF rows if:
			## (1) characters other than "+" and "-" are in the first column
			## (2) names different from node names (or integer numbers larger
			## than network size) are in the second column
			if(is.list(DF) == FALSE & length(DF) > 1){
				if(is.matrix(DF) == TRUE & is.character(DF) == TRUE & ncol(DF) == 2){
					remo1 <- which(is.na(match(DF[,1],c("+","-"))) == TRUE)
					remo2 <- which(is.na(match(DF[,2],names_n)) == TRUE)
					##
					remo_n1 <- which(is.na(suppressWarnings(as.numeric(DF[,2]))) == FALSE)
					{
					if(length(remo_n1) != 0){
						num_to_nod <- which(is.na(match(c(1:n_nodes),remo_n1)) == FALSE)
						if(length(num_to_nod) != 0){
							num_to_nod_sele <- as.numeric(DF[num_to_nod,2])
							nod_from_num <- names_n[num_to_nod_sele]
							remo_n2 <- which(is.na(match(c(DF[-num_to_nod,2],nod_from_num),names_n)) == TRUE)
							remo <- unique(c(remo1,remo_n2))
							DF[num_to_nod,2] <- nod_from_num
							}
						}
					else remo <- unique(c(remo1,remo2))
					}
					{
					if(length(remo) == nrow(DF)){
						DF <- NA
						}
					else{
						if(length(remo) > 0 & length(remo) < nrow(DF)){
							DF <- DF[-remo,]
							}
						}
					}
				}
			}
			##
			if(is.vector(DF) == TRUE & length(DF) == 2)DF <- matrix(DF, nrow = 1)
			##
			{
			if(is.list(DF) == FALSE & length(DF) > 1){
				if(is.matrix(DF) == TRUE & is.character(DF) == TRUE & ncol(DF) == 2){
					in_sgn <- rep(1,nrow(DF))
					l_sgn <- length(in_sgn)
					nega <- which(DF[,1] == "-")
					if(length(nega) != 0)in_sgn[which(DF[,1] == "-")] <- -1
					in_var <- DF[,2]
					##
					{
					if(l_sgn == 1){
						conc_name <- paste0("\npredictions for input on node ", in_var[1])
						}
						else if(l_sgn == 2){
							conc_name <- paste0("\nconcurrent predictions for inputs on nodes ", in_var[1], " and ", in_var[2])
							}
							else{
								conc_name1 <- paste(in_var[c(1:(l_sgn-1))], collapse = ", ")
								conc_name2 <- paste0(conc_name1, " and ", in_var[l_sgn])
								conc_name <- paste0("\nconcurrent predictions for inputs on nodes ", conc_name2)
								}
					}
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
					cat(conc_name)
					cat("\n\n")
					print(levins_conc_list)
					cat("\n")
					##
					## error message in case of some non-existing target nodes in DF (i.e. DFin)
					if(l_sgn != 0){
						if((nrow(levins_conc_list)-2) < nrow(DFin)){
							cat("error: inputs with signs different from \"+\" and \"-\" or not targeting network nodes cannot be considered")
							cat("\n\n")
						}
					}
					##
					## the results are saved in a text file
					{
					if(is.character(fname) == TRUE){
						file_name_txt <- paste(fname, "_levins_concurrent_ns.txt", sep = "")
						sink(file_name_txt)
						cat(conc_name)
						cat("\n\n")
						print(levins_conc_list)
						sink()
						cat(paste("results are in the file \"", file_name_txt, "\"", sep = ""))
						cat("\n\n")
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
					}
					##
					## invisible object for the results
					invisible(levins_conc_list)
				}
				else{
					levins_conc_list <- NA
					cat("\nerror: the input for concurrent predictions must be a two-column matrix of type \"character\"\n")
					cat("first column with input signs (i.e. \"+\" or \"-\") and second column with node names")
					cat("\n\n")
					##
					## invisible object for the results
					invisible(NULL)
					}
				}
			else{
				levins_conc_list <- NA
				cat("\nerror: the input for concurrent predictions must be a two-column matrix of type \"character\"\n")
				cat("first column with input signs (i.e. \"+\" or \"-\") and second column with node names")
				cat("\n\n")
				##
				## invisible object for the results
				invisible(NULL)
				}
			}
		}
	}
	else{
		cat("\nerror: OUT_LP must be a list generated by the function \"levins_predictions\"")
		cat("\n\n")
		invisible(NULL)
		}
	}
}
