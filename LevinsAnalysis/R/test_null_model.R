test_null_model <- function(OUT_LP1 = NA, DF_list1 = NA, MV_target1 = NA, r_model = "r_graph",
n_graphs = 999, NTS = 1000, fname1 = "community"){
	##
	## model validation
	results_s11 <- model_validation(OUT_LP = OUT_LP1, DF_list = DF_list1, MV_target = MV_target1)
	##
	{
	if(length(results_s11) != 0){
		##
		## community matrix
		JJ <- OUT_LP1[[4]]
		##
		## number of rows
		n_rows <- nrow(JJ)
		ll_1 <- length(results_s11)
		results_s1 <- results_s11[[ll_1]]
		results_s11 <- list.remove(results_s11, range = ll_1)
		ll_2 <- length(results_s11)
		##
		## target for qualitative predictions
		MV_target2 <- results_s11[[ll_2]]
		##
		## scenarios' inputs
		DF_list2 <- list.remove(results_s11, range = ll_2)
		##
		n_prediz <- nrow(results_s1)
		results_s1[which(is.na(results_s1)==TRUE)] <- 0
		tot_pred_model <- apply(results_s1,2,sum)[2]
		max_pred_model <- apply(results_s1,2,max)[2]
		tot_correct_pred_model <- apply(results_s1,2,sum)[1]
		V_spa <- 0
		##
		## type of random graph to be used
		names_r_graphs <- c("r_graph", "r_nd_graph", "r_symm_graph", "r_foodweb_graph")
		names_r_graphs_m <- length(which(is.na(match(r_model, names_r_graphs))==FALSE))
		if(names_r_graphs_m != 1){
			cat("\nwarning: the type of random graph must be indicated using one of the following options:\n")
			cat("\"r_graph\", \"r_nd_graph\", \"r_symm_graph\" or \"r_foodweb_graph\"\n")
			cat("r_model was reset to default option (i.e. \"r_graph\")\n\n")
			r_model <- "r_graph"
			V_spa <- 1
		}
		##
		## check for the consistency concerning the type of random graph generation
		JJs <- JJ
		JJs[which(JJ < 0)] <- 1
		if(r_model == "r_symm_graph"){
			if(isSymmetric(JJs) == FALSE){
				cat("\nwarning: the community matrix is asymmetric and the \"r_symm_graph\" option cannot be used\n")
				cat("r_model was reset to default option (i.e. \"r_graph\")\n\n")
				r_model <- "r_graph"
				V_spa <- 1
			}
		}
		JSND <- JJ
		diag(JSND) <- rep(0, n_rows)
		sum_tJ <- JSND + t(JSND)
		l_sum_tJ <- length(which(sum_tJ!=0))
		if(r_model == "r_foodweb_graph"){
			if(l_sum_tJ != 0){
				cat("\nwarning: the community matrix does not represent predator-prey relations and the \"r_foodweb_graph\" option cannot be used\n")
				cat("r_model was reset to default option (i.e. \"r_graph\")\n\n")
				r_model <- "r_graph"
				V_spa <- 1
			}
		}
		##
		## number of random graphs generated
		number_r_graphs <- c(9, 99, 999, 9999)
		number_r_graphs_m <- length(which(is.na(match(n_graphs, number_r_graphs))==FALSE))
		if(number_r_graphs_m != 1){
			if(V_spa == 1)cat("warning: the number of random graphs generated must be one of the following: 99, 999 or 9999\n")
			else cat("\nwarning: the number of random graphs generated must be one of the following: 99, 999 or 9999\n")
			cat("the number was reset to default option (i.e. 999)\n\n")
			n_graphs <- 999
			V_spa <- 1
		}
		##
		## number of simulation runs with levins_predictions_notext
		number_NTS_graphs <- c(10, 100, 1000)
		number_NTS_graphs_m <- length(which(is.na(match(NTS, number_NTS_graphs))==FALSE))
		if(number_NTS_graphs_m != 1){
			if(V_spa == 1)cat("warning: the number of runs for levins_predictions must be one of the following: 10, 100 or 1000\n")
			else cat("\nwarning: the number of runs for levins_predictions must be one of the following: 10, 100 or 1000\n")
			cat("NTS was reset to default option (i.e. 1000)\n\n")
			NTS <- 1000
			V_spa <- 1
		}
		##
		M_rnd_test <- matrix(rep(0,(n_graphs * n_prediz)), nrow = n_graphs)
		colnames(M_rnd_test) <- rownames(results_s1)
		##
		## test with random graphs
		for(i in 1:n_graphs){
			{
			if(r_model == "r_graph")J_rnd <- random_graph_notext(JJ)
			else if(r_model == "r_nd_graph")J_rnd <- random_nd_graph_notext(JJ)
				else if(r_model == "r_symm_graph")J_rnd <- random_symm_graph_notext(JJ)
					else if(r_model == "r_foodweb_graph") J_rnd <- random_foodweb_graph_notext(JJ)
			}
			##
			{
			if(length(J_rnd) == 0){
				if(V_spa == 1)cat("error: the generation of a connected random graph failed after 1000 attempts\n")
				else cat("\nerror: the generation of a connected random graph failed after 1000 attempts\n")
				cat("please, verify the numeric community matrix used to get OUT_LP1 is connected\n\n")
				return(invisible(NULL))
				break
				}
			else{
				OUT_LP2 <- levins_predictions_notext(J_rnd, NT = NTS)
				{
				if(length(OUT_LP2) != 0){
					results_s22 <- model_validation_notext(OUT_LP = OUT_LP2,
					DF_list = DF_list2, MV_target = MV_target2)
					##
					tt_1 <- length(results_s22)
					results_s2 <- results_s22[[tt_1]]
					for(j in 1:n_prediz)M_rnd_test[i,j] <- results_s2[j,1]
					}
				else{
					if(V_spa == 1)cat("error: none of the the random graphs constructed with r_model = \"", r_model, "\" was stable\n", sep = "")
					else cat("\nerror: none of the the random graphs constructed with r_model = \"", r_model, "\" was stable\n", sep = "")
					cat("please, consider using a different option for r_model (i.e. to generate random graphs)\n\n")
					return(invisible(NULL))
					break
					}
				}
				##
				}
			}
		}
		M_rnd_test_F <- rbind(M_rnd_test,results_s1[,1])
		tot_pred_random <- apply(M_rnd_test_F,1,function(x)sum(x, na.rm = TRUE))
		V_ref <- c(0:tot_pred_model)
		V_max <- c(0:max_pred_model)
		tot_pred_random_hist <- unlist(lapply(V_ref,function(x)length(which(tot_pred_random == x))))
		names(tot_pred_random_hist) <- V_ref
		for(j in 1:n_prediz){
			if(j == 1)MAT_pred_hist <- matrix(unlist(lapply(V_max,function(x)length(
			which(M_rnd_test_F[,j] == x)))), nrow = 1)
			else MAT_pred_hist <- rbind(MAT_pred_hist,
			unlist(lapply(V_max,function(x)length(which(M_rnd_test_F[,j] == x)))))
		}
		scenarios_n <- unlist(lapply(rownames(results_s1),function(x)strsplit(x,"DF_list_")[[1]][2]))
		rownames(MAT_pred_hist) <- scenarios_n
		colnames(MAT_pred_hist) <- c(0:max_pred_model)
		##
		results_s1M <- matrix(results_s1[,2], nrow = 1)
		colnames(results_s1M) <- scenarios_n
		rownames(results_s1M) <- ""
		results_s2M <- matrix(results_s1[,1], nrow = 1)
		colnames(results_s2M) <- scenarios_n
		rownames(results_s2M) <- ""
		##
		## list where results are saved
		results_list <- as.list(rep(NA,9))
		results_list[[1]] <- JJ
		results_list[[2]] <- r_model
		results_list[[3]] <- n_graphs
		names(tot_pred_model) <- NULL
		results_list[[4]] <- tot_pred_model
		names(tot_correct_pred_model) <- NULL
		results_list[[5]] <- tot_correct_pred_model
		results_list[[6]] <- tot_pred_random_hist
		results_list[[7]] <- results_s1M
		results_list[[8]] <- results_s2M
		results_list[[9]] <- MAT_pred_hist 
		nameslist <- c("CM", "RGT", "NRG", "TPE", "TCPRE", "TCPRG", "SPE", "SCPRE", "SCPRG")
		names(results_list) <- nameslist
		##
		cat("\ncommunity matrix")
		cat("\n\n")
		print(results_list[[1]])
		cat("\n")
		cat("random graph type: \"", results_list[[2]], "\"", sep = "")
		cat("\n\n")
		cat("number of random graphs generated: ", results_list[[3]], sep = "")
		cat("\n\n")
		cat("total number of predictions evaluated: ", results_list[[4]], sep = "")
		cat("\n\n")
		cat("total number of correct predictions obtained with the reference model: ", results_list[[5]], sep = "")
		cat("\n\n")
		cat("total number of correct predictions obtained with the random graphs:\n")
		print(results_list[[6]])
		cat("\n")
		cat("number of predictions evaluated for each scenario:\n")
		print(results_list[[7]], rownames = NULL)
		cat("\n")
		cat("number of correct predictions obtained for each scenario with the reference model:\n")
		print(results_list[[8]], rownames = NULL)
		cat("\n")
		cat("number of correct predictions obtained for each scenario with the random graphs:\n")
		print(results_list[[9]])
		cat("\n")
		##
		{
		if(length(fname1) == 1){
			if(is.character(fname1) == TRUE){
				file_name1_txt <- paste(fname1, "_test_null_model.txt", sep = "")
				sink(file_name1_txt)
				cat("community matrix")
				cat("\n")
				print(results_list[[1]])
				cat("\n")
				cat("random graph type: ", results_list[[2]], sep = "")
				cat("\n\n")
				cat("number of random graphs generated: ", results_list[[3]], sep = "")
				cat("\n\n")
				cat("total number of predictions evaluated: ", results_list[[4]], sep = "")
				cat("\n\n")
				cat("total number of correct predictions obtained with the reference model: ", results_list[[5]], sep = "")
				cat("\n\n")
				cat("total number of correct predictions obtained with the random graphs:\n")
				print(results_list[[6]])
				cat("\n")
				cat("number of predictions evaluated for each scenario:\n")
				print(results_list[[7]])
				cat("\n")
				cat("number of correct predictions obtained for each scenario with the reference model:\n")
				print(results_list[[8]])
				cat("\n")
				cat("number of correct predictions obtained for each scenario with the random graphs:\n")
				print(results_list[[9]])
				sink()
				cat("results are in the file \"", file_name1_txt, "\"", sep = "")
				cat("\n\n")
				}
			else{
				fname1 <- "community"
				file_name1_txt <- paste(fname1, "_test_null_model.txt", sep = "")
				sink(file_name1_txt)
				cat("community matrix")
				cat("\n")
				print(results_list[[1]])
				cat("\n")
				cat("random graph type: ", results_list[[2]], sep = "")
				cat("\n\n")
				cat("number of random graphs generated: ", results_list[[3]], sep = "")
				cat("\n\n")
				cat("total number of predictions evaluated: ", results_list[[4]], sep = "")
				cat("\n\n")
				cat("total number of correct predictions obtained with the reference model: ", results_list[[5]], sep = "")
				cat("\n\n")
				cat("total number of correct predictions obtained with the random graphs:\n")
				print(results_list[[6]])
				cat("\n")
				cat("number of predictions evaluated for each scenario:\n")
				print(results_list[[7]])
				cat("\n")
				cat("number of correct predictions obtained for each scenario with the reference model:\n")
				print(results_list[[8]])
				cat("\n")
				cat("number of correct predictions obtained for each scenario with the random graphs:\n")
				print(results_list[[9]])
				sink()
				cat("warning: the filename must a single element be of type \"character\"\n")
				cat("the filname was reset to default and results saved in the file \"", file_name1_txt, "\"", sep = "")
				cat("\n\n")
				}
			}
		else{
			fname1 <- "community"
			file_name1_txt <- paste(fname1, "_test_null_model.txt", sep = "")
			sink(file_name1_txt)
			cat("community matrix")
			cat("\n")
			print(results_list[[1]])
			cat("\n")
			cat("random graph type: ", results_list[[2]], sep = "")
			cat("\n\n")
			cat("number of random graphs generated: ", results_list[[3]], sep = "")
			cat("\n\n")
			cat("total number of predictions evaluated: ", results_list[[4]], sep = "")
			cat("\n\n")
			cat("total number of correct predictions obtained with the reference model: ", results_list[[5]], sep = "")
			cat("\n\n")
			cat("total number of correct predictions obtained with the random graphs:\n")
			print(results_list[[6]])
			cat("\n")
			cat("number of predictions evaluated for each scenario:\n")
			print(results_list[[7]])
			cat("\n")
			cat("number of correct predictions obtained for each scenario with the reference model:\n")
			print(results_list[[8]])
			cat("\n")
			cat("number of correct predictions obtained for each scenario with the random graphs:\n")
			print(results_list[[9]])
			sink()
			cat("warning: the filename must a single element be of type \"character\"\n")
			cat("the filname was reset to default and results saved in the file \"", file_name1_txt, "\"", sep = "")
			cat("\n\n")
			}
		}
		##
		return(invisible(results_list))
	}
	else{
		cat("error: there are no results from \"model_validation\" that can be tested with random graphs")
		cat("\n\n")
		return(invisible(NULL))
		}
	}
}
