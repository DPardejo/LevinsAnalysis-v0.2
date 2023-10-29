model_validation_notext <- function(OUT_LP = NA, DF_list = NA, MV_target = NA){
	names_predi <- colnames(OUT_LP[[8]])
	n_nodes <- length(names_predi)
	names(DF_list) <- c(1:length(DF_list))
	##
	for(m in 1:length(DF_list)){
		if(length(DF_list[[m]]) == 2 & is.vector(DF_list[[m]]) == TRUE){
			DF_list[[m]] <- matrix(DF_list[[m]], nrow = 1)
		}
		remo_n1 <- which(is.na(suppressWarnings(as.numeric(DF_list[[m]][,2]))) == FALSE)
		if(length(remo_n1) != 0){
			num_to_nod <- which(is.na(match(c(1:n_nodes),remo_n1)) == FALSE)
			if(length(num_to_nod) != 0){
				num_to_nod_sele <- as.numeric(DF_list[[m]][num_to_nod,2])
				nod_from_num <- names_predi[num_to_nod_sele]
				remo_n2 <- which(is.na(match(c(DF_list[[m]][-num_to_nod,2],nod_from_num),
				names_predi)) == TRUE)
				DF_list[[m]][num_to_nod,2] <- nod_from_num
			}
		}
	}
	##
	if(is.vector(MV_target) == TRUE){
		nomi_target <- names(MV_target)
		MV_target <- matrix(MV_target, nrow = 1)
		colnames(MV_target) <- nomi_target
	}
	##
	if(is.matrix(MV_target) == TRUE){
		nomi_target <- colnames(MV_target)
		rownames(MV_target) <- c(1:nrow(MV_target))
		nomi_match <- match(names_predi, nomi_target)
		nomi_target_ele <- which(is.na(nomi_match)==TRUE)
		nomi_target_l <- length(nomi_target_ele)
	}
	##
	ll_list <- length(DF_list)
	matches_pe <- matrix(rep(0,(ll_list * 2)), nrow = ll_list)
	##
	for(p in 1:length(DF_list)){
		ris_s1 <- levins_concurrent_ns_notext(OUT_LP, DF = DF_list[[p]])
		ris_s1_sign <- ris_s1[nrow(ris_s1),]
		##
		ris_s1_sign_p <- which(ris_s1_sign == "?+")
		ris_s1_sign_n <- which(ris_s1_sign == "?-")
		ris_s1_sign_z <- which(ris_s1_sign == "0*")
		##
		if(length(ris_s1_sign_p) != 0)ris_s1_sign[ris_s1_sign_p] <- "+"
		##
		if(length(ris_s1_sign_n) != 0)ris_s1_sign[ris_s1_sign_n] <- "-"
		##
		if(length(ris_s1_sign_z) != 0)ris_s1_sign[ris_s1_sign_z] <- "0"
		##					
		counting <- 0
		ris_s1_sign_sel <- ris_s1_sign[names_predi[which(is.na(nomi_match)==FALSE)]]
		con_s1_sign_sel <- MV_target[p,nomi_target[nomi_match[which(is.na(nomi_match)==FALSE)]]]
		for(v in 1:length(ris_s1_sign_sel)){
			if(is.na(con_s1_sign_sel[v]) == FALSE)
				if(ris_s1_sign_sel[v] == con_s1_sign_sel[v])counting <- counting + 1
		}
		matches_pe[p,1] <- counting
		matches_pe[p,2] <- length(which(is.na(con_s1_sign_sel)==FALSE))
	}
	##
	if(is.vector(matches_pe) == TRUE)matches_pe <- matrix(matches_pe, nrow = 1)
	colnames(matches_pe) <- c("correct", "total")
	rownames(matches_pe) <- names(DF_list)
	nn_r <- nrow(matches_pe)
	nn_n <- rownames(matches_pe)
	nn_r_LISTA <- as.list(rep(NA,nn_r))
	numbero_lett_V <- rep(NA,nn_r)
	##
	for(x in 1:nn_r){
		numbero <- as.integer(nn_n[x])
		numbero_lett <- as.character(numbero)
		levins_concurrent_ns_notext(OUT_LP, DF = DF_list[[numbero_lett]])
		nn_r_LISTA[[x]] <- DF_list[[numbero_lett]]
		numbero_lett_V[x] <- numbero_lett
	}
	##
	if(is.vector(MV_target)==FALSE){
		sleli <- which(is.na(match(rownames(MV_target), numbero_lett_V))==FALSE)
		MV_target <- MV_target[sleli,]
	}
	nn_r_LISTA[[(length(nn_r_LISTA)+1)]] <- MV_target
	nn_r_LISTA[[(length(nn_r_LISTA)+1)]] <- matches_pe
	names(nn_r_LISTA) <- c(rownames(matches_pe), "target", "validation")
	return(invisible(nn_r_LISTA))
}
