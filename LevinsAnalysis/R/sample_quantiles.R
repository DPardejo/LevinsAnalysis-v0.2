sample_quantiles <- function(DF_TNM = NA, alternative = "two-sided", conf_level = 0.95, 
dataset_ws = "whole"){
	##
	## error message if the input is not a list object with the following nine elements:
	## [1] CM = community matrix
	## [2] RGT = random graph type
	## [3] NRG = number of random graphs generated
	## [4] TPE = total number of predictions evaluated
	## [5] TCPRE = total number of correct predictions obtained with the reference model
	## [6] TCPRG = total number of correct predictions obtained with the random graphs
	## [7] SPE = number of predictions evaluated for each scenario
	## [8] SCPRE = number of correct predictions obtained for each scenario with the reference model
	## [9] SCPRG = number of correct predictions obtained for each scenario with the random graphs
	##
	nomes_TNM <- c("CM", "RGT", "NRG", "TPE", "TCPRE", "TCPRG", "SPE", "SCPRE", "SCPRG")
	##
	a_allowed <- c("two-sided", "less", "greater")
	cl_allowed <- c(0.90, 0.95, 0.99)
	ds_allowed <- c("whole", "scenario")
	##
	V_spa <- 0
	##
	{
	if(length(names(DF_TNM)) == 0) l_ml <- 1
	else l_ml <- length(which(is.na(match(names(DF_TNM),nomes_TNM))==TRUE))
	}
	##
	{
	if(l_ml == 0){
		##
		## alternative hypothesis must be one of "two-sided" (default), "greater" or "less"
		a_length <- length(which(is.na(match(alternative, a_allowed))==FALSE))
		if(a_length != 1){
			cat("\nwarning: alternative hypothesis must be \"two-sided\" (default), \"greater\" or \"less\"\n")
			cat("alternative was reset to default option (i.e. \"two-sided\")\n\n")
			alternative <- "two-sided"
			V_spa <- 1
		}
		##
		## confidence level of the interval must be 0.90, 0.95 or 0.99
		cl_length <- length(which(is.na(match(conf_level, cl_allowed))==FALSE))
		if(cl_length != 1){
			if(V_spa == 1)cat("warning: the confidence level of the interval must be one of the following: 0.90, 0.95 or 0.99\n")
			else cat("\nwarning: the confidence level of the interval must be one of the following: 0.90, 0.95 or 0.99\n")
			cat("conf_level was reset to default option (i.e. 0.95)\n\n")
			conf_level = 0.95
			V_spa <- 1
		}
		##
		## sample quantiles must be computed either for the whole results or for each scenario
		ds_length <- length(which(is.na(match(dataset_ws, ds_allowed))==FALSE))
		if(ds_length != 1){
			if(V_spa == 1)cat("warning: the sample quantiles must be defined for either the \"whole\" results or for each \"scenario\"\n")
			else cat("\nwarning: the sample quantiles must be defined for either the \"whole\" results or for each \"scenario\"\n")
			cat("dataset_ws was reset to default option (i.e. \"whole\")\n\n")
			dataset_ws = "whole"
		}
		##
		## sample quantiles corresponding to the probability defined by the user
		{
		if(dataset_ws == "whole"){
			v_d <- DF_TNM$TCPRG
			if(length(which(DF_TNM$TCPRG==0))!=0)v_d <- v_d[-which(DF_TNM$TCPRG==0)]
			v_n <- as.integer(names(v_d))
			for(i in 1:length(v_d)){
				if(i == 1)v_h <- rep(v_n[i], v_d[i])
				else v_h <- c(v_h, rep(v_n[i], v_d[i]))
			}
			{
			if(V_spa == 1){
				cat("sample quantiles of the whole results at ", conf_level,
				" confidence interval (alternative hypothesis: ", alternative, "):\n\n", sep = "")
				}
			else {
				cat("\nsample quantiles of the whole results at ", conf_level,
				" confidence interval (alternative hypothesis: ", alternative, "):\n\n", sep = "")
				}
			}
			{
			if(alternative == "two-sided"){
				prob_lo <- (1 - conf_level)/2
				prob_hi <- conf_level + prob_lo
				v_q <- quantile(v_h, probs = c(prob_lo, prob_hi), type = 1)
				names(v_q) <- paste(c(prob_lo, prob_hi) * 100, "%", sep = "")
				print(v_q)
				cat("\n")
				return(invisible(v_q))
				}
			else{
				if(alternative == "less"){
					prob_lo <- (1 - conf_level)
					v_q <- quantile(v_h, probs = c(prob_lo), type = 1)
					names(v_q) <- paste(prob_lo * 100, "%", sep = "")
					print(v_q)
					cat("\n")
					return(invisible(v_q))
					}
				else if(alternative == "greater"){
					prob_hi <- conf_level
					v_q <- quantile(v_h, probs = c(prob_hi), type = 1)
					names(v_q) <- paste(prob_hi * 100, "%", sep = "")
					print(v_q)
					cat("\n")
					return(invisible(v_q))
					}			
				}
			}
			##
			}
		else{
			lunh <- nrow(DF_TNM$SCPRG)
			colu <- 0
			for(j in 1:lunh){
				row_sel <- DF_TNM$SCPRG[j,]
				v_d <- row_sel
				if(length(which(row_sel==0))!=0)v_d <- v_d[-which(row_sel==0)]
				v_n <- as.integer(names(v_d))
				for(i in 1:length(v_d)){
					if(i == 1)v_h <- rep(v_n[i], v_d[i])
					else v_h <- c(v_h, rep(v_n[i], v_d[i]))
				}
				{
				if(V_spa == 1 & j == 1){
					cat("sample quantiles of the scenario(s) results at ", conf_level,
					" confidence interval (alternative hypothesis: ", alternative, "):\n\n", sep = "")
					}
				else if(j == 1){
					cat("\nsample quantiles of the scenario(s) results at ", conf_level,
					" confidence interval (alternative hypothesis: ", alternative, "):\n\n", sep = "")
					}
				}
				{
				if(alternative == "two-sided"){
					prob_lo <- (1 - conf_level)/2
					prob_hi <- conf_level + prob_lo
					proba <- c(prob_lo, prob_hi)
					if(j == 1)v_q <- quantile(v_h, probs = proba, type = 1)
					else v_q <- rbind(v_q, quantile(v_h, probs = proba, type = 1))
					}
				else{
					if(alternative == "less"){
						prob_lo <- (1 - conf_level)
						proba <- prob_lo
						if(j == 1)v_q <- quantile(v_h, probs = proba, type = 1)
						else v_q <- c(v_q, quantile(v_h, probs = proba, type = 1))
						colu <- 1
						}
					else if(alternative == "greater"){
						prob_hi <- conf_level
						proba <- c(prob_hi)
						if(j == 1)v_q <- quantile(v_h, probs = proba, type = 1)
						else v_q <- c(v_q, quantile(v_h, probs = proba, type = 1))
						colu <- 1
						}			
					}
				}
			}
			if(j == 1)v_q <- matrix(v_q, nrow = 1)
			else if(colu == 1)v_q <- matrix(v_q, nrow = lunh, byrow = TRUE)
			colnames(v_q) <- paste(proba * 100, "%", sep = "")
			rownames(v_q) <- rownames(DF_TNM$SCPRG)
			print(v_q)
			cat("\n")
			return(invisible(v_q))
			}
		}	
	}
	else{
		cat("\nerror: DF_TNM must be a list generated by the function \"test_null_model\"")
		cat("\n\n")
		return(invisible(NULL))
		}
	}
}