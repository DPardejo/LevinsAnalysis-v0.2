#'Levins Concurrent Prediction
#'
#'Computes a table of predictions of the direction of change of the state of the
#'Community Network matrix \code{J} variables due to concurrent effects of
#'multiple inputs in variables of the same matrix.
#'
#'@param J a square matrix describing the qualitative direct casual effect
#'  between the variables of a community network. See \code{Note} for more
#'  details.
#'@param DF vector containing the sign and nodes upon which the effect of
#'  concurrent input over the Community Network matrix \code{J} is to be
#'  investigated. See \code{Details} for more information on use.
#'@param RM the method to be used for random distribution for the interaction
#'  strength values: \code{"uniform"} (default), \code{"normal"},
#'  \code{"pareto"}. See \code{Details} for more information on use.
#'@param INT_MIN matrix defining the minimum value to be used by random
#'  distribution assignment of interaction strength values of each variable
#'  pair. The default value \code{NA} defines the minimum interaction strength
#'  value as \code{1e-6}. See \code{Details} for more information on use..
#'@param INT_MAX matrix defining the maximum value to be used for the random
#'  distribution assignment of interaction strength values of each variable
#'  pair. The default value \code{NA} defines the maximum interaction strength
#'  value as \code{1}. See \code{Details} for more information on use.
#'@param NT defines the number of random interaction strength matrices to be
#'  randomly generated by the function. The default value is \code{1000}. See
#'  \code{Details} for more details.
#'@param fname the name for the exporting file. By default is \code{NA} and no
#'  file is exported. See \code{Details} for more information on use.
#'
#'@export
#'@importFrom MASS ginv
#'@importFrom msm rtnorm
#'@importFrom stats runif
#'
#'@details \itemize{ \item Through the argument \code{DF} the function allows
#'  the user to indicate the direction and input nodes of a perturbation which
#'  will be used to compute the table of predictions of the direction of change
#'  in the state of the Community Network matrix \code{J} variables. The
#'  function allows this through the use of a vector as follows: \itemize{ \item
#'  odd elements indicate the sign of the input. Should be indicated as \code{1}
#'  for positive input and \code{-1} for negative input; \item even elements
#'  indicate the variable of the Community Network matrix \code{J} upon which
#'  the input is acting; \item e.g. \code{DF = c(1,"A", -1, "E")}; \item see
#'  \code{Examples} for pratical use.} \item Through the argument \code{RM} the
#'  user can define the type of random distribution method to be used on the
#'  random assignment of values of interaction strength for the variable pairs
#'  of the Community Matrix \code{J}, by choosing the appropriate option as
#'  follows: \itemize{ \item \code{"uniform"} random assignment method follows a
#'  Uniform probability distribution (default); \item \code{"normal"} random
#'  assignment method follows a Normal probability distribution; \item
#'  \code{"pareto"} random assignment method follows a Pareto probability
#'  distribution; \item see \code{Note} for more information and \code{Examples}
#'  for practical use.} \item Through the use of matrices for the arguments
#'  \code{INT_MIN} and \code{INT_MAX} the user can change the default interval
#'  values for the random assignment ([1e-6, 1]). This can be done by use of
#'  square matrices of the same length and writing format as the Community
#'  Network matrix \code{J}, as follows: \itemize{ \item \code{INT_MIN} the
#'  user defined values for each variable pair will be used by the process as
#'  the minimum value in the interval for random assignment ([x; 1]); \item
#'  \code{INT_MAX} the user defined values for each variable pair will be used
#'  by the process as the maximum value in the interval for random assignment
#'  ([1e-6; x]); \item the user defined values in \code{INT_MIN} and
#'  \code{INT_MAX} should be within the interval (0, 1]; \item to define an
#'  exact interaction strength value for any given variable pair (and not an
#'  interval), the user should indicate same value in both matrices; \item for
#'  variable pairs which the user does not pretend to change the default
#'  interval values, fill with value \code{0}; \item see \code{Note} for more
#'  information and \code{Examples} for practical use.} \item To compute a
#'  \code{interaction Strength Matrix} the program will compute \code{x} number
#'  of random matrices which will then be screened for their stability, with
#'  those accepeted being used to compute an average interaction strength
#'  matrix. Through the use of the argument \code{NT} the user can define the
#'  number of runs of the process. The user defined number of runs follows the
#'  equation: \code{Number of Runs = (Nº of variables of J * NT)}. See
#'  \code{Examples} for pratical use; \item \code{fname} should be used giving
#'  the name intended for the exporting files noting that a suffix
#'  \code{"_levins_conc_pred.txt"} will be added to the file name. See
#'  \code{Examples} for pratical use.}
#'
#'@note \itemize{ \item The \code{J} matrix should be formatted as follows:
#'  \itemize{ \item should be written by indicating the qualitative direct
#'  casual effect from j (rows) variables on i (columns) variables; \item the
#'  signs of qualitative direct casual effect allowed are: \code{1} for
#'  positive, \code{-1} for negative and \code{0} for null.} \item Regarding the
#'  \code{INT_MIN} and \code{INT_MAX} matrices: \itemize{ \item the matrices
#'  should be written by indicating the interaction strength of the direct
#'  casual effect from j (rows) variables on i (columns) variables; \item it
#'  should be noted that these arguments only work with the argument \code{ RM =
#'  "uniform"}. If a user defined matrix is passed on either \code{INT_MIN} or
#'  \code{INT_MAX} and a different option other than the default is chosen for
#'  the argument \code{RM}, the interval considered for the value assignment
#'  will be the default [1e-6; 1].}}
#'
#'@return It returns the total number of random matrices generated, number of
#'  stable matrices used for the calculations and Levins Table of Concurrent
#'  Predictions.
#'
#'  If \code{fname} is indicated, the results will be exported to a
#'  \code{fname.txt} file.
#'
#'  An \code{invisible} list with the results is generated. It can be assigned
#'  to an object.
#'
#'
#'@keywords arith, array, distribution, math, print
#'@concept concurrent, interaction, levins, loop, prediction, network, sign,
#'  strength
#'
#'@seealso \code{\link{cm_strength}}, \code{\link{cm_stability}},
#'  \code{\link{levins_prediction}}, \code{\link{M}}, \code{\link{pred_graph}}
#'
#' @examples ## To compute the table of predictions for two positive concurrent effects on "A" and "E", without printing the results to a .txt file
#' levins_concurrent(M, DF = c(1,"A", 1, "E"), RM = "uniform", INT_MIN = NA, INT_MAX = NA, NT = 100, fname = NA)
#' ## or (to assign the invisible results list)
#' results <- levins_concurrent(M, DF = c(1,"A", 1, "E"), RM = "uniform", INT_MIN = NA, INT_MAX = NA, NT = 1000, fname = NA)
#'
#' ## To compute the table of predictions for two positive concurrent effects on "A" and "E", printing the results to a .txt file
#' levins_concurrent(M, DF = c(1,"A", 1, "E"), RM = "uniform", INT_MIN = NA, INT_MAX = NA, NT = 100, fname = "M")
#'
#' ## To compute the table of predictions for three positive concurrent effects on "A", "C" and "E", printing the results to a .txt file
#' levins_concurrent(M, DF = c(1,"A", 1, "C", 1, "E"), RM = "uniform", INT_MIN = NA, INT_MAX = NA, NT = 100, fname = "M")
#'
#' ## To compute the table of predictions for two opposite signed concurrent effects on "A" and "E"
#' levins_concurrent(M, DF = c(1,"A", -1, "E"), RM = "uniform", INT_MIN = NA, INT_MAX = NA, NT = 100, fname = NA)
#'
#' ## To compute the table of predictions for two opposite signed concurrent effects on "A" and "E", using a Normal distribution
#' levins_concurrent(M, DF = c(1,"A", -1, "E"), RM = "normal", INT_MIN = NA, INT_MAX = NA, NT = 100, fname = NA)
#'
#' ## To compute the table of predictions for two opposite signed concurrent effects on "A" and "E", using a Pareto distribution
#' levins_concurrent(M, DF = c(1,"A", -1, "E"), RM = "pareto", INT_MIN = NA, INT_MAX = NA, NT = 100, fname = NA)
#'
#' ## To compute the table of predictions for two opposite signed concurrent effects on "A" and "E", using a user defined minimum strength interval value
#' matrix_min<-matrix(rep(0,(length(M))), nrow = sqrt(length(M)))
#' matrix_min[1, 5]<- 0.3
#' levins_concurrent(M, DF = c(1,"A", -1, "E"), RM = "pareto", INT_MIN = matrix_min, INT_MAX = NA, NT = 100, fname = NA)
#'
#' ## To compute the table of predictions for two opposite signed concurrent effects on "A" and "E", using a user defined maximum strength interval value
#' matrix_max<-matrix(rep(0,(length(M))), nrow = sqrt(length(M)))
#' matrix_max[1, 5]<- 0.7
#' levins_concurrent(M, DF = c(1,"A", -1, "E"), RM = "pareto", INT_MIN = NA, INT_MAX = matrix_max, NT = 100, fname = NA)
#'
#' ## To compute the table of predictions for two opposite signed concurrent effects on "A" and "E", where [1, 5] pair has user defined fixed strength value
#' matrix_min<-matrix(rep(0,(length(M))), nrow = sqrt(length(M)))
#' matrix_min[1, 5]<- 0.4
#' matrix_max<-matrix(rep(0,(length(M))), nrow = sqrt(length(M)))
#' matrix_max[1, 5]<- 0.4
#' levins_concurrent(M, DF = c(1,"A", -1, "E"), RM = "pareto", INT_MIN = matrix_min, INT_MAX = matrix_max, NT = 100, fname = NA)
#'
#' ## To compute the table of predictions for two opposite signed concurrent effects on "A" and "E" with user defined number of runs
#' levins_concurrent(M, DF = c(1,"A", -1, "E"), RM = "uniform", INT_MIN = NA, INT_MAX = NA, NT = 300, fname = NA)
#'
#'@author Daniel Pereira \email{dsldlf@@unife.it}, Marta Rocchi, Antonio Bodini,
#'  Marco Scotti \email{marcoscot@@gmail.com}
#'
#'@references  \itemize{ \item Our manual; \item and something from Loop
#'  analysis}
#'
#'
#'
#'
#'
levins_concurrent <- function(J, DF = NA, RM = "uniform", INT_MIN = NA, INT_MAX = NA, NT = 1000, fname = NA){
	##
	## error message if the input is not a square numeric matrix
	if(is.matrix(J) == FALSE | is.numeric(J) == FALSE){
		cat("\nerror: the input must be a square numeric matrix")
		cat("\n\n")
		invisible(NULL)
	}
	else{
		## number of rows and columns
		n_rows <- nrow(J)
		n_cols <- ncol(J)
		##
		## square matrix check
		{
		if(n_rows == n_cols){
			##
			## initial data for concurrent inputs
			DF_ini <- DF
			##
			## names of the nodes
			names_n <- colnames(J)
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
			## total number of runs = number of nodes * NT
			{
			if(is.numeric(NT) == FALSE){
				NT <- 1000
				cat("\nwarning: NT must be strictly positive and integer and was reset to NT = 1000")
				cat("\n")
				}
			else{
				if(is.list(NT) == TRUE | is.matrix(NT) == TRUE | length(NT) > 1){
				NT <- 1000
				cat("\nwarning: NT must be strictly positive and integer and was reset to NT = 1000")
				cat("\n")
				}
				else if(check_positive_integers(NT) == FALSE){
					NT <- 1000
					cat("\nwarning: NT must be strictly positive and integer and was reset to NT = 1000")
					cat("\n")
					}
				}
			}
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
			## check of the type of probability distribution
			all_distr <- c("uniform","normal","pareto")
			##
			## empty vectors to count positive, negative or null responses
			n_plus <- as.vector(rep(0,n_nodes2), mode = "integer")
			n_min <- as.vector(rep(0,n_nodes2), mode = "integer")
			n_o <- as.vector(rep(0,n_nodes2), mode = "integer")
			##
			{
			if(length(RM) != 1){
				RM <- "uniform"
				cat("\nwarning: RM must be one of the three following probability distributions: \"uniform\", \"normal\" or \"pareto\"\n")
				cat("RM was reset to default option (i.e. RM = \"uniform\")")
				cat("\n")
				}
			else{
				if(length(which(all_distr == RM)) != 1){
					RM <- "uniform"
					cat("\nwarning: RM must be one of the three following probability distributions: \"uniform\", \"normal\" or \"pareto\"\n")
					cat("RM was reset to default option (i.e. RM = \"uniform\")")
					cat("\n")
					}
				}
			}
			##
			## check of the format of the matrices that define maximum and minimum thresholds for random sampling
			conteggio <- 0
			{
			if(is.matrix(INT_MIN) == TRUE & is.numeric(INT_MIN) == TRUE){
				if(nrow(INT_MIN) == n_nodes & ncol(INT_MIN) == n_nodes){
					if(max(INT_MIN) <= 1 & min(INT_MIN) >= 0){
						t_user_d_int_min <- t(INT_MIN)
						if(RM == "normal" | RM == "pareto"){
							conteggio <- 1
							cat("\nwarning: user-defined thresholds for random sampling only apply to \"uniform\" distribution\n")
							cat("RM was reset to default option (i.e. RM = \"uniform\")")
							cat("\n\n")
							RM <- "uniform"
							}
						}
					else{
						for(i in 1:n_nodes){
							for(j in 1:n_nodes){
								if(INT_MIN[i,j] > 1 | INT_MIN[i,j] < 0)INT_MIN[i,j] <- 0
								}
							}
							t_user_d_int_min <- t(INT_MIN)
							if(RM == "normal" | RM == "pareto"){
								conteggio <- 1
								cat("\nwarning: user-defined thresholds for random sampling only apply to \"uniform\" distribution\n")
								cat("RM was reset to default option (i.e. RM = \"uniform\")")
								cat("\n")
								RM <- "uniform"
							}
							cat("\nwarning: all values of INT_MIN must be in the interval (0, 1]\n")
							cat("values outside such interval were reset to default option (i.e. 1e-6)")
							cat("\n")
						}
					}
				else{
					t_user_d_int_min <- matrix(rep(0,n_nodes2), nrow = n_nodes)
					if(RM == "normal" | RM == "pareto"){
						conteggio <- 1
						cat("\nwarning: user-defined thresholds for random sampling only apply to \"uniform\" distribution\n")
						cat("RM was reset to default option (i.e. RM = \"uniform\")")
						cat("\n")
						RM <- "uniform"
					}
					cat("\nwarning: INT_MIN must be a square matrix with same size as J\n")
					cat("INT_MIN was reset to default size with all values = 1e-6")
					cat("\n\n")
					}
				}
			else{
				if(is.na(INT_MIN[1]) == TRUE){
					t_user_d_int_min <- matrix(rep(0,n_nodes2), nrow = n_nodes)
					}
				else{
					t_user_d_int_min <- matrix(rep(0,n_nodes2), nrow = n_nodes)
					if(RM == "normal" | RM == "pareto"){
						conteggio <- 1
						cat("\nwarning: user-defined thresholds for random sampling only apply to \"uniform\" distribution\n")
						cat("RM was reset to default option (i.e. RM = \"uniform\")")
						cat("\n")
						RM <- "uniform"
					}
					cat("\nwarning: INT_MIN must be a square matrix with same size as J\n")
					cat("INT_MIN was reset to default size with all values = 1e-6")
					cat("\n\n")
					}
				}
			}
			##
			{
			if(is.matrix(INT_MAX) == TRUE & is.numeric(INT_MAX) == TRUE){
				if(nrow(INT_MAX) == n_nodes & ncol(INT_MAX) == n_nodes){
					if(max(INT_MAX) <= 1 & min(INT_MAX) >= 0){
						t_user_d_int_max <- t(INT_MAX)
						if(RM == "normal" | RM == "pareto"){
							if(conteggio == 0){
								cat("\nwarning: user-defined thresholds for random sampling only apply to \"uniform\" distribution\n")
								cat("RM was reset to default option (i.e. RM = \"uniform\")")
								cat("\n")
								RM <- "uniform"
								}
							}
						}
					else{
						for(i in 1:n_nodes){
							for(j in 1:n_nodes){
								if(INT_MAX[i,j] > 1 | INT_MAX[i,j] < 0)INT_MAX[i,j] <- 1
								}
							}
							t_user_d_int_max <- t(INT_MAX)
							if(RM == "normal" | RM == "pareto"){
								if(conteggio == 0){
									cat("\nwarning: user-defined thresholds for random sampling only apply to \"uniform\" distribution\n")
									cat("RM was reset to default option (i.e. RM = \"uniform\")")
									cat("\n")
									RM <- "uniform"
								}
							}
							cat("\nwarning: all values of INT_MAX must be in the interval (0, 1]\n")
							cat("values outside such interval were reset to default option (i.e. 1)")
							cat("\n")
						}
					}
				else{
					t_user_d_int_max <- matrix(rep(0,n_nodes2), nrow = n_nodes)
					if(RM == "normal" | RM == "pareto"){
						if(conteggio == 0){
							cat("\nwarning: user-defined thresholds for random sampling only apply to \"uniform\" distribution\n")
							cat("RM was reset to default option (i.e. RM = \"uniform\")")
							cat("\n")
							RM <- "uniform"
						}
					}
					cat("\nwarning: INT_MAX must be a square matrix with same size as J\n")
					cat("INT_MAX was reset to default size with all values = 1")
					cat("\n\n")
					}
				}
			else{
				if(is.na(INT_MAX[1]) == TRUE){
					t_user_d_int_max <- matrix(rep(0,n_nodes2), nrow = n_nodes)
					}
				else{
					t_user_d_int_max <- matrix(rep(0,n_nodes2), nrow = n_nodes)
					if(RM == "normal" | RM == "pareto"){
						if(conteggio == 0){
							cat("\nwarning: user-defined thresholds for random sampling only apply to \"uniform\" distribution\n")
							cat("RM was reset to default option (i.e. RM = \"uniform\")")
							cat("\n")
							RM <- "uniform"
						}
					}
					cat("\nwarning: INT_MAX must be a square matrix with same size as J\n")
					cat("INT_MAX was reset to default size with all values = 1")
					cat("\n\n")
					}
				}
			}
			##
			matrice_somma_user <- t_user_d_int_min + t_user_d_int_max
			indicatore <- sum(matrice_somma_user)
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
			## check there are no entries in the matrix of minimum thresholds that exceed those in the matrix of maximum thresholds
			coun <- 0
			for(i in 1:n_nodes){
				for(j in 1:n_nodes){
					if(t_user_d_int_min_f[i,j] > t_user_d_int_max_f[i,j]){
						t_user_d_int_min_f[i,j] <- 1e-6
						t_user_d_int_max_f[i,j] <- 1
						{
						if(coun == 0){
							cat("\nwarning: each entry in INT_MIN must not exceed the corresponding one in INT_MAX\n")
							cat("the following entries were reset to default values (INT_MIN = 1e-6, INT_MAX = 1):\n")
							cat(rownames(J)[j], colnames(J)[i], "\n", sep = " ")
							coun <- 1
							}
						else{
							cat(rownames(J)[j], colnames(J)[i], "\n", sep = " ")
							}
						}
					}
				}
			}
			##
			##
			k <- 1
			##
			##
			for(k in 1:ntent){
				##
				casuale <- matrix(rep(0,n_nodes2), nrow = n_nodes)
				##
				## uniform distribution
				if(RM == "uniform"){
					if(indicatore == 0){
						for(i in 1:n_nodes){
							for(j in 1:n_nodes){
								casuale[i,j] <- runif(n = 1, min = 1e-6, max = 1)
							}
						}
					}
					##
					else{
						for(i in 1:n_nodes){
							for(j in 1:n_nodes){
								casuale[i,j] <- runif(n = 1, min = t_user_d_int_min_f[i,j],
								max = t_user_d_int_max_f[i,j])
							}
						}
					}
				}
				##
				## normal distribution
				if(RM == "normal"){
					for(i in 1:n_nodes){
						for(j in 1:n_nodes){
							casuale[i,j] <- rtnorm(1, mean = 0.5, sd = 0.2,
							lower = 1e-6, upper = 1)
						}
					}
				}
				##
				## Pareto distribution
				if(RM == "pareto"){
					qpareto <- function(u, shape = 1, location = 1)
					location/(1 - u)^(1/shape);
					rpareto <- function(n, shape = 1, location = 1)
					qpareto(runif(n), shape, location);
					for(i in 1:n_nodes){
						for(j in 1:n_nodes){
							repeat{
								vu <- rpareto(1, shape = 8.8) - 1
								if(vu >= 1e-6 & vu <= 1){
									break
								}
							}
							casuale[i,j] <- vu
						}
					}
				}
				##
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
				## concurrent predictions
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
				if(length(DF) > 1){
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
				if(length(DF) > 1){
					if(is.matrix(DF) == TRUE & is.character(DF) == TRUE & ncol(DF) == 2){
						in_sgn <- rep(1,nrow(DF))
						l_sgn <- length(in_sgn)
						nega <- which(DF[,1] == "-")
						if(length(nega) != 0)in_sgn[which(DF[,1] == "-")] <- -1
						in_var <- DF[,2]
						##
						{
						if(l_sgn == 1){
							conc_name <- paste0("predictions for input on node ", in_var[1])
							}
							else if(l_sgn == 2){
								conc_name <- paste0("concurrent predictions for inputs on nodes ", in_var[1], " and ", in_var[2])
								}
								else{
									conc_name1 <- paste(in_var[c(1:(l_sgn-1))], collapse = ", ")
									conc_name2 <- paste0(conc_name1, " and ", in_var[l_sgn])
									conc_name <- paste0("concurrent predictions for inputs on nodes ", conc_name2)
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
							if(in_sgn[i] > 0)levins_conc[i,] <- round(per_p[in_var[i],] - per_m[in_var[i],], 3)
							else levins_conc[i,] <- round(per_m[in_var[i],] - per_p[in_var[i],], 3)
							##
							levins_zeros[i,] <- round(per_o[in_var[i],], 3)
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
					}
					else levins_conc_list <- NA
					}
				else levins_conc_list <- NA
				}
				##
				## list where results are saved
				list_results <- as.list(rep(NA,9))
				list_results[[1]] <- ntent
				list_results[[2]] <- n_stable
				list_results[[3]] <- RM
				list_results[[4]] <- J
				list_results[[5]] <- per_p
				list_results[[6]] <- per_m
				list_results[[7]] <- per_o
				list_results[[8]] <- levins_pred
				list_results[[9]] <- levins_conc_list
				{
				if(l_sgn == 0){
					nameslist <- c("total number of matrices generated", "total number of stable matrices",
					"probability distribution", "community matrix", "(%) +", "(%) -", "(%) 0",
					"table of predictions", "no concurrent predictions")
					}
					else{
						if(l_sgn == 1){
							nameslist <- c("total number of matrices generated", "total number of stable matrices",
							"probability distribution", "community matrix", "(%) +", "(%) -", "(%) 0",
							"table of predictions", conc_name) ## "predictions generated with a single input")
						}
						else{
							nameslist <- c("total number of matrices generated", "total number of stable matrices",
							"probability distribution", "community matrix", "(%) +", "(%) -", "(%) 0",
							"table of predictions", conc_name) ## "concurrent predictions")
						}
					}
				}
				names(list_results) <- nameslist
				##
				cat(paste("\ntotal number of matrices generated: ", list_results[[1]], sep = ""))
				cat("\n\n")
				cat(paste("total number of stable matrices: ", list_results[[2]], sep = ""))
				cat("\n\n")
				cat(paste("the probability distribution used is:", RM, sep = " "))
				cat("\n\n\n")
				cat("community matrix")
				cat("\n\n")
				print(list_results[[4]])
				cat("\n\n")
				cat("(%) + ")
				cat("\n\n")
				print(list_results[[5]])
				cat("\n\n")
				cat("(%) - ")
				cat("\n\n")
				print(list_results[[6]])
				cat("\n\n")
				cat("(%) 0 ")
				cat("\n\n")
				print(list_results[[7]])
				cat("\n\n")
				cat("table of predictions")
				cat("\n\n")
				print(list_results[[8]])
				{
				if(l_sgn != 0)cat("\n\n")
				else cat("\n")
				}
				{
				if(l_sgn != 0){
					if(l_sgn == 1){
						## cat("predictions generated with a single input")
						cat(conc_name)
						cat("\n\n")
						print(list_results[[9]])
						cat("\n")
						}
					else{
						## cat("concurrent predictions")
						cat(conc_name)
						cat("\n\n")
						print(list_results[[9]])
						cat("\n")
						}
					}
				else{
					if(length(DF) == 1){
						if(length(DF_ini) == 1){
							if(is.na(DF_ini) == FALSE){
								cat("error: the input for concurrent predictions must be a two-column matrix of type \"character\"\n")
								cat("first column with input signs (i.e. \"+\" or \"-\") and second column with node names")
								cat("\n\n")
								}
							else{
								cat("error: no concurrent predictions because no inputs were specified by the user")
								cat("\n\n")
								}
							}
						else{
							cat("error: the input for concurrent predictions must be a two-column matrix of type \"character\"\n")
							cat("first column with input signs (i.e. \"+\" or \"-\") and second column with node names")
							cat("\n\n")
							}
						}
					else{
						cat("error: the input for concurrent predictions must be a two-column matrix of type \"character\"\n")
						cat("first column with input signs (i.e. \"+\" or \"-\") and second column with node names")
						cat("\n\n")
						}
					}
				}
				##
				## error message in case of some non-existing target nodes in DF (i.e. DFin)
				if(l_sgn != 0){
					if((nrow(list_results[[9]])-2) < nrow(DFin)){
						cat("error: inputs with signs different from \"+\" and \"-\" or not targeting network nodes cannot be considered")
						cat("\n\n")
					}
				}
			}
			else{
				cat("\nerror: none of the matrices generated was stable")
				cat("\n\n")
				list_results <- NULL
				}
			}
			##
			## the list of results is saved in a text file
			{
			if(is.character(fname) == TRUE){
				file_name_txt <- paste(fname, "_levins_concurrent.txt", sep = "")
				if(length(list_results) != 0){
					sink(file_name_txt)
					print(list_results)
					sink()
					cat(paste("results are in the file \"", file_name_txt, "\"", sep = ""))
					cat("\n\n")
					}
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
			##
			}
			##
			## invisible object for the list of results
			invisible(list_results)
		}
		else{
			cat("\nerror: the input must be a square numeric matrix")
			cat("\n\n")
			invisible(NULL)
			}
		}
	}
}
