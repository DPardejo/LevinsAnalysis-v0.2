#' Community Network Stability
#'
#' Tests the stability of the Community Network.
#'
#' @param J a square matrix describing the qualitative direct casual effect
#'   between the variables of a Community Ntwork. See \code{Note} for more
#'   details.
#' @export
#' @details The function tests if the Community Network matrix  \code{J} is stable by
#'   verifying if all eigen values are negative. See \code{Note} for detail.
#'
#' @note \itemize{ \item The \code{J} matrix should be formatted as follows:
#'   \itemize{ \item should be written by indicating the qualitative direct
#'   casual effect from j (rows) variables on i (columns) variables; \item the
#'   signs of qualitative direct casual effect allowed are: \code{1} for
#'   positive, \code{-1} for negative and \code{0} for null.} \item Regarding
#'   Community Network matrix stability: \itemize{ \item one requisite for  Community Network matrix stability is
#'   the sign of the the determinant. For matrices wih even number of variables,
#'   the \code{determinant} must be positive, while for matrices with odd number of
#'   variables the \code{determinant} must be negative; \item a second requisite for
#'   stability is for all \code{eigen values} to be negative; for this to be
#'   true in case of even numbered matrices the \code{determinant} must be
#'   positive and in case of odd numbered matrices the \code{determinant} must be negative.}}
#' @return A print indicating the result of the test: if the Community Network matrix \code{J} is
#'   stable or unstable.
#' @keywords array, math
#' @concept network, stability
#' @seealso \code{\link{cm_structure}}, \code{\link{J}}
#' @examples ## If M is stable:
#' ## cm_stability(M)
#' ## [1] The community matrix is stable
#'
#' ## If M is not stable
#' ## cm_stability(M)
#' ## [1] The Community Network is unstable
#'
#' @author Stefania Favilla (check with Antonio if its just this), Antonio Bodini
#' @references \itemize{ \item Our manual; \item Puccia, C. J. and Levins, R. (1986) Qualitative Modeling of
#'   Complex Systems: An Introduction to Loop Analysis and Time Averaging.
#'   Cambridge: Harvard University Press.}
#'
#'
la_stability <- function(J){
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
			##
			## the criterion for stability is that the real part
			## of every eigenvalue of J must be negative
			eig_CM_vre <- round(Re(eigen(J)$values), digits = 6)
			##
			## the imaginary part of eigenvalues of J is not considered
			## eig_CM_vim <- round(Im(eigen(J)$values), digits = 6)
			##
			nn <- length(eig_CM_vre)
			cnt <- length(which(eig_CM_vre < 0))
			{
			if(cnt < nn){
				cat("\nthe community matrix is unstable")
				cat("\n\n")
				invisible(FALSE)
				}
			else{
				cat("\nthe community matrix is stable")
				cat("\n\n")
				invisible(TRUE)
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
