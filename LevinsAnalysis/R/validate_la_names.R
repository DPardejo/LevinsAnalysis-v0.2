validate_la_names <- function(J){
	##
	if(identical(rownames(J), NULL)){
		CM_Name_Val <- c(1)
		}
    {
	if("" %in% rownames(J) | NA %in% rownames(J)){
		CM_Name_Val <- c(2)
		}
    else{
		CM_Name_Val <- c(3)
		}
	}
    if(identical(colnames(J), NULL)){
		CM_Name_Val <- c(CM_Name_Val, 1)
		}
    {
	if("" %in% colnames(J) | NA %in% colnames(J)){
		CM_Name_Val <- c(CM_Name_Val, 2)
		}
	else{
		CM_Name_Val <- c(CM_Name_Val, 3)
		}
	}
    if(identical(CM_Name_Val[1], 3) & identical(CM_Name_Val[2],3) & !identical(rownames(J), colnames(J))){
		cat("\nwarning: variable names are different for rows and columns\n")
		cat("rownames are used as variable names")
		cat("\n\n")
		}
	return(CM_Name_Val)
}
