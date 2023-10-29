check_positive_integers_char <- function(x){
	if(is.numeric(x) == FALSE)return(FALSE)
	else{
		ve <- 0
		{
		if(length(x) > 1){
			for(i in 1:length(x))if(x[i] != round(x[i]) | x[i] <= 0)ve <- 1
			}
		else{
			if(x != round(x) | x <= 0)ve <- 1
			}
		}
		{
		if(ve == 0)return(TRUE)
		else return(FALSE)
		}
	}
}
