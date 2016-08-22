

convert_to_stanreg <- function(object, m, ...) {

	named_elements <- c("coefficients",
	                    "ses",
						"fitted.values",
						"linear.predictors",
						"residuals",
						"y",
						"x",
						"weights",
						"formula",
						"glmod")
	out <- sapply(named_elements, function(nm) object[[nm]][[m]], simplify = FALSE)
	
	out$na.action 	<- object$na.action
	out$call 		    <- object$call
	out$algorithm 	<- object$algorithm
	out$stanfit 	  <- object$stanfit
		
	structure(out, class = c("stanreg", "lmerMod"))
}