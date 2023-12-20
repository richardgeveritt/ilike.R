evaluate_log_factor <- function(factor,
                                parameters,
                                data)
{
  if ("evaluate_log_prior" %in% names(factor))
  {
    return(factor$evaluate_log_prior(parameters=parameters))
  }

  if ("evaluate_log_likelihood" %in% names(factor))
  {
    return(factor$evaluate_log_likelihood(parameters=parameters,
                                          data=data))
  }

  return(0)
}

evaluate_gradient_log_factor <- function(factor,
                                         variable,
                                         parameters,
                                         data)
{
  if ("evaluate_gradient_log_prior" %in% names(factor))
  {
    return(as.matrix(factor$evaluate_gradient_log_prior(variable=variable,
                                                        parameters=parameters)))
  }

  if ("evaluate_gradient_log_likelihood" %in% names(factor))
  {
    return(as.matrix(factor$evaluate_gradient_log_likelihood(variable-variable,
                                                             parameters=parameters,
                                                             data=data)))
  }

  param = as.matrix(parameters[[variable]])
  return(matrix(0,nrow(param),ncol(param)))
}

evaluate_gradient_log_likelihood <- function(factor,
                                             variable,
                                             parameters,
                                             data)
{
  if ("evaluate_gradient_log_likelihood" %in% names(factor))
  {
    return(as.matrix(factor$evaluate_gradient_log_likelihood(variable=variable,
                                                             parameters=parameters,
                                                             data=data)))
  }

  param = as.matrix(parameters[[variable]])
  return(matrix(0,nrow(param),ncol(param)))
}

#' Evaluate the log of the factors in the model.
#'
#' @param model The output of calling ilike::compile on a file.
#' @param parameters A list containing the parameters.
#' @param data A list containing the data.
#' @param factor_index (optional) The indices giving which factors to evaluate (default is to evaluate all factors in the model).
#' @return A number giving the log probability,
#' @export
evaluate_log_factors <- function(model,
                                 parameters,
                                 data,
                                 factor_index=NULL)
{

  factors = model[["factor"]]

  if (is.null(factor_index))
  {
    factor_index = 1:length(factors)
  }

  return(sum(sapply(1:length(factor_index),FUN = function(i) { return(evaluate_log_factor(factors[[factor_index[i]]],parameters,data) ) } )))
}


#' Evaluate the log of the factors in the model.
#'
#' @param model The output of calling ilike::parse_ilike_model on a file.
#' @param variable The name of the variable to take the derivative with respect to.
#' @param parameters A list containing the parameters.
#' @param data A list containing the data.
#' @param factor_index (optional) The indices giving which factors to evaluate (default is to evaluate all factors in the model).
#' @return A matrix giving the gradient.
#' @export
evaluate_gradient_log_factors <- function(model,
                                          variable,
                                          parameters,
                                          data,
                                          factor_index=NULL)
{

  factors = model[["factor"]]

  if (is.null(factor_index))
  {
    factor_index = 1:length(factors)
  }

  all_gradients = sapply(1:length(factor_index),FUN = function(i) { return(as.matrix(evaluate_gradient_log_factor(factors[[factor_index[i]]],variable,parameters,data) ) ) } )

  if (is.matrix(all_gradients))
  {
    return(as.matrix(rowSums(all_gradients)))
  }
  else
  {
    return(as.matrix(sum(all_gradients)))
  }
}

#' Evaluate the log of the factors in the model.
#'
#' @param model The output of calling ilike::parse_ilike_model on a file.
#' @param variable The name of the variable to take the derivative with respect to.
#' @param parameters A list containing the parameters.
#' @param data A list containing the data.
#' @param factor_index (optional) The indices giving which factors to evaluate (default is to evaluate all factors in the model).
#' @return A matrix giving the gradient.
#' @export
evaluate_gradient_log_likelihoods <- function(model,
                                              variable,
                                              parameters,
                                              data,
                                              factor_index=NULL)
{

  factors = model[["factor"]]

  if (is.null(factor_index))
  {
    factor_index = 1:length(factors)
  }

  all_gradients = sapply(1:length(factor_index),FUN = function(i) { return(as.matrix(evaluate_gradient_log_likelihood(factors[[factor_index[i]]],variable,parameters,data) ) ) } )

  if (is.matrix(all_gradients))
  {
    return(as.matrix(rowSums(all_gradients)))
  }
  else
  {
    return(as.matrix(sum(all_gradients)))
  }
}

#' Evaluate the log of the factors in the model.
#'
#' @param model The output of calling ilike::parse_ilike_model on a file.
#' @param variables The name of the variables to take the derivative with respect to.
#' @param parameters A list containing the parameters.
#' @param data A list containing the data.
#' @param factor_index (optional) The indices giving which factors to evaluate (default is to evaluate all factors in the model).
#' @return A matrix giving the gradient.
#' @export
evaluate_gradient_log_likelihoods_all_parameters <- function(model,
                                                             variables,
                                                             parameters,
                                                             data,
                                                             factor_index=NULL)
{

  factors = model[["factor"]]

  if (is.null(factor_index))
  {
    factor_index = 1:length(factors)
  }

  return(as.matrix(unlist(lapply(1:length(variables),FUN = function(j) { all_gradients = sapply(1:length(factor_index),FUN = function(i) { return(evaluate_gradient_log_likelihood(factors[[factor_index[i]]],variables[j],parameters,data) ) } ); if (is.matrix(all_gradients)) { return(rowSums(all_gradients)) } else { return(sum(all_gradients))}  } ) ) ) )
}
