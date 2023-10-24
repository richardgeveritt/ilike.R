#' Load the data from the model.
#'
#' @param sampler_output The output of sampler giving the points we will use to estimate the active subspace, formatted using the ilike::load_xxxx_output functions.
#' @param model The output of calling ilike::parse_ilike_model on a file, giving the model for which we wish to find the active subspace.
#' @param use_weights (optional) If FALSE, will ignore particle weights in the line graph. If TRUE, will use the particle weights. (defaults to TRUE)
#' @return An eigendecomposition from which we can find the active subspace.
#' @export
active_subspace_eigen_decomposition <- function(sampler_output,
                                                model,
                                                use_weights=TRUE)
{
  data = get_data(model)

  output_to_use = sampler_output
  if ("ExternalTarget" %in% names(sampler_output))
  {
    output_to_use = dplyr::filter(sampler_output,ExternalTarget==max(output_to_use$ExternalTarget))
  }

  if ("Target" %in% names(sampler_output))
  {
    output_to_use = dplyr::filter(sampler_output,Target==max(output_to_use$Target))
  }

  parameters = sampler_output_to_list_of_matrices(output_to_use)
  if (use_weights==TRUE)
  {
    if ("LogWeight" %in% names(sampler_output))
    {
      log_weights = sampler_output_log_weights_to_list_of_matrices(output_to_use)
    }
    else
    {
      stop('"use_weights" set to TRUE, but no LogWeight column found in "sampler_output".')
    }
  }

  parameter_names = names(parameters)

  # For each variable, calculate gradient and stack them up.

  # single_point = get_single_point_from_list_of_matrices(parameters,parameter_names,1)
  #
  # evaluate_gradient_log_likelihoods_all_parameters(model,
  #                                                  parameter_names,
  #                                                  single_point,
  #                                                  data)

  number_of_points = nrow(parameters[[parameter_names[1]]])

  all_gradients = as.matrix(sapply(1:number_of_points,FUN=function(i){ single_point = get_single_point_from_list_of_matrices(parameters,parameter_names,i); return(evaluate_gradient_log_likelihoods_all_parameters(model,
                                                                                                                                                                                                                                parameter_names,
                                                                                                                                                                                                                                single_point,
                                                                                                                                                                                                                                data)) } ))




  if (use_weights==TRUE)
  {
    for_eigen = lapply(1:number_of_points,FUN = function(i) { exp(log_weights[i]) * all_gradients[,i] %*% t(all_gradients[,i]) } )
  }
  else
  {
    for_eigen = lapply(1:number_of_points,FUN = function(i) { (1/number_of_points) * all_gradients[,i] %*% t(all_gradients[,i]) } )
  }

  return(eigen(Reduce('+', for_eigen)))
}
