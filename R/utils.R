extract_single_parameter <- function(single_point_output,
                                     parameter_name)
{
  return(as.matrix(dplyr::filter(single_point_output,ParameterName==parameter_name)$Value))
}

extract_single_parameter_to_matrix <- function(output,parameter_name)
{
  new_variable_names = mapply(FUN = function(a,b) { paste(a,"_",b,sep="") },output$ParameterName,output$Dimension)
  output = subset(output,select = -c(ParameterName,Dimension))
  output = subset(output,select = -c(Time,NormalisingConstant,ISESS))

  if ("AncestorIndex" %in% names(output))
  {
    output = subset(output,select = -c(AncestorIndex))
  }

  if ("LogWeight" %in% names(output))
  {
    output = subset(output,select = -c(LogWeight))
  }

  output$Parameter = new_variable_names
  output = dplyr::distinct(output)

  output_to_use = tidyr::pivot_wider(output,names_from=Parameter,values_from=Value)

  output_to_use = as.matrix(output_to_use[,3:length(names(output_to_use))])
  colnames(output_to_use)<-NULL
  return(output_to_use)
}

extract_parameters <- function(single_point_output)
{
  parameter_names = unique(single_point_output$ParameterName)
  parameters = lapply(1:length(parameter_names),FUN=function(i){return( extract_single_parameter(single_point_output,parameter_names[i]) )})
  names(parameters) = parameter_names
  return(parameters)
}

extract_chain <- function(chain_output)
{
  return(lapply(1:unique(chain_output$Iteration),FUN=function(i){return( extract_parameters(dplyr::filter(chain_output,Iteration==i)) )}))
}

extract_list_of_matrices <- function(chain_output)
{
  if ("ExternalIndex" %in% names(chain_output))
  {
    chain_output = subset(chain_output,select = -c(ExternalIndex))
  }

  if ("ExternalTarget" %in% names(chain_output))
  {
    chain_output = subset(chain_output,select = -c(ExternalTarget))
  }

  if ("Target" %in% names(chain_output))
  {
    chain_output = subset(chain_output,select = -c(Target))
  }

  parameter_names = unique(chain_output$ParameterName)
  parameter_matrices = lapply(1:length(parameter_names),FUN=function(i){return( extract_single_parameter_to_matrix(dplyr::filter(chain_output,ParameterName==parameter_names[i]),parameter_names[i]) )})
  names(parameter_matrices) = parameter_names
  return(parameter_matrices)
}

extract_log_weights <- function(chain_output)
{
  for_output = as.matrix(chain_output$LogWeight)
  colnames(for_output)<-NULL
  return(for_output)
}

extract_target <- function(target_output)
{
  if ("Particle" %in% names(target_output))
  {
    parameter_list = lapply(unique(target_output$Particle),FUN=function(i){return(extract_chain(dplyr::filter(target_output,Particle==i)))})
  }
  else if ("Chain" %in% names(target_output))
  {
    parameter_list = lapply(unique(target_output$Chain),FUN=function(i){return(extract_chain(dplyr::filter(target_output,Chain==i)))})
  }

  return(unlist(parameter_list,recursive = FALSE))
}

extract_target_to_list_of_matrices <- function(target_output)
{
  return(extract_list_of_matrices(target_output))
}

extract_log_weights_target_to_list_of_matrices <- function(target_output)
{
  return(extract_log_weights(target_output))
}

extract_external_target <- function(external_target_output)
{
  if ("Target" %in% names(external_target_output))
  {
    parameter_list = lapply(unique(external_target_output$Target),FUN=function(i){return(extract_target(dplyr::filter(external_target_output,Target==i)))})
    if (length(unique(sampler_output$Target))==1)
    {
      parameter_list = parameter_list[[1]]
    }
  }
  else
  {
    if ("Particle" %in% names(external_target_output))
    {
      parameter_list = lapply(unique(external_target_output$Particle),FUN=function(i){return(extract_chain(dplyr::filter(external_target_output,Particle==i)))})
    }
    else if ("Chain" %in% names(external_target_output))
    {
      parameter_list = lapply(unique(external_target_output$Chain),FUN=function(i){return(extract_chain(dplyr::filter(external_target_output,Chain==i)))})
    }
    return(unlist(parameter_list,recursive = FALSE))
  }
}

extract_external_target_to_list_of_matrices <- function(external_target_output)
{
  if ("Target" %in% names(external_target_output))
  {
    parameter_list = lapply(unique(external_target_output$Target),FUN=function(i){return(extract_target_to_list_of_matrices(dplyr::filter(external_target_output,Target==i)))})
    if (length(unique(sampler_output$Target))==1)
    {
      parameter_list = parameter_list[[1]]
    }
  }
  else
  {
    return(extract_list_of_matrices(external_target_output))
  }
}

extract_log_weights_external_target_to_list_of_matrices <- function(external_target_output)
{
  if ("Target" %in% names(external_target_output))
  {
    parameter_list = lapply(unique(external_target_output$Target),FUN=function(i){return(extract_log_weights_target_to_list_of_matrices(dplyr::filter(external_target_output,Target==i)))})
    if (length(unique(sampler_output$Target))==1)
    {
      parameter_list = parameter_list[[1]]
    }
  }
  else
  {
    return(extract_log_weights(external_target_output))
  }
}

sampler_output_to_nested_list <- function(sampler_output)
{
  if ( ("ExternalIndex" %in% names(sampler_output)) && (length(unique(sampler_output$ExternalIndex))>1) )
  {
    stop('sampler_output contains an ExternalIndex with more than one value: not yet supported in sampler_output_to_parameter_list.')
  }

  if ("ExternalTarget" %in% names(sampler_output))
  {
    parameter_list = lapply(unique(sampler_output$ExternalTarget),FUN=function(i){return(extract_external_target(dplyr::filter(sampler_output,ExternalTarget==i)))})
    if (length(unique(sampler_output$ExternalTarget))==1)
    {
      parameter_list = parameter_list[[1]]
    }
  }
  else if ("Target" %in% names(sampler_output))
  {
    parameter_list = lapply(unique(sampler_output$Target),FUN=function(i){return(extract_target(dplyr::filter(sampler_output,Target==i)))})
    if (length(unique(sampler_output$Target))==1)
    {
      parameter_list = parameter_list[[1]]
    }
  }
  else
  {
    if ("Particle" %in% names(external_target_output))
    {
      parameter_list = lapply(unique(external_target_output$Particle),FUN=function(i){return(extract_chain(dplyr::filter(external_target_output,Particle==i)))})
    }
    else if ("Chain" %in% names(external_target_output))
    {
      parameter_list = lapply(unique(external_target_output$Chain),FUN=function(i){return(extract_chain(dplyr::filter(external_target_output,Chain==i)))})
    }
    parameter_list = unlist(parameter_list,recursive = FALSE)
  }

  return(parameter_list)
}

sampler_output_to_list_of_matrices <- function(sampler_output)
{
  if ( ("ExternalIndex" %in% names(sampler_output)) && (length(unique(sampler_output$ExternalIndex))>1) )
  {
    stop('sampler_output contains an ExternalIndex with more than one value: not yet supported in sampler_output_to_parameter_list.')
  }

  if ("ExternalTarget" %in% names(sampler_output))
  {
    parameter_list = lapply(unique(sampler_output$ExternalTarget),FUN=function(i){return(extract_external_target_to_list_of_matrices(dplyr::filter(sampler_output,ExternalTarget==i)))})
    if (length(unique(sampler_output$ExternalTarget))==1)
    {
      parameter_list = parameter_list[[1]]
    }
  }
  else if ("Target" %in% names(sampler_output))
  {
    parameter_list = lapply(unique(sampler_output$Target),FUN=function(i){return(extract_target_to_list_of_matrices(dplyr::filter(sampler_output,Target==i)))})
    if (length(unique(sampler_output$Target))==1)
    {
      parameter_list = parameter_list[[1]]
    }
  }
  else
  {
    return(extract_list_of_matrices(sampler_output))
  }

  return(parameter_list)
}

sampler_output_log_weights_to_list_of_matrices <- function(sampler_output)
{
  sampler_output = dplyr::filter(sampler_output,Dimension==1)
  if ( ("ExternalIndex" %in% names(sampler_output)) && (length(unique(sampler_output$ExternalIndex))>1) )
  {
    stop('sampler_output contains an ExternalIndex with more than one value: not yet supported in sampler_output_to_parameter_list.')
  }

  if (!("LogWeight" %in% names(sampler_output)))
  {
    stop('"LogWeight must be in output in order to use this function.')
  }

  if ("ExternalTarget" %in% names(sampler_output))
  {
    parameter_list = lapply(unique(sampler_output$ExternalTarget),FUN=function(i){return(extract_log_weights_external_target_to_list_of_matrices(dplyr::filter(sampler_output,ExternalTarget==i)))})
    if (length(unique(sampler_output$ExternalTarget))==1)
    {
      parameter_list = parameter_list[[1]]
    }
  }
  else if ("Target" %in% names(sampler_output))
  {
    parameter_list = lapply(unique(sampler_output$Target),FUN=function(i){return(extract_log_weights_target_to_list_of_matrices(dplyr::filter(sampler_output,Target==i)))})
    if (length(unique(sampler_output$Target))==1)
    {
      parameter_list = parameter_list[[1]]
    }
  }
  else
  {
    return(extract_log_weights(sampler_output))
  }

  return(parameter_list)
}

get_single_point_from_list_of_matrices <- function(parameters,parameter_names,i)
{
  parameter_list = lapply(parameter_names,FUN = function(n) { return(parameters[[n]][i,]) } )
  names(parameter_list) = parameter_names
  return(parameter_list)
}
