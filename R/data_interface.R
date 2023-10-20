#' Load the data from the model.
#'
#' @param model The output of calling ilike::parse_ilike_model on a file.
#' @return A list containing the data.
#' @export
get_data <- function(model)
{
  return(model[["data"]][[1]]$data())
}
