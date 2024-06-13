new_regression_model <- function(design, outcome, model_name) {
  model <- list(design = design, outcome = outcome, name = model_name)
  class(model) <- paste(model_name, "model", sep = "_")
  return(model)
}

get_logit_outcome_pair <- function(model) {
  outcome <- model$outcome
  if (is.list(outcome)) {
    n_success <- outcome$n_success
    n_trial <- outcome$n_trial
  } else {
    n_success <- outcome
    n_trial <- 1 # Assume binary outcome for all observations
  }
  return(list(n_success = n_success, n_trial = n_trial))
}

matvec_by_design <- function(model, v, subset_ind = NULL) {
  if (is.null(subset_ind)) {
    return(as.vector(model$design %*% v))  
  }
  else {
    return(row_subset_matvec_via_transpose(t(model$design), v, subset_ind))
  }
}

matvec_by_design_transp <- function(model, w, subset_ind = NULL) {
  if (is.null(subset_ind)) {
    return(as.vector(t(w) %*% model$design))  
  }
  else {
    return(t(model$design[subset_ind,]) %*% w)
    #return(transpose_row_subset_matvec(model$design, w, subset_ind))
  }
}
