new_regression_model <- function(design, outcome, model_name, via_transp) {
  n_obs <- nrow(design)
  n_pred <- ncol(design)
  if (via_transp) {
    model <- list(design_transpose = t(design), outcome = outcome, name = model_name, n_obs = n_obs, n_pred = n_pred)   
  }
  else {
    model <- list(design = design, outcome = outcome, name = model_name, n_obs = n_obs, n_pred = n_pred)
  }
  class(model) <- paste(model_name, "model", sep = "_")
  return(model)
}

get_logit_outcome_pair <- function(model, subset_ind = NULL) {
  outcome <- get_outcome(model, subset_ind)
  if (is.list(outcome)) {
    n_success <- outcome$n_success
    n_trial <- outcome$n_trial
  } else {
    n_success <- outcome
    n_trial <- 1 # Assume binary outcome for all observations
  }
  return(list(n_success = n_success, n_trial = n_trial))
}

matvec_by_design <- function(model, v, subset_ind = NULL, via_transp, use_rcpp) {
  if (is.null(subset_ind)) {
    return(as.vector(model$design %*% v))  
  }
  else {
    if (via_transp) {
      if (use_rcpp) {
        return(row_subset_matvec_via_transpose(model$design_transpose, v, subset_ind))
      }
      else {
        return(as.vector(t(v) %*% model$design_transpose[, subset_ind]))
      }
    }
    else {
      if (use_rcpp) {
        return(row_subset_matvec(model$design, v, subset_ind))
      }
      else {
        return(model$design[subset_ind, ] %*% v)
      }
    }
  }
}

matvec_by_design_transp <- function(model, w, subset_ind = NULL, via_transp, use_rcpp) {
  if (is.null(subset_ind)) {
    return(as.vector(t(w) %*% model$design))  
  }
  else {
    if (via_transp) {
      if (use_rcpp) {
        return(col_subset_matvec(model$design_transpose, w, subset_ind))
      }
      else {
        return(model$design_transpose[, subset_ind] %*% w)
      }  
    }
    else {
      design_transpose <- t(model$design)
      if (use_rcpp) {
        return(col_subset_matvec(design_transpose, w, subset_ind))
      }
      else {
        return(design_transpose[, subset_ind] %*% w)
      }
    }
  }
}

get_outcome <- function(model, subset_ind) {
  if (is.null(subset_ind)) {
    return(model$outcome)  
  }
  else {
    return(model$outcome[subset_ind])
  }
}

initialize_option <- function(option) {
  if(is.null(option$mle_solver)) {
    option$mle_solver = "newton"
  }
  if (option$mle_solver == "SGD") {
    if(is.null(option$use_matvec_via_transp)) {
      option$use_matvec_via_transp = TRUE
    }
    if(is.null(option$use_cpp_matvec)) {
      option$use_cpp_matvec = TRUE
    }
  }
  else {
    option$use_matvec_via_transp = FALSE
    option$use_cpp_matvec = TRUE
  }
  return(option)
}