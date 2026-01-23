pegs <- function(Y, X, maxit = 100, logtol = -4.0, NNC = FALSE, 
                 covbend = 1.1, covMinEv = 10e-4, XFA = -1) {
  
  # --- 1. Input Validation and Preparation for Y ---
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (typeof(Y) == "integer") storage.mode(Y) <- "numeric"
  
  # Store Y column names for later use in output
  trait_names <- colnames(Y)
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait", seq_len(ncol(Y)))
    colnames(Y) <- trait_names
  }
  
  # --- 2. Input Validation and Preparation for X ---
  # If the user provides a single matrix for X, wrap it in a list for compatibility.
  if (!is.list(X)) {
    X <- list(Effect1 = X)
  }
  
  # Ensure all elements in the list are numeric matrices
  X <- lapply(X, function(mat) {
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    if (typeof(mat) == "integer") storage.mode(mat) <- "numeric"
    return(mat)
  })
  
  # Get names for the output lists. Create default names if missing.
  effect_names <- names(X)
  if (is.null(effect_names) || any(effect_names == "")) {
    effect_names <- paste0("Effect", seq_along(X))
    names(X) <- effect_names
  }
  
  # --- 3. Call the C++ function ---
  res <- .Call("PEGS_lapack", # Call the new multi-effect function
               Y,
               X,               # Pass the list of matrices
               as.integer(maxit),
               as.double(logtol),
               as.logical(NNC),
               as.double(covbend),
               as.double(covMinEv),
               as.integer(XFA)
  )
  
  # --- 4. Process and Name the Output for Clarity ---
  # Rename list elements returned from C++ to be more user-friendly
  names(res)[names(res) == "b_list"] <- "b"
  names(res)[names(res) == "GC_list"] <- "GC"
  names(res)[names(res) == "bend_list"] <- "bend"
  
  # Add names to the top-level list elements
  names(res$b) <- effect_names
  names(res$GC) <- effect_names
  names(res$bend) <- effect_names
  
  # Add dimension names to the single-item outputs
  dimnames(res$hat) <- dimnames(Y)
  names(res$h2) <- trait_names
  names(res$mu) <- trait_names
  
  # Loop through the list-based outputs to name their dimensions
  for (i in seq_along(effect_names)) {
    effect <- effect_names[i]
    
    # Get column names for the current X matrix
    marker_names <- colnames(X[[effect]])
    if (is.null(marker_names)) {
      marker_names <- paste0("M", seq_len(ncol(X[[effect]])))
    }
    
    # Set dimnames for the coefficient matrix b[[i]]
    dimnames(res$b[[effect]]) <- list(marker_names, trait_names)
    
    # Set dimnames for the genetic correlation matrix GC[[i]]
    dimnames(res$GC[[effect]]) <- list(trait_names, trait_names)
  }
  
  
  # If there is a single random effect, unlist the result for convenience
  if(length(res$b) == 1){
    res$b <- res$b[[1]]
    res$GC <- res$GC[[1]]
    res$bend <- res$bend[[1]]
  }
  
  
  return(res)
}