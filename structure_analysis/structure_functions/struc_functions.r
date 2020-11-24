# Functions used during analysis of STRUCTURE data

order_struc_samps <- function(clumpp_result_df, pop_order, zero_cut = 0.1){
  ####################
  # Function for ordering samples based on group membership and specific
  #   order of populations
  # INPUTS
  # clumpp_result_df = dataframe with results from CLUMPP
  # pop_order = vector of the order that populations should be displayed in
  #               the figure
  # zero_cut = the value at which membership is set to 0; makes it easier
  #              to "clump" and sort "pure" samples
  # OUTPUT
  # vector of indices for the order of samples in clumpp_result_df
  ###################
  test_results <- clumpp_result_df
  # make low values = 0 to make clumping "pure" samples easier
  test_results[test_results < zero_cut] <- 0
  #
  # get the rank-order of membership for each sample
  order_mat <- apply(test_results[, c(2:ncol(test_results))], 1,
    order, decreasing = T)
  #
  # while-statement to cycle through all the population comparisons
  comp_interval = 1
  tot_ord_inds <- c()
  pop_comp_vec <- c()
  while(length(pop_order) >= (2*comp_interval)){
    for(i in seq(comp_interval)){
      tmp_pop_order <- pop_order[seq(from = i, to = length(pop_order),
        by = comp_interval)]
      for(j in seq(length(tmp_pop_order))){
        pop_1 <- tmp_pop_order[j]
        if(j == length(tmp_pop_order)){
          pop_2 <- tmp_pop_order[1]
        } else{
          pop_2 <- tmp_pop_order[j+1]
        }
        # need to know if the first and last pops have been compared before
        p12_comp_name <- paste(pop_1, pop_2, sep = ':')
        p21_comp_name <- paste(pop_2, pop_1, sep = ':')
        p1_top_inds <- apply(order_mat, 2, function(x) x[1] == pop_1)
        p2_top_inds <- apply(order_mat, 2, function(x) x[1] == pop_2)
        p1_2nd_inds <- apply(order_mat, 2, function(x) x[2] == pop_1)
        p2_2nd_inds <- apply(order_mat, 2, function(x) x[2] == pop_2)
        #
        p12_inds <- intersect(which(p1_top_inds), which(p2_2nd_inds))
        # if have already compared first and last pops, then skip adding the
        #    comparison to prevent duplicating samples
        if((p12_comp_name %in% pop_comp_vec) | 
             (p21_comp_name %in% pop_comp_vec)){
          p12_inds <- c()
        }
        if(length(p12_inds) > 0){
          p12_ord <- p12_inds[order(test_results[p12_inds, (pop_1 + 1)],
            decreasing = T)]
          tot_ord_inds <- c(tot_ord_inds, p12_ord)
        }
        #
        p21_inds <- intersect(which(p2_top_inds), which(p1_2nd_inds))
        if((p12_comp_name %in% pop_comp_vec) | 
             (p21_comp_name %in% pop_comp_vec)){
          p21_inds <- c()
        }
        if(length(p21_inds) > 0){
          p21_ord <- p21_inds[order(test_results[p21_inds, (pop_2 + 1)],
            decreasing = F)]
          tot_ord_inds <- c(tot_ord_inds, p21_ord)
        }
        pop_comp_vec <- c(pop_comp_vec, p12_comp_name, p21_comp_name)
      }
    }
    comp_interval <- comp_interval + 1
  }
  return(tot_ord_inds)
}


