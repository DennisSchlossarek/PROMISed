# Functions for Shiny App: "PROMIS Data Processing"
#
# Content:
# Main Functions directly used for Processing App:
#   1. data.array.shiny()         - Organizes data in 4-dim-Array and applies FilterPeaksSpanFractions and 
#                                    optionally Threshold, Max-Normalization, Loess-Regression. Includes Progressbar
#   2. combine.reps.shiny()       - Combines Replicates based chosen Pearson-Correlation coefficent. Includes Progressbar
#   3. data.unarray()             - Reorganizes 4-dim-Array in 2-dim dataframe to enable table-output with write.table()
#   4. manhattan.anova.shiny()    - Calculates "dis-elution-score" (?): Compares manhattan-distances between Replicates of 
#                                    treatments with distances of replicates of the same treatment using ANOVA and TukeyHSD
#   4a.manhattan.anova.row()      - Calculates "dis-elution-score" (?): Compares manhattan-distances between Replicates of 
#                                    treatments with distances of replicates of the same treatment using ANOVA and TukeyHSD,
#                                    here onyl for selected row. Used to create Boxplots
#
#   5. deconvolute.new()          - Deconvolutes Peaks based on 3 criteria. Used in deconvolute.array(). Output also used
#                                    for reactive plotting with plotDeconvoluted. Original Script by Michal Gorka.
#   6. deconvolute.array()        - Applies deconvolute.new to an array created by data.array.shiny() 

# Plotting Functions
#   7. PlotMyArray.row()          - Creates lineplots from chosen row of an array created by data.array.shiny(), 
#                                    used for interactive plotting
#   8. PlotMyArray.shiny()        - Creates lineplots from an array created by data.array.shiny(),
#                                    used to create PDF of all lineplots from an array
#   9. PlotDeconvolited()         - Creates lineplots of deconvoluted profiles. Used as reactive plot to choose
#                                    deconvolution critera
#   10.manhattan.anova.boxplot()        <- Boxplot created form manhattan.anova.row results, shows Distance-vectors in
#                                    boxplot and the results of TukesHSD
#
# "Secondary" Functions used by Functions listed above:
#   10.namesCombn()               - combines Names of two Treatments using combn
#   11.maxNormalizeML()           - Normalizes Profiles to maximum Intensity. Original Script by Marcin.
#   12.loessDS()                  - wraps up loess-regression and post-smoothing normalization, removal of negative values
#   13.FilterPeaksSpanFraction_ML - Sets Single-Peaks to 0. Orginal Script by Marcin

# Functions for Data-Integration APP: 
#   b1. deconvoluted.list()       - rearanges input data-frame into saparate peak lists for all treatments
#   b2. scale_coreness()          - creates vector of scaled coreness and a vecotr of l = 20 containing a heatmap color scheme

colours <- as.character(brewer.pal(n = 10, name = "Paired"))
colours2 <- c(colours[1],colours[3],colours[5],colours[7],colours[9])
colours1 <- c(colours[2],colours[4],colours[6],colours[8],colours[10])
colours1 <- c(colours1, colours1, colours1)
colours2 <- c(colours2, colours2, colours2)

####################################
# Extract meta data information from colnames
####################################

colnames2metadata <- function(my_data, filename){
  
  my_extracted_colnames <- matrix(nrow = ncol(my_data), ncol = 4)
  
  i <- 1
  while(i <= ncol(my_data)){
    
    my_extracted_colnames[i,] <- strsplit(colnames(my_data)[i], split = "_")[[1]]
    
    i <- i + 1
    
  }
  
  # Extracted colnames are passed as arguments for data.array.shiny()-function AND OTHERS!
  
  extracted_metadata <- list(
    names_treatments  <- unique(my_extracted_colnames[,1]),
    nr_treatments     <- length(unique(my_extracted_colnames[,1])),
    names_rep         <- unique(my_extracted_colnames[,2]),
    nr_rep            <- max(as.numeric(unique(my_extracted_colnames[,3]))),
    names_column      <- unique(my_extracted_colnames[,4]),
    nr_column         <- length(unique(my_extracted_colnames[,4])),
    my_file_name <- paste0(strsplit(filename, ".txt")[[1]])
  )
  
}

####################################
# Arrange data in 4 dim Array, includes normalisation and smoothing
####################################

data.array.shiny <- function(x, filterpeak, normalization, lower_bound, smoothing, span_value,
                       nr_treatments, names_treatments, nr_rep, names_rep, nr_column, names_column){
  
  the_list <- list()
  rep_i <- list()
  repi_treatj <- data.frame()
  
#  withProgress(message = "Rearranging data", min = 0, max = nr_rep*nr_treatments, value = 0, {
    
  i <- 1
  while(i <= nr_rep){
    
    
    j <- 1
    while(j <= nr_treatments){
      
      repi_treatj <- data.frame(x[, ((grep(pattern = paste0(names_treatments[j], "_", names_rep, "_", i,"_", names_column[1]), colnames(x)))[1] : 
                                             (grep(pattern = paste0(names_treatments[j], "_", names_rep, "_", i,"_", names_column[1]), colnames(x)) + nr_column -1)[1])])
      
      # GREPs data of treatment "i" and replicate "j" to one data.frame based on the names and values set outside the function
      
      repi_treatj[is.na(repi_treatj)] <- 0  # replaces NAs with 0 (zero)
      
      if(filterpeak == TRUE){
      repi_treatj <- filterPeaksSpanFractions_ML(repi_treatj) # removes peaks consisting of single values
      }
      if(normalization == TRUE){
      
      repi_treatj <- maxNormalizeML(repi_treatj)  # normalizes elution profile for each row to max intensity 
      repi_treatj[is.na(repi_treatj)] <- 0  # replaces NAs with 0 (zero)
      repi_treatj[repi_treatj <= lower_bound] <- 0  # "cuts off" peaks resulting from data noise which are lower than "lower_bound" (recommended as 0.1)
        
      }
      
      if(smoothing == TRUE){
        
       repi_treatj <- loessDS(y = repi_treatj, span = span_value, normalization = normalization)
        
      }
      
      rownames(repi_treatj) <- rownames(x)
      colnames(repi_treatj) <- colnames(names_column)
      
      rep_i[[j]] <- repi_treatj
      
 #     incProgress(amount = 1)
      
      j <- j + 1
    }
  
    the_list[[i]] <- rep_i
    
    i <- i + 1
  }
#  })
  the_array <- array(unlist(the_list), dim = c(nrow(x), length(names_column), nr_treatments, nr_rep)
                     ,dimnames = list(rownames(the_list[[1]][[1]]), names_column, names_treatments, paste0(names_rep, c(1:nr_rep))))
  
  # Converts nested list to 4 dimensional array with dim = c(nrow(x), ncol(x), nr_treatments, nr_rep)
  
}


###################################
# Combine reps, including trigger to use single-occurence replicates
###################################

combine.reps.shiny <- function(x, PCC_thresh, names_treatments, rep_box, cor_method){
  
  reps_combined <- array(dim = c(dim(x)[1], dim(x)[2], dim(x)[3], 1))
  dimnames(reps_combined) <- list(rownames(x), colnames(x), names_treatments)
  
#  withProgress(message = "Combining Replicates", min = 0, max = dim(x)[1], value = 0, {
  
  comb_ind <- combn(dim(x)[4], 2)
  
  if(dim(x)[4] > 1){
    
    i <- 1
    while(i <= nrow(x)){                # 1st loop: Rows
      
      j <- 1
      while(j <= dim(x)[3]){            # 2nd loop: Treatments
        
        treat_ij <- t(x[i,,j, ])        # subsets values of x: x[row_i, (all columns), treatmant_j, (all replicates)]
        
        reps_ij <- data.frame(matrix(nrow = 0, ncol = dim(x)[2]))
        
        k <- 1
        while(k <= dim(comb_ind)[2]){
          
          n <- comb_ind[1,k]   
          m <- comb_ind[2,k]
          
          if(sum(treat_ij[n,]) > 0 & sum(treat_ij[m,]) > 0){
            
            pcc_nm <- cor(treat_ij[n,], treat_ij[m,], method = cor_method)
            pcc_nm[is.na(pcc_nm)] <- 0
            
            if(pcc_nm >= PCC_thresh){
              
              reps_ij <- rbind(reps_ij, treat_ij[n,], treat_ij[m,])
              
            }
            
          } else if(sum(treat_ij[n,]) > 0 & sum(treat_ij[m,]) == 0 & rep_box == TRUE){
            
            reps_ij <- rbind(reps_ij, treat_ij[n,])
            
          } else if(sum(treat_ij[m,]) > 0 & sum(treat_ij[n,]) == 0 & rep_box == TRUE){
            
            reps_ij <- rbind(reps_ij, treat_ij[m,])
            
          }
          
          k <- k + 1 # close loop 3
        }
        
        rep_ij <- maxNormalizeML(t(data.frame(colSums(unique(reps_ij)))))
        rep_ij[is.na(rep_ij)] <- 0
        
        reps_combined[i,,j,1] <- t(rep_ij)
        
        j <- j + 1 # close loop 2
        
      }
      
#      incProgress(amount = 1)
      
      i <- i + 1 # close loop 1
    }
    
    reps_combined[is.nan(reps_combined)] <- 0
    
    reps_combined
  }
#  })
} # close function


###################################
# Re-arrange your 4-dim array into a boring, 2-dim table 
###################################

data.unarray <- function(x, names_treatments, names_rep, names_column){    
  
  data_table <- data.frame(matrix(nrow = nrow(x), ncol = 0))
  rownames(data_table) <- rownames(x)
  
  i <- 1
  while(i <= dim(x)[3]){
    
    j <- 1
    while(j <= dim(x)[4]){
      
      data_table_ij <- x[,,i,j]
      colnames(data_table_ij) <- paste0(names_treatments[i], "_", names_rep, "_", j, "_", names_column)
      
      data_table <- cbind(data_table, data_table_ij)
      
      j <- j + 1
    }
    
    i <- i + 1
  }
  
  data_table   
}

#######################################################################################################################################
# Gives the Pr>F of a one-way Anova between Manhattan distance between two treatments and Manhattan distances within those treatments #
#######################################################################################################################################

####################################
# Manhattan-Anova for Shiny: Row- and pair-wise # <<<----------------<<<----------------<<<----------------<<<----------------<<<----------------
####################################

manhattan.anova.shiny <- function(x, names_treatments, pvalue){
  
 # withProgress(message = "Calculating Dis-Elution-Score", min = 0, max = nrow(x), value = 0, {
  
  comb_treat <- combn(dim(x)[3],2) 
  
  man_anova <- data.frame(matrix(nrow = nrow(x), ncol = ncol(combn(dim(x)[3], 2))))   # Data.frame to "collect"  Pr>F values
  man_anova_cor <- data.frame(matrix(nrow = nrow(x), ncol = ncol(combn(dim(x)[3], 2))))   # Data.frame to "collect" FDR corrected Pr>F values
  
  man_anova_check <- data.frame(matrix(nrow = nrow(x), ncol = ncol(combn(dim(x)[3], 2))))   # Data.frame to "collect" FDR corrected Pf>F values that passed "Tukey-Test"
  rownames(man_anova_check) <- rownames(x)
  colnames(man_anova_check) <- namesCombn(names_treatments)  
  
  tukey_df <- data.frame(matrix(nrow = 0, ncol = (ncol(combn(dim(x)[3], 2))*3) ))
  
  i <- 1
  while(i <= nrow(x)){                # 1st loop: Rows
  
    tukey_comb <- c()
    
    j <- 1                         # j == Nr of Combination
    while(j <= ncol(comb_treat)){     # 2nd loop: Combinations of treatments
      
      manhattan <- data.frame(matrix(nrow = dim(x)[4], ncol = dim(x)[4]))   # creates temporal data.frame to store manhattan distances of a-b
      rownames(manhattan) <- paste0("a",c(1:dim(x)[4]))
      colnames(manhattan) <- paste0("b",c(1:dim(x)[4]))
      
      auto_a <- data.frame(matrix(nrow = dim(x)[4], ncol = dim(x)[4]))      # creates temporal data.frame to store "auto"-manhattan distances of a-a
      auto_b <- data.frame(matrix(nrow = dim(x)[4], ncol = dim(x)[4]))      # creates temporal data.frame to store "auto"-manhattan distances of b-b
      
      treat_a <- x[i,,comb_treat[1,j], ]                                    # subsets values of x: x[row, (all columns), treatmant-a, (all replicates)]
      treat_b <- x[i,,comb_treat[2,j], ]                                    # subsets values of x: x[row, (all columns), treatmant-b, (all replicates)]
      # treatment-a and -b are defined by comb_treat[c(1,2),j]
      
      k <- 1
      while(k <= dim(x)[4]){  # 3rd loop: loop through rows of combinations of a and b
        
        l <- 1
        while(l <= ncol(treat_a)){ # 4th loop: loop trough columns of combinations of a and b
          
          manhattan[k,l]  <- MESS::auc(c(1:nrow(treat_a)), abs(treat_a[,k] - treat_b[,l]))    # calculates manhattan-distance between replicate k of a and replicate l of b
          auto_a[k,l]     <- MESS::auc(c(1:nrow(treat_a)), abs(treat_a[,k] - treat_a[,l]))    # calculates manhattan-distance between replicate k and l of a
          auto_b[k,l]     <- MESS::auc(c(1:nrow(treat_b)), abs(treat_b[,k] - treat_b[,l]))    # calculates manhattan-distance between replicate k and l of b
          
          l <- l + 1 # close loop 4
        }
        
        k <- k + 1 # close loop 3
      }
      
      manhattan[is.na(manhattan)] <- 0
      auto_a[is.na(auto_a)] <- 0
      auto_b[is.na(auto_b)] <- 0
      
        manhattan_v <- as.vector(t(manhattan))                   # transforms manhattan matrix to vector, used as input for anova
        auto_a_v    <- as.vector(t(auto_a[lower.tri(auto_a)]))   # transforms lower half of auto_a matrix to vector, used as input for anova
        auto_b_v    <- as.vector(t(auto_b[lower.tri(auto_b)]))   # transforms lower half of auto_b matrix to vector, used as input for anova
        
        anova_input <- data.frame(Y = c(manhattan_v, auto_a_v, auto_b_v),                                # creates input data for anova
                                  Treatment = factor(rep(c("manhattan", "auto_a", "auto_b"), 
                                                         times = c(length(manhattan_v), length(auto_a_v), length(auto_b_v)))))
        
        aov_x <- aov(Y~Treatment, data = anova_input)    # Anova test part 1
        
        anova_x <- anova(aov_x) # Anova test part 2
        anova_x[is.na(anova_x)] <- 1    # Sets NAs to 1 (maximum not-significant) to prevent loop from arresting
        
        tukey_x <- TukeyHSD(aov_x)[[1]] # Post-hoc TukeyHSD 
        tukey_x[is.na(tukey_x)] <- 1    # Sets NAs to 1 (maximum not-significant) to prevent loop from arresting
        
        tukey_comb <- c(tukey_comb,tukey_x[1,4], tukey_x[2,4], sum(c(tukey_x[1,4], tukey_x[2,4]) <= pvalue))
        
        # CHECK HERE
  #      if(tukey_x[1,4] <= pvalue & tukey_x[2,4] <= pvalue){ # & tukey_x[3,4] > pvalue){   # Checks if there is a significant difference between Manhattan and a, and Manhattan and b
  #        
          man_anova_check[i,j] <- TRUE                 
  #        
  #      } else {
  #        
  #        man_anova_check[i,j] <- FALSE                
  #        
  #      }
        
        man_anova[i,j] <- anova_x[1,5]  # Collects the uncorrected P.values for each combination
   
      j <- j + 1 # close loop 2
   }
     incProgress(amount = 1, detail = paste0(format((i / nrow(x)*100), digits = 3), "%"))
    
    tukey_df[i,] <- tukey_comb
    
    i <- i + 1 # close loop 1
  }
  
  colnames(tukey_df) <- unlist(lapply(X = namesCombn(names_treatments), FUN = function(X){paste(c("TukeyA", "TukeyB",  paste0("Tukey<=", pvalue)),X)}))
  
  man_anova_check[is.na(man_anova_check)] <- FALSE
  
  m <- 1
  while(m <= ncol(man_anova)){
    
    man_anova_cor[,m] <- p.adjust(man_anova[,m], method = "fdr")   # Corrects p.values of each combination using FDR corrections
    man_anova_cor[is.na(man_anova_cor)] <- 1
    
    man_anova_check[(man_anova_check[,m] == TRUE),m] <- man_anova_cor[(man_anova_check[,m] == TRUE),m] # passes corrected p.values for combination which passed "Tukey-Test"
    
    m <- m + 1
  }
  
  colnames(man_anova_check) <- apply(X = combn(names_treatments, 2, simplify = TRUE), MARGIN = 2, FUN = paste, collapse = "-")
  
  man_anova_check <- cbind(man_anova_check, tukey_df)
  
  man_anova_check  # lets output data.frame escape the function
#  })

} # close function

# 4a manhattan.row:

manhattan.row <- function(x, names_treatments, select_rows, j){
  
  comb_treat <- combn(dim(x)[3],2) 
  
  man_anova <- data.frame(matrix(nrow = nrow(x), ncol = ncol(combn(dim(x)[3], 2))))         # Data.frame to "collect" Pr>F values
  man_anova_cor <- data.frame(matrix(nrow = nrow(x), ncol = ncol(combn(dim(x)[3], 2))))     # Data.frame to "collect" FDR corrected Pr>F values
  
  man_anova_check <- data.frame(matrix(nrow = nrow(x), ncol = ncol(combn(dim(x)[3], 2))))   # Data.frame to "collect" FDR corrected Pf>F values that passed "Tukey-Test"
  rownames(man_anova_check) <- rownames(x)
  colnames(man_anova_check) <- namesCombn(names_treatments)  
  
  i <- which(rownames(x) == select_rows)                          # Select single row to analyse from array 
  
  manhattan <- data.frame(matrix(nrow = dim(x)[4], ncol = dim(x)[4]))   # creates temporal data.frame to store manhattan distances of a-b
  rownames(manhattan) <- paste0("a",c(1:dim(x)[4]))
  colnames(manhattan) <- paste0("b",c(1:dim(x)[4]))
  
  auto_a <- data.frame(matrix(nrow = dim(x)[4], ncol = dim(x)[4]))      # creates temporal data.frame to store "auto"-manhattan distances of a-a
  auto_b <- data.frame(matrix(nrow = dim(x)[4], ncol = dim(x)[4]))      # creates temporal data.frame to store "auto"-manhattan distances of b-b
  
  treat_a <- x[i,,comb_treat[1,j], ]                                    # subsets values of x: x[row, (all columns), treatmant-a, (all replicates)]
  treat_b <- x[i,,comb_treat[2,j], ]                                    # subsets values of x: x[row, (all columns), treatmant-b, (all replicates)]
  # treatment-a and -b are defined by comb_treat[c(1,2),j]
  
  k <- 1
  while(k <= dim(x)[4]){  # 3rd loop: loop through rows of combinations of a and b
    
    l <- 1
    while(l <= ncol(treat_a)){ # 4th loop: loop trough columns of combinations of a and b
      
      manhattan[k,l]  <- MESS::auc(c(1:nrow(treat_a)), abs(treat_a[,k] - treat_b[,l]))    # calculates manhattan-distance between replicate k of a and replicate l of b
      auto_a[k,l]     <- MESS::auc(c(1:nrow(treat_a)), abs(treat_a[,k] - treat_a[,l]))    # calculates manhattan-distance between replicate k and l of a
      auto_b[k,l]     <- MESS::auc(c(1:nrow(treat_b)), abs(treat_b[,k] - treat_b[,l]))    # calculates manhattan-distance between replicate k and l of b
      
      l <- l + 1 # close loop 4
    }
    
    k <- k + 1 # close loop 3
  }
  
  manhattan[is.na(manhattan)] <- 0
  auto_a[is.na(auto_a)] <- 0
  auto_b[is.na(auto_b)] <- 0
  
  
  if(sum(manhattan) != 0 & sum(auto_a) != 0 & sum(auto_b) != 0 ){   # Bypasses the actual ANOVA testing, if all Inputs to ANOVA ar 0 (or NAs) 
    
    manhattan_v <- as.vector(t(manhattan))                   # transforms manhattan matrix to vector, used as input for anova
    auto_a_v    <- as.vector(t(auto_a[lower.tri(auto_a)]))   # transforms lower half of auto_a matrix to vector, used as input for anova
    auto_b_v    <- as.vector(t(auto_b[lower.tri(auto_b)]))   # transforms lower half of auto_b matrix to vector, used as input for anova
    
    anova_input_list <- list(manhattan_v, auto_a_v, auto_b_v)
  }
  
  anova_input_list
  
}

####################################
# Deconvolution
####################################

deconvolute.new <- function(my_data, var_min_peak, var_limit_small, var_limit_large){
  
  # 1. Finding local maxima ###################################################################
  
  my_data <- cbind(0 ,my_data, 0)
  peak_position <- matrix(nrow = dim(my_data)[1], ncol = dim(my_data)[2], 0)
  
  i <- 1
  for(i in 1:dim(my_data)[1]){
    
    max_i       <- max(my_data[i,])
    min_peak_i  <- max_i * var_min_peak 
    my_data_v  <- as.vector(t(my_data[i,]))
    peaks_i     <- extract(turnpoints(my_data_v), dim(my_data)[2], peak = 1, pit = 0)
    
    peaks_i[which(my_data_v < min_peak_i)] <- 0
    
    peak_position[i,] <- peaks_i
    
  }
  
  peak_position <- peak_position[rowSums(peak_position) != 0,, drop = FALSE]
  
  # 2. Deconvolution ###########################################################################
  
  my_data_decon  <- numeric()
  peak_width      <- numeric()
  
  i <- 1
  for(i in 1:dim(my_data)[1]){ # loop over rows
    
    
    var_record  <- 0  # will be set to 1 when intensity relative to max-intensity is greater than var_limit_small
    peak_count  <- 1
    rowname_i   <- rownames(my_data)[i]
    random3     <- 3
    
    my_data_i       <- my_data[i,, drop = FALSE] # subsetting my_data 
    peak_position_i <- peak_position[i,, drop = FALSE] # subsetting peak_position
    peak_record_pos <- rep(0, dim(my_data_i)[2]) # create vector to record position?
    peak_record_int <- rep(0, dim(my_data_i)[2]) # create vector to record peak height?
    
    #  j <- 1
    for(j in 1:dim(my_data_i)[2]){ # loop over columns
      
      max_peak <- max(peak_record_int)
      my_data_relative <-  my_data_i[,j]/max(my_data_i)
      
      if(var_record == 0 & my_data_relative > var_limit_small){ # checks, if intensity relative to max-intensity is greater than var_limit_small -> Increasing slope?
        
        var_record <- 1
        
        if(j > 1){
          peak_record_int[(j-1)] <- my_data_i[,(j-1)] # records last data-point before threshold (var_limit_small) was reached 
        }
      }
      
      if(var_record == 1 & my_data_relative > var_limit_small){ 
        
        peak_record_pos[j] <- peak_position_i[,j]
        peak_record_int[j] <- my_data_i[,j]
        
        if(my_data_i[,j] < max_peak*var_limit_large){ # checks if intensity is lower than previous recorded max-intensity * var_limit_large -> Peak-maximum reached?
          
          var_record <- 2
        }
      }
      
      if(var_record == 1 & my_data_relative <= var_limit_small){  # checks, if intensity relative to max-intensity is lower than var_limit_small -> End of peak?
        # If reached, "saves" recorded peak to separate row in my_data_decon
        peak_record_pos[j] <- peak_position_i[,j]
        peak_record_int[j] <- my_data_i[,j]
        
        if(sum(peak_record_pos) > 0){ #cheks, if a maximum has been "hit" already
          
          my_data_decon <- rbind(my_data_decon, peak_record_int) # save recorded peak in own row
          rownames(my_data_decon)[dim(my_data_decon)[1]] <- paste0(rowname_i, "_K_", which.max(peak_record_int)) # label peak with rowname + separator + peak_count
          peak_count <- peak_count + 1
        }
        
        var_record  <- 0  # will be set to 1 when intensity relative to max-intensity is greater than var_limit_small
        peak_record_int  <- rep(0, dim(my_data_i)[2]) # create NEW vector to record peak height?
        peak_record_pos  <- rep(0, dim(my_data_i)[2]) # create NEW vector to record position?
        
      }
      
      if(var_record == 2){
        
        if(my_data_relative > var_limit_small &  my_data_i[,j]*var_limit_large <= my_data_i[,(j-1)]){ # is intensity still falling?
          
          peak_record_pos[j] <- peak_position_i[,j]
          peak_record_int[j] <- my_data_i[,j]
          
          if(my_data_i[,j] < max_peak*var_limit_large){ # checks if intensity is lower than previous recorded max-intensity * var_limit_large -> Peak-maximum reached?
            
            var_record <- 2
          }
        } 
        if(my_data_relative <= var_limit_small){  # checks, if intensity relative to max-intensity is lower than var_limit_small -> End of peak?
          # If reached, "saves" recorded peak to separate row in my_data_decon
          peak_record_pos[j] <- peak_position_i[,j]
          peak_record_int[j] <- my_data_i[,j]
          
          if(sum(peak_record_pos) > 0){ #cheks, if a maximum has been "hit" already
            
            my_data_decon <- rbind(my_data_decon, peak_record_int) # save recorded peak in own row
            rownames(my_data_decon)[dim(my_data_decon)[1]] <- paste0(rowname_i, "_K_", which.max(peak_record_int)) # label peak with rowname + separator + peak_count
            peak_count <- peak_count + 1
          }
          
          var_record  <- 0  # will be set to 1 when intensity relative to max-intensity is greater than var_limit_small
          peak_record_int  <- rep(0, dim(my_data_i)[2]) # create NEW vector to record peak height?
          peak_record_pos  <- rep(0, dim(my_data_i)[2]) # create NEW vector to record position?
          
        }
        
        if(my_data_i[,j] * var_limit_large > my_data_i[,(j-1)]){  # is intesity increasing again -> start of new peak
          # stop recording peaks
          if(sum(peak_record_pos) > 0){ #cheks, if a maximum has been "hit" already
            
            my_data_decon <- rbind(my_data_decon, peak_record_int) # save recorded peak in own row
            rownames(my_data_decon)[dim(my_data_decon)[1]] <- paste0(rowname_i, "_K_", which.max(peak_record_int)) # label peak with rowname + separator + peak_count
            peak_count <- peak_count + 1
          }
          
          var_record  <- 1  # will be set to 1 when intensity relative to max-intensity is greater than var_limit_small
          peak_record_int  <- rep(0, dim(my_data_i)[2]) # create NEW vector to record peak height?
          peak_record_pos  <- rep(0, dim(my_data_i)[2]) # create NEW vector to record position?
          peak_record_int[j] <- my_data_i[,j]
          peak_record_int[j-1] <- my_data_i[,j-1]
          peak_record_pos[j] <- peak_position_i[,j]
        }
      }
    }
  }
  
  my_data_decon <- my_data_decon[,-(dim(my_data_decon)[2]), drop = FALSE] # remove last column of zeroes
  my_data_decon <- my_data_decon[,-1, drop = FALSE] # remove first column of zeroes
  
  my_data_decon
  
} # End of Function


#######################
# Deconvolute Array-structured data-sets
#######################

deconvolute.array <- function(x, var_min_peak, var_limit_small, var_limit_large, names_treatment, name_col, nr_col){

  my_data_deconvo <- data.frame()
  
  k <- 1
  while(k <= dim(x)[3]){
    
    g <- x[,,k,]
    g[is.na(g)] <- 0
    g <- g[(rowSums(g) > 0),]
    
    f <- data.frame(deconvolute.new(g, var_min_peak = var_min_peak, var_limit_small = var_limit_small, var_limit_large = var_limit_large))
    rownames(f) <- paste0(rownames(f), "_", names_treatment[k])
    f <- cbind(treatment = names_treatment[k], f)
    
    my_data_deconvo <- rbind(my_data_deconvo, f)
    
    k <- k + 1
  }
  
  colnames(my_data_deconvo) <- c("Treatment", name_col)
  
  my_data_deconvo
  
}

####################################
# Plot all replicateas from Array
####################################

# Plot single Metabolites based on rowname

PlotMyArray.row <- function(x, my_row_name, names_treatments, titel_name, plot_what){
  
  par(mar = c(12,4,4,2))
  
  lwidth <- 1.6

  i <- which(rownames(x) == my_row_name)
  
    plot(c(1:dim(x)[2]), cex = 0, lwd = lwidth, xlab = "fraction", ylab = "relative intensity", ylim = c(0,max(x[i,,,])))
    title(titel_name, cex = 0.5)    # sets plot title to rowname " t " <- Metabolite name
    legend(x = 1, 
           y = -(max(x[i,,,])*0.15), 
           legend = names_treatments[c(1:dim(x)[3])], 
           #col = brewer.pal("Set2", n = 8)[c(1:dim(x)[3])], 
           col = pal_startrek()(7)[c(1:dim(x)[3])],
           lwd = lwidth, 
           xpd = TRUE, 
           cex = 1, 
           title = "Treatment")
    
    #    legend(x = 0, y = -0.12, legend = c(1:nrow(x)), col = colours1[1:nrow(x)], lwd = lwidth, xpd = TRUE, title = "Peaks") 
    
    j <- 1
    while(j <= dim(x)[3]){ # 2nd loop throudh treatments
      
      k <- 1
      while(k <= dim(x)[4]){ # 3rd loop through replicates
        
        if(plot_what[j,k] == TRUE){
        
        lines(c(1:dim(x)[2]), 
              x[i, ,j,k], 
              cex = 2, 
              lwd = lwidth,
              type = "l", 
           #   col = brewer.pal("Set2", n = 8)[j], 
              col = pal_startrek()(7)[j],
              ylim = c(0, max(x[i,,,])))
          #   ,ylim = c(0,1) ,
          #   xlab = "fraction",
          #   ylab = "relative intensity")
        
        }
        
        k <- k + 1
      }
      
      j <- j + 1 
    }
 
}

#####################

PlotMyArray.shiny <- function(x,  names_treatments){
  
  par(mar = c(12,4,4,2))
  lwidth <- 1.3
  
  i <- 1
  while(i <= dim(x)[1]){    # 1st loop through rows
    
    plot(c(1:dim(x)[2]), 
         cex = 0, 
         lwd = lwidth, 
         ylim = c(0,1), 
         xlab = "fraction", 
         ylab = "relative intensity")
    title(paste0("#", i ,"_",rownames(x)[i]), cex = 0.5)    # sets plot titel to rowname " t " <- Metabolite name
    
    legend(x = 1, 
           y = -(max(x)*0.15), 
           legend = names_treatments[c(1:dim(x)[3])], 
           col = pal_startrek()(7)[c(1:dim(x)[3])], 
           lwd = lwidth, 
           xpd = TRUE, 
           cex = 1, 
           title = "Treatment")
    
    j <- 1
    while(j <= dim(x)[3]){ # 2nd loop throudh treatments
      
      k <- 1
      while(k <= dim(x)[4]){ # 3rd loop through replicates
        
        lines(c(1:dim(x)[2]), 
              x[i, ,j,k], 
              cex = 2, 
              lwd = lwidth, 
              ylim = c(0,1), 
              type = "l", 
              col = pal_startrek()(7)[j] ) 
        #, xlab = "fraction", 
        #ylab = "relative intensity")
        
        k <- k + 1
      }
      
      j <- j + 1 
    }
    
    i <- i + 1  
  }
}

# Plot Deconvoluted Profiles 

PlotDeconvoluted <- function(rep_comb, deconvoluted, my_row_name, treatment){
  
  inputrow <- which(rownames(rep_comb) == my_row_name)
  
  my_rep_comb  <- rep_comb[inputrow,,treatment,]
  
  deconvoluted_select <- deconvoluted[grep(pattern = my_row_name, x = rownames(deconvoluted), fixed = TRUE),]
  
  if(nrow(deconvoluted_select) > 0){
    my_legend <- c(my_row_name, paste0("Peak_", c(1 : nrow(deconvoluted_select))))
  } else {
    my_legend <- c(my_row_name)
  }
  
  # Plotting
  par(mar = c(12,4,4,2))
  lwidth <- 2
  
  plot(my_rep_comb, type = "l", cex = 0, lwd = 1, ylim = c(0,1), xlab = "fraction", ylab = "relative intensity")
  legend(x = 1, y = -(max(my_rep_comb)*0.15), legend = my_legend, col = c(1:length(my_legend)), lwd = lwidth, xpd = TRUE, cex = 1, title = "Peaks")
  
 if(nrow(deconvoluted_select) > 0) {
    
   i <- 1
    while(i <= nrow(deconvoluted_select)){
    
      lines(c(1:ncol(deconvoluted_select)), deconvoluted_select[i,], type = "l", col = i + 1, lwd = lwidth)
    
    i <- i + 1
   }
 }
    
}


# 10 manhattan.anova.boxplot

manhattan.anova.boxplot <- function(manhattan_distances, name_panel, names_treatments, j){
  
  distance    <- manhattan_distances
  manhattan_v <- distance[[1]]
  auto_a_v    <- distance[[2]]
  auto_b_v    <- distance[[3]]
  

  combination_names <- c(name_panel[j], paste0("|", combn(names_treatments, 2)[1,j], "|") , paste0("|", combn(names_treatments, 2)[2,j], "|"))
  
  anova_input <-  data.frame(man_dist = c(manhattan_v, auto_a_v, auto_b_v),
                             Combinations = factor(rep(combination_names,
                                                       times = c(length(manhattan_v), length(auto_a_v), length(auto_b_v)))))
  
  yheight <- max(anova_input$man_dist)
  yheight1 <- yheight + (0.12 * yheight)
  yheight2 <- yheight + (0.20 * yheight)
  yheight3 <- yheight + (0.28 * yheight)
  
  model = lm(man_dist ~ Combinations, data = anova_input)
  aov_x = aov(model)
  TUKEY <- TukeyHSD(aov_x)
  
  pvalues <- format(TUKEY[[1]][,4], digits = 2, scientific = TRUE)
  
  p <- ggplot(data = anova_input, aes(x = Combinations , y = man_dist, fill = Combinations)) +  
    geom_boxplot(outlier.shape = NA, color = "black", fill = "white", varwidth = FALSE) +
    ylim(0, (yheight + 0.33*yheight)) +
    ggtitle("Comparing Manhattan Distances between Treatments and Replicates") +
    xlab("Comparison") + 
    ylab("Manhattan Distance") +
    
    theme_set(theme_bw()) +
    theme(legend.position = "none") +
    theme(text = element_text(size = 15)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") +
    
    geom_segment(aes(x = 1, y = (yheight1) , xend = 2, yend = (yheight1) )) +
    geom_segment(aes(x = 1,  xend = 1, y = (yheight1 - 0.03*yheight1) , yend = (yheight1) )) +
    geom_segment(aes(x = 2,  xend = 2, y = (yheight1 - 0.03*yheight1) ,  yend = (yheight1) )) +
    geom_text(label = pvalues[1], aes(x = 1.5, y = yheight1 + yheight1 * 0.03), size = 5) +
    
    geom_segment(aes(x = 2, y = (yheight2) , xend = 3, yend = (yheight2) )) +
    geom_segment(aes(x = 2,  xend = 2, y = (yheight2 - 0.03*yheight2) , yend = (yheight2) )) +
    geom_segment(aes(x = 3,  xend = 3, y = (yheight2 - 0.03*yheight2) ,  yend = (yheight2) )) +
    geom_text(label = pvalues[3], aes(x = 2.5, y = yheight2 + yheight2 * 0.03), size = 5) +
    
    geom_segment(aes(x = 1, y = (yheight3) , xend = 3, yend = (yheight3) )) +
    geom_segment(aes(x = 1,  xend = 1, y = (yheight3 - 0.03*yheight3) , yend = (yheight3) )) +
    geom_segment(aes(x = 3,  xend = 3, y = (yheight3 - 0.03*yheight3) ,  yend = (yheight3) )) +
    geom_text(label = pvalues[2], aes(x = 2, y = yheight3 + yheight3 * 0.03), size = 5) 
  
  p
  
  tt <- ttheme_default()
  
  tbl_anova <- format(anova(aov_x), scientific = TRUE, digits = 2)
  tbl_tukey <- format(TUKEY[[1]], scientific = TRUE, digits = 4)
  
  blank <- grid.rect(gp=gpar(col="white"))
  
  tbl_anova <- tableGrob(tbl_anova, theme = tt)
  tbl_tukey <- tableGrob(tbl_tukey, theme = tt)
  tbls <- grid.arrange(blank, tbl_anova, tbl_tukey, blank, blank, blank, blank, nrow = 7, as.table = TRUE )#, heights = c(2,2)) #, widths = unit(10, "cm"))
  
  grid.arrange(p, tbls, nrow = 1, ncol = 2, widths = c(3,2), as.table = TRUE)#, heights = c(1,1)) #, width = c(3,1)) # 3,2
  
}




############################
# Smaller functions used: 
############################

# Gives combination of all treatment names

namesCombn <- function(x) {
  
  pair_names <- data.frame(combn(x, 2, simplify = TRUE ))
  comb_names <- c()
  i <- 1
  while(i <= length(pair_names)){
    
    comb_names[i] <- paste("[" ,pair_names[1,i], " - ", pair_names[2,i], "]", sep = "")
    
    i <- i +1
  }
  comb_names
}

# Max normalization

maxNormalizeML <- function(x){      
  
  dataNew <- x
  
  #x[x == 0] <- NA
  
  maxFactor <- apply(x, 1, max, na.rm = TRUE)
  
  dataNew <- sweep(x, 1, apply(x, 1, max, na.rm = TRUE), FUN = "/")
  
  dataNew <- as.data.frame(dataNew)
}

#### Looping for loess-regression

loessDS <- function (y, span, normalization) {
  
  y_loess <- data.frame(matrix(nrow = nrow(y), ncol = ncol(y), NA))
  
  i <- 1
  while(i <= nrow(y)){
    
    y_loess[i,] <- loess(t(y[i,]) ~ c(1:ncol(y)), span = span)$fitted
    y_loess[i,(y_loess[i,] <= 0)]  <- 0

    if(normalization == TRUE){y_loess[i,] <- maxNormalizeML(y_loess[i,])}
    
    y_loess[i, is.na(y_loess[i,])]  <- 0
    
    i <- i + 1
    
  } 
  
  y_loess  
}

# Removes single-value-peaks

filterPeaksSpanFractions_ML <- function(z){
  dataConsPeaks <- as.data.frame(NULL)
  w <- 1
  while (w <= nrow(z)){
    tmpNos <- z[w,] # copy intensities from row w of data x
    tmpCount <- which(tmpNos>0) # check which "cells" are higher than 0
    tmp<-c() # create empty vector
    tmpr<-c() # create empty vector
    
    for (p in 1:length(tmpCount)){
      tmp[p]<- tmpCount[p+1]-(tmpCount[p]) # check what is the distance between "non-zero" values - 1st way
    }
    
    for (q in 1:length(tmpCount)){
      tmpr[q]<- tmpCount[q]-(tmpCount[q-1]) # check what is the distance between "non-zero" values - 2nd way
    }
    
    tmpb <- NULL
    for(i in 1:length(tmp)){
      tmpb[i] <- tmp[i] == 1 | tmpr[i] == 1 # check if at least one method gives 1 (positive value)
    }
    
    tmpb[is.na(tmpb)] <- FALSE # remove NA values from tmpb and replace by FALSE
    
    tmpCountReal <- tmpCount[tmpb] # check which cells are higher than 0 AND the distance between non zero values is 1
    #tmpNosReal <- tmpNos[,tmpCountReal] # probably no needed
    
    tmpm <- as.data.frame(matrix(0, nrow = nrow(tmpNos), ncol = ncol(tmpNos)))
    x <- 1
    y <- 1
    
    while(x <= length(tmpNos)){
      while(y <= length(tmpCountReal)){
        if(x == tmpCountReal[y]){
          tmpm[,x] <- tmpNos[x]
          y <- y+1
        }else{
          y <- y+1
        }
      }
      y <- 1
      x <- x+1
    }
    
    colnames(tmpm) <- colnames(tmpNos)
    rownames(tmpm) <- rownames(tmpNos)
    
    dataConsPeaks <- rbind(dataConsPeaks, tmpm)
    w <- w+1
  }
  dataConsPeaks # give the result of the function to the "clipboard"
}

#################################
# Function for Data-Integration #
#################################

deconvoluted.list <- function(x){
  
  names_treatments <- unique(x[,1])
  nr_treatments <- length(names_treatments)
  
  decon_list <- list()
  
  i <- 1
  while(i <= nr_treatments){
    
    y <- x[which(x[,1] == names_treatments[i]),]
    
    decon_list[[i]] <- y[,-1]
    
    i <- i + 1
  }
  decon_list
}

scale_coreness <- function(x){
  
  coreness_scale <- c()
  coreness_rel <- x/(max(x, na.rm = FALSE))
  
  i <- 1
  while(i <= 20){
    
    coreness_scale[which(coreness_rel >= ((i-1)*0.05) & coreness_rel <= (i+1)*0.05)] <- i
    
    i <- i + 1
  }
  
#  coreness_color   <- rev(rainbow(max(coreness_scale, na.rm = FALSE)*1.5)[0:max(coreness_scale, na.rm = FALSE)])
  coreness_color <- c("#440154", "#481568", "#482677", "#453781", "#3F4788",
                      "#39558C", "#32648E", "#2D718E", "#287D8E", "#238A8D", 
                      "#1F968B", "#20A386", "#29AF7F", "#3CBC75", "#56C667", 
                      "#74D055", "#94D840", "#B8DE29", "#DCE318", "#FDE725")
  
  list(coreness_scale, coreness_color)
}


get_peak_names <- function(x){
  
  names_decon1 <- data.frame(rownames(x))
  colnames(names_decon1) <- c("Data_names1")
  names_decon1_split <- c()
  
  k <- 1
  while(k <= nrow(x)){
    
    names_decon1_k <- as.character(names_decon1[k,])
    
    names_decon1_split[k] <- unlist(strsplit(names_decon1_k, split = c("_K_")))[1]
    
    k <- k + 1
  }
  
  names_decon1_split <- data.frame(names_decon1_split)
  my_names1 <- unique(names_decon1_split)
  
}

cor_reactive <- function(data1, data2, pattern, method){
  
    my_x <- as.matrix(data1)
    my_y <- as.matrix(data2)
    
    my_x <- my_x[grep(pattern = pattern, rownames(my_x), fixed = TRUE),, drop = FALSE]
    my_correl_table <- cor(x = t(my_y), y = t(my_x), method = method)
    
    my_correl_table
}


Plot2selections <- function(data1, data2, selector1, selector2){
  
  x <- data1[grep(pattern = selector1, rownames(data1), fixed = TRUE),, drop = FALSE]
  x <- x / max(x)
  
  y <- data2[grep(pattern = selector2, rownames(data2), fixed = TRUE),, drop = FALSE]
  y <- y / max(y)
  
  par(mar = c(12,4,4,2))
  lwidth <- 2
  
  plot(data.frame(c(1:dim(x)[2])), cex = 0, lwd = lwidth, ylim = c(0,1),xlab = "", ylab = "Relative Intensity", xaxt = "n")
  axis(1, seq(0, ncol(x), 1), tick = TRUE, at = c(0:dim(data1)[2]), label = as.character(c(0:ncol(x))))
  text(x = 38, y = -(max(x)*0.2), pos = 2, labels = c("Fraction"), xpd = TRUE, cex  = 0)
  
  
  title(paste0( selector1," vs ", selector2), cex = 0.5)
  legend(x = 0, y = -0.15, legend = rownames(x), col = c("black"), lwd = lwidth, xpd = TRUE, title = "Selection 1") 
  legend(x = 10, y = -0.15, legend = rownames(y), col =  c("red"), lwd = lwidth, xpd = TRUE, title = "Selection 2") 
  
  l <- 1
  while(l <= dim(x)[2]){
    
    lines(c(1:dim(x)[2]), x[l,], cex = 2, lwd = lwidth, type = "l", col = c("black") )
    lines(c(1:dim(y)[2]), y[l,], cex = 2, lwd = lwidth, type = "l", col = c("red") )
    
    l <- l + 1
  }  
  
}  

PlotBestHit <- function(decon_data1, decon_data2, selector, pcc_threshold, method){
  
  my_correl_table1 <- reactive({
    my_x <- as.matrix(decon_data1)
    my_y <- as.matrix(decon_data2)
    
    my_x <- my_x[grep(pattern = selector, rownames(my_x), fixed = TRUE),, drop = FALSE]
    my_correl_table <- cor(x = t(my_y), y = t(my_x), method = method)
    
    my_correl_table
  })
  
  pcc_row <- data.frame(my_correl_table1())
  
  profiles_hit <- data.frame(matrix(nrow = 0, ncol = ncol(decon_data1)))
  
  i<- 1
  while(i <= ncol(pcc_row)){
    
    if(max(pcc_row[,i]) >= as.numeric(pcc_threshold)){
      profiles_hit[i,] <- decon_data2[grep(pattern = rownames(pcc_row[which.max(pcc_row[,i]),, drop = FALSE]), 
                                           rownames(pcc_row), fixed = TRUE), ]
    } else {
      profiles_hit[i,] <-  data.frame(matrix(nrow = 1, ncol = ncol(decon_data1), 0))
    }
    i <- i + 1
  }
  
  x <- decon_data1[grep(pattern = selector, rownames(decon_data1), fixed = TRUE),, drop = FALSE]
  x <- x/ max(x)
  
  y <- profiles_hit
  
  par(mar = c(12,4,4,2))
  lwidth <- 2
  
  plot(c(1:dim(x)[2]), cex = 0, lwd = lwidth, ylim = c(0,1), xlab = "Fraction", ylab = "Relative Intensity", xaxt = "n")
  axis(1, seq(0, ncol(x), 1), tick = TRUE, at = c(0:dim(decon_data1)[2]), label = as.character(c(0:ncol(x))))
  
  title(selector, cex = 0.5)
  legend(x = 0, y = -0.15, legend = c(1:nrow(x)), col = colours1[1:nrow(x)], lwd = lwidth, xpd = TRUE, title = "Peaks") 
  legend(x = 10, y = -0.15, legend = rownames(y), col = colours2[1:nrow(y)], lwd = lwidth, xpd = TRUE, title = "Co-Eluting Peaks") 
  
  i <- 1
  while(i <= dim(x)[2]){
    
    lines(c(1:dim(x)[2]), x[i,], cex = 2, lwd = lwidth, type = "l", col = colours1[i] )
    lines(c(1:dim(y)[2]), maxNormalizeML(y[i,]) * max(x[i,]), cex = 2, lwd = lwidth, type = "l", col = colours2[i] )
    
    i <- i + 1
  }
}


ListCoelutions <- function(data1, data2, selector, pcc_threshold, method, names_treatments, nr_treatments){
  
  diagram_venn1 <- reactive({    # Create reactive correlation_list of input$select_front_1
    
    coeluting_list <- list()  
    
    my_color <-  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    listnames <- names_treatments
    drop <- c()
    mycol <- c()
    
    i <- 1
    while(i <= nr_treatments){
      
      x <- as.matrix(data1[[i]])
      y <- as.matrix(data2[[i]])
      
      x <- x[grep(pattern = selector, rownames(x), fixed = TRUE),, drop = FALSE]
      my_correl_table <- cor(x = t(y), y = t(x), method = method)
      
      my_correl_table <- my_correl_table >= pcc_threshold
      
      coeluting <- rownames(my_correl_table)[rowSums(my_correl_table) > 0]
      
      if(length(coeluting) > 0){
        
        coeluting_names <- c()
        
        k <- 1
        while(k <= length(coeluting)){
          
          coeluting_names[k] <- strsplit(coeluting[k], "_K_")[[1]][1] 
          
          k <- k + 1  
        }
        
        coeluting_list[[i]] <- unique(coeluting_names)
        mycol <- c(mycol, my_color[i])
        
      } else {
        
        drop <- c(drop, i)
        
      }
      
      i <- i + 1
    }
    
    if(length(drop) > 0){
      coeluting_list <- coeluting_list[-drop]
      listnames <- listnames[-drop]
    }
    names(coeluting_list) <- listnames
    coeluting_list 
  })
}

GetIntersections <- function(x){
  
  overlap <- calculate.overlap(x)
  
  if(length(x) == 2){
    names(overlap) <- c("a1", "a2", "a12")
  } else if(length(x) == 3){
    names(overlap) <- c("a123", "a12", "a13", "a23", "a1", "a2", "a3")
  } else if(length(x) == 4){
    names(overlap) <- c("a1234", "a123", "a124", "a134", "a234", "a12", "a13", "a14", "a23", "a24", "a34", "a1", "a2", "a3", "a4")
  }
  
  maxrow <- c()
  k <- 1
  while(k <= length(overlap)){
    
    maxrow[k] <- length(overlap[[k]])
    
    k <- k + 1
  }
  
  overlap_matrix <- matrix(nrow = max(maxrow), ncol = length(overlap), NA)
  colnames(overlap_matrix) <- names(overlap)
  
  k <- 1
  while(k <= length(overlap)){
    
    if(length(c(overlap[[k]])) > 0){
      overlap_matrix[(1:maxrow[k]), k] <- c(overlap[[k]])  
    }
    k <- k + 1
  }
  
  overlap_matrix

}

custom_igraph <- function(cor_table, pcc){
  
  if(nrow(cor_table) == 0 | ncol(cor_table) == 0 | sum(cor_table) == 0 ){
    
    cor_table <- matrix(nrow = 2, ncol = 2, 1)
    rownames(cor_table) <- c("Warning:", "No Interactions Found")
    colnames(cor_table) <- c("Warning:", "No Interactions Found")
  
  }
  
  if(nrow(cor_table) == ncol(cor_table)){
  
    if(all(colnames(cor_table) == rownames(cor_table))){
    
      network_igraph <- graph_from_adjacency_matrix(cor_table, 
                                                    mode = "undirected", 
                                                    weighted = TRUE, 
                                                    diag = TRUE)
    
    
    } else if(all(colnames(cor_table) != rownames(cor_table))){
    
      my_links <- data.frame(matrix(nrow = 0, ncol = 3))
      colnames(my_links) <- c("from", "to", "weight")
    
      i <- 1         
      while(i <= ncol(cor_table)){
      
        my_links_sub <- cbind( from = as.character(colnames(cor_table)[i]),
                               to   = as.character(rownames(cor_table)),
                               weight = as.numeric(cor_table[,i]))
      
        my_links <- rbind(my_links, my_links_sub)
      
        i <- i + 1
      }  
    
    my_links <- unique(my_links[as.numeric(as.character(my_links[,3])) >= pcc,])
    
    network_igraph <- graph_from_data_frame(my_links, directed = FALSE)
  }
  } else {
    
    my_links <- data.frame(matrix(nrow = 0, ncol = 3))
    colnames(my_links) <- c("from", "to", "weight")
    
    i <- 1         
    while(i <= ncol(cor_table)){
      
      my_links_sub <- cbind( from = as.character(colnames(cor_table)[i]),
                             to   = as.character(rownames(cor_table)),
                             weight = as.numeric(cor_table[,i]))
      
      my_links <- rbind(my_links, my_links_sub)
      
      i <- i + 1
    }  
    
    my_links <- unique(my_links[as.numeric(as.character(my_links[,3])) >= pcc,])
    
    network_igraph <- graph_from_data_frame(my_links, directed = FALSE)
    
  }
  network_igraph  
}
