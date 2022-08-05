# modified code from the 'rifi' R package

fragment_delay <- function(probe, pen, pen_out) {
  stranded <- 1

  num_args <- list(pen, pen_out)
  names(num_args) <- c("pen", "pen_out")

  # the dataframe is sorted by strand and position.
  probe <- probe[with(probe, order(-xtfrm(probe$strand), probe$position)), ]

  probe[probe$strand == "-", ] <- probe[probe$strand == "-", ][
    order(probe[probe$strand == "-", ]$position, decreasing = TRUE), ]

  probe[, "delay_fragment"] <- NA
  probe[, "velocity_fragment"] <- NA
  probe[, "intercept"] <- NA
  probe[, "slope"] <- NA

  # tmp_df takes all relevant variables from probe
  tmp_df <-
    data.frame(
      ID = probe$ID,
      val = probe$delay,
      position = probe$position,
      seg = probe$position_segment
    )
  # if strand information is given and stranded is TRUE, all "-" are inverted...
  # ...and strand is added to tmp_df.
  if (stranded == TRUE) {
    tmp_df$strand <- probe$strand
    # the positions are inverted for the "-" strand
    # so that the slope is seen as increasing by the scoring function.
    tmp_df[tmp_df$strand == "-", "position"] <-
      ((tmp_df[tmp_df$strand == "-", "position"]) -
      (tmp_df[tmp_df$strand == "-", ][1, "position"])) * -1
  }
  # All lines with any NA are ignored at this point and are taken care of later
  # in the function.
  tmp_df <- na.omit(tmp_df)

  # makes a vector of all position segments (S_1,S_2,...)
  unique_seg <- unlist(unique(tmp_df$seg))

  # count is needed to name the fragments.
  count <- 1

  # II. Dynamic Programming: the scoring function is interpreted

  # the foreach loop iterates over each unique segment
  frags <-list()
  
  for(k in seq_along(unique_seg)) {
    # only the part of the tmp_df that responds to the respective segment is
    # picked
    section <- tmp_df[which(tmp_df$seg == unique_seg[k]), ]

    # best_frags collects all scores that the dp is referring to
    best_frags <- c()
    # best names collects the names, and its last element is returned as result
    best_names <- c()

    if (nrow(section) > 2) {
      # only segments with more than two values are grouped...*

      for (i in 3:nrow(section)) {
        # the loop iterates over each value in the segment
        # this part always goes from position 1 to the referred position
        # 1:3,1:4...
        tmp_score <-
          score_fun_linear(section[seq_len(i), "val"],
                           section[seq_len(i), "position"],
                           section[seq_len(i), "ID"], pen_out, stranded)
        tmp_name <- names(tmp_score)
        # in this loop all smaller parts are scored e.g (i = 6) 6:6,5:6,4:6...
        # they are then combined with the former score e.g 1,2,3,4,5|6,
        # 1,2,3,4|5,6...
        if (i > 5) {
          # only parts bigger than 6 are accepted as three is the smallest
          # possible fragment size
          for (j in (i - 2):4) {
            tmp_val <- section[j:i, "val"]
            tmp_position <- section[j:i, "position"]
            tmp_ID <- section[j:i, "ID"]
            # penalty for a new fragment and former scores are added
            tmp <-
              score_fun_linear(tmp_val, tmp_position, tmp_ID, pen_out,
                               stranded) + pen + best_frags[j - 3]
            tmp_score <- c(tmp_score, tmp)
            # the new fragment is pasted to its corresponding former fragment
            tmp_n <- paste0(best_names[j - 3], "|", names(tmp))
            tmp_name <- c(tmp_name, tmp_n) # the names is cached
          }
        }
        # from the first score eg 1:6 and the smaller scores from the loop
        # 1,2,3,4,5|6, 1,2,3,4|5,6... the smallest is chosen and passed to
        # best_frags and best_names for the next iteration
        pos <-
          which(tmp_score == min(tmp_score))[1] # lowest score is collected
        tmp_score <- tmp_score[pos]
        tmp_name <- tmp_name[pos]
        best_frags <- c(best_frags, tmp_score)
        best_names <- c(best_names, tmp_name)
      }
    } else {
      #* ...all segments with less than three values are grouped automatically
      tmp_score <-
        score_fun_linear(section[, "val"], section[, "position"],
                         section[, "ID"], pen_out, stranded)
      tmp_name <- names(tmp_score)
      best_names <- c(best_names, tmp_name)
    }
    # the final result put into a list called frags
    frags[[k]]<-best_names[length(best_names)]
  }

  # III. Fill the dataframe

  # this loop iterates over the segments to fill the dataframe
  for (k in seq_along(frags)) {
    # the single fragments are split by |
    na <- strsplit(frags[[k]], "\\|")[[1]]

    # the loop goes over each fragment
    for (i in seq_along(na)) {
      # trgt are all IDs in the fragment
      tmp_trgt <- strsplit(na[i], "_")[[1]][1]
      trgt <- strsplit(tmp_trgt, ",")[[1]]
      # This statement checks if there are outliers
      if (length(strsplit(na[i], "_")[[1]]) == 4) {
        # outl are the IDs of all outliers
        tmp_outl <- strsplit(na[i], "_")[[1]][4]
        outl <- strsplit(tmp_outl, ",")[[1]]
        # while the original trgt gets stripped of the outliers, so trgt are
        # now all valid IDs
        trgt <- trgt[-which(trgt %in% outl)]
        # now the real outliers are named
        rows <- match(outl, probe[, "ID"])
        nam <- paste0("D_", count, "_O")
        probe[rows, "delay_fragment"] <- nam
        # and the values are assigned (normal outliers are allowed to share...
        # ...the same values, while terminal outliers are not)
        probe[rows, "slope"] <-
          as.numeric(strsplit(na[i], "_")[[1]][2])
        probe[rows, "intercept"] <-
          as.numeric(strsplit(na[i], "_")[[1]][3])
        probe[rows, "velocity_fragment"] <-
          1 / (as.numeric(strsplit(na[i], "_")[[1]][2]))
      }
      # last thing is that the real fragment is named and the values are...
      # ...assigned the same way as for the outliers
      rows <- match(trgt, probe[, "ID"])
      nam <- paste0("D_", count)
      probe[rows, "delay_fragment"] <- nam
      probe[rows, "slope"] <- as.numeric(strsplit(na[i], "_")[[1]][2])
      probe[rows, "intercept"] <-
        as.numeric(strsplit(na[i], "_")[[1]][3])
      probe[rows, "velocity_fragment"] <-
        1 / (as.numeric(strsplit(na[i], "_")[[1]][2]))
      # the count increases by one for the next iteration
      count <- count + 1
    }
  }
  # in this part the NAs are dealt with initial check if the first lines are NA
  if (sum(cumprod(is.na(probe$delay_fragment))) > 0) {
    probe$delay_fragment[seq_len(sum(cumprod(is.na(probe$delay_fragment))))] <-
      "D_0.5_NA"
  }
  # all rows with NA in delay fragment are collected
  row_NA <- which(is.na(probe$delay_fragment))
  # check if there are any NAs
  if (length(row_NA) > 0) {
    # group is initiated with the first NA
    group <- c(row_NA[1])
    # the loop iterates over all NAs
    for (i in seq_along(row_NA)) {
      # if there are more NAs following the first NA, they will be put
      # together in group
      if (is.na(probe$delay_fragment[row_NA[i] + 1])) {
        group <- c(group, row_NA[i] + 1)
      }
      # if there are no more NAs to add, it will be checked if the fragment
      # above and below is the same or not. Normal outliers are ignored in
      # this case
      else if (gsub("_O", "", probe$delay_fragment[group[1] - 1]) ==
               gsub("_O", "", probe$delay_fragment[group[length(group)] + 1])) {
        # if they are the same, so the NAs are in the middle of one fragment
        # they get the name marked with an _NA and the values of the fragment
        probe$delay_fragment[group] <-
          paste0(gsub("_O", "", probe$delay_fragment[group[1] - 1]), "_NA")
        probe$velocity_fragment[group] <-
          probe$velocity_fragment[group[1] - 1]
        probe$intercept[group] <- probe$intercept[group[1] - 1]
        probe$slope[group] <- probe$slope[group[1] - 1]
        # group is initiated with the next NA after the ones that were just
        # dealt with
        group <- row_NA[i + 1]
      } else if (probe$delay_fragment[group[1] - 1] !=
                 probe$delay_fragment[group[length(group)] + 1]) {
        # if the NAs are between two fragments, it makes a new fragment
        # between the two
        probe$delay_fragment[group] <-
          paste0("D_", as.numeric(gsub(
            "D_|_O", "", probe$delay_fragment[group[1] - 1]
          )) + 0.5, "_NA")
        # group is initiated with the next NA after the ones that were just
        # dealt with
        group <- row_NA[i + 1]
      }
    }
    # this deals with the case, that NAs might be ate the very end of the
    # dataframe
    probe$delay_fragment[is.na(probe$delay_fragment)] <-
      paste0("D_", as.numeric(gsub("D_|_O", "",
               probe$delay_fragment[!is.na(probe$delay_fragment)]
               [length(probe$delay_fragment[!is.na(probe$delay_fragment)])])) +
               0.5, "_NA")
  }
  # the function ends with the return of the dataframe
  probe
}


score_fun_linear <-
  function(y,
           x,
           z = x,
           pen,
           stran,
           n_out = min(10, max(1, 0.2 * length(x)))) {
    if (length(unique(x)) > 1) {
      # a linear regression can not be used on just one position
      n_out <-
        min(n_out, length(x) - 3) # the number of allowed outliers is decided
      mo <- lm(y ~ x) # the linear fit is performed
      mo_save <- mo
      tmp <- abs(residuals(mo)) # the residuals are cached
      tmp_velo <- coef(mo)[[2]] # velocity...
      tmp_inter <- coef(mo)[[1]] # .. and intercept are cached
      if (coef(mo_save)[[2]] < 0 &
        stran == TRUE) {
        # if the stranded option is active and the slope is negative, the
        # residuals are calculated as if the slope was 0.
        tmp <- abs(y - (mean(y))) # the score is overwritten
        tmp_velo <- 0 # velocity...
        tmp_inter <- mean(y) # .. and intercept are overwritten
      }
      if (coef(mo_save)[[2]] < (- (1 / 60)) &
        stran == FALSE) {
        # if the stranded option is inactive and ...
        # ...the slope is smaller - 1/60 the residuals are calculated as if
        # the slope was -1/60
        mo <- lm(y - I((-1 / 60)) * x ~ 1)
        tmp <- abs(residuals(mo))
        tmp_velo <- (-1 / 60) # velocity...
        tmp_inter <- coef(mo)[[1]] # .. and intercept are overwritten
      }
      if (coef(mo_save)[[2]] > (1 / 60)) {
        # if the slope is bigger 1/60, the residuals are calculated as
        # if the slope was 1/60
        mo <- lm(y - I(1 / 60) * x ~ 1)
        tmp <- abs(residuals(mo))
        tmp_velo <- 1 / 60 # velocity...
        tmp_inter <- coef(mo)[[1]] # .. and intercept are overwritten
      }
      tmp_score <- sum(tmp) # the sum of residuals is the temporary score
      names(tmp) <-
        z # tmp, y and x are named by x to decide for outliers
      names(y) <- z
      names(x) <- z
      out <-
        sort(tmp, decreasing = TRUE)[1] # out is the sorted vector of residuals
      if (n_out >= 1) {
        # checks if more than 0 outliers are allowed
        for (i in seq_len(n_out)) {
          # the loop iterates over the number of allowed outliers
          tmp_n <- names(out) # all but the i worst are selected
          tmp_y <-
            y[!names(y) %in% tmp_n] # new y is chosen (without outliers)
          tmp_x <-
            x[!names(y) %in% tmp_n] # new x is chosen (without outliers)
          mo <- lm(tmp_y ~ tmp_x) # new linear fit
          mo_save <- mo
          tmp_velo_o <- coef(mo)[[2]] # velocity is cached
          if (coef(mo_save)[[2]] < (- (1 / 60)) &
            stran == FALSE) {
            # if the stranded option is inactive and ...
            # ...the slope is smaller - 1/60 the residuals are calculated as
            # if the slope was -1/60
            mo <- lm(tmp_y - I((-1 / 60)) * tmp_x ~ 1)
            tmp_velo_o <- (-1 / 60) # velocity is overwritten
          }
          if (coef(mo_save)[[2]] > (1 / 60)) {
            # if the slope is bigger 1/60, the residuals are calculated as
            # if the slope was 1/60
            mo <- lm(tmp_y - I(1 / 60) * tmp_x ~ 1)
            tmp_velo_o <- 1 / 60 # velocity is overwritten
          }
          tmp_inter_o <- coef(mo)[[1]] # intercept is cached
          # the sum of residuals is the temporary score with outlier
          # penalty times i
          tmp_score_o <- sum(abs(residuals(mo))) + pen * i
          tmp_tmp <- abs(residuals(mo))
          if (coef(mo_save)[[2]] < 0 &
            stran == TRUE) {
            # if the stranded option is active and ...
            # ...the slope is negative, the residuals are calculated as
            # if the slope was 0
            tmp_score_o <-
              sum(abs(tmp_y - (mean(tmp_y)))) + pen * i # the score is
            #overwritten
            tmp_velo_o <- 0 # velocity...
            tmp_inter_o <- mean(tmp_y) # .. and intercept are overwritten
            tmp_tmp <- abs(tmp_y - (mean(tmp_y)))
          }
          tmp_velo <-
            c(tmp_velo, tmp_velo_o) # velocity, intercept and score...
          tmp_inter <-
            c(tmp_inter, tmp_inter_o) # ...are put into one long vector
          tmp_score <- c(tmp_score, tmp_score_o)
          out <- c(out, sort(tmp_tmp, decreasing = TRUE)[1])
        }
      }
      mi <-
        which(tmp_score == min(tmp_score))[1] # the lowest score is chosen
      nam <- paste0(z, collapse = ",") # the name is all IDs divided by ","
      nam <- paste0(nam, "_", tmp_velo[mi]) # the velocity is pasted behind that
      #by "_"
      nam <- paste0(nam, "_", tmp_inter[mi]) # same with the intercept
      if (mi > 1) {
        # if multiple tmp_scores exist, outliers were found
        outlier <- paste(names(out)[seq_len(mi - 1)], collapse = ",")
        nam <- paste0(nam, "_", outlier) # the outliers are pasted behind the
        # name with "_"
      }
      res <- tmp_score[mi] # the final result is cached
      names(res) <- nam # and the name is added
    } else {
      # if only one value is given..
      nam <- paste0(z, collapse = ",") # ... the names are pasted together..
      nam <- paste0(nam, "_", "0") # ...the velocity is 0..
      nam <-
        paste0(nam, "_", mean(y)) # ... intercept is the mean of the input
      #values
      res <- 0 # and the score is 0 (perfect score)
      names(res) <- nam
    }
    res # the score is returned with additional info in the name
  }

