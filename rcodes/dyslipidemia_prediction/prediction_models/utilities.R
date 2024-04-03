## for binary variables, show counts and percentage
get_N_percentage <- function(x, digits = 1) {
  if (any(is.na(x))) {
    warning(sprintf("%s observation(s) missing", sum(is.na(x))))
  }
  paste0(sum(x, na.rm = T), " (", format(round(mean(x, na.rm = T) * 100, digits), nsmall = digits), "%)")
}

get_percentage_N <- function(x, digits = 1) {
  if (any(is.na(x))) {
    warning(sprintf("%s observation(s) missing", sum(is.na(x))))
  }
  paste0(format(round(mean(x, na.rm = T) * 100, digits), nsmall = digits),
         "% (", sum(x, na.rm = T), ")")
}

get_mean_sd <- function(x, digits = 2) {
  if (any(is.na(x))) {
    warning(sprintf("%s observation(s) missing", sum(is.na(x))))
  }
  paste0(format(round(mean(x, na.rm = T), digits), nsmall = digits), " (",
         format(round(sd(x, na.rm = T), digits), nsmall = digits), ")")
}

get_median_range <- function(x, digits = 2) {
  if (any(is.na(x))) {
    warning(sprintf("%s observation(s) missing", sum(is.na(x))))
  }
  paste0(format(round(median(x, na.rm = T), digits), nsmall = digits), " (",
         format(round(range(x, na.rm = T)[1], digits), nsmall = digits), ", ",
         format(round(range(x, na.rm = T)[2], digits), nsmall = digits), ")")
}


## get the summary statistics for a regression model
get_summary_reg <- function(reg, digits = 3, pval_digits = 3,
                            transformation = NULL, ...) {
  reg_tab <- data.frame(est = coef(reg),
                        ci_lb = confint.default(reg)[, 1],
                        ci_ub = confint.default(reg)[, 2],
                        pval = summary(reg)$coefficients[, "Pr(>|z|)"])
  
  if (!is.null(transformation)) {
    reg_tab <- reg_tab %>%
      mutate(est = transformation(est),
             ci_lb = transformation(ci_lb),
             ci_ub = transformation(ci_ub))
  }
  
  result <- with(reg_tab, data.frame(
    ESTIMATE =  format(round(est, digits), digits),
    CI = paste0(format(round(ci_lb, digits), digits), ", ",
                format(round(ci_ub, digits), digits)),
    PVAL = format(round(pval, pval_digits), pval_digits)
  ))
  
  rownames(result) <- names(coef(reg))
  return(result)
}
## calculate the year difference between startdate and stopdate
yrdif <- function(startdate, stopdate, type = 'actual'){
  ndate <- length(startdate)
  diff <- rep(NA, ndate)
  for (i in 1:ndate) {
    if (type == 'actual'){
      ys <- year(startdate[i])
      ye <- year(stopdate[i])
      
      if (is.na(ys) | is.na(ye)) {
        diff[i] <- NA
      } else {
        yslast <- as.Date(paste0(year(startdate[i]) + 1,'-01-', '01'))
        yelast <- as.Date(paste0(year(stopdate[i]) + 1,'-01-', '01'))
        
        if (ys == ye){
          
          if (leap_year(ys)) {
            diff[i] <- (stopdate[i] - startdate[i]) / 366
          }
          
          if (!leap_year(ys)){
            diff[i] <- (stopdate[i] - startdate[i]) / 365
          }
        }
        
        if (ys != ye){
          
          denom <- rep(1, length(ys:ye))
          
          if (leap_year(ys)) {denom[1] <- (yslast - startdate[i]) / 366}
          if (!leap_year(ys)) {denom[1] <- (yslast - startdate[i]) / 365}
          
          if (leap_year(ye)) {denom[length(ys:ye)] <- abs(1 - (yelast - stopdate[i]) / 366)}
          if (!leap_year(ye)) {denom[length(ys:ye)] <- abs(1 - (yelast - stopdate[i]) / 365)}
          
          diff[i] <- sum(denom)
          
        }
      }
    }
    
    if (type == '365'){
      diff[i] <- as.numeric(as.Date(stopdate[i]) - as.Date(startdate[i])) / 365
    }
  }
  return(diff)
}

## get race/ethnicity
get_race_ethnicity <- function(racegrp, hispanic) {
  race_ethnicity <- rep(NA, length(racegrp))
  
  for (i in 1:length(racegrp)) {
    racegrp_i <- racegrp[i]
    hispanic_i <- hispanic[i]
    if (racegrp_i == "White" & hispanic_i == "Non Hispanic/Latino") {
      race_ethnicity[i] <- "Non-Hispanic White"
    } else if (racegrp_i == "Black" & hispanic_i == "Non Hispanic/Latino") {
      race_ethnicity[i] <- "Non-Hispanic Black"
    } else if (hispanic_i == "Hispanic/Latino") {
      race_ethnicity[i] <- "Hispanic"
    } else if (racegrp_i == "Other") {
      race_ethnicity[i] <- "Other"
    } 
  }
  return(race_ethnicity)
}


## among a list of visits, match an mrn and a date with a specific visit
## method can be one of "closest", "closest_prior", "closest_after",
## "within", "within_prior", "within_after"
## action_missing_date_a can be one of "NA", "earliest" or "latest"
match_visit <- function(id = NULL, id_a = NULL, id_b = NULL,
                        data_a, data_b,
                        date_a, date_b,
                        date_name = date_a,
                        var_b = NULL,
                        method = "closest", 
                        action_missing_date_a = "earliest",
                        within = 180){
  
  
  if (!is.null(id)) {
    id_a <- id_b <- id
  }
  
  
  if (is.null(var_b)) {
    data_b1 <- data_b
  } else {
    data_b1 <- data_b %>% 
      filter(c(!is.na(data_b[, var_b]))) 
  }
  
  data_b_names <- setdiff(names(data_b1), c(id_b, date_b))
  
  data_out <- data.frame(person_name = data_a[, id_a],
                         person_visit_date = data_a[, date_a]) 
  
  
  for (i in 1:nrow(data_a)) {
    id_ai <- unlist(data_a[i, id_a])
    date_ai <- unlist(data_a[i, date_a])
    data_bi <- data_b1[data_b1[, id_b] == id_ai, ]
    
    if (is.na(date_ai)) {
      if (nrow(data_bi) > 0) {
        if (action_missing_date_a == "earliest") {
          id_select <- which.min(data_bi[, date_b])
          data_out[i, data_b_names] <- data_bi[id_select, data_b_names]
        } else if (action_missing_date_a == "lastest") {
          id_select <- which.max(data_bi[, date_b])
          data_out[i, data_b_names] <- data_bi[id_select, data_b_names]
        } else if (action_missing_date_a == "na") {
          data_out[i, data_b_names] <- NA
        }
      }
    } else if (nrow(data_bi) > 0) {
      diff_date <- as.Date(date_ai) - 
        as.Date(unlist(data_bi[, date_b]))
      if (method == "closest") {
        id_select <- which.min(abs(diff_date))[1]
        data_out[i, data_b_names] <- data_bi[id_select, data_b_names]
      }
      
      if (method == "within") {
        elig_dates_n <- sum(abs(diff_date) <= within)
        if (elig_dates_n >= 1) {
          id_select <- which.min(abs(diff_date))[1]
          data_out[i, data_b_names] <- data_bi[id_select, data_b_names]
        }
      }
    }
  }
  return(data_out)
}



###################################################
## Derived Variables
###################################################
## get participants' SJLIFE status
get_status <- function(x) {
  y <- rep(NA, length(x))
  
  for (i in seq_along(x)) {
    if (x[i] == 3) {
      y[i] <- 1 ## Eligible survivor participants
    } else if (x[i] %in% c(8, 22)) {
      y[i] <- 2 ## SJLIFE ineligible
    } else if (x[i] %in% c(15, 18)) {
      y[i] <- 3 ## Not recruited
    } else if (x[i] %in% c(1, 2, 13, 21, 11)) {
      y[i] <- 4 ## No visited yet
    } else if (x[i] %in% c(17, 20)) {
      y[i] <- 5 ## died prior to visit
    } else if (x[i] %in% c(4, 5, 24)) {
      y[i] <- 6 ## survey only
    } else if (x[i] %in% c(9, 12, 19, 23)) {
      y[i] <- 7 ## Refused
    } else if (x[i] %in% c(6, 7, 10, 14, 99)) {
      y[i] <- 8 ## Nonresponse
    }
  }
  
  return(y)
}



get_hsct <- function(sjlife_path = "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/",
                     mrn, dates, datesname = "date") {
  transplant <- read_sas(paste0(sjlife_path, "Clinical data/transplant.sas7bdat"))%>%
    rename_with(tolower) %>%
    arrange(mrn)
  
  data_out <- data.frame(mrn = mrn, dates = dates)
  names(data_out)[2] <- datesname
  data_out$hsc_2g <- NA
  
  hsc_tran <- transplant %>%
    dplyr::select(mrn, tpdt, tptype)
  
  for (i in 1:nrow(data_out)) {
    hsc_tran_i <- filter(hsc_tran,
                         mrn == data_out$mrn[i] &
                           tpdt <= data_out[, datesname][i])
    if (nrow(hsc_tran_i) > 1) {
      data_out$hsc_2g[i] <- 1
    } else {
      data_out$hsc_2g[i] <- 0
    }
  }
  
  return(data_out)
}


## adolescent physical activity from health habits survey
get_adolescent_pa <- function(sjlife_path = "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/") {
  adolescent_healthhabits <- read_sas(paste0(sjlife_path, "Survey data/adolescent_healthhabits.sas7bdat")) %>%
    rename_with(tolower)
  
  ## compute MET for each adolescent with MRN in the study population
  pa_adoles <- adolescent_healthhabits %>%
    dplyr::select(mrn, datecomp, percomp,vpa10, vpadays, vpamin,
           mpa10, mpadays, mpamin, lpa10, lpadays, lpamin,
           activity, play60, pa20days, patone,
           screen, relation, survey) %>%
    mutate(wvpa = NA, wmpa = NA, mvpawk = NA, wmspa = NA, mets = NA)
  
  
  for (i in 1:nrow(pa_adoles)) {
    if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 1 & is.na(pa_adoles$vpadays[i])) {
      pa_adoles$vpadays[i] <- 1
    }
    
    if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 1 & is.na(pa_adoles$vpamin[i])) {
      pa_adoles$vpamin[i] <- 10
    }
    
    if (!is.na(pa_adoles$vpamin[i]) & pa_adoles$vpamin[i] > 360) {
      pa_adoles$vpamin[i] <- 360
    }
    
    if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 1) {
      pa_adoles$wvpa[i] = pa_adoles$vpadays[i] * pa_adoles$vpamin[i]
    }
    
    if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2) {
      pa_adoles$wvpa[i] = 0
    }
    
    if (is.na(pa_adoles$wvpa[i])) {
      if ((is.na(pa_adoles$vpa10[i]) | (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2)) &
          (!is.na(pa_adoles$activity[i]) & pa_adoles$activity[i] == 2)) {
        pa_adoles$wvpa[i] <- 0
      } 
      
      if ((is.na(pa_adoles$vpa10[i])|(!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2)) &
          !(is.na(pa_adoles$play60[i])|(!is.na(pa_adoles$play60[i]) & pa_adoles$play60[i] == 1))) {
        pa_adoles$wvpa[i] <- 0
      }
      
      if ((is.na(pa_adoles$vpa10[i]) | (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2)) &
          !(is.na(pa_adoles$pa20days[i]) | (!is.na(pa_adoles$pa20days[i]) & pa_adoles$pa20days[i] == 1))) {
        pa_adoles$wvpa[i] <- (pa_adoles$pa20days[i] - 1) * 20
      }
    }
    
    if (is.na(pa_adoles$mpa10[i]) & 
        ((!is.na(pa_adoles$mpadays[i]))|(!is.na(pa_adoles$mpamin[i])))) {
      pa_adoles$mpa10[i] <- 1
    }
    
    if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 1 & is.na(pa_adoles$mpadays[i])) {
      pa_adoles$mpadays[i] <- 1
    }
    
    if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 1 & is.na(pa_adoles$mpamin[i])) {
      pa_adoles$mpamin[i] <- 10
    }
    
    if (!is.na(pa_adoles$mpamin[i]) & pa_adoles$mpamin[i] > 360) pa_adoles$mpamin[i] <- 360
    if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 1) {
      pa_adoles$wmpa[i] <- 
        pa_adoles$mpadays[i] * pa_adoles$mpamin[i]
    }
    if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 2) {
      pa_adoles$wmpa[i] <-0
    }
    
    if (is.na(pa_adoles$wmpa[i])) {
      if ((is.na(pa_adoles$mpa10[i]) | pa_adoles$mpa10[i] == 2) &
          (!is.na(pa_adoles$activity[i]) & pa_adoles$activity[i] == 2)) {
        pa_adoles$wmpa[i] <- 0
      } 
      
      if ((is.na(pa_adoles$mpa10[i]) | pa_adoles$mpa10[i] == 2) &
          !(is.na(pa_adoles$play60[i]) | pa_adoles$play60[i] == 1)) {
        pa_adoles$wmpa[i] <- (pa_adoles$play60[i] - 1) * 60
      } 
      
      if (!is.na(pa_adoles$wvpa[i])) pa_adoles$wmpa[i] <- 0
    }
    
    
    if (is.na(pa_adoles$wmpa[i]) & is.na(pa_adoles$wvpa[i])) {
      if (!is.na(pa_adoles$play60[i]) & pa_adoles$play60[i] == 1) {
        pa_adoles$wvpa[i] <- pa_adoles$wmpa[i] <- 0
      }
    }
    
    if (!is.na(pa_adoles$patone[i]) & pa_adoles$patone[i] > 1) {
      pa_adoles$wmspa[i] <- (pa_adoles$patone[i] - 1) * 20
    } else {
      pa_adoles$wmspa[i] <- 0
    }
    
    pa_adoles$mvpawk[i] <- pa_adoles$wmpa[i] + 2 * pa_adoles$wvpa[i] +
      2 * pa_adoles$wmspa[i]
    
    if (is.na(pa_adoles$wvpa[i])) pa_adoles$mvpawk[i] <- pa_adoles$wmpa[i]
    
    if (!is.na(pa_adoles$mvpawk[i]) & pa_adoles$mvpawk[i] > 2520) pa_adoles$mvpawk[i] <- 2520
    
    pa_adoles$mets[i] <- 3 * pa_adoles$mvpawk[i]
    
  }
  
  pa_adoles_final <- dplyr::select(pa_adoles, mrn, datecomp, relation, percomp, 
                            survey, screen, wvpa, wmpa, mvpawk, mets)
  return(pa_adoles_final)
}

## adult pa from health habits survey
get_adult_pa <- function(sjlife_path = "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/") {
  adult_healthhabits <- read_sas(paste0(sjlife_path, "Survey data/adult_healthhabits.sas7bdat")) %>%
    rename_with(tolower)
  
  pa_adults <- adult_healthhabits %>%
    dplyr::select(mrn, datecomp, relation, nopa,
           vpa10, vpadays, vpamin, mpa10, mpadays, mpamin,
           lpa10, lpadays, lpamin, ltpaw, wtlt, yoga, pa20) %>%
    mutate(wvpa = NA, wmpa = NA, mvpawk = NA, mets = NA)
  
  
  for (i in 1:nrow(pa_adults)) {
    if (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 1 & is.na(pa_adults$vpadays[i])) {
      pa_adults$vpadays[i] <- 1
    }
    
    if (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 1 & is.na(pa_adults$vpamin[i])) {
      pa_adults$vpamin[i] <- 10
    }
    
    if (!is.na(pa_adults$vpamin[i]) & pa_adults$vpamin[i] > 360) {
      pa_adults$vpamin[i] <- 360
    }
    
    if (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 1) {
      pa_adults$wvpa[i] <- pa_adults$vpadays[i] * pa_adults$vpamin[i]
    }
    
    if (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 2) {
      pa_adults$wvpa[i] <- 0
    }
    
    if (is.na(pa_adults$wvpa[i])) {
      if ((is.na(pa_adults$vpa10[i]) | (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 2)) &
          (!is.na(pa_adults$nopa[i]) & pa_adults$nopa[i] == 1) &
          !(is.na(pa_adults$pa20[i]) | (!is.na(pa_adults$pa20[i]) & pa_adults$pa20[i] == 1))) {
        pa_adults$wvpa[i] <- 0
      } 
      
      if ((is.na(pa_adults$vpa10[i])|(!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 2)) &
          !(is.na(pa_adults$pa20[i])|(!is.na(pa_adults$pa20[i]) & pa_adults$pa20[i] == 1))) {
        pa_adults$wvpa[i] <- (pa_adults$pa20[i] - 1) * 20
      }
      
      if ((is.na(pa_adults$vpa10[i]) | (!is.na(pa_adults$vpa10[i]) & pa_adults$vpa10[i] == 2)) &
          (!is.na(pa_adults$nopa[i]) & pa_adults$nopa[i] == 2) &
          (is.na(pa_adults$pa20[i]) | (!is.na(pa_adults$pa20[i]) & pa_adults$pa20[i] == 1))) {
        pa_adults$wvpa[i] <- 0
      }
    }
    
    if (is.na(pa_adults$mpa10[i]) & 
        ((!is.na(pa_adults$mpadays[i]))|(!is.na(pa_adults$mpamin[i])))) {
      pa_adults$mpa10[i] <- 1
    }
    
    if (!is.na(pa_adults$mpa10[i]) & pa_adults$mpa10[i] == 1 & is.na(pa_adults$mpadays[i])) {
      pa_adults$mpadays[i] <- 1
    }
    
    if (!is.na(pa_adults$mpa10[i]) & pa_adults$mpa10[i] == 1 & is.na(pa_adults$mpamin[i])) {
      pa_adults$mpamin[i] <- 10
    }
    
    if (!is.na(pa_adults$mpamin[i]) & pa_adults$mpamin[i] > 360) pa_adults$mpamin[i] <- 360
    if (!is.na(pa_adults$mpa10[i]) & pa_adults$mpa10[i] == 1) {
      pa_adults$wmpa[i] <- 
        pa_adults$mpadays[i] * pa_adults$mpamin[i]
    }
    if (!is.na(pa_adults$mpa10[i]) & pa_adults$mpa10[i] == 2) {
      pa_adults$wmpa[i] <-0
    }
    
    if (is.na(pa_adults$wmpa[i])) {
      if ((is.na(pa_adults$mpa10[i]) | 
           (!is.na(pa_adults$mpa10[i]) & pa_adults$mpa10[i] == 2)) &
          (!is.na(pa_adults$nopa[i]) & pa_adults$nopa[i] == 1)) {
        pa_adults$wmpa[i] <- 0
      } 
      
      if (!is.na(pa_adults$wvpa[i])) {
        pa_adults$wmpa[i] <- 0
      } 
      
    }
    
    
    
    
    pa_adults$mvpawk[i] <- pa_adults$wmpa[i] + 2 * pa_adults$wvpa[i] 
    
    if (!is.na(pa_adults$mvpawk[i]) & pa_adults$mvpawk[i] > 2520) pa_adults$mvpawk[i] <- 2520
    
    pa_adults$mets[i] <- 3 * pa_adults$mvpawk[i]
    
  }
  
  pa_adults_final <- dplyr::select(pa_adults, mrn, datecomp, relation, wmpa, wvpa, mvpawk, mets)
  return(pa_adults_final)
}

get_offspring <- function(stlife_path = "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/",
                          mrn_list) {
  adult_home <- read_sas(paste0(sjlife_path, "Survey data/adult_home.sas7bdat")) %>%
    rename_with(tolower) %>%
    filter(mrn %in% mrn_list)
  
  everpreg <- adult_home$everpreg
  n_preg <- adult_home$n_preg
  n_offspring <- rep(NA, length(everpreg))
  pregout1 <- adult_home$pregout1
  pregout2 <- adult_home$pregout2
  pregout3 <- adult_home$pregout3
  pregout4 <- adult_home$pregout4
  pregout5 <- adult_home$pregout5
  pregout6 <- adult_home$pregout6
  pregout7 <- adult_home$pregout7
  pregout8 <- adult_home$pregout8
  
  for (i in 1:length(n_offspring)) {
    if (everpreg[i] %in% 2) {
      n_preg[i] <- 0
      n_offspring[i] <- 0
    } else if (n_preg[i] %in% 1) {
      n_offspring[i] <- sum(pregout1[i] == 1)
    } else if (n_preg[i] %in% 2) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) 
    } else if (n_preg[i] %in% 3) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) +
        as.numeric(pregout3[i] == 1)
    } else if (n_preg[i] %in% 4) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) +
        as.numeric(pregout3[i] == 1) +
        as.numeric(pregout4[i] == 1)
    } else if (n_preg[i] %in% 5) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) +
        as.numeric(pregout3[i] == 1) +
        as.numeric(pregout4[i] == 1) + 
        as.numeric(pregout5[i] == 1)
    } else if (n_preg[i] %in% 6) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) +
        as.numeric(pregout3[i] == 1) +
        as.numeric(pregout4[i] == 1) + 
        as.numeric(pregout5[i] == 1) +
        as.numeric(pregout6[i] == 1)
    } else if (n_preg[i] %in% 7) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) +
        as.numeric(pregout3[i] == 1) +
        as.numeric(pregout4[i] == 1) + 
        as.numeric(pregout5[i] == 1) +
        as.numeric(pregout6[i] == 1) +
        as.numeric(pregout7[i] == 1)
    } else if (is.na(n_preg[i]) | n_preg[i] >= 8) {
      n_offspring[i] <- sum(pregout1[i] == 1, pregout2[i] == 1,
                            pregout3[i] == 1, pregout4[i] == 1,
                            pregout5[i] == 1, pregout6[i] == 1,
                            pregout7[i] == 1, pregout8[i] == 1,
                            na.rm = T)
    }
  }
  
  return(data.frame(mrn = adult_home$mrn,
                    datecomp = adult_home$datecomp,
                    relation = adult_home$relation,
                    everpreg = everpreg,
                    n_preg = n_preg,
                    n_offspring = n_offspring))
}
## method: closest, most recent, or past max

get_ctcae_grade <- function(stlife_path = "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/",
                            mrn, dates, conditions, 
                            method = "most recent", within = 180) {
  ctcaegrades <- read_sas(paste0(sjlife_path, "Event data/ctcaegrades.sas7bdat")) %>% 
    rename_with(tolower) %>%
    arrange(mrn) %>%
    filter(!grade %in% c(NA, -9))
  
  data_out <- data.frame(mrn = mrn, 
                         dates = dates) %>%
    cbind(matrix(NA, nrow = length(mrn), ncol = length(conditions)))
  
  names(data_out)[-c(1, 2)] <- conditions %>%
    tolower() %>%
    stringr::str_replace_all(" ", "_")
  
  for (cc in 1:length(conditions)) {
    sub_grades <- ctcaegrades %>% 
      filter(condition == conditions[cc]) %>%
      dplyr::select(mrn, gradedt, grade)
    
    for (i in 1:nrow(data_out)) {
      date_i <- data_out$dates[i]
      sub_grades_i <- filter(sub_grades, mrn == data_out$mrn[i])
      
      if (nrow(sub_grades_i) > 0) {
        diff_date <- date_i - sub_grades_i$gradedt
        if (method == "closest") {
          id_select <- which.min(abs(diff_date))[1]
          data_out[i, 2 + cc] <- sub_grades_i[id_select, "grade"]
        }
        
        if (method == "most recent") {
          elig_dates_n <- sum(diff_date >= -within)
          if (elig_dates_n >= 1) {
            id_select <- which(diff_date == min(diff_date[diff_date >= -within]))[1]
            data_out[i, 2 + cc] <- sub_grades_i[id_select, "grade"]
          }
        }
      }
    }
  }
  
  return(data_out)
}
