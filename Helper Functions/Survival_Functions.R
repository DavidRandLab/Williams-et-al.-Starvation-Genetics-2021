# removes the t in txxhrs from count sheets 
qf<- function(x){
  sub('t','',x)
}

# depreciated
# splits the levels in survival output, gets rid of the 'levelname=' level 
surv.delevel<- function(x){
  sub('.*=','',x,perl=T)
}

# pipable version of linear model
pipeable.lm <- function(data,...){
  lm(..., data = data )
}

# function to return the granualar survival innfo for ALL levels in the model
# provides the input for log mortality function
# requires survival object to have been run, input is summary of fit
get_survival.time_info<- function(surv.sum, strata = surv.sum$strata[1]){
  age.main = surv.sum$time[surv.sum$strata == strata]
  n.event.main = surv.sum$n.event[surv.sum$strata == strata]
  n.risk.main = surv.sum$n.risk[surv.sum$strata == strata]
  tibble(strata = as.character(strata), age= age.main,
         n.event.main = n.event.main, n.risk.main = n.risk.main) %>%
    return()
}

# surv.sum<- summary(x)
# 
# surv_complete <- tibble()
# for (n in surv.sum$strata){
#   surv_complete <- surv_complete %>%
#     bind_rows(
#       get_survival.time_info(surv.sum, strata = n))
# }

log.mort <- function(data, strata){
  ### log mortality 
  surv.sum<- summary(data)
  
  surv_complete <- tibble()
  for (n in surv.sum$strata){
    surv_complete <- surv_complete %>%
      bind_rows(
        get_survival.time_info(surv.sum, strata = n))
  }
  
  surv_complete %>%
    mutate(qx.main = n.event.main/n.risk.main) %>%
    mutate(px.main = 1 - qx.main) %>%
    mutate(logMortality = -log10(px.main)) %>%
    filter(!is.infinite(logMortality))%>%
    separate(strata, strata ,', ') %>%
    mutate_at(.vars = strata, .funs = ~str_replace_all(.,'.*=','')) %>%
    #mutate_at(.vars= vars(age:logMortality), .funs = as.numeric) %>%
    
    return()
}
  
  
# function to plot survival from survival fit object 


survival_df <- function(x){
  x.sum <- summary(x) # summarize fit 
  Condition <- x.sum$strata %>%
    as.character() 
  
  # create data.frame with points to plot 
  df <- data.frame(Condition, time = x$time, survival = x$surv) 
  df <- mutate(df, Condition = factor(Condition))
  # add 0 time columns t df 
  for ( levs in levels(df$Condition)){
    df <- df %>%
      add_row(Condition = levs, time = 0, survival = 1 )
  }
  return(df)
}

