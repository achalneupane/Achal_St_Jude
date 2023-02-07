# smoker_former_or_never_yn >> Current_smoker_yn
# NOT_RiskyHeavyDrink_yn > RiskyHeavyDrink_yn
# Not_obese_yn >> Obese_yn

edit_lifestyle <- function(ALL.LIFESTYLE){
ALL.LIFESTYLE$Current_smoker_yn[ALL.LIFESTYLE$smoker_former_or_never_yn == 0] <- "Yes"
ALL.LIFESTYLE$Current_smoker_yn[ALL.LIFESTYLE$smoker_former_or_never_yn == 1] <- "No"
ALL.LIFESTYLE$Current_smoker_yn[ALL.LIFESTYLE$smoker_former_or_never_yn == "Unknown"] <- "Unknown"
ALL.LIFESTYLE$Current_smoker_yn <- factor(ALL.LIFESTYLE$Current_smoker_yn, level = c("No", "Yes", "Unknown")) 


ALL.LIFESTYLE$RiskyHeavyDrink_yn[ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn == 0] <- "Yes"
ALL.LIFESTYLE$RiskyHeavyDrink_yn[ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn == 1] <- "No"
ALL.LIFESTYLE$RiskyHeavyDrink_yn[ALL.LIFESTYLE$NOT_RiskyHeavyDrink_yn == "Unknown"] <- "Unknown"
ALL.LIFESTYLE$RiskyHeavyDrink_yn <- factor(ALL.LIFESTYLE$RiskyHeavyDrink_yn, level = c("No", "Yes", "Unknown")) 

ALL.LIFESTYLE$Obese_yn[ALL.LIFESTYLE$Not_obese_yn == 0] <- "Yes"
ALL.LIFESTYLE$Obese_yn[ALL.LIFESTYLE$Not_obese_yn == 1] <- "No"
ALL.LIFESTYLE$Obese_yn[ALL.LIFESTYLE$Not_obese_yn == "Unknown"] <- "Unknown"
ALL.LIFESTYLE$Obese_yn <- factor(ALL.LIFESTYLE$Obese_yn, level = c("No", "Yes", "Unknown")) 

ALL.LIFESTYLE$PhysicalActivity_yn <- as.character(ALL.LIFESTYLE$PhysicalActivity_yn)
ALL.LIFESTYLE$PhysicalActivity_yn[ALL.LIFESTYLE$PhysicalActivity_yn == 1] <- "Yes"
ALL.LIFESTYLE$PhysicalActivity_yn[which(ALL.LIFESTYLE$PhysicalActivity_yn == 0)] <- "No"
ALL.LIFESTYLE$PhysicalActivity_yn[ALL.LIFESTYLE$PhysicalActivity_yn == "Unknown"] <- "Unknown"
ALL.LIFESTYLE$PhysicalActivity_yn <- factor(ALL.LIFESTYLE$PhysicalActivity_yn, level = c("Yes", "No", "Unknown")) 

ALL.LIFESTYLE$HEALTHY_Diet_yn <- as.character(ALL.LIFESTYLE$HEALTHY_Diet_yn)
ALL.LIFESTYLE$HEALTHY_Diet_yn[ALL.LIFESTYLE$HEALTHY_Diet_yn == 1] <- "Yes"
ALL.LIFESTYLE$HEALTHY_Diet_yn[which(ALL.LIFESTYLE$HEALTHY_Diet_yn == 0)] <- "No"
ALL.LIFESTYLE$HEALTHY_Diet_yn[ALL.LIFESTYLE$HEALTHY_Diet_yn == "Unknown"] <- "Unknown"
ALL.LIFESTYLE$HEALTHY_Diet_yn <- factor(ALL.LIFESTYLE$HEALTHY_Diet_yn, level = c("Yes", "No", "Unknown")) 

return(ALL.LIFESTYLE)
}


# table(ALL.LIFESTYLE$Obese_yn)
# table(ALL.LIFESTYLE$Not_obese_yn)
