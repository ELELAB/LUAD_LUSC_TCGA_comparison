
#function to get patient code
get_patient <- function(list_barcodes){
  IDs <- strsplit(c(as.character(list_barcodes)), "-")
  IDs <- ldply(IDs, rbind)
  IDs$patient <- apply(IDs[,c(1,2,3)],1,paste,collapse = "-" )
  return(IDs$patient)
}

# get survival time: if dead patient time=days_to_death, 
# if it is alive time:days_to_last_followup
get_survival_time <- function(clinical){
  
  for(i in 1:nrow(clinical)){
    if(clinical$vital_status[i]=="dead")
      clinical$time[i] <- clinical$days_to_death[i]
    else
      clinical$time[i] <- clinical$days_to_last_follow_up[i]
  }
  return(clinical)
}

# add the label "low" and "high"
get_type_sample <-function(x,patient_up,patient_down){
  if(x %in% patient_down){v <- "low"}
  else 
    if(x %in% patient_up)
    {v <-"high"}
}

get_age_range <- function(x,median){
  
  if(x > median){v <- paste0("> ",median)}
  else
  {v <- paste0("< ",median)}
}


#--------------------------------------------------------------------------------------

# method: "cox" or "KM" (Kaplan-Meier)

get_survival_table <- function(dataFilt, clinical, gene, threshDown, threshUp, method){
  
  #get only tumor samples
  dataSmTP <- TCGAquery_SampleTypes(colnames(dataFilt),"TP")
  
  #calculate quantile to be used as threshold to define two groups of samples to compare
  quantile_down <- quantile(dataFilt[gene,dataSmTP],threshDown)
  quantile_up <- quantile(dataFilt[gene,dataSmTP],threshUp)
  #print(paste0("the lower percentile is ",quantile_down," and the higher percentile is ",quantile_up))
  samples_down <- colnames(dataFilt[,dataSmTP])[which(dataFilt[gene,dataSmTP]< quantile_down)] 
  samples_up <- colnames(dataFilt[,dataSmTP])[which(dataFilt[gene,dataSmTP]> quantile_up)]
  patient_down <- get_patient(samples_down)
  patient_up <- get_patient(samples_up)
  #print(paste0("there are ",length(intersect(patient_up,patient_down))," patients in common between two groups:",intersect(patient_up,patient_down)))
  
  clinical <- subset(clinical, clinical$bcr_patient_barcode %in% union(patient_up,patient_down),
                     select = c("bcr_patient_barcode","vital_status","days_to_death","days_to_last_follow_up","gender","tumor_stage","age_at_diagnosis","cigarettes_per_day"))
  #add survival time in clinical data
  clinical <- get_survival_time(clinical)
  # get time in years
  clinical$time <- clinical$time/365
  clinical$group <- as.character(lapply(clinical$bcr_patient_barcode,function(x) get_type_sample(x,patient_up,patient_down)))
  clinical$age_at_diagnosis <- (clinical$age_at_diagnosis)/365
  
  #include age-info and cigarettes-info only if we decide to perform cox regression
  if(method=="cox"){
  clinical <- subset(clinical,!is.na(clinical$age_at_diagnosis))
  median_age <- floor(median(clinical$age_at_diagnosis))
  clinical$age <- as.character(lapply(clinical$age_at_diagnosis,function(x) get_age_range(x,median_age)))
  #clinical <- subset(clinical,!is.na(clinical$cigarettes_per_day))
  #median_cigarettes <- floor(median(clinical$cigarettes_per_day))
  #clinical$cigarettes <- as.character(lapply(clinical$cigarettes_per_day,function(x) get_age_range(x,median_cigarettes)))
  }
  # get vital_status: dead=1, alive=0
  clinical$status <- as.numeric(lapply(clinical$vital_status, FUN=function(x) if(x=="dead") v<-1 else v<-0))
  clinical <- subset(clinical, !is.na(clinical$time))
  #print(paste0("there are ",length(which(clinical$group=="low"))," patients in the low group"))
  #print(paste0("and ",length(which(clinical$group=="high"))," patients in the high group"))
  
  return(clinical)
}

survival_plot <- function(clinical,gene,cancer_type,percentile){
  
  surv_object<-survfit(Surv(time,status)~group, data=clinical)
  # Log-rang test
  sdf <- survdiff(Surv(time,status)~group, data = clinical)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  
  #pdf(paste0("survival_plot_",gene,"_",cancer_type,"_",percentile,".pdf"), width=9, height=7)
  myPlot <- autoplot(surv_object, conf.int = FALSE, ylab = "% surviving",xlab = "time (years)", surv.size = 1.5)+
    ggtitle(paste0("Logrank p-value=",signif(p.val)))+  
    theme(axis.title=element_text(size=20),axis.text=element_text(size=20),
          legend.text=element_text(size=20),legend.title=element_blank(),
          plot.title = element_text(hjust = 0.5,size=20), plot.subtitle = element_text(size=18),
          panel.background = element_rect(fill = "white",colour = "black",color = "black"))+
    scale_color_manual(values = c("red","blue"))
    #scale_fill_manual(values = c("red","blue"), breaks=c("high","low"),
                      #labels=c(paste0("high (n=",length(which(clinical$group=="high")),")"),
                               #paste0("low (n=",length(which(clinical$group=="low")),")")))
  ggsave(filename = paste0("survival_plot_",gene,"_",cancer_type,"_",percentile,".pdf"),plot = myPlot)
}

survival_plot1 <- function(clinical,gene,cancer_type,percentile){
  
  surv_object<-survfit(Surv(time,status)~group, data=clinical)
  # Log-rang test
  sdf <- survdiff(Surv(time,status)~group, data = clinical)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  plot(surv_object,col=c("red","blue"))
  legend("topright", c("high","low"), lty=c(1,1),col=c("black","grey"))
}


get_signif_list <- function(dataFilt,clinical,gene_name,threshDown, threshUp,method,list){
  clinical <- get_survival_table(dataFilt,clinical,gene_name,threshDown, threshUp,method)
  # including tumor stage
  clinical <- subset(clinical, !clinical$tumor_stage=="not reported")
  a <- unlist(lapply(as.character(clinical$tumor_stage), function(x) strsplit(x," ")[[1]][[2]]))
  #remove "a", "b" and "c" labels from stages
  s <- gsub("a","",a)
  s <- gsub("b","",s)
  s <- gsub("c","",s)
  clinical$tumor_stage <- paste0("stage"," ",s)
  cox.model <- coxph(Surv(time,status)~group+tumor_stage+age+gender,data=clinical)
  #test the proportionality
  time_test <- cox.zph(cox.model)
  test_result <- as.vector(time_test$table[,"p"]>0.05)
  if(length(which(test_result=="TRUE"))==length(test_result)){
    list <- append(list,gene_list[i])
  }
  return(list)
}



########################################################

cox_regression <- function(clinical,dataFilt,gene_name, threshDown, threshUp,method,final_table,num){
  
  clinical <- get_survival_table(dataFilt,clinical,gene_name,threshDown, threshUp,method)
  # including tumor stage
  clinical <- subset(clinical, !clinical$tumor_stage=="not reported")
  a <- unlist(lapply(as.character(clinical$tumor_stage), function(x) strsplit(x," ")[[1]][[2]]))
  #remove "a", "b" and "c" labels from stages
  s <- gsub("a","",a)
  s <- gsub("b","",s)
  s <- gsub("c","",s)
  clinical$tumor_stage <- paste0("stage"," ",s)
  
  cox.model <- coxph(Surv(time,status)~group+tumor_stage+age+gender,data=clinical)
    sum_cox <- summary(cox.model)$coefficients
    sum_cox <- sum_cox[,-c(3,4)]
    colnames(sum_cox) <- c("coef","exp(coef)","p-value")
    sum_cox <- signif(sum_cox,digits = 3)
    fdr <- signif(p.adjust(sum_cox[,"p-value"],method = "BH",n = num),digits = 3)
    sum_cox <- cbind(sum_cox,fdr)
    #vector <- c(sum_cox[1:nrow(sum_cox),])
    vector <- c(sum_cox[1,],sum_cox[2,],sum_cox[3,],sum_cox[4,],sum_cox[5,],sum_cox[6,])
    names(vector) <- c("coef(group_low)","exp(coef)(group_low)","p-value(group_low)","fdr(group_low)",
                       "coef(stageii)","exp(coef)(stageii)","p-value(stageii)","fdr(stageii)",
                       "coef(stageiii)","exp(coef)(stageiii)","p-value(stageiii)","fdr(stageiii)",
                       "coef(stageiv)","exp(coef)(stageiv)","p-value(stageiv)","fdr(stageiv)",
                       "coef(age>65)","exp(coef)(age>65)","p-value(age>65)","fdr(age>65)",
                       "coef(gendermale)","exp(coef)(gendermale)","p-value(gendermale)","fdr(gendermale)")
    final_table <- rbind(final_table,vector)#}
  
    rownames(final_table)[nrow(final_table)] <- gene_name
    
  return(final_table)
  
}