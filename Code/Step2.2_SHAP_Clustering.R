
#-------------------------------------------------------------------------------------------------
# Identify subgroups of recreational activity into subgroups
# ------------------------------------------------------------------------------------------------

in.path<-"./SHAP_Result/"

f.lst <- list.files(path= in.path, pattern = "*_rf_training.n.mean.csv", all.files = TRUE, full.names = TRUE, recursive = TRUE)

shap.dt.tr <- read.csv(file=paste0(in.path, "../Data/ph.acts.binary_tr_29act.csv"))
shap.dt.tt <- read.csv(file=paste0(in.path, "../Data/ph.acts.binary_tt_29act.csv"))

#----------------------------------------------------------------------
# dt1, dt2 contain original binary sai.p data
# fn.shap1, fn.shap2 contain filenames for the calculated SHAP values

# use all available data
extract.shap.test<-function(dt1=shap.dt.1, dt2=shap.dt.2, fn.shap1, fn.shap2, out.dir, out.fn, type.outcome="cbcl"){
  shap.tr<-shap.tt<-shap.tr0<-shap.tt0<-NULL
  
  # this reads in data with the SHAP values
  dt.shap1 <- read.csv(paste0(fn.shap1))
  dt.shap2 <- read.csv(paste0(fn.shap2))
  sai.nm <- colnames(dt.shap1)[grep("sai.p",colnames(dt.shap1))]
  
  tt12.mean<-tt22.mean<-NULL
  tt.count<-matrix(0,ncol=4, nrow=length(sai.nm))
  tt.minmax<-matrix(0,ncol=8, nrow=length(sai.nm))
  for(i in sai.nm){
    sel.id0<- intersect(as.character(dt.shap1$src_subject_id), as.character(dt1$src_subject_id[which(dt1[[i]]==0)]))
    sel.id20<- intersect(as.character(dt.shap2$src_subject_id), as.character(dt2$src_subject_id[which(dt2[,i]==0)]))
    tt12.mean<-c(tt12.mean, dt.shap1[[i]][which(dt.shap1$src_subject_id %in% sel.id0)])
    tt22.mean<-c(tt22.mean, dt.shap2[[i]][which(dt.shap2$src_subject_id %in% sel.id20)])
  }
  tt12.mean <- mean(tt12.mean)
  tt22.mean <- mean(tt22.mean)
  
  # print(c(tt12.mean,tt22.mean))
  
  p.cut <- 0.05/length(sai.nm)
  k<-0
  for(i in sai.nm){
    k<-k+1
    print(i)
    # calculate the proportion of participants endorsing the activity
    
    # tt11 has length of 20 around the proportion of participant endorsing this activity from top
    # among users, take shap values across full ranges
    sel.id1<- intersect(as.character(dt.shap1$src_subject_id), as.character(dt1$src_subject_id[which(dt1[,i]==1)]))
    sel.id0<- intersect(as.character(dt.shap1$src_subject_id), as.character(dt1$src_subject_id[which(dt1[,i]==0)]))
    sel.id21<- intersect(as.character(dt.shap2$src_subject_id), as.character(dt2$src_subject_id[which(dt2[,i]==1)]))
    sel.id20<- intersect(as.character(dt.shap2$src_subject_id), as.character(dt2$src_subject_id[which(dt2[,i]==0)]))
    
    tt11<-dt.shap1[[i]][which(dt.shap1$src_subject_id %in% sel.id1)]
    tt12<-dt.shap1[[i]][which(dt.shap1$src_subject_id %in% sel.id0)]
    
    tt21<-dt.shap2[[i]][which(dt.shap2$src_subject_id %in% sel.id21)]
    tt22<-dt.shap2[[i]][which(dt.shap2$src_subject_id %in% sel.id20)]
    
    tt.count[k,1]<- max(length(which(tt11>0)), length(which(tt11<0)))
    tt.count[k,2]<- max(length(which(tt12>0)), length(which(tt12<0)))
    tt.count[k,3]<- max(length(which(tt21>0)), length(which(tt21<0)))
    tt.count[k,4]<- max(length(which(tt22>0)), length(which(tt22<0)))
    
    tt.minmax[k,] <-c(range(tt12),range(tt11),range(tt22),range(tt21))
    # test whether there are significant differences in SHAP values between these with and without engaged in the activity
    tryCatch({
      #----- modify 05/23/2024 -----#
      temp <- merge(dt.shap1[c("src_subject_id", i)], 
                    dt.grp.tr[c("src_subject_id", i, "rel_family_id", "site_id_l")] %>%
                      setNames(c("src_subject_id", "grp.index", "rel_family_id", "site_id_l")))
      formula <- paste0(i, " ~ grp.index + (1|rel_family_id)")
      mod <- lmer(formula, data=temp, REML = FALSE)
      tt11.pvalue <- summary(mod)$coefficients["grp.index", "Pr(>|t|)"]
      print(tt11.pvalue)
      #----- END -----#
    }, error = function(e) {
      # In case of error, set p-value to 1
      tt11.pvalue <- 1
    }, finally = {
      tt11.pvalue <- tt11.pvalue
    }
    )
    
    tryCatch({
      #----- modify 05/23/2024 -----#
      temp <- merge(dt.shap2[c("src_subject_id", i)], 
                    dt.grp.tt[c("src_subject_id", i, "rel_family_id", "site_id_l")] %>%
                      setNames(c("src_subject_id", "grp.index", "rel_family_id", "site_id_l")))
      formula <- paste0(i, " ~ grp.index + (1|rel_family_id)")
      mod <- lmer(formula, data=temp, REML = FALSE)
      tt21.pvalue <- summary(mod)$coefficients["grp.index", "Pr(>|t|)"]
      print(tt21.pvalue)
      #----- END -----#
      
    }, error = function(e) {
      # In case of error, set p-value to 1
      tt21.pvalue <- 1
    }, finally = {
      tt21.pvalue <- tt21.pvalue
    }
    )
    
    # due to algorithm limitation, allow a small number of participants with outlying values
    # change this on 03/30/2024
    if(min(length(tt11),length(tt21))>=100){
      perc<-1
    }else{
      perc<-1
    }
    
    # added on 03/31/2024 to zero out SHAP values with small deviation from zero
    # zero.cut was selected based on 1% of min(mean(int), mean(ext)) and min(mean(picvocab), mean(pic)): 0.04, 0.08
    # zero.cut was selected based on 1% of max(median(int), median(ext)) and max(median(picvocab), median(pic)): 0.03, 0.10
    if(type.outcome=="cbcl"){
      zero.cut<-0.02
    }else{
      zero.cut<-0.09
    }
    
    len11.eff<-length(tt11[which(abs(tt11)>zero.cut)])
    # need to make sure there are at least 0.5*len participants with the same sign
    # and the mean SHAP will be larger for these who engaged in the activity than these who don't
    tt11.same.sign<-max(sum(sign(tt11)==1), sum(sign(tt11)==-1))
    
    # new_scoring_05 used this abs(mean(tt11)-mean(tt12))>=abs(mean(tt12))
    if ((min(tt11)<min(tt12)) & (max(tt11)>max(tt12))) {
      shap.tr<-rbind(shap.tr, tt1.comb<-c(0,0))
    } else {
      if(tt11.pvalue < p.cut & tt11.same.sign >= max(ceiling(0.5*length(sel.id1)), floor(perc*len11.eff)) & abs(mean(tt11))>=max(abs(tt12.mean),zero.cut)){
        if((max(tt11)<=min(tt12)+zero.cut & any(tt12>= -zero.cut)) || (min(tt11)>=max(tt12)-zero.cut & any(tt12<= zero.cut))){
          tt1.comb<-c(mean(tt12),mean(tt11))
          shap.tr<-rbind(shap.tr, tt1.comb)
        }else{
          shap.tr<-rbind(shap.tr, tt1.comb<-c(0,0))
        }
      }else{
        shap.tr<-rbind(shap.tr, tt1.comb<-c(0,0))
      }
    }
    
    len21.eff<-length(tt21[which(abs(tt21)>zero.cut)])
    tt21.same.sign<-max(sum(sign(tt21)==1), sum(sign(tt21)==-1))
    
    # new_scoring_05 used this abs(mean(tt21)-mean(tt22))>=abs(mean(tt22))
    if ((min(tt21)<min(tt22)) & (max(tt21)>max(tt22))) {
      shap.tt<-rbind(shap.tt, tt2.comb<-c(0,0))
    } else {
      if(tt21.pvalue < p.cut & tt21.same.sign >= max(ceiling(0.5*length(sel.id21)), floor(perc*len21.eff)) & abs(mean(tt21))>=max(abs(tt22.mean),zero.cut)){
        # add a threshold to mark activity with no effect?? TBD 03/29/2024 --> hold on (marked 04/01/2024)
        if((max(tt21)<=min(tt22)+zero.cut & any(tt22>= -zero.cut)) || (min(tt21)>=max(tt22)-zero.cut & any(tt22<= zero.cut))){
          tt2.comb<-c(mean(tt22),mean(tt21))
          shap.tt<-rbind(shap.tt, tt2.comb)
        }else{
          shap.tt<-rbind(shap.tt, tt2.comb<-c(0,0))
        }
      }else{
        shap.tt<-rbind(shap.tt, tt2.comb<-c(0,0))
      }
    }
    
  }

  colnames(tt.count)<-c("tr1.count","tr0.count","tt1.count","tt0.count")
  colnames(tt.minmax)<-c("tr0.min","tr0.max","tr.min","tr.max","tt0.min","tt0.max","tt.min","tt.max")
  
  shap.mean<-data.frame(cbind(sai.nm, shap.tr, shap.tt))
  colnames(shap.mean)<-c("sai.nm","mean.tr0","mean.tr","mean.tt0","mean.tt")
  
  if(type.outcome=="cbcl"){
    shap.mean$tr.idx<-ifelse(shap.mean$mean.tr < 0, "pos", ifelse(shap.mean$mean.tr > 0, "neg","mixed"))
    shap.mean$tt.idx<-ifelse(shap.mean$mean.tt < 0, "pos", ifelse(shap.mean$mean.tt > 0, "neg","mixed"))
  }else{
    shap.mean$tr.idx<-ifelse(shap.mean$mean.tr > 0, "pos", ifelse(shap.mean$mean.tr < 0, "neg","mixed"))
    shap.mean$tt.idx<-ifelse(shap.mean$mean.tt > 0, "pos", ifelse(shap.mean$mean.tt < 0, "neg","mixed")) 
  }
  shap.mean$comb.idx<-ifelse(shap.mean$tt.idx != shap.mean$tr.idx, "mixed", shap.mean$tr.idx)
  shap.mean$sai.nm <- sai.nm
  
  write.csv(cbind(shap.mean,tt.minmax), file=paste0(out.dir, out.fn), row.names = F)
  
  return(list(count=data.frame(tt.count), shap.mean=data.frame(cbind(shap.mean,tt.minmax))))
}

out.dir<-"./Table/"

# main analysis
nsc.cbcl.total<-extract.shap.test(dt1=shap.dt.tr, dt2=shap.dt.tt, fn.shap1=f.lst[grep("cbcl_discovery",f.lst)], 
                           fn.shap2=f.lst[grep("cbcl_holdout",f.lst)], 
                           out.dir=out.dir, out.fn="index_cbcl.total.csv", type.outcome = "cbcl")

nsc.totalcomp<-extract.shap.test(dt1=shap.dt.tr, dt2=shap.dt.tt, fn.shap1=f.lst[grep("cog_discovery",f.lst)], 
                                fn.shap2=f.lst[grep("cog_holdout",f.lst)], 
                                out.dir=out.dir, out.fn="index_cog_total.csv", type.outcome = "cog")

table(nsc.cbcl.total$shap.mean$tr.idx,nsc.cbcl.total$shap.mean$tt.idx)
table(nsc.totalcomp$shap.mean$tr.idx,nsc.totalcomp$shap.mean$tt.idx)

nsc.shap.mean.all<-data.frame(cbind(nsc.cbcl.total$shap.mean[,c("sai.nm","tr.idx","tt.idx","comb.idx","mean.tr","mean.tt")],
  nsc.totalcomp$shap.mean[,c("tr.idx","tt.idx","comb.idx","mean.tr","mean.tt")]))
colnames(nsc.shap.mean.all)<-c("sai.nm",
                               paste0("cbcltotal.",c("tr.idx","tt.idx","comb.idx","tr","tt")),
                               paste0("totalcomp.",c("tr.idx","tt.idx","comb.idx","tr","tt")))     
      
write.csv(nsc.shap.mean.all, file= paste0(out.dir, "index_lifetime_total_r.csv"), row.names = F)

