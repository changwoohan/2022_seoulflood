
# create function for gsynth analysis 
# with 95 CI

han_gsynth_95 = function(dat.use, val.name, date, fol, fig_name) {
  
  fig.list=NULL
  in.list=NULL
  out.list=NULL
  pe <- c(0, 14, 0, 21, 0, 28)
  date2 = paste0("figures/", date)
  path <-here::here(date2)
  dir.create(path)
  basic <- COUNT_A_1_d0 ~ int
  
  nint <- dat.use %>% 
    filter(int==1) %>% 
    distinct(REGION_EMD) %>% nrow()
  
  nam <- dat.use %>% 
    filter(int==1) %>% 
    distinct(REGION_EMD) 
  
  for (i in (1:length(val.name))){ 

    val.name[i]
    fig_tit <- d_name[d_name$xy ==  val.name[i],]$z
    tab_tit <- d_name[d_name$xy ==  val.name[i],]$y
    uf  <- update.formula(basic, as.formula(paste(val.name[i], "~ .")))
    out <- gsynth(uf, data = dat.use, index = c("REGION_EMD","time"), EM= F,
                  CV = TRUE, r = c(0, 5), force = "two-way",
                  nboots = 1000, inference = "parametric", se = TRUE, parallel = TRUE, seed=1234)
    
    bef <- dat.use %>% ungroup() %>%
      filter(REGION_EMD %in% nam$REGION_EMD & eq==0) %>%
      filter(time>-15 & time<0) %>% 
      select(val.name[i]) 
    
    bef_sum<- round(colSums(bef)/nint,1)
    
    # figure generation
    
    b1 <-as.data.frame(out$Y.co)
    b2 <-apply(b1,1,quantile,c(0.05,0.95))
    b3 <-as.data.frame(t((b2)))
    
    b4 <-b3 %>% mutate(time=as.numeric(row.names(b3))+1)
    colnames(b4) <-c("lower", "upper", "time")
    
    k1 <-as.data.frame(out$est.att) %>% mutate(time=as.numeric(row.names(out$est.att)))
    k2 <-as.data.frame(out$Y.bar)  %>% mutate(time=as.numeric(row.names(out$Y.bar))+1)
    
    
  k3 <-k2 %>% left_join(k1, by=c("time")) %>%
    left_join(b4, by=c("time")) %>%
    mutate(Exposed=Y.ct.bar+ATT,
           Controls=Y.ct.bar,
           exp_up=Y.ct.bar+CI.upper,
           exp_low=Y.ct.bar+CI.lower)
  
  k4 <- k3 %>% 
    pivot_longer(
      cols =c("Exposed", "Controls"),
      names_to="value",
      values_to="rate") %>%  
    mutate(value=as.factor(value))
  
  
  k5<-k4 %>% 
    filter(time>-40 & time<40) %>% 
    select(rate)
  
  k4$value <- relevel(k4$value,"Exposed")
    
    fig.list[[i]] <- ggplot(k4, aes_string(x="time", y="rate", group=as.factor(k4$value))) + 
      geom_line(aes(lty = as.factor(value)), size=1.0) +
      labs(title=fig_tit,
           x="Days from the flood", 
           y="Number of hospital visits") +
      coord_cartesian(ylim=c(0, max(k5$rate)*2), xlim=c(-40,40)) +
      geom_vline(xintercept =0, linetype="longdash", size=1.0) +
      geom_vline(xintercept =28, linetype="longdash", size=1.0) +
      scale_x_continuous(breaks = c(-20, 0, 20)) +
      geom_ribbon(aes(x = time, ymin = exp_low, ymax = exp_up),
                  fill = "grey", alpha = 0.4) +
      theme(
        plot.title = element_text(hjust=0.5, size=20, face="bold"),  
        panel.background = element_rect(fill = "white"),         # Set plot background to white
        legend.key  = element_rect(fill = "white"),              # Set legend item backgrounds to white
        axis.line.x = element_line(colour = "black", size = 1),  # Add line to x axis
        axis.line.y = element_line(colour = "black", size = 1),   # Add line to y axis
        legend.position = c(.95, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_blank(),
        legend.text=element_text(size=20),
        legend.key.height = unit(0.4,"cm"),
        legend.key.width = unit(3,"cm"),
        axis.title=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15)) 
    
    
    tiff(file=paste0(path,"/",fig_name,"_",val.name[i],".tiff"), 
         res=300, width=12, height=8, units="in",  pointsize=60, compression = 'lzw')
    
    plot(fig.list[[i]])
    
    dev.off()
    
    
    for (k in (1:3)) {
      
      cumu <- cumuEff(out, cumu = TRUE, id = NULL, period = c(pe[2*k -1],pe[2*k]))
      cumu <- as.data.frame(cumu$est.catt)  
      cumu <- cumu %>% mutate (dd=out$Y,
                               CATT=round(CATT,1),
                               low=round(CI.lower,1),
                               high=round(CI.upper,1), 
                               NO=row.names(cumu)) %>% 
        select(c("dd","CATT","low","high","p.value","NO")) %>% 
        filter(row_number()==n())
      
      cumu_1<- cbind(cumu, bef_sum) %>% 
        mutate(pe   = round(CATT/bef_sum *100,1),
               l_pe = round(low/bef_sum *100,1),
               h_pe = round(high/bef_sum *100,1),
               DIS  = tab_tit) %>% 
        mutate(pre_event = bef_sum) %>% 
        mutate(diff_out = paste0(CATT," (",low,", ",high,")")) %>% 
        mutate(diff_per = paste0(pe," (",l_pe,", ",h_pe,")")) %>% 
        mutate(DIS2 = val.name[i]) %>% 
        mutate(pvalue =p.value) %>% 
        select(DIS2, pvalue, DIS, pre_event, diff_out, diff_per)  
      
      in.list[[k]]  <- cumu_1  
      out.list[[i]] <- as.data.frame(do.call(cbind, in.list))
      
    } 
    
  }
  
  final_out   <- as.data.frame(do.call(rbind, out.list)) 


  date3 = paste0("output/",date)
  path <-here::here(date3)
  dir.create(path)
  write.csv(final_out, paste0(path,"/", fol,"flood.csv"))
  
}





han_gsynth_ALL_fig = function(dat.use, folder, fig_name) {
  
  
  date2 = paste0("figures/", folder)
  path <-here::here(date2)
  dir.create(path)
  
  fig.list=NULL
  in.list=NULL
  out.list=NULL
  
  x1 <-c("ALL_A_T","ALL_E_T",
         "ALL_I_T","ALL_J_T",
         "ALL_O_T","ALL_S_T")
  
  title <-c("Certain infectious and parasitic diseases\n (A00-B99)",
            "Endocrine, nutritional and metabolic diseases\n (E00-E90)",
            "Diseases of the circulatory system\n (I00-I99)",
            "Diseases of the respiratory system\n (J00-J99)",
            "Pregnancy, childbirth and the puerperium\n (O00-O99)",
            "Injury, poisoning and certain other consequences\n of external causes (S00-T98)")
  
  ymax <-c(200,200,200,250,5,200)
  
  basic <- RATE_A_INC_d0 ~ int
  
  
  for (i in (1:length(x1))){ 
    x1[i]
    uf <-update.formula(basic, as.formula(paste(x1[i], "~ .")))
    out <- gsynth(uf, data = dat.use, index = c("REGION_EMD","time"), EM= F,
                  CV = TRUE, r = c(0, 5), force = "two-way",
                  nboots = 1000, inference = "parametric", se = TRUE, parallel = TRUE, seed=1234)
    
    b1 <-as.data.frame(out$Y.co)
    b2 <-apply(b1,1,quantile,c(0.05,0.95))
    b3 <-as.data.frame(t((b2)))
    
    b4 <-b3 %>% mutate(time=as.numeric(row.names(b3))+1)
    colnames(b4) <-c("lower", "upper", "time")
    
    k1 <-as.data.frame(out$est.att) %>% mutate(time=as.numeric(row.names(out$est.att)))
    k2 <-as.data.frame(out$Y.bar)  %>% mutate(time=as.numeric(row.names(out$Y.bar))+1)
    
    k3 <-k2 %>% left_join(k1, by=c("time")) %>%
      left_join(b4, by=c("time")) %>%
      mutate(int=Y.ct.bar+ATT) %>% 
      select("time", "Y.tr.bar", "Y.ct.bar","upper","lower") %>% 
      mutate(Flooded=Y.tr.bar,
             Controls=Y.ct.bar)
    
    k4 <- k3 %>% 
      pivot_longer(
        cols =c("Flooded", "Controls"),
        names_to="value",
        values_to="rate") %>%  
      mutate(value=as.factor(value))
    
    k4$value <- relevel(k4$value,"Flooded")
    
    fig.list[[i]] <- ggplot(k4, aes_string(x="time", y="rate", group=as.factor(k4$value))) + 
      geom_line(aes(lty = as.factor(value)), size=1.0) +
      labs(title=title[i],
           x="Days from the flood", 
           y="Daily number of hospital visits") +
      xlim(-30,30) +
      coord_cartesian(ylim=c(0, ymax[i])) +
      geom_vline(xintercept =0, linetype="longdash", size=1.0) +
      geom_vline(xintercept =13, linetype="longdash", size=1.0) +
      
      theme(
        plot.title = element_text(hjust=0.5, size=35, face="bold"),  
        panel.background = element_rect(fill = "white"),         # Set plot background to white
        legend.key  = element_rect(fill = "white"),              # Set legend item backgrounds to white
        axis.line.x = element_line(colour = "black", size = 1),  # Add line to x axis
        axis.line.y = element_line(colour = "black", size = 1),   # Add line to y axis
        legend.position = c(.95, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_blank(),
        legend.text=element_text(size=25),
        legend.key.height = unit(0.4,"cm"),
        legend.key.width = unit(3,"cm"),
        axis.title=element_text(size=25),
        axis.title.x=element_text(size=25),
        axis.title.y=element_text(size=25),
        axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25)) 
    plot(fig.list[[i]])
    
  }
  
  tiff(file=paste0(path,"/", fig_name,".tiff"), 
       res=300, width=35, height=30, units="in",  pointsize=60, compression = 'lzw')
  
  grid.arrange(fig.list[[1]], fig.list[[2]],
               fig.list[[3]], fig.list[[4]],
               fig.list[[5]], fig.list[[6]],
               nrow=3, ncol=2)
  dev.off()
  
  
  
}






han_gsynth_SPE_fig = function(dat.use, folder, fig_name) {
  
  date2 = paste0("figures/", folder)
  path <-here::here(date2)
  dir.create(path)
  
  fig.list=NULL
  in.list=NULL
  out.list=NULL
  
  x1 <-c("B_4_T","J_1_T","J_5_T",
         "S_6_T","S_7_T","S_10_T")
  
  title <-c("Other viral diseases\n (B25-B34)",
            "Acute upper respiratory infections\n (J00-J06)",
            "Chronic lower respiratory diseases\n (J40-J47)",
            "Injuries to the elbow and forearm\n (S50-S59)",
            "Injuries to the wrist and hand\n (S70-799)",
            "Injuries to the ankle and foot\n (S90-S99)")
  
  ymax <-c(15,120,25,20,30,30)
  
  basic <- RATE_A_INC_d0 ~ int
  
  
  for (i in (1:length(x1))){ 
    x1[i]
    uf <-update.formula(basic, as.formula(paste(x1[i], "~ .")))
    out <- gsynth(uf, data = dat.use, index = c("REGION_EMD","time"), EM= F,
                  CV = TRUE, r = c(0, 5), force = "two-way",
                  nboots = 1000, inference = "parametric", se = TRUE, parallel = TRUE, seed=1234)
    
    b1 <-as.data.frame(out$Y.co)
    b2 <-apply(b1,1,quantile,c(0.05,0.95))
    b3 <-as.data.frame(t((b2)))
    
    b4 <-b3 %>% mutate(time=as.numeric(row.names(b3))+1)
    colnames(b4) <-c("lower", "upper", "time")
    
    k1 <-as.data.frame(out$est.att) %>% mutate(time=as.numeric(row.names(out$est.att)))
    k2 <-as.data.frame(out$Y.bar)  %>% mutate(time=as.numeric(row.names(out$Y.bar))+1)
    
    k3 <-k2 %>% left_join(k1, by=c("time")) %>%
      left_join(b4, by=c("time")) %>%
      mutate(int=Y.ct.bar+ATT) %>% 
      select("time", "Y.tr.bar", "Y.ct.bar","upper","lower") %>% 
      mutate(Flooded=Y.tr.bar,
             Controls=Y.ct.bar)
    
    k4 <- k3 %>% 
      pivot_longer(
        cols =c("Flooded", "Controls"),
        names_to="value",
        values_to="rate") %>%  
      mutate(value=as.factor(value))
    
    k4$value <- relevel(k4$value,"Flooded")
    
    fig.list[[i]] <- ggplot(k4, aes_string(x="time", y="rate", group=as.factor(k4$value))) + 
      geom_line(aes(lty = as.factor(value)), size=1.0) +
      labs(title=title[i],
           x="Days from the flood", 
           y="Daily number of hospital visits") +
      xlim(-30,30) +
      coord_cartesian(ylim=c(0, ymax[i])) +
      geom_vline(xintercept =0, linetype="longdash", size=1.0) +
      geom_vline(xintercept =13, linetype="longdash", size=1.0) +
      
      theme(
        plot.title = element_text(hjust=0.5, size=35, face="bold"),  
        panel.background = element_rect(fill = "white"),         # Set plot background to white
        legend.key  = element_rect(fill = "white"),              # Set legend item backgrounds to white
        axis.line.x = element_line(colour = "black", size = 1),  # Add line to x axis
        axis.line.y = element_line(colour = "black", size = 1),   # Add line to y axis
        legend.position = c(.95, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_blank(),
        legend.text=element_text(size=25),
        legend.key.height = unit(0.4,"cm"),
        legend.key.width = unit(3,"cm"),
        axis.title=element_text(size=25),
        axis.title.x=element_text(size=25),
        axis.title.y=element_text(size=25),
        axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25)) 
    plot(fig.list[[i]])
    
  }
  
  
  tiff(file=paste0(path,"/", fig_name,".tiff"), 
       res=300, width=35, height=30, units="in",  pointsize=60, compression = 'lzw')
  
  grid.arrange(fig.list[[1]], fig.list[[2]],
               fig.list[[3]], fig.list[[4]],
               fig.list[[5]], fig.list[[6]],nrow=3, ncol=2)
  dev.off()
  
}