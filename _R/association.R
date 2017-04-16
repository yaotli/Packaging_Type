# aim: to test the association of substitution and subtype
#

library(seqinr)
library(stringr)

  # read-in

gsgd_fasta <- read.fasta("~/Desktop/RNA_GsGD.fasta")

  seq_name0 = attributes(gsgd_fasta)$names
       seq0 = getSequence(gsgd_fasta)
      
  # eliminate unsubtyped virus     
      
  nonsubtype <- which(is.na(str_match(pattern = "_(H[0-9]+N[0-9]+)_", seq_name0)[,1]) == TRUE)

  seq_name <- seq_name0[-nonsubtype]
       seq <- seq0[-nonsubtype]

  # input data        
       
  seq_matrix <- do.call(rbind, seq)

  # extract subtype data  
 
     n1_id <- which( str_match(pattern = "_(H[0-9]+N[0-9]+)_", seq_name)[,2] == "H5N1" )
  nonN1_id <- seq(1:length(seq_name))[-h5_id]
 
  # out_n1
  out_Nx <- rep(0, length(seq_name))
  out_Nx[nonN1_id] = 1
 
  # out_nt
  out_sub <- rep(0, length(seq_name))
 
# HA gene-wise searching

  # 2x2 Chi-Q     
  
HA_Chi_p <- 
  apply(seq_matrix, 2, function(x){
    
    # remove "-"
    
    out_sub <- out_sub[-which(x == "-")]
     out_Nx <- out_Nx[-which(x == "-")]
          x <- x[-which(x == "-")]
    
    # dist of all nt    
    
          nt = 
            as.data.frame(
              table(
                x[which(out_Nx == 1)] ), stringsAsFactors = FALSE)$Var1[
                  
                  # consensus of N1      
                  which.max(
                    as.data.frame(
                      table(
                        x[which(out_Nx == 1)] )
                    )$Freq
                  ) 
                  ]
   
    out_sub[which(x == nt)] = 1
          
    if (
      
      length(which(out_sub == 0) ) == 0
      
    ){ return(1) }else{
    
    totest <- data.frame(Nx = out_Nx, Sub = out_sub)
      
    p <- -log10((chisq.test(totest$Nx, totest$Sub, correct = TRUE)$p.value))
    
return(p)
    }
    
  })
 
HA_Chi_value <- 
  apply(seq_matrix, 2, function(x){
    
    # remove "-"
    
    out_sub <- out_sub[-which(x == "-")]
    out_Nx <- out_Nx[-which(x == "-")]
    x <- x[-which(x == "-")]
    
    # dist of all nt    
    
    nt = 
      as.data.frame(
        table(
          x[which(out_Nx == 1)] ), stringsAsFactors = FALSE)$Var1[
            
            # consensus of N1      
            which.max(
              as.data.frame(
                table(
                  x[which(out_Nx == 1)] )
              )$Freq
            ) 
            ]
    
    out_sub[which(x == nt)] = 1
    
    if (
      
      length(which(out_sub == 0) ) == 0
      
    ){ return(1) }else{
      
      totest <- data.frame(Nx = out_Nx, Sub = out_sub)
      
      p <- chisq.test(totest$Nx, totest$Sub, correct = TRUE)$statistic
      
      return(p)
    }
    
  })

df_HA_Chi_value <- data.frame(pos = seq(1, length(HA_Chi_value)), ChiQ = HA_Chi_value)


  # OR

HA_LR_p <- 
  apply(seq_matrix, 2, function(x){
    
    # remove "-"
    
    out_sub <- out_sub[-which(x == "-")]
    out_Nx <- out_Nx[-which(x == "-")]
    x <- x[-which(x == "-")]
    
    # dist of all nt    
    
    nt = 
      as.data.frame(
        table(
          x[which(out_Nx == 1)] ), stringsAsFactors = FALSE)$Var1[
            
            # consensus of N1      
            which.max(
              as.data.frame(
                table(
                  x[which(out_Nx == 1)] )
              )$Freq
            ) 
            ]
    
    out_sub[which(x == nt)] = 1
    
    if (
      
      length(which(out_sub == 0) ) == 0
      
    ){ return(0) }else{
      
      totest <- data.frame(Nx = out_Nx, Sub = out_sub)
      
      LR <- glm(Nx ~ Sub, data = totest, family = binomial(link="logit") )
      
       p = -log10(summary(LR)$coefficients[2,4])
      
      return(p)
    }
    
  })

HA_LR_OR <- 
  apply(seq_matrix, 2, function(x){
    
    # remove "-"
    
    out_sub <- out_sub[-which(x == "-")]
    out_Nx <- out_Nx[-which(x == "-")]
    x <- x[-which(x == "-")]
    
    # dist of all nt    
    
    nt = 
      as.data.frame(
        table(
          x[which(out_Nx == 1)] ), stringsAsFactors = FALSE)$Var1[
            
            # consensus of N1      
            which.max(
              as.data.frame(
                table(
                  x[which(out_Nx == 1)] )
              )$Freq
            ) 
            ]
    
    out_sub[which(x == nt)] = 1
    
    if (
      
      length(which(out_sub == 0) ) == 0
      
    ){ return(0) }else{
      
      totest <- data.frame(Nx = out_Nx, Sub = out_sub)
      
      LR <- glm(Nx ~ Sub, data = totest, family = binomial(link="logit") )
      
      p = exp(summary(LR)$coefficients[2,1])
      
      return(p)
    }
    
  })


df_HA_LR <- data.frame(pos = seq(1, length(HA_LR_p)), OR = HA_LR_OR, pvalue = HA_LR_p)




# figure 

  library(ggplot2)

# fig1
  p0 = 
  ggplot(df_HA_Chi_value, aes(x = pos)) + 
    geom_line(aes(y = ChiQ)) + 
    theme_bw() + 
    geom_text(x = 1651, y = 2500, label = "1651 (108)", size = 2) + 
    geom_text(x = 1724, y = 2300, label = "1724 (35)", size = 2) + 
    
    geom_text(x = 38, y = 5500, label = "38(3')", color = "RED", size = 2) + 
    geom_text(x = 1634, y = 5500, label = "124(5')", color = "RED", size = 2) + 
    
    geom_vline(xintercept = 38, color = "RED", linetype="dashed") + 
    geom_vline(xintercept = 1634, color = "RED", linetype="dashed") + 
    scale_x_continuous(breaks = seq(0, 1700, by = 250)) + 
    xlab('')
  
  
# fig2 
  p1 =
  ggplot(df_HA_LR, aes(x = pos)) + 
    geom_line(aes(y = OR)) + 
    theme_bw() + 
    geom_vline(xintercept = 1651,  linetype="dashed") + 
    geom_vline(xintercept = 1724,  linetype="dashed") +
    
    geom_text(x = 1651, y = 5000000, label = "69", size = 2) + 
    geom_text(x = 1724, y = 5500000, label = "1844", size = 2) + 
    
    geom_vline(xintercept = 38, color = "RED", linetype="dashed") + 
    geom_vline(xintercept = 1634, color = "RED", linetype="dashed") + 
    
    scale_x_continuous(breaks = seq(0, 1700, by = 250)) + 
    xlab("") 
  
  p2 = 
  ggplot(df_HA_LR, aes(x = pos)) + 
    geom_line(aes(y = pvalue)) + 
    theme_bw() + 
    geom_vline(xintercept = 1651,  linetype="dashed") + 
    geom_vline(xintercept = 1724,  linetype="dashed") + 
    
    geom_vline(xintercept = 38, color = "RED", linetype="dashed") + 
    geom_vline(xintercept = 1634, color = "RED", linetype="dashed") + 
    
    scale_x_continuous(breaks = seq(0, 1700, by = 250)) + 
    ylab("p-value (-log10)")
  
  
multiplot(p0, p1, p2, ncol=1)

  
  
  

 
 
 
 
  
