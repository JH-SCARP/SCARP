#1-day risk analysis for progression to severe disease or death in patients hospitalized with COVID-19  

library(tidyverse)
library(magrittr)
library(dplyr)
library(pROC)
library(rfSLAM)
library(splines)
library(grid)
library(gridExtra)


### 1-day predictions
# rf.df.1 is the dataframe with data in the 6-hour interval format
node.size <- round((dim(rf.df.1)[1])*0.1)
set.seed(321)
start_time <- Sys.time()
mymodel.1.full <- rfSLAM::rfsrc(i.sev_died.1 ~., data = rf.df.1, nodesize = node.size, ntree = 100, na.action = "na.impute", splitrule="poisson.split1", risk.time = rf.cpiu.1$rt_1p, stratifier=rf.df.1$q6, membership = TRUE, bootstrap = "by.user", samp = rfslam.samp.1, var.used = "by.tree")
end_time <- Sys.time()

# obtain predicted event rates
p.cpiu.be.ppl <- risk.adjust.be(rf = mymodel.1.full, status = rf.df.1$i.sev_died.1, rt = rf.cpiu.1$rt_1p, alpha.tm = 0)


# variable importance
vars.tree <- mymodel.1.full$var.used

for(m in 1:100){
  for(n in 1:105){
    vars.tree[m,n] <- ifelse(vars.tree[m,n]>1,1,vars.tree[m,n])
  }
}

# number of trees using each variable
var.tree.sumC <- colSums(vars.tree) 
# number of unique variables used in each tree
var.tree.sumR <- rowSums(vars.tree)


var.used <- data.frame(Variable = names(var.tree.sumC), perc_tree = unname(var.tree.sumC)) 
var.used <- left_join(var.used,var_key, by = "Variable") 

var.used <- var.used %>% filter(perc_tree >=5) 
var.used$Type <- ifelse(var.used$Type %in% c("CBC", "Other lab", "Inflammatory marker", "Lab ratio", "BMP"), "Lab", var.used$Type)

library(ggsci) 

p1 <- ggplot(data = var.used, aes(x=reorder(Key,perc_tree), y = perc_tree, fill = Type)) + geom_bar(stat="identity", color = "black") + theme_bw() + coord_flip() 

p1 + scale_fill_npg() + scale_y_continuous(limits = c(0,70), expand = c(0,0), breaks = seq(0, 70, by = 10)) + xlab("") + 
  
  ylab("Percent of Trees Using Variable") + 
  
  theme_minimal() + 
  
  theme(plot.title = element_text(hjust = 0.5,face = "bold"), axis.text=element_text(size=6)) + 
  
  ggtitle("Variable Importance") 

### Calibration

## 5-fold cross-validation for performance metric

set.seed(321)
rf.cpiu.1$p.hat <- 1-exp(-p.cpiu.be.ppl) # predicted event probabilities 
rf.cpiu.1$ni.sca <- as.numeric(rf.cpiu.1$i.sev_died.1)-1
shuffled <- rf.cpiu.1[sample(nrow(rf.cpiu.1)),] 

# create 5 equally sized folds 

folds <- cut(seq(1,nrow(shuffled)), breaks = 5, labels = FALSE) 

shuffled$fold <- folds 

shuffled$risk <- NA 

for(i in 1:5){ 
  
  testing <- shuffled[shuffled$fold == i,] 
  
  training <- shuffled[shuffled$fold != i,] 
  
  lr_model = glm(as.factor(ni.sca) ~ ns(p.hat,2)*ns(q6,2),data = training, family = binomial)
  
  p.hat.df <- data.frame(p.hat = testing$p.hat, q6 = testing$q6)
  
  lr_probs = predict(lr_model,  
                     newdata = p.hat.df, 
                     type = "response")
  
  shuffled[shuffled$fold == i,"risk"] <- lr_probs
  
} 


## Performance

# AUC
db2 <- shuffled

auc.df <- data.frame(matrix(NA, nrow = 2, ncol = 4))

names(auc.df) <- c("time", "lci", "auc", "uci")

auc.df$time <- c(1, 2)

# week 1
t0 <- db2 %>% filter(q6 <= 27)
auc.df[1,2:4] <- unname(ci.auc(t0$sev_died_within_1_days, t0$risk))


# week 2
t0 <- db2 %>% filter(q6 > 27)
auc.df[2,2:4] <- unname(ci.auc(t0$sev_died_within_1_days, t0$risk))


## calibration with smooth curve

g1 <- mutate(db2, bin = ntile(risk, 10)) %>% 
  # Bin prediction into 10ths
  group_by(bin) %>%
  mutate(n = n(), # Get ests and CIs
         bin_pred = mean(risk), 
         bin_prob = mean(sev_died_within_1_days), 
         se = sqrt((bin_prob * (1 - bin_prob)) / n), 
         ul = bin_prob + 1.96 * se, 
         ll = bin_prob - 1.96 * se) %>%
  ungroup() %>%
  ggplot(aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
  geom_pointrange(size = 0.5, color = "black") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  geom_abline() + # 45 degree line indicating perfect calibration
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", 
              color = "black", formula = y~-1 + x) + 
  geom_smooth(aes(x = risk, y = sev_died_within_1_days), 
              color = "gray", se = FALSE, method = "loess") + 
  # loess fit through estimates
  xlab("") +
  ylab("Observed Probability") + theme_classic() + ggtitle("1-Day Risk Predictions") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=12))


g2 <- ggplot(db2, aes(x = risk)) +
  geom_histogram(fill = "black", bins = 200) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  xlab("Predicted Probability") +
  ylab("") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 8000),breaks = seq(0, 8000, by = 2000)) +
  theme(panel.grid.minor = element_blank(), text = element_text(size=15))

g <- arrangeGrob(g1, g2, respect = TRUE, heights = c(1, 0.25), ncol = 1)
grid.newpage()
grid.draw(g)
