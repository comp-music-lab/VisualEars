library(ggplot2)
library(ggpubr)
library(effectsize)
library(pwr)

#load data
sh1 <- read.csv("sh1.csv")

#assign unique id to pairs
uniqueval <- unique(sh1[, c("Participant", "Music.country", "Solo.group")])
sh1$groupid <- 0
for (i in 1:dim(uniqueval)[1]) {
  idx <- sh1$Participant == uniqueval$Participant[i] & sh1$Music.country == uniqueval$Music.country[i] & sh1$Solo.group == uniqueval$Solo.group[i]
  sh1$groupid[idx] <- i
}

#plot - density
pd <- position_dodge(0.3)
g_list <- vector(mode = "list", length = 2)

sh1_solo <- sh1[sh1$Solo.group == "Solo", ]
g <- ggplot(sh1_solo, aes(x = Tempo, y = density.tempo, color=Participant.country, shape=Music.country)) +
  geom_violin(aes(group = Tempo), draw_quantiles = 0.5) +
  geom_point(position = pd, size = 3) +
  geom_line(aes(group = groupid), size = 0.5, linetype = "dashed", position = pd, alpha = 0.5) +
  labs(y= "Density", x = "Tempo", title = "Solo") +
  scale_x_discrete(limits = c("high", "low")) + 
  scale_size_discrete(breaks = c("Solo", "Group")) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
g_list[[1]] <- g

sh1_group <- sh1[sh1$Solo.group == "Group", ]
g <- ggplot(sh1_group, aes(x = Tempo, y = density.tempo, color=Participant.country, shape=Music.country)) +
  geom_violin(aes(group = Tempo), draw_quantiles = 0.5) +
  geom_point(position = pd, size = 3) +
  geom_line(aes(group = groupid), size = 0.5, linetype = "dashed", position = pd, alpha = 0.5) +
  labs(y = "", x = "Tempo", title = "Group") +
  scale_x_discrete(limits = c("high", "low")) + 
  scale_size_discrete(breaks = c("Solo", "Group")) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
g_list[[2]] <- g

g <- ggarrange(plotlist = g_list, ncol = 2, nrow = 1, common.legend = TRUE)

plot(g)

#estimate overal effect size:
cohens_d(density.tempo ~ Tempo, data = sh1, paired=TRUE) #Cohen's d = 0.83

#estimate effect sizes for each participant group:
cohens_d(density.tempo ~ Tempo, data = subset(sh1,Participant.country=="Japan"), paired=TRUE) #Japan: Cohen's d = 0.95
cohens_d(density.tempo ~ Tempo, data = subset(sh1,Participant.country=="Canada"), paired=TRUE) #Canada: Cohen's d = 1.08
cohens_d(density.tempo ~ Tempo, data = subset(sh1,Participant.country=="Iran"), paired=TRUE) #Iran: Cohen's d = 0.59


#plot - arousal
pd <- position_dodge(0.3)
g_list <- vector(mode = "list", length = 2)

sh1_solo <- sh1[sh1$Solo.group == "Solo", ]
g <- ggplot(sh1_solo, aes(x = Tempo, y = arousal.tempo, color=Participant.country, shape=Music.country)) +
  geom_violin(aes(group = Tempo), draw_quantiles = 0.5) +
  geom_point(position = pd, size = 3) +
  geom_line(aes(group = groupid), size = 0.5, linetype = "dashed", position = pd, alpha = 0.5) +
  labs(y= "Arousal", x = "Tempo", title = "Solo") +
  scale_x_discrete(limits = c("high", "low")) + 
  scale_size_discrete(breaks = c("Solo", "Group")) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
  ylim(c(1, 5))
g_list[[1]] <- g

sh1_group <- sh1[sh1$Solo.group == "Group", ]
g <- ggplot(sh1_group, aes(x = Tempo, y = arousal.tempo, color=Participant.country, shape=Music.country)) +
  geom_violin(aes(group = Tempo), draw_quantiles = 0.5) +
  geom_point(position = pd, size = 3) +
  geom_line(aes(group = groupid), size = 0.5, linetype = "dashed", position = pd, alpha = 0.5) +
  labs(y = "", x = "Tempo", title = "Group") +
  scale_x_discrete(limits = c("high", "low")) + 
  scale_size_discrete(breaks = c("Solo", "Group")) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
  ylim(c(1, 5))
g_list[[2]] <- g

g <- ggarrange(plotlist = g_list, ncol = 2, nrow = 1, common.legend = TRUE)

plot(g)

#estimate overal effect size:
cohens_d(arousal.tempo ~ Tempo, data = sh1, paired=TRUE) #Cohen's d = 0.80

#estimate effect sizes for each participant group:
cohens_d(arousal.tempo ~ Tempo, data = subset(sh1,Participant.country=="Japan"), paired=TRUE) #Japan: Cohen's d = 0.82
cohens_d(arousal.tempo ~ Tempo, data = subset(sh1,Participant.country=="Canada"), paired=TRUE) #Canada: Cohen's d = 1.06
cohens_d(arousal.tempo ~ Tempo, data = subset(sh1,Participant.country=="Iran"), paired=TRUE) #Iran: Cohen's d = 0.58


#power analysis
pwr.t.test(d=0.4,power=0.95,sig.level=0.05/2,type="paired",alternative="greater") #n = 83.16425 pairs