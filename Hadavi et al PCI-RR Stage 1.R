library(ggplot2)
library(readxl)
library(effectsize)
library(pwr)

sh1 <- read_excel("sh1.xlsx")

#draw a plot for arousal based on tempo changes
ggplot(sh1, aes(x = Tempo, y = arousal.tempo, color=Participant.country, shape=Music.country, size=Solo.group))   +    geom_jitter(width=0.3, height =0.15) + labs(y= "Arousal", x = "Tempo")

#estimate overal effect size:
cohens_d(arousal.tempo ~ Tempo, data = sh1, paired=TRUE) #Cohen's d = 0.80

#estimate effect sizes for each participant group:
cohens_d(arousal.tempo ~ Tempo, data = subset(sh1,Participant.country=="Japan"), paired=TRUE) #Japan: Cohen's d = 0.82
cohens_d(arousal.tempo ~ Tempo, data = subset(sh1,Participant.country=="Canada"), paired=TRUE) #Canada: Cohen's d = 1.06
cohens_d(arousal.tempo ~ Tempo, data = subset(sh1,Participant.country=="Iran"), paired=TRUE) #Iran: Cohen's d = 0.58


#same for tempo-density:
ggplot(sh1, aes(x = Tempo, y = density.tempo, color=Participant.country, shape=Music.country, size=Solo.group))   +    geom_jitter(width=0.3, height =0.15) + labs(y= "Density", x = "Tempo")
cohens_d(density.tempo ~ Tempo, data = sh1, paired=TRUE)

#estimate overal effect size:
cohens_d(density.tempo ~ Tempo, data = sh1, paired=TRUE) #Cohen's d = 0.83

#estimate effect sizes for each participant group:
cohens_d(density.tempo ~ Tempo, data = subset(sh1,Participant.country=="Japan"), paired=TRUE) #Japan: Cohen's d = 0.95
cohens_d(density.tempo ~ Tempo, data = subset(sh1,Participant.country=="Canada"), paired=TRUE) #Canada: Cohen's d = 1.08
cohens_d(density.tempo ~ Tempo, data = subset(sh1,Participant.country=="Iran"), paired=TRUE) #Iran: Cohen's d = 0.59


#power analysis
pwr.t.test(d=0.4,power=0.95,sig.level=0.05/2,type="paired",alternative="greater") #n = 83.16425 pairs