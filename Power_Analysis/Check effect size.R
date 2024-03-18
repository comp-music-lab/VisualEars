library(tidyverse)
library(ggpubr)
library(emmeans)
library(ordinal)
library(kableExtra)
library(performance)
library(scales)
library(effectsize)

source("clmm-power-library.R", local = knitr::knit_global())

pval.lev <- function(pvals) {
  ifelse(pvals < 0.0001,
         "\\textbf{< 0.0001}",
         ifelse(pvals < 0.001,
                "\\textbf{< 0.001}",
                ifelse(pvals < 0.05,
                       paste0("\\textbf{", round(pvals, 4), "}"),
                       round(pvals, 2))))
}

summary.sig <- function(mod, custom_caption) {
  modTab <- data.frame(summary(mod)$coefficients) |>
    rownames_to_column() |>
    mutate_at("rowname", str_replace_all, ":", " × ") |>
    mutate_at("rowname", str_replace_all, "`", "") |>
    mutate("rowname" = str_replace_all(rowname,
                                       c("Tempo1" = "Tempo [Low]",
                                         "Solo.group1" = "Instrumentation [Group]",
                                         "Participant.country1" = "Participant country [Canada]",
                                         "Participant.country2" = "Participant country [Iran]",
                                         "Music.country1" = "Music country [Canada]",
                                         "Music.country2" = "Music country [Iran]"))) |> 
    mutate(Pr...z.. = pval.lev(Pr...z..)) |> 
    kable(digits = 2,
          booktabs = TRUE,
          align = c("l", rep("c", 4)),
          linesep = "",
          caption = custom_caption,
          col.names = c("Effects",
                        "Estimate", 
                        "Std. Error", 
                        "$z$", 
                        "$p$"),
          escape = FALSE) |> 
    kable_styling(latex_options = c("HOLD_position", "scale_down")) |>
    pack_rows(group_label = "Thresholds",
              start_row = 1,
              end_row = 4,
              hline_after = TRUE,
              bold = TRUE) |>
    pack_rows(group_label = "Terms",
              start_row = 5,
              end_row = 39,
              hline_before = TRUE,
              hline_after = TRUE,
              bold = TRUE)
  return(modTab)
}

contr.stars <- function(emms){
  require(emmeans)
  x <- tibble(as.data.frame(emms$contrasts))
  x <- separate(x,
                col = 1, 
                into = c("group1", "group2"), 
                sep = " - ", 
                remove = TRUE)
  x$p.signif <- ifelse(x$p.value < 0.0001, "****",
                       ifelse(x$p.value < 0.001, "***",
                              ifelse(x$p.value < 0.01, "**",
                                     ifelse(x$p.value < 0.05, "*", NA))))
  x <- x |>
    mutate_at("group1", str_replace_all, "[()]", "") |> 
    mutate_at("group2", str_replace_all, "[()]", "")
  return(x)
}

tempo.contr <- function(mod){
  emms.mod <- emmeans(mod, pairwise~Tempo)
  emms.mod.tab <- tibble(data.frame(emms.mod$emmeans))
  t.mod <- contr.stars(emms.mod) |>
    mutate(p.value = pval.lev(p.value))
  emm.contr.tab <- merge(emms.mod.tab, t.mod, by = 0, all = TRUE) |>
    select(-c(1,5,12,15)) |>
    unite(Contrast, group1, group2, sep = " - ") |>
    mutate_at("Contrast", str_replace_all, "NA - NA", " ") |>
    kable(digits = 2,
          booktabs = TRUE,
          align = c("l", rep("c", 5), "l", rep("c", 5)),
          linesep = "",
          caption = "Estimated marginal means and contrasts between Tempo conditions",
          col.names = c("Tempo",
                        "EMM",
                        "$SE$",
                        "$2.5\\% CI$",
                        "$97.5\\% CI$",
                        "Contrast",
                        "Difference",
                        "$SE$",
                        "$z$",
                        "$p$"),
          escape = FALSE) |>
    add_header_above(c(" " = 5, "Contrasts" = 5)) |>
    kable_styling(latex_options = c("HOLD_position")) |>
    footnote(general = "EMM = estimated marginal mean.",
             threeparttable = TRUE,
             footnote_as_chunk = TRUE,
             escape = FALSE)
  return(emm.contr.tab)
}


# Load data
sh1 <- read.csv("../sh1.csv") |> 
  mutate(Tempo = str_replace(Tempo, "low", "Low")) |> 
  mutate(Tempo = str_replace(Tempo, "high", "High")) |> 
  mutate_if(is.character,as.factor) |> 
  mutate(Tempo = fct_relevel(Tempo, c("Low", "High")))

# Assign unique id to pairs
uniqueval <- unique(sh1[, c("Participant", "Music.country", "Solo.group")])
sh1$groupid <- 0
for (i in 1:dim(uniqueval)[1]) {
  idx <- sh1$Participant == uniqueval$Participant[i] & 
    sh1$Music.country == uniqueval$Music.country[i] & 
    sh1$Solo.group == uniqueval$Solo.group[i]
  sh1$groupid[idx] <- i
}

# Select only relevant variables
data <- sh1  |> 
  select(1:9,20)

# Plot distribution of results for density ratings
p1a <- ggplot(data, aes(x = density.tempo, y = after_stat(density),
                        fill = Tempo, color = Tempo)) +
  scale_fill_hue(direction = -1) + scale_colour_hue(direction = -1) +
  geom_histogram(alpha = 0.3, position = "identity", binwidth = 1) +
  labs(y= "Probability", x = "Rating", title = "Density") +
  geom_text(aes(label = format(after_stat(density), digits = 1), y= after_stat(density)), 
            stat= "bin", binwidth = 1, 
            vjust = -0.2,
            show.legend = FALSE) +
  facet_wrap(~Solo.group)

# Plot distribution of results for arousal ratings
p1b <- ggplot(data, aes(x = arousal.tempo, y = after_stat(density),
                        fill = Tempo, color = Tempo)) +
  scale_fill_hue(direction = -1) + scale_colour_hue(direction = -1) +
  geom_histogram(alpha = 0.3, position = "identity", binwidth = 1) +
  labs(y= NULL, x = "Rating", title = "Arousal") +
  geom_text(aes(label = format(after_stat(density), digits = 1), y= after_stat(density)), 
            stat= "bin", binwidth = 1, 
            vjust = -0.2,
            show.legend = FALSE) +
  facet_wrap(~Solo.group)

# Arrange plots
p1 <- ggarrange(p1a, p1b,
                common.legend = TRUE,
                legend = "bottom",
                labels = "AUTO")
p1

# Density ratings for group music
dat.den.gr.PILOT <- clmm_generate_data(n_participants = 100,
                                       n_trials = 3,
                                       control_distribution = c(.04, .11, .59, .22, .04),
                                       effect = 1,
                                       participant_variation = 1,
                                       within_subject = TRUE,
                                       control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Density.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 200)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 100), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 200)) |>
  mutate(Solo.group = rep("Group", TIMES = 600)) |>
  select(c(1,6,5,2,7,3))

# Density ratings for solo music
dat.den.so.PILOT <- clmm_generate_data(n_participants = 100,
                                       n_trials = 3,
                                       control_distribution = c(.55, .30, .04, .10, .01),
                                       effect = 1,
                                       participant_variation = 1,
                                       within_subject = TRUE,
                                       control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Density.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 200)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 100), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 200)) |>
  mutate(Solo.group = rep("Solo", TIMES = 600)) |>
  select(c(1,6,5,2,7,3))

# Arousal ratings for group music
dat.aro.gr.PILOT <- clmm_generate_data(n_participants = 100,
                                       n_trials = 3,
                                       control_distribution = c(.11, .11, .52, .19, .07),
                                       effect = 1,
                                       participant_variation = 1,
                                       within_subject = TRUE,
                                       control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Arousal.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 200)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 100), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 200)) |>
  mutate(Solo.group = rep("Group", TIMES = 600)) |>
  select(c(1,6,5,2,7,3))

# Arousal ratings for solo music
dat.aro.so.PILOT <- clmm_generate_data(n_participants = 100,
                                       n_trials = 3,
                                       control_distribution = c(.43, .41, .14, .01, .01),
                                       effect = 1,
                                       participant_variation = 1,
                                       within_subject = TRUE,
                                       control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Arousal.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 200)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 100), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 200)) |>
  mutate(Solo.group = rep("Solo", TIMES = 600)) |>
  select(c(1,6,5,2,7,3))

dat.sim.PILOT <- left_join(bind_rows(dat.den.gr.PILOT, dat.den.so.PILOT),
                           bind_rows(dat.aro.gr.PILOT, dat.aro.so.PILOT),
                           relationship = "many-to-many") |> 
  mutate(Tempo = fct_relevel(Tempo, c("Low", "High"))) |>
  mutate_if(is.character,as.factor) |> 
  mutate(Density.tempo = as.factor(Density.tempo)) |>
  mutate(Arousal.tempo = as.factor(Arousal.tempo))

glimpse(dat.sim.PILOT)

options(contrasts = c("contr.sum","contr.poly"))
model.den <- clmm(Density.tempo ~ Tempo * Solo.group +
                    (1 + Tempo | Participant),
                  data = dat.sim.PILOT)

#summary.sig(model.den, "Visual density ratings model summary")

coefs.den <- data.frame(summary(model.den)$coefficients)
confs.den <- data.frame(confint(model.den))

emms <- emmeans(model.den, pairwise ~ Tempo | Solo.group)
data.frame(emms) |> 
  mutate(pred = ifelse(emmean <= coefs.den$Estimate[1], 1,
                       ifelse(emmean <= coefs.den$Estimate[2], 2,
                              ifelse(emmean <= coefs.den$Estimate[3], 3,
                                            ifelse(emmean <= coefs.den$Estimate[4], 4, 5)))))

emmip(model.den, ~ Tempo | Solo.group, CIs = TRUE) +
  geom_hline(yintercept = coefs.den$Estimate[1:4])

emmip(model.den, ~ Tempo , CIs = TRUE) +
  geom_hline(yintercept = coefs.den$Estimate[1:4])

p2a <- ggplot(data, aes(x = density.tempo, y = after_stat(density),
                        fill = Tempo, color = Tempo)) +
  scale_fill_hue(direction = -1) + scale_colour_hue(direction = -1) +
  geom_histogram(alpha = 0.3, position = "identity", binwidth = 1) +
  labs(y= "Probability", x = "Rating", title = "Pilot data") +
  geom_text(aes(label = format(after_stat(density), digits = 1), y= after_stat(density)), 
            stat= "bin", binwidth = 1, 
            vjust = -0.2,
            show.legend = FALSE) +
  facet_wrap(~Solo.group)

p2b <- ggplot(dat.sim.PILOT, aes(x = as.integer(Density.tempo), y = after_stat(density),
                        fill = Tempo, color = Tempo)) +
  scale_fill_hue(direction = -1) + scale_colour_hue(direction = -1) +
  geom_histogram(alpha = 0.3, position = "identity", binwidth = 1) +
  labs(y= "Probability", x = "Rating", title = "Simulated data") +
  geom_text(aes(label = format(after_stat(density), digits = 1), y= after_stat(density)), 
            stat= "bin", binwidth = 1, 
            vjust = -0.2,
            show.legend = FALSE) +
  facet_wrap(~Solo.group)

ggarrange(p2a, p2b)


emmip(model.den, ~ Tempo | Solo.group, CIs = TRUE) +
  geom_hline(yintercept = coefs.den$Estimate[1:4],
             color = "darkred") +
  annotate('rect', xmin=-Inf, xmax=Inf, 
           ymin = confs.den$X2.5..[1:4], 
           ymax = confs.den$X97.5..[1:4], 
           alpha=.2, fill='darkred') +
  scale_y_continuous("Linear prediction", 
    sec.axis = sec_axis(~ .,
                        breaks = coefs.den$Estimate[1:4],
                        labels = row.names(coefs.den)[1:4], 
                        name = "Visual density")) +
  theme(axis.title.y.right = element_text(color="darkred"),
        axis.text.y.right  = element_text(color="darkred"),
        axis.ticks.y.right = element_line(color="darkred"))


emmip(model.den, ~ Tempo | Solo.group, CIs = TRUE) +
  geom_hline(yintercept = coefs.den$Estimate[1:4],
             color = "darkred") +
  annotate('rect', xmin=-Inf, xmax=Inf, 
           ymin = confs.den$X2.5..[1:4], 
           ymax = confs.den$X97.5..[1:4], 
           alpha = 0.1, fill='darkred') +
  scale_y_continuous("Linear prediction (latent, continuous variable)", 
                     sec.axis = sec_axis(~ .,
                                         breaks = c(-2.5, 
                                                    mean(coefs.den$Estimate[1:2]), 
                                                    mean(coefs.den$Estimate[2:3]), 
                                                    mean(coefs.den$Estimate[3:4]), 
                                                    3.8),
                                         labels = 1:5, 
                                         name = "Visual density (original Likert scale)")) +
  theme(axis.title.y.right = element_text(color="darkred"),
        axis.text.y.right  = element_text(color="darkred"),
        axis.ticks.y.right = element_line(color="darkred"))
  


