---
title: 'Cross-cultural relationships between music, emotion, and visual imagery: A comparative study of Iran, Canada, and Japan [Stage 1 Registered Report]'
subtitle: 'Code and analyses: Simulation-based power analysis'
author:
  - name: Shafagh Hadavi \orcidlink{0009-0008-1184-7238}
    correspondence: false
    institute: keiko
  - name: Junji Kuroda
    correspondence: false
    institute: yamaha
  - name: Taiki Shimozono
    correspondence: false
    institute: yamaha
  - name: Juan David Leongómez \orcidlink{0000-0002-0092-6298}
    correspondence: false
    institute: ueb
  - name: Patrick E. Savage \orcidlink{0000-0001-6996-7496}
    correspondence: false
    institute: keiko
institute:
  - keiko: Graduate School of Media and Policy / Faculty of Environment and Information Studies, Keio University, Fujisawa, Japan.
  - yamaha: YAMAHA Corporation, Hamamatsu, Japan.
  - ueb: Human Behaviour and Evolution Lab, Faculty of Psychology, Universidad El Bosque, Bogota, Colombia.
date: "18 March, 2024"
output:
  bookdown::pdf_document2:
    citation_package: biblatex
    highlight: zenburn
    number_sections: yes
    keep_tex:  true
    toc: no
    pandoc_args:
      - '--lua-filter=lua/scholarly-metadata.lua'
      - '--lua-filter=lua/author-info-blocks.lua'
classoption: 
      - bookmarksnumbered
editor_options:
  chunk_output_type: console
geometry: margin=2cm
header-includes: 
  \usepackage{caption} 
  \usepackage{float} 
  \floatplacement{figure}{H} 
  \usepackage[utf8]{inputenc} 
  \usepackage{fancyhdr}
  \pagestyle{fancy} 
  \usepackage{hanging}
  \usepackage{csquotes}
  \usepackage[style=apa,backend=biber]{biblatex}
  \lhead{Hadavi et al.} 
  \rhead{\textit{Power analysis}} 
  \renewcommand{\abstractname}{Description} 
  \usepackage{multicol}
  \usepackage{orcidlink}
  \newcommand{\opensupplement}{\setcounter{table}{0}
    \renewcommand{\thetable}{S\arabic{table}} \setcounter{figure}{0}
    \renewcommand{\thefigure}{S\arabic{figure}}}
  \newcommand{\closesupplement}{\setcounter{table}{0}
    \renewcommand{\thetable}{\arabic{table}} \setcounter{figure}{0}
    \renewcommand{\thefigure}{\arabic{figure}}}
  \usepackage{multirow,booktabs,setspace}
  \DeclareCaptionLabelSeparator{point}{. }
  \DeclareCaptionLabelSeparator{point}{. }
  \captionsetup[table]{labelfont=bf,
    textfont=it,
    format=plain,
    labelsep=point,
    skip=5pt}
  \captionsetup[figure]{labelfont=bf,
    format=plain,
    justification=justified,
    singlelinecheck=false,
    labelsep=point,
    skip=5pt}
always_allow_html: yes
csl: https://raw.githubusercontent.com/citation-style-language/styles/master/apa.csl
bibliography: Bib/Bibliography.bib
urlcolor: blue
citecolor: gray
linkcolot: gray
link-citations: true
---

```{=tex}
\begin{center}
Correspondence to:
SH: \href{mailto:shafagh@keio.jp}{shafagh@keio.jp}. 
PES: \href{mailto:psavage@sfc.keio.ac.jp}{psavage@sfc.keio.ac.jp}.
```

------------------------------------------------------------------------

```{=tex}
\textbf{Description}
\end{center}

\par
\begingroup
\leftskip3em
\rightskip\leftskip
```

This document contains all code, and step by step explanations for the simulation-based power analysis (including supplementary figures and tables) for:

```{=latex}
\begin{hangparas}{.25in}{1}
Hadavi, S., Kuroda, J., Shimozono, T., Leongómez, J. D., \& Savage, P. E. (in prep). \textit{Cross-cultural relationships between music, emotion, and visual imagery: A comparative study of Iran, Canada, and Japan} [Stage 1 Registered Report]. \url{https://doi.org/10.31234/osf.io/26yg5}
\end{hangparas}
```

Data available from the *VisualEars* repository in GitHub: <https://github.com/comp-music-lab/VisualEars>. This document and its underlying code were created in `R Markdown` by Juan David Leongómez using \LaTeX.

------------------------------------------------------------------------

```{=latex}
\par
\endgroup

{\hypersetup{hidelinks}
\setcounter{tocdepth}{6}
\tableofcontents
}
\opensupplement
```



------------------------------------------------------------------------
 
# Power analysis strategy

Our power analysis utilizes Monte Carlo simulations to account for the characteristics of our study data and the statistical model we are using. We employ Cumulative Link Mixed Models (CLMM) due to the nature of our dependent variables, which are based on 5-point Likert scales rather than normally distributed continuous data. Additionally, the paired responses from each participant are not independent of one another.

While analytical power analysis can be conducted for certain statistical tests involving a single independent variable, such as correlation, simple regression, t-tests, and one-way ANOVA, more complex models and factorial designs require a simulation-based approach [e.g., @lakensSimulationBasedPowerAnalysis2021].

In this approach, a statistical model is defined with a specific effect size of interest. Multiple synthetic datasets are generated based on this model, and the proportion of datasets resulting in significant results (i.e., rejecting the null hypothesis with $p < \alpha$) is computed ($\alpha = 0.025$ after applying Bonferroni correction). This proportion represents the Monte Carlo estimate of the test's power to detect the desired effect size.

To perform our analysis, we simulated a population of 10,000 participants, each providing 6 paired responses, resulting in a total of 120,000 ratings for both visual density and arousal. From this simulated population, we randomly selected 1,000 samples. The sample size for each was adjusted until we achieved the desired statistical power. For each sample, we fitted models to assess both arousal and visual density ratings, incorporating fixed effects for Tempo (Low, High), Instrumentation (Group, Solo), Participant country (Iran, Canada, Japan), and music country (Iran, Canada, Japan), along with their interactions. Random intercepts and slopes between Tempo conditions were included for each participant to ensure a fair representation of the experimental design [@barrRandomEffectsStructure2013; @debruineUnderstandingMixedEffectsModels2021].

However, since our hypotheses focus solely on the effects of manipulating tempo, the power analysis primarily examines the planned contrast between Tempo conditions, independent of instrumentation, participant country, and music country levels.

# Preliminaries

## Load packages

This file was created using `knitr` [@knitrcit], mostly using `tidyverse` [@tidyversecit] syntax. As such, data wrangling was mainly done using packages such as `dplyr` [@dplyrcit], and most figures were created using `ggplot2` [@ggplotcit]. Tables were created using `knitr::kable` and `kableExtra` [@kableExtracit].

Cumulative Link Mixed Models (CLMM) were fitted using `ordinal` [@ordinalcit], and contrasts and interactions were explored using `emmeans` [@emmeanscit]. In all cases, threshold between levels of the dependent (ordinal) variable, were set as flexible.

All packages used in this file can be directly installed from the Comprehensive R Archive Network ([CRAN](https://cran.r-project.org/)). For a complete list of packages used to create this file, and their versions, see section \@ref(session), at the end of the document.


```r
library(tidyverse)
library(ggpubr)
library(emmeans)
library(ordinal)
library(kableExtra)
library(performance)
library(scales)
```

## Custom functions

### `clmm_generate_data`

To simulate ordinal data with a specific distribution, we used the function `clmm_generate_data`, created by @bordersPowerAnalysisOrdinal2022 as part of their tutorial (https://osf.io/e6usd/). To run this function, a number of other custom functions need to be defined and implemented as well. The original code can be found in the `clmm-power-library.R` file (https://osf.io/tjpkf). Here, it is run from source, but in a slightly modified version of the file in which we changed the `mc.cores` parameter to 1, to avoid errors.


```r
source("clmm-power-library.R", local = knitr::knit_global())
```

### `pval.lev`

The function `pval.lev` takes p-values and formats them in \LaTeX, highlighting significant results in bold.


```r
pval.lev <- function(pvals) {
  ifelse(pvals < 0.0001,
         "\\textbf{< 0.0001}",
         ifelse(pvals < 0.001,
                "\\textbf{< 0.001}",
                ifelse(pvals < 0.05,
                       paste0("\\textbf{", round(pvals, 4), "}"),
                       round(pvals, 2))))
}
```

### `summary.sig`

The function `summary.sig` was created to bold significant *p* values from summary model tables. It highlights significant $p$ values, and formats the output in \LaTeX, ready to be used with `kable`.


```r
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
```

### `contr.stars`

Function to create a data frame of model contrasts, representing significance levels from an `emmeans::emmeans` output. These data frames are formatted to be called by the `ggpubr::stat_pvalue_manual` function used in model figures.


```r
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
```

### `tempo.contr`

This functions creates a table combining both estimated marginal means and contrasts between levels of Tempo, and formats the output in \LaTeX, ready to be used with `kable`. 


```r
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
```

# Data

## Pilot data

To simulate data, first we looked at the distribution of the pilot data. 


```r
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
```

![(\#fig:unnamed-chunk-7)Distribution of results from pilot data for ratings of both Density (**A**) and Arousal (**B**), by Instrumentation (Group, Solo), and Tempo (High, Low).](Power_analysis_TELESCOPE_files/figure-latex/unnamed-chunk-7-1.pdf) 

## Data simulation

While basing a power analysis on pilot data is problematic and can bias the sample size estimation [see e.g., @albersWhenPowerAnalyses2018], we used this distribution as a starting point.

The `clmm_generate_data` function [@bordersPowerAnalysisOrdinal2022] allows to simulate data with random and fixed effects with an ordinal outcome. In this case, we simulated data for Tempo, the within-subject variable for which there was an *a-priori* prediction, replicating the distribution of the first level (Low) using a specified probability of each rating score (using the argument `control_distribution`). The second level of that within-subject variable (Tempo = High) is generated based on a hypothesized effect (using the argument `effect`). Other important factors such as participant country and music country were randomly assigned, as there are no *a-priori* predictions regarding these factors.

### Replication of the pilot data distribution

To obtain distributions similar to the probabilities observed in the pilot data, we first simulated individual data for Density and Arousal ratings, separated by Instrumentation type (Group, Solo). In all cases, we simulated data from 10000 participants, and changed the value assigned to the `effect` argument until the distribution of the second level of the within-subject variable (Tempo = High) resembled the pilot data.


```r
# Density ratings for group music
dat.den.gr.PILOT <- clmm_generate_data(n_participants = 10000,
                           n_trials = 3,
                           control_distribution = c(.04, .11, .59, .22, .04),
                           effect = 3.2,
                           participant_variation = 1,
                           within_subject = TRUE,
                           control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Density.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 20000)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 10000), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 20000)) |>
  mutate(Solo.group = rep("Group", TIMES = 60000)) |>
  select(c(1,6,5,2,7,3))

# Density ratings for solo music
dat.den.so.PILOT <- clmm_generate_data(n_participants = 10000,
                                 n_trials = 3,
                                 control_distribution = c(.55, .30, .04, .10, .01),
                                 effect = 1.2,
                                 participant_variation = 1,
                                 within_subject = TRUE,
                                 control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Density.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 20000)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 10000), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 20000)) |>
  mutate(Solo.group = rep("Solo", TIMES = 60000)) |>
  select(c(1,6,5,2,7,3))

# Arousal ratings for group music
dat.aro.gr.PILOT <- clmm_generate_data(n_participants = 10000,
                                 n_trials = 3,
                                 control_distribution = c(.11, .11, .52, .19, .07),
                                 effect = 3.1,
                                 participant_variation = 1,
                                 within_subject = TRUE,
                                 control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Arousal.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 20000)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 10000), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 20000)) |>
  mutate(Solo.group = rep("Group", TIMES = 60000)) |>
  select(c(1,6,5,2,7,3))

# Arousal ratings for solo music
dat.aro.so.PILOT <- clmm_generate_data(n_participants = 10000,
                                 n_trials = 3,
                                 control_distribution = c(.43, .41, .14, .01, .01),
                                 effect = 1.6,
                                 participant_variation = 1,
                                 within_subject = TRUE,
                                 control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Arousal.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 20000)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 10000), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 20000)) |>
  mutate(Solo.group = rep("Solo", TIMES = 60000)) |>
  select(c(1,6,5,2,7,3))
```

We then merged these four data frames to obtain a simulated data frame.


```r
dat.sim.PILOT <- left_join(bind_rows(dat.den.gr.PILOT, dat.den.so.PILOT),
                           bind_rows(dat.aro.gr.PILOT, dat.aro.so.PILOT),
                           relationship = "many-to-many") |> 
  mutate(Tempo = fct_relevel(Tempo, c("Low", "High")))
```

The distribution of this replication of the pilot data distribution is represented in \@ref(fig:pilot-rep)


```r
# Plot distribution of results for density ratings
p2a <- ggplot(dat.sim.PILOT, aes(x = Density.tempo, y = after_stat(density),
                                 fill = Tempo, color = Tempo)) +
  scale_fill_hue(direction = -1) + scale_colour_hue(direction = -1) +
  geom_histogram(alpha = 0.3, position = "identity", binwidth = 1) +
  labs(y= "Probability", x = "Rating") +
  geom_text(aes(label = format(after_stat(density), digits = 1), y= after_stat(density)), 
            stat= "bin", binwidth = 1, 
            vjust = -0.2,
            show.legend = FALSE) +
  facet_wrap(~Solo.group)

# Plot distribution of results for arousal ratings
p2b <- ggplot(dat.sim.PILOT, aes(x = Arousal.tempo, y = after_stat(density),
                                 fill = Tempo, color = Tempo)) +
  scale_fill_hue(direction = -1) + scale_colour_hue(direction = -1) +
  geom_histogram(alpha = 0.3, position = "identity", binwidth = 1) +
  labs(y= NULL, x = "Rating") +
  geom_text(aes(label = format(after_stat(density), digits = 1), y= after_stat(density)), 
            stat= "bin", binwidth = 1, 
            vjust = -0.2,
            show.legend = FALSE) +
  facet_wrap(~Solo.group)

# Arrange plots
p2 <- ggarrange(p2a +
                  labs(subtitle = "Density"), 
                p2b  +
                  labs(subtitle = "Arousal"),
                common.legend = TRUE,
                legend = "bottom",
                labels = "AUTO")
p2
```

![(\#fig:pilot-rep)Distribution of results from a simulated replication of the pilot data distribution for ratings of both Density (**A**) and Arousal (**B**), by Instrumentation (Group, Solo), and Tempo (High, Low).](Power_analysis_TELESCOPE_files/figure-latex/pilot-rep-1.pdf) 

### Simulation of data with a more conservative effect

While the data simulated in the previous section resembled the distribution of the pilot data, using these data would be problematic for at least two reasons:

First, the effects of Tempo were extremely large (3.2, 1.2, 3.1, 1.6), particularly for ratings of Group instrumentation (3.2, and 3.1). With such differences, a very high statistical power of about $1 - \beta = 0.95$  can be achieved with 2 or 3 participants. And second, as mentioned before, estimating effects from pilot data is problematic, as it biases the results [@albersWhenPowerAnalyses2018]. 

For this reason, we decided to maintain the distributions of probabilities for the first Tempo level (Low), as it is a good starting point to estimate the distribution of ratings that participants would assign, but simulated the second level (Tempo = High) with a much more conservative effect of 1 in all cases. 


```r
# Density ratings for group music
dat.den.gr <- clmm_generate_data(n_participants = 10000,
                                 n_trials = 3,
                                 control_distribution = c(.04, .11, .59, .22, .04),
                                 effect = 0.4,
                                 participant_variation = 1,
                                 within_subject = TRUE,
                                 control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Density.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 20000)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 10000), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 20000)) |>
  mutate(Solo.group = rep("Group", TIMES = 60000)) |>
  select(c(1,6,5,2,7,3))

# Density ratings for solo music
dat.den.so <- clmm_generate_data(n_participants = 10000,
                                 n_trials = 3,
                                 control_distribution = c(.55, .30, .04, .10, .01),
                                 effect = 0.4,
                                 participant_variation = 1,
                                 within_subject = TRUE,
                                 control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Density.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 20000)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 10000), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 20000)) |>
  mutate(Solo.group = rep("Solo", TIMES = 60000)) |>
  select(c(1,6,5,2,7,3))

# Arousal ratings for group music
dat.aro.gr <- clmm_generate_data(n_participants = 10000,
                                 n_trials = 3,
                                 control_distribution = c(.11, .11, .52, .19, .07),
                                 effect = 0.4,
                                 participant_variation = 1,
                                 within_subject = TRUE,
                                 control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Arousal.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 20000)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 10000), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 20000)) |>
  mutate(Solo.group = rep("Group", TIMES = 60000)) |>
  select(c(1,6,5,2,7,3))

# Arousal ratings for solo music
dat.aro.so <- clmm_generate_data(n_participants = 10000,
                                 n_trials = 3,
                                 control_distribution = c(.43, .41, .14, .01, .01),
                                 effect = 0.4,
                                 participant_variation = 1,
                                 within_subject = TRUE,
                                 control_weight = .5) |>
  mutate(group = ifelse(group == 0, "Low", "High")) |>
  dplyr::rename(Arousal.tempo = pas) |>
  dplyr::rename(Tempo = group) |>
  dplyr::rename(Participant = id) |>
  mutate(Participant = paste0("p", Participant)) |>
  select(1:3) |>
  mutate(pair = rep(1:3, times = 20000)) |>
  mutate(Participant.country = rep(rep(c("Iran", "Canada", "Japan"), each = 10000), 2)) |>
  mutate(Music.country = rep(c("Iran", "Canada", "Japan"), times = 20000)) |>
  mutate(Solo.group = rep("Solo", TIMES = 60000)) |>
  select(c(1,6,5,2,7,3))
```

Again, we then merged these four data frames to obtain a final simulated data frame.


```r
dat.sim <- left_join(bind_rows(dat.den.gr, dat.den.so),
                     bind_rows(dat.aro.gr, dat.aro.so),
                     relationship = "many-to-many") |> 
  mutate(Tempo = fct_relevel(Tempo, c("Low", "High")))

# Assign unique id to pairs
uniqueval <- unique(dat.sim[, c("Participant", "Music.country", "Solo.group")])
dat.sim$groupid <- 0
for (i in 1:dim(uniqueval)[1]) {
  idx <- dat.sim$Participant == uniqueval$Participant[i] & 
    dat.sim$Music.country == uniqueval$Music.country[i] & 
    dat.sim$Solo.group == uniqueval$Solo.group[i]
  dat.sim$groupid[idx] <- i
}
```

The distribution of these simulated data is represented in \@ref(fig:final-sim)


```r
# Plot distribution of results for density ratings
p3a <- ggplot(dat.sim, aes(x = Density.tempo, y = after_stat(density),
                           fill = Tempo, color = Tempo)) +
  scale_fill_hue(direction = -1) + scale_colour_hue(direction = -1) +
  geom_histogram(alpha = 0.3, position = "identity", binwidth = 1) +
  labs(y= "Probability", x = "Rating") +
  geom_text(aes(label = format(after_stat(density), digits = 1), y= after_stat(density)), 
            stat= "bin", binwidth = 1, 
            vjust = -0.2,
            show.legend = FALSE) +
  facet_wrap(~Solo.group)

# Plot distribution of results for arousal ratings
p3b <- ggplot(dat.sim, aes(x = Arousal.tempo, y = after_stat(density),
                           fill = Tempo, color = Tempo)) +
  scale_fill_hue(direction = -1) + scale_colour_hue(direction = -1) +
  geom_histogram(alpha = 0.3, position = "identity", binwidth = 1) +
  labs(y= NULL, x = "Rating") +
  geom_text(aes(label = format(after_stat(density), digits = 1), y= after_stat(density)), 
            stat= "bin", binwidth = 1, 
            vjust = -0.2,
            show.legend = FALSE) +
  facet_wrap(~Solo.group)

# Arrange plots
p3 <- ggarrange(p3a +
                  labs(subtitle = "Density"), 
                p3b  +
                  labs(subtitle = "Arousal"),
                common.legend = TRUE,
                legend = "bottom",
                labels = "AUTO")
p3
```

![(\#fig:final-sim)Distribution of results from a conservative simulated distribution for ratings of both Density (**A**) and Arousal (**B**), by Instrumentation (Group, Solo), and Tempo (High, Low).](Power_analysis_TELESCOPE_files/figure-latex/final-sim-1.pdf) 

## Comparisons of the distribuition of results between pilot and simulated data


```r
p4 <- ggarrange(ggarrange(p1a + labs(title = "Pilot data"), 
                          p1b + labs(title = " "),
                          legend = "none"), 
                ggarrange(p2a + labs(title = "Simulation (pilot data replication)"), 
                          p2b + labs(title = " "),
                          legend = "none"),
                ggarrange(p3a + labs(title = "Final simulation (conservative effect)"), 
                          p3b + labs(title = " "),
                          common.legend = TRUE, legend = "bottom"),
                common.legend = TRUE,
                labels = "AUTO",
                nrow = 3,
                heights = c(1,1,1.1))
p4
```

![(\#fig:unnamed-chunk-12)Distribution of results of ratings of both Density and Arousal (**B**), by Instrumentation (Group, Solo), and Tempo (High, Low). **A** pilot data; **B** simulation replicating pilot data distribution; **C** simulation with a more conservative effect.](Power_analysis_TELESCOPE_files/figure-latex/unnamed-chunk-12-1.pdf) 

# Power analysis {#power-section}

Using the final, simulated population (*N* = 10,000 participants, with 6 paired observations per participant), we conducted a simulation-based power analysis. To do this, we first defined the desired number of simulations, the sample size (per country) of each simulation, and the $\alpha$ value (significance level) for all statistical tests:


```r
# Number of simulations to run
num_sims.clmm = 100 

# Sample size for each simulation (number of participants per country)
sample_size.clmm = 80 

# Significance level
alpha.clmm = 0.05
```

Then, models were fitted from the defined 100 random samples extracted from the simulated population. 

From each sample, two Cumulative Link Mixed Models (CLMM) were fitted [see e.g., @taylorRatingNormsShould2022]: one for Density ratings, and one for arousal ratings. All models were fitted with the following call, simply changing the dependent variable (DV): `clmm(DV ~ Tempo * Solo.group * Participant.country * Music.country + (1 + Tempo | Participant)`. 

Thus, all models had the same structure, which included the main effects and all possible interactions between Tempo (Low, High), Instrumentation (Group, Solo), Participant country (Iran, Canada, Japan), and Music country (Iran, Canada, Japan) as fixed effects, as well as random intercepts and random slopes between Tempo conditions for each participant. 

For each random sample, we adjusted the number of participants per country (in this case, *n* = 80) to achieve the target statistical power of at least $1 - \beta = 0.95$ for both the Density and Arousal models, for the contrast between Tempo levels (i.e., averaged over the levels of instrumentation, participant country and music country).

## Select random samples


```r
# Get participant IDs
part <- dat.sim  |> 
  group_by(Participant.country) |>
  distinct(Participant) |>
  ungroup()

# Select samples of random participants 
samples.clmm <- map_dfr(seq_len(num_sims.clmm), ~part |>
                        sample_n(sample_size.clmm) |>
                        mutate(sample.clmm = as.factor(.x)))

# Final data base of data for each participant selected for each sample
samples.clmm.long <- left_join(samples.clmm, dat.sim, by = c("Participant", 
                                                             "Participant.country"),
                                relationship = "many-to-many") |>
  mutate(Density.tempo = as.factor(Density.tempo)) |>
  mutate(Arousal.tempo = as.factor(Arousal.tempo)) |>
  mutate(Tempo = as.factor(Tempo))
```

## Visual density power analysis

Please be aware that this part of the simulation can take many ours (or even days) to run.


```r
# Create empty data frame
clmm.comps.den <- data.frame(Sample = 1:num_sims.clmm)

# Loop to fit a model for each random sample and extract the p-value of the Tempo contrasts
for (i in 1:num_sims.clmm){
  samp <- samples.clmm.long |>
    filter(sample.clmm == i)
  tryCatch(mod <- clmm(Density.tempo ~ 
                         Tempo * Solo.group * Participant.country * Music.country + 
                         (1 + Tempo | Participant), data = samp),
           error = function(e) {})
  tryCatch(comp <- emmeans(mod, pairwise~Tempo),
           error = function(e) {})
  clmm.comps.den$p[i] <- round(data.frame(comp$contrast)$p.value, 3)
  clmm.comps.den$`Statistical signicance`[i] <- ifelse(clmm.comps.den$p[i] <= alpha.clmm, 
                                                       "Significant", "Non-significant")
}

# Calculate power (proportion of simulations in which the Tempo contrast was significant)
power.clmm.den <- sum(clmm.comps.den$`Statistical signicance` == "Significant")/num_sims.clmm
```

## Arousal power analysis

Again, please keep in mind that this part of the simulation can take many ours (or even days) to run.


```r
# Create empty data frame
clmm.comps.aro <- data.frame(Sample = 1:num_sims.clmm)

# Loop to fit a model for each random sample and extract the p-value of the Tempo contrast
for (i in 1:num_sims.clmm){
  samp <- samples.clmm.long |>
    filter(sample.clmm == i)
  tryCatch(mod <- clmm(Arousal.tempo ~ 
                         Tempo * Solo.group * Participant.country * Music.country + 
                         (1 + Tempo | Participant), data = samp),
           error = function(e) {})
  tryCatch(comp <- emmeans(mod, pairwise~Tempo),
           error = function(e) {})
  clmm.comps.aro$p[i] <- round(data.frame(comp$contrast)$p.value, 3)
  clmm.comps.aro$`Statistical signicance`[i] <- ifelse(clmm.comps.aro$p[i] <= alpha.clmm, 
                                                       "Significant", "Non-significant")
}

# Calculate power (proportion of simulations in which the Tempo contrast was significant)
power.clmm.aro <- sum(clmm.comps.aro$`Statistical signicance` == "Significant")/num_sims.clmm
```

## Achieved power plot


```r
p5a <- ggplot(clmm.comps.den, aes(x = p, fill = `Statistical signicance`)) +
  scale_fill_hue(direction = -1) +
  geom_histogram(bins = 1/alpha.clmm, breaks = seq(0, 1, alpha.clmm)) +
  labs(y = "Count", x = "p-value", title = "Density") +
  annotate("text", x = 0.5, 
           y = sum(clmm.comps.den$`Statistical signicance` == "Significant") * 0.9,
           label = paste0("Power = ", round(power.clmm.den, 3),
                          "\nSample size = ", sample_size.clmm,
                          " participants per country\n(6 paired responses per participant)"))

p5b <- ggplot(clmm.comps.aro, aes(x = p, fill = `Statistical signicance`)) +
  scale_fill_hue(direction = -1) +
  geom_histogram(bins = 1/alpha.clmm, breaks = seq(0, 1, alpha.clmm)) +
  labs(y = "Count", x = "p-value", title = "Arousal") +
  annotate("text", x = 0.5, 
           y = sum(clmm.comps.aro$`Statistical signicance` == "Significant") * 0.9, 
           label = paste0("Power = ", round(power.clmm.aro, 3), 
                          "\nSample size = ", sample_size.clmm, 
                          " participants per country\n(6 paired responses per participant)"))

p5 <- ggarrange(p5a, p5b,
                common.legend = TRUE,
                legend = "bottom",
                labels = "AUTO")
p5
```

![(\#fig:unnamed-chunk-15)Histogram of *p*-values for the contrasts between Tempo conditions (Low, High) for all simulations. The estatistical power is the proportion of significant results. Significant results are in red.](Power_analysis_TELESCOPE_files/figure-latex/unnamed-chunk-15-1.pdf) 

# Analysis script with example data  {#example-data}

Here, using a randomly selected sample, we fit one model as if it was the data collected from the experiment. This is not only to make sure that the data are sufficient and the structure is adequate to fit the model, but to plan in advance all the code.


```r
set.seed(1)
# Select a random sample
ex.data <- samples.clmm.long  |> 
  filter(sample.clmm == sample(1:num_sims.clmm, 1, replace = TRUE)) |> 
  mutate_if(is.character,as.factor)
```

Because only the effect of Tempo was predicted and simulated *a-priori*, all other effects and interactions are completely random.


```r
p6a <- ggplot(ex.data, aes(x = Tempo, y = Density.tempo, 
                           color = Participant.country, shape = Music.country)) +
  geom_violin(aes(group = Tempo), draw_quantiles = 0.5) +
  geom_point(position = position_dodge(0.3), size = 3) +
  geom_line(aes(group = factor(groupid)), 
            linewidth = 0.5, linetype = "dashed", position = position_dodge(0.3), alpha = 0.5) +
  labs(y= "Density", x = "Tempo", title = "Density") +
  scale_x_discrete(limits = c("Low", "High")) + 
  scale_size_discrete(breaks = c("Solo", "Group")) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
  facet_wrap(~Solo.group)

p6b <- ggplot(ex.data, aes(x = Tempo, y = Arousal.tempo, 
                           color = Participant.country, shape = Music.country)) +
   geom_violin(aes(group = Tempo), draw_quantiles = 0.5) +
  geom_point(position = position_dodge(0.3), size = 3) +
  geom_line(aes(group = factor(groupid)), 
            linewidth = 0.5, linetype = "dashed", position = position_dodge(0.3), alpha = 0.5) +
  labs(y= "Arousal", x = "Tempo", title = "Arousal") +
  scale_x_discrete(limits = c("Low", "High")) + 
  scale_size_discrete(breaks = c("Solo", "Group")) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
  facet_wrap(~Solo.group)

p6 <- ggarrange(p6a, p6b,
                common.legend = TRUE,
                legend = "bottom")
p6
```

![(\#fig:unnamed-chunk-17)Data for Density and Arousal ratings from a ramdomly selected, simulated sample. Dashed linesrepresent the within-subject effect of the Tempo manipulation.](Power_analysis_TELESCOPE_files/figure-latex/unnamed-chunk-17-1.pdf) 

We fitted the models with the same structure, simply changing the dependent variable (DV): `clmm(DV ~ Tempo * Solo.group * Participant.country * Music.country + (1 + Tempo | Participant)`. To obtain *p*-values that represent main effects and interactions in an ANOVA-type manner (i.e. the intercept is the grand mean of all cells, and estimates are differences between each category mean and the mean of all categories), we used *sum-to-zero* contrasts[^1] [see e.g., @kaufmanContrastCodingLeast1974; @keppelDataAnalysisResearch1989].

[^1]: We did not do this for the models fitted for the power analysis, as it would not change the contrasts between Tempo conditions, but it is useful to display the full model results as we used `summary` (regression-type tables of estimates) instead on ANOVA-type tables (which are not available for models of `clmm` class).

In all cases, results are for main effects and all possible interactions between Tempo (Low, High), Instrumentation (Group, Solo), Participant country (Iran, Canada, Japan), and Music country (Iran, Canada, Japan) as fixed effects, as well as random intercepts and random slopes between Tempo conditions for each participant. *High* was used as reference category for Tempo, and *Japan* for both participant and music country. Contrasted levels are in square brackets. Significant effects are in bold.

As pointed out in a comment by @532079, author of the `emmeans` package [-@emmeanscit], it is crucial to consider that ordinal models rely on estimates of a dependent variable (*y*) following a logistic distribution. However, *y* itself cannot be directly observed; rather, only the interval in which it lies is observable. Consequently, *y* is an unobserved latent variable. Thus, in order to interpret the estimates, the model also calculates the thresholds (cut points) for the ordinal levels of the dependent variable.

## Visual density example model


```r
# Fit Density model
options(contrasts = c("contr.sum","contr.poly"))
model.den <- clmm(Density.tempo ~ Tempo * Solo.group * Participant.country * Music.country + 
                    (1 + Tempo | Participant),
                  data = ex.data)

# Summary table
summary.sig(model.den, "Visual density ratings model summary")
```

\begin{table}[H]
\centering
\caption{(\#tab:mod-ex-den)Visual density ratings model summary}
\centering
\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{
\begin{tabular}[t]{lcccc}
\toprule
Effects & Estimate & Std. Error & $z$ & $p$\\
\midrule
\addlinespace[0.3em]
\multicolumn{5}{l}{\textbf{Thresholds}}\\
\hline
\hspace{1em}1|2 & -1.58 & 0.12 & -13.02 & \textbf{< 0.0001}\\
\hspace{1em}2|3 & -0.04 & 0.11 & -0.40 & 0.69\\
\hspace{1em}3|4 & 1.90 & 0.12 & 15.24 & \textbf{< 0.0001}\\
\hspace{1em}4|5 & 4.29 & 0.22 & 19.70 & \textbf{< 0.0001}\\
\addlinespace[0.3em]
\hline
\multicolumn{5}{l}{\textbf{Terms}}\\
\hline
\hspace{1em}Tempo [Low] & -0.13 & 0.06 & -2.06 & \textbf{0.0392}\\
\hspace{1em}Instrumentation [Group] & 1.27 & 0.08 & 16.67 & \textbf{< 0.0001}\\
\hspace{1em}Participant country [Canada] & -0.03 & 0.15 & -0.21 & 0.83\\
\hspace{1em}Participant country [Iran] & -0.18 & 0.13 & -1.44 & 0.15\\
\hspace{1em}Music country [Canada] & -0.05 & 0.09 & -0.51 & 0.61\\
\hspace{1em}Music country [Iran] & -0.07 & 0.09 & -0.76 & 0.45\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] & 0.11 & 0.06 & 1.72 & 0.09\\
\hspace{1em}Tempo [Low] × Participant country [Canada] & 0.08 & 0.10 & 0.85 & 0.39\\
\hspace{1em}Tempo [Low] × Participant country [Iran] & -0.20 & 0.08 & -2.43 & \textbf{0.015}\\
\hspace{1em}Instrumentation [Group] × Participant country [Canada] & 0.06 & 0.10 & 0.67 & 0.5\\
\hspace{1em}Instrumentation [Group] × Participant country [Iran] & 0.07 & 0.08 & 0.81 & 0.42\\
\hspace{1em}Tempo [Low] × Music country [Canada] & -0.08 & 0.09 & -0.89 & 0.37\\
\hspace{1em}Tempo [Low] × Music country [Iran] & 0.05 & 0.09 & 0.60 & 0.55\\
\hspace{1em}Instrumentation [Group] × Music country [Canada] & 0.07 & 0.09 & 0.80 & 0.42\\
\hspace{1em}Instrumentation [Group] × Music country [Iran] & 0.01 & 0.09 & 0.09 & 0.93\\
\hspace{1em}Participant country [Canada] × Music country [Canada] & -0.08 & 0.14 & -0.57 & 0.57\\
\hspace{1em}Participant country [Iran] × Music country [Canada] & -0.03 & 0.12 & -0.26 & 0.8\\
\hspace{1em}Participant country [Canada] × Music country [Iran] & 0.06 & 0.13 & 0.47 & 0.64\\
\hspace{1em}Participant country [Iran] × Music country [Iran] & 0.04 & 0.11 & 0.35 & 0.73\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Canada] & -0.03 & 0.10 & -0.33 & 0.74\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Iran] & 0.04 & 0.08 & 0.49 & 0.62\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Music country [Canada] & 0.04 & 0.09 & 0.47 & 0.64\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Music country [Iran] & -0.05 & 0.09 & -0.57 & 0.57\\
\hspace{1em}Tempo [Low] × Participant country [Canada] × Music country [Canada] & 0.02 & 0.14 & 0.17 & 0.87\\
\hspace{1em}Tempo [Low] × Participant country [Iran] × Music country [Canada] & -0.06 & 0.12 & -0.54 & 0.59\\
\hspace{1em}Tempo [Low] × Participant country [Canada] × Music country [Iran] & -0.06 & 0.13 & -0.48 & 0.63\\
\hspace{1em}Tempo [Low] × Participant country [Iran] × Music country [Iran] & 0.10 & 0.12 & 0.91 & 0.37\\
\hspace{1em}Instrumentation [Group] × Participant country [Canada] × Music country [Canada] & -0.06 & 0.14 & -0.43 & 0.67\\
\hspace{1em}Instrumentation [Group] × Participant country [Iran] × Music country [Canada] & 0.08 & 0.12 & 0.67 & 0.5\\
\hspace{1em}Instrumentation [Group] × Participant country [Canada] × Music country [Iran] & 0.08 & 0.13 & 0.59 & 0.56\\
\hspace{1em}Instrumentation [Group] × Participant country [Iran] × Music country [Iran] & -0.13 & 0.11 & -1.11 & 0.27\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Canada] × Music country [Canada] & 0.11 & 0.14 & 0.82 & 0.41\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Iran] × Music country [Canada] & -0.06 & 0.12 & -0.54 & 0.59\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Canada] × Music country [Iran] & -0.18 & 0.13 & -1.31 & 0.19\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Iran] × Music country [Iran] & 0.04 & 0.11 & 0.31 & 0.76\\
\bottomrule
\end{tabular}}
\end{table}

## Arousal example model


```r
# Fit Arousal model
options(contrasts = c("contr.sum","contr.poly"))
model.aro <- clmm(Arousal.tempo ~ Tempo * Solo.group * Participant.country * Music.country + 
                    (1 + Tempo | Participant),
                  data = ex.data)

# Summary table
summary.sig(model.aro, "Arousal ratings model summary")
```

\begin{table}[H]
\centering
\caption{(\#tab:mod-ex-aro)Arousal ratings model summary}
\centering
\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{
\begin{tabular}[t]{lcccc}
\toprule
Effects & Estimate & Std. Error & $z$ & $p$\\
\midrule
\addlinespace[0.3em]
\multicolumn{5}{l}{\textbf{Thresholds}}\\
\hline
\hspace{1em}1|2 & -1.68 & 0.12 & -14.34 & \textbf{< 0.0001}\\
\hspace{1em}2|3 & 0.03 & 0.10 & 0.32 & 0.75\\
\hspace{1em}3|4 & 2.30 & 0.13 & 17.78 & \textbf{< 0.0001}\\
\hspace{1em}4|5 & 4.13 & 0.20 & 21.15 & \textbf{< 0.0001}\\
\addlinespace[0.3em]
\hline
\multicolumn{5}{l}{\textbf{Terms}}\\
\hline
\hspace{1em}Tempo [Low] & -0.11 & 0.06 & -1.69 & 0.09\\
\hspace{1em}Instrumentation [Group] & 1.46 & 0.08 & 18.22 & \textbf{< 0.0001}\\
\hspace{1em}Participant country [Canada] & -0.13 & 0.14 & -0.93 & 0.35\\
\hspace{1em}Participant country [Iran] & -0.05 & 0.12 & -0.42 & 0.67\\
\hspace{1em}Music country [Canada] & 0.01 & 0.09 & 0.15 & 0.88\\
\hspace{1em}Music country [Iran] & -0.19 & 0.09 & -2.14 & \textbf{0.0324}\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] & 0.05 & 0.06 & 0.78 & 0.44\\
\hspace{1em}Tempo [Low] × Participant country [Canada] & 0.06 & 0.10 & 0.65 & 0.51\\
\hspace{1em}Tempo [Low] × Participant country [Iran] & -0.04 & 0.08 & -0.52 & 0.6\\
\hspace{1em}Instrumentation [Group] × Participant country [Canada] & -0.02 & 0.10 & -0.25 & 0.81\\
\hspace{1em}Instrumentation [Group] × Participant country [Iran] & 0.01 & 0.08 & 0.13 & 0.9\\
\hspace{1em}Tempo [Low] × Music country [Canada] & -0.19 & 0.09 & -2.11 & \textbf{0.0346}\\
\hspace{1em}Tempo [Low] × Music country [Iran] & 0.01 & 0.09 & 0.08 & 0.93\\
\hspace{1em}Instrumentation [Group] × Music country [Canada] & -0.11 & 0.09 & -1.21 & 0.23\\
\hspace{1em}Instrumentation [Group] × Music country [Iran] & 0.11 & 0.09 & 1.21 & 0.22\\
\hspace{1em}Participant country [Canada] × Music country [Canada] & -0.20 & 0.13 & -1.51 & 0.13\\
\hspace{1em}Participant country [Iran] × Music country [Canada] & 0.02 & 0.12 & 0.20 & 0.84\\
\hspace{1em}Participant country [Canada] × Music country [Iran] & 0.01 & 0.14 & 0.04 & 0.96\\
\hspace{1em}Participant country [Iran] × Music country [Iran] & 0.19 & 0.12 & 1.65 & 0.1\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Canada] & 0.00 & 0.10 & 0.02 & 0.98\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Iran] & -0.05 & 0.08 & -0.56 & 0.58\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Music country [Canada] & 0.01 & 0.09 & 0.17 & 0.87\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Music country [Iran] & 0.10 & 0.09 & 1.12 & 0.26\\
\hspace{1em}Tempo [Low] × Participant country [Canada] × Music country [Canada] & -0.16 & 0.13 & -1.19 & 0.23\\
\hspace{1em}Tempo [Low] × Participant country [Iran] × Music country [Canada] & 0.02 & 0.12 & 0.21 & 0.83\\
\hspace{1em}Tempo [Low] × Participant country [Canada] × Music country [Iran] & 0.09 & 0.14 & 0.68 & 0.5\\
\hspace{1em}Tempo [Low] × Participant country [Iran] × Music country [Iran] & -0.01 & 0.12 & -0.06 & 0.95\\
\hspace{1em}Instrumentation [Group] × Participant country [Canada] × Music country [Canada] & -0.07 & 0.13 & -0.53 & 0.6\\
\hspace{1em}Instrumentation [Group] × Participant country [Iran] × Music country [Canada] & 0.32 & 0.12 & 2.80 & \textbf{0.0051}\\
\hspace{1em}Instrumentation [Group] × Participant country [Canada] × Music country [Iran] & -0.09 & 0.14 & -0.65 & 0.52\\
\hspace{1em}Instrumentation [Group] × Participant country [Iran] × Music country [Iran] & -0.21 & 0.12 & -1.77 & 0.08\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Canada] × Music country [Canada] & -0.16 & 0.13 & -1.17 & 0.24\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Iran] × Music country [Canada] & 0.06 & 0.12 & 0.52 & 0.6\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Canada] × Music country [Iran] & 0.30 & 0.14 & 2.15 & \textbf{0.0312}\\
\hspace{1em}Tempo [Low] × Instrumentation [Group] × Participant country [Iran] × Music country [Iran] & -0.15 & 0.12 & -1.26 & 0.21\\
\bottomrule
\end{tabular}}
\end{table}

## Estimated marginal means and confidence intervals for both models

Estimated marginal means and confidence intervals for this simulation are represented in \@ref(fig:plot-mod-ex).


```r
p7a <- emmip(model.den, Solo.group ~ Tempo | Participant.country + Music.country,
      CIs = TRUE) +
  labs(title = "Density")

p7b <- emmip(model.aro, Solo.group ~ Tempo | Participant.country + Music.country,
      CIs = TRUE) +
  labs(title = "Arousal")

p7 <- ggarrange(p7a, p7b,
                common.legend = TRUE,
                legend = "bottom")
p7
```

![(\#fig:plot-mod-ex)Estimated marginal means and confidence intervals for the within-subject effects of Tempo from a ramdomly selected, simulated sample. Columns represent data from simulated participants from each country, while rows represent the simulated response to music from each country.](Power_analysis_TELESCOPE_files/figure-latex/plot-mod-ex-1.pdf) 

## Effect of Tempo

Estimated marginal mean and contrast between the levels of Tempo (Low, High) averaged over the levels of instrumentation (`Solo.group`), participant country (`Participant.country`) and music country (`Music.country`). To obtain the estimated distributions of the ordinal levels of the dependent variable, the argument `mode = "mean.class"` can be added to the `emmeans` function [see @532079].

### Visual density model

Table of estimated marginal means and contrasts between tempo conditions. All estimated marginal means and contrasts were calculated using the `emmeans` function from the `emmeans` package [@emmeanscit].


```r
tempo.contr(model.den)
```

\begin{table}[H]
\centering
\caption{(\#tab:tab-den-emms)Estimated marginal means and contrasts between Tempo conditions}
\centering
\begin{threeparttable}
\begin{tabular}[t]{lccccclccc}
\toprule
\multicolumn{5}{c}{ } & \multicolumn{5}{c}{Contrasts} \\
\cmidrule(l{3pt}r{3pt}){6-10}
Tempo & EMM & $SE$ & $2.5\% CI$ & $97.5\% CI$ & Contrast & Difference & $SE$ & $z$ & $p$\\
\midrule
Low & -1.27 & 0.13 & -1.52 & -1.02 & Low - High & -0.26 & 0.13 & -2.06 & \textbf{0.0392}\\
High & -1.01 & 0.12 & -1.25 & -0.77 &  &  &  &  & \\
\bottomrule
\end{tabular}
\begin{tablenotes}[para]
\item \textit{Note: } 
\item EMM = estimated marginal mean.
\end{tablenotes}
\end{threeparttable}
\end{table}

Figure of estimated marginal means and contrasts between tempo conditions.


```r
emmip(model.den, ~Tempo, CIs = TRUE)
```

![(\#fig:unnamed-chunk-18)Estimated marginal means and confidence intervals for the within-subject effects of Tempo on Visual density ratings, from a ramdomly selected, simulated sample, averaged across the levels of instrumentation, participant country and music country.](Power_analysis_TELESCOPE_files/figure-latex/unnamed-chunk-18-1.pdf) 

### Arousal model

Table of estimated marginal means and contrasts between tempo conditions. All estimated marginal means and contrasts were calculated using the `emmeans` function from the `emmeans` package [@emmeanscit].


```r
tempo.contr(model.aro)
```

\begin{table}[H]
\centering
\caption{(\#tab:tab-aro-emms)Estimated marginal means and contrasts between Tempo conditions}
\centering
\begin{threeparttable}
\begin{tabular}[t]{lccccclccc}
\toprule
\multicolumn{5}{c}{ } & \multicolumn{5}{c}{Contrasts} \\
\cmidrule(l{3pt}r{3pt}){6-10}
Tempo & EMM & $SE$ & $2.5\% CI$ & $97.5\% CI$ & Contrast & Difference & $SE$ & $z$ & $p$\\
\midrule
Low & -1.30 & 0.12 & -1.54 & -1.07 & Low - High & -0.21 & 0.13 & -1.69 & 0.09\\
High & -1.09 & 0.11 & -1.31 & -0.87 &  &  &  &  & \\
\bottomrule
\end{tabular}
\begin{tablenotes}[para]
\item \textit{Note: } 
\item EMM = estimated marginal mean.
\end{tablenotes}
\end{threeparttable}
\end{table}

Figure of estimated marginal means and contrasts between tempo conditions.


```r
emmip(model.aro, ~Tempo, CIs = TRUE)
```

![(\#fig:unnamed-chunk-19)Estimated marginal means and confidence intervals for the within-subject effects of Tempo on Arousal ratings, from a ramdomly selected, simulated sample, averaged across the levels of instrumentation, participant country and music country.](Power_analysis_TELESCOPE_files/figure-latex/unnamed-chunk-19-1.pdf) 



# Conclusion

Based on the conducted power analysis (section \@ref(power-section)), it can be inferred that a sample size of 80 participants per country (240 participants in total) would yield a statistical power ($1 - \beta$) greater than 0.95 for the planned contrasts between Tempo conditions when the effect size of this Tempo manipulation is 1. While the pilot data indicated a larger effect size, a conservative estimate was utilized to ensure sufficient statistical power even if the true effect size is significantly smaller than observed in the pilot study. This approach guarantees a high likelihood of correctly rejecting the null hypothesis in the presence of a true effect, providing confidence in the study's ability to detect meaningful findings.

In addition, as shown when fitting and reporting the models with a randomly selected sample (section \@ref(example-data)), the sample size and data structure are adequate to fit the planned models.

# References {.unnumbered #refs}

\begin{multicols}{2}
\AtNextBibliography{\footnotesize}
\printbibliography[heading=none]
\normalsize
\end{multicols}

\def\printbibliography{}

# Session info (for reproducibility) {.unnumbered #session}


```r
library(pander)
pander(sessionInfo(), locale = FALSE)
```

**R version 4.3.3 (2024-02-29 ucrt)**

**Platform:** x86_64-w64-mingw32/x64 (64-bit) 


**attached base packages:** 
_stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_pander(v.0.6.5)_, _beepr(v.1.3)_, _scales(v.1.3.0)_, _performance(v.0.10.9)_, _kableExtra(v.1.4.0)_, _ordinal(v.2023.12-4)_, _emmeans(v.1.10.0)_, _ggpubr(v.0.6.0)_, _lubridate(v.1.9.3)_, _forcats(v.1.0.0)_, _stringr(v.1.5.1)_, _dplyr(v.1.1.4)_, _purrr(v.1.0.2)_, _readr(v.2.1.5)_, _tidyr(v.1.3.1)_, _tibble(v.3.2.1)_, _ggplot2(v.3.5.0)_, _tidyverse(v.2.0.0)_, _knitr(v.1.45)_ and _tictoc(v.1.2)_

**loaded via a namespace (and not attached):** 
_tidyselect(v.1.2.0)_, _viridisLite(v.0.4.2)_, _farver(v.2.1.1)_, _fastmap(v.1.1.1)_, _TH.data(v.1.1-2)_, _digest(v.0.6.34)_, _estimability(v.1.5)_, _timechange(v.0.3.0)_, _lifecycle(v.1.0.4)_, _survival(v.3.5-8)_, _magrittr(v.2.0.3)_, _compiler(v.4.3.3)_, _rlang(v.1.1.3)_, _tools(v.4.3.3)_, _utf8(v.1.2.4)_, _yaml(v.2.3.8)_, _ggsignif(v.0.6.4)_, _labeling(v.0.4.3)_, _xml2(v.1.3.6)_, _multcomp(v.1.4-25)_, _abind(v.1.4-5)_, _rsconnect(v.1.2.1)_, _withr(v.3.0.0)_, _numDeriv(v.2016.8-1.1)_, _grid(v.4.3.3)_, _fansi(v.1.0.6)_, _xtable(v.1.8-4)_, _colorspace(v.2.1-0)_, _MASS(v.7.3-60.0.1)_, _tinytex(v.0.49)_, _insight(v.0.19.8)_, _cli(v.3.6.2)_, _mvtnorm(v.1.2-4)_, _rmarkdown(v.2.25)_, _generics(v.0.1.3)_, _rstudioapi(v.0.15.0)_, _tzdb(v.0.4.0)_, _audio(v.0.1-11)_, _splines(v.4.3.3)_, _vctrs(v.0.6.5)_, _Matrix(v.1.6-5)_, _sandwich(v.3.1-0)_, _carData(v.3.0-5)_, _bookdown(v.0.37)_, _car(v.3.1-2)_, _hms(v.1.1.3)_, _rstatix(v.0.7.2)_, _systemfonts(v.1.0.5)_, _glue(v.1.7.0)_, _codetools(v.0.2-19)_, _cowplot(v.1.1.3)_, _stringi(v.1.8.3)_, _gtable(v.0.3.4)_, _munsell(v.0.5.0)_, _pillar(v.1.9.0)_, _htmltools(v.0.5.7)_, _R6(v.2.5.1)_, _ucminf(v.1.2.1)_, _evaluate(v.0.23)_, _lattice(v.0.22-5)_, _highr(v.0.10)_, _backports(v.1.4.1)_, _broom(v.1.0.5)_, _Rcpp(v.1.0.12)_, _svglite(v.2.1.3)_, _coda(v.0.19-4.1)_, _gridExtra(v.2.3)_, _nlme(v.3.1-164)_, _xfun(v.0.42)_, _zoo(v.1.8-12)_ and _pkgconfig(v.2.0.3)_
