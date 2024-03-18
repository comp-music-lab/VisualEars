library(tictoc)
tic()
rmarkdown::render("Power_analysis_TELESCOPE.Rmd", clean = FALSE)
#Sys.sleep(1)
toc(log = TRUE)
tic.log(format = TRUE)

13736/60
6475
