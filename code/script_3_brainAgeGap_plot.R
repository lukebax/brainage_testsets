library(dabestr)

df <- read.csv("../data/data_dabestrDF_leuven.csv")
dfROB <- read.csv("../data/data_dabestrDF_ROB_leuven.csv")

shared.control <- dabest(df, Group, BrainAgeGap, idx = c("Normal", "Mild", "Severe"))
shared.controlROB <- dabest(dfROB, Group, BrainAgeGap, idx = c("Normal", "Mild", "Severe"))

shared.control.mean_diff <- shared.control %>% mean_diff()
shared.controlROB.mean_diff <- shared.controlROB %>% mean_diff()

plot(shared.control.mean_diff, rawplot.type = "swarmplot", rawplot.ylim = c(0, 2.6), effsize.ylim = c(-0.6, 1.5))
plot(shared.controlROB.mean_diff, rawplot.type = "swarmplot", rawplot.ylim = c(0, 2.6), effsize.ylim = c(-0.6, 1.5))
