library(dabestr)

df <- read.csv("../data/data_dabestrDF_leuven.csv")
dfROB <- read.csv("../data/data_dabestrDF_ROB_leuven.csv")

shared_control <- load(df, x = Group, y = BrainAgeGap, idx = c("Normal", "Mild", "Severe"))
shared_controlROB <- load(dfROB, x = Group, y = BrainAgeGap, idx = c("Normal", "Mild", "Severe"))

shared_control.mean_diff <- mean_diff(shared_control)
shared_controlROB.mean_diff <- mean_diff(shared_controlROB)

dabest_plot(shared_control.mean_diff,
            swarm_label = "Brain Age Gap",
            swarm_x_text = 13,
            contrast_x_text = 13,
            swarm_ylim = c(0, 2.7),
            contrast_ylim = c(-0.6, 1.5),
            custom_palette = "d3")

dabest_plot(shared_controlROB.mean_diff,
            swarm_label = "Brain Age Gap",
            swarm_x_text = 13,
            contrast_x_text = 13,
            swarm_ylim = c(0, 2.7),
            contrast_ylim = c(-0.6, 1.5),
            custom_palette = "d3")
