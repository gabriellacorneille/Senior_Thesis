
setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

library(cowplot)

#focal plots
fig_2 <- plot_grid(fig_2_pt1, fig_2_pt2,
                   labels = LETTERS[1:2], nrow = 2) ; fig_2

ggsave("TablesFigures/focal ID -> density/2016_south_combined_plot.png",
       width = 8,
       height = 10,
       dpi = 600)

#to do this you need to have both figures in your environment, so run plots from length vs density (female) and female length vs position
#female length plots
fig_4 <- plot_grid(fig_4_pt1, fig_4_pt2,
                   labels = LETTERS[1:2], nrow = 2) ; fig_4

ggsave("TablesFigures/female length vs harem position/female length plots.png",
       width = 8,
       height = 10,
       units = "in",
       dpi = 600)

#run the male length plots now (male length plots)
fig_6 <- plot_grid(fig_6_pt1, fig_6_pt2,
                   labels = LETTERS[1:2], nrow = 2); fig_6

ggsave(
  "TablesFigures/bigger male bigger harem?/male length plots.png",
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)
