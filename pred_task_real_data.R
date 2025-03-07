# library(evd)
# library(extRemes)
library(mvtnorm)
library(tidyverse)
library(patchwork)
library(grafify)
library(qs)

data_type <- "redstone"
data <- qread("data/redstone_expo.qs")$cloud_tib

data |> as_tibble() |> ggplot(aes(x = x, y = y)) + 
  # geom_point(size = 0.5) +
  geom_point() +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0,12), expand = expansion(mult = c(0,0.01))) +
  scale_y_continuous(limits = c(0, 15), expand = expansion(mult = c(0,0.01))) +
  # annotate("rect", xmin = 7, xmax = 9, ymin = 12,  ymax = 14, fill = get_graf_colours("ok_orange"), color = "black", alpha = 0.75) +
  # annotate("rect", xmin = 4, xmax = 6, ymin = 8,  ymax = 10, fill = get_graf_colours("ok_bluegreen"), color = "black", alpha = 0.75) +
  # annotate("rect", xmin = 8, xmax = 10, ymin = 2,  ymax = 4, fill = get_graf_colours("ok_redpurple"), color = "black", alpha = 0.75) +
  xlab(expression("ERC")) + ylab(expression("FWI"))

ggsave("~/Desktop/research/posters-presentations/NSA_talk/redstone_pred_task_b3.pdf",
       dpi = 320,
       width = 4,
       height = 4,
       bg = 'transparent')
knitr::plot_crop("~/Desktop/research/posters-presentations/NSA_talk/redstone_pred_task_b3.pdf")
