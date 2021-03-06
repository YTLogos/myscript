library(tidyverse)
# devtools::install_github("EmilHvitfeldt/paletteer")
library(paletteer)
library(lubridate)
barchart_theme <-  theme(axis.text.y   = element_text(size=13, face="bold", colour = "black"),
                  axis.text.x   = element_text(size=13, face="bold", colour = "black"),
                  axis.title.x  = element_text(size=13, face="bold"),
                  axis.title.y  = element_text(size=13, face = "bold"),
                  panel.background = element_rect(fill = "white"),
                  axis.ticks.length = unit(.25, "cm"),
                  plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5),
                  plot.caption = element_text(size = 10),
                  plot.subtitle = element_text(hjust = 0.5, size = 10),
                  strip.background = element_rect(fill = "white"),
                  strip.text = element_text(size = 18, hjust = 0, colour = "black", face ="bold"),
                  legend.position = "none",
                  legend.title = element_text(size = 13, face = "bold"),
                  legend.text = element_text(size = 13),
                  panel.border = element_rect(color = "black", fill = NA, size = 2),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed"))