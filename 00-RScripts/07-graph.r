# Just a simple script to graph data. No really
# documentation here. Mainly because I did this quickly

library("ggplot2");
library("ggpubr");
library("data.table");
library("viridis");

fileStr="07-artic-1st-bench/07-artic-1st-bench-stats-filter.tsv";
dataDF = read.csv(fileStr,sep = "\t",header = TRUE);

graphObj = NULL;

# For some odd reason R is ignoring the first column
graphObj = 
   ggplot(
      data = dataDF[
         dataDF$lenRef > 20000 &
         dataDF$index < 1200
         ,
      ],
      aes(
         x = noReads / 100, # 1% = 100 otherwise
         y = matches,
         color = as.character(index),
         fill = as.character(index),
            # R is ignoring the first column
         shape = as.character(index),
     ) # aes
); # make graph

# These are the 0.2% and 0.1% lines
graphObj = graphObj+geom_hline(yintercept = 0.002 * 29803);
graphObj =
   graphObj +
   geom_label(
      data = NULL,
      mapping = NULL,
      label = "0.2% difference",
      x = 1,
      y = 0.002 * 29803,
      inherit.aes = FALSE
   );

graphObj =
   graphObj + 
   geom_hline(yintercept = 0.001 * 29903, color="PURPLE");
graphObj =
   graphObj +
   geom_label(
      data = NULL,
      mapping = NULL,
      label = "0.1% difference",
      x = 1,
      y = 0.001 * 29803,
      inherit.aes = FALSE,
      color = "PURPLE"
   );

graphObj =
   graphObj + 
   geom_hline(yintercept = 0.0005 * 29903, color="PURPLE");
graphObj =
   graphObj +
   geom_label(
      data = NULL,
      mapping = NULL,
      label = "0.05% difference",
      x = 1,
      y = 0.0005 * 29803,
      inherit.aes = FALSE,
      color = "MAGENTA"
   );

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 0.3)
); # Plot points

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = -1,
      option = "D",
      name = "Depth"
); # Add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = -1,
      option = "D",
      name = "Depth"
); # Add in color scheme


graphObj = graphObj + scale_shape_discrete(name = "Depth");

graphObj =
   graphObj +
   xlab("Reference % difference from sequence");
graphObj = graphObj + ylab("Number of false snps");

graphObj = graphObj + theme_pubr();

ggsave("artic-initial-bench.svg",device = "svg",dpi = 300);

