# This should be called as "Rscript 00-Rscripts/10-graph.r"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-00: 
#  - Functions and libraries
#  o sec-00 sub-01:
#    - Libraries used
#  o sec-00 sub-02:
#    - saveggplot: save a ggplot graph (wrapper for ggsave)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-00 Sub-01:
#  - Libraries used
#**********************************************************

#args = commandArgs(); # User command line input
   # openbsd or vim does not supply file path to script

library("ggplot2");
library("ggpubr");
library("data.table");
library("viridis");

#**********************************************************
# Sec-00 Sub-02:
#  - save ggplot
#**********************************************************

#----------------------------------------------------------
# Name: saveggplot
# Use:
#  - Is a wrapper for ggsave (saves ggplot graph)
# Input:
#  - nameStr:
#    o Name of the output graph (no extension)
#  - svgBl:
#    o 1: Ouput graph as an svg file
#    o 0: Output graph as an png file
# Output:
#  - graph saved as an svg or png (depends on svgBl)
#----------------------------------------------------------
saveggplot = function (
   nameStr,  # Name of gaph to save
   svgBl = 1 # 1: Save as svg file; 0 save as png
) { # Save a ggplot graph (wrapper for ggsave)
    if(svgBl > 0)
    { # If I am saving this as a svg file
       ggsave(paste(nameStr, ".svg", sep = ""), 
              device = "svg", # save as tiff file
              dpi = 300,
       ); # ggsave (save the graph)
    } # If I am saving this as a svg file

    else
    { # If I am saving this as a png file
       ggsave(paste(nameStr, ".png", sep = ""), 
              device = "png", # save as tiff file
              dpi = 300,
       ); # ggsave (save the graph)
    } # If I am saving this as a png file
} # saveggplot

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01: 
#  - Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

fileStr="10-new-bench/10-bench-artic-ivar-stats.tsv";
dataDF = read.table(fileStr, sep = "\t", header = TRUE);
asSvgBl = 1; # 1: save files as svg; 0: save as png

graphObj = NULL;
lastCompleteRepI = 8;

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02: 
#  - Prepare data
#  o sec-02 sub-01:
#    - Clean up invalid entries (errored out)
#  o sec-02 sub-02:
#    - Set up test datatype (program-ivar-scaffoldMethod)
#  o sec-02 sub-03:
#    - Re-factor data (output order)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-02 Sub-01:
#  - Clean up invalid entries (errored out)
#**********************************************************

dataDF = dataDF[dataDF$snp != "ERR",];

#**********************************************************
# Sec-02 Sub-02:
#  - Set up test datatype (program-ivar-scaffoldMethod)
#**********************************************************

# Add in if used stich/scaffold builder
dataDF$test =
   paste(
      dataDF$program,
      dataDF$conScaffoldBuilder,
      dataDF$conStich,
      sep = "-"
); # paste test together

dataDF$test = sub("-false-false", "-NA", dataDF$test);
dataDF$test = sub("-true-false", "-scaffold", dataDF$test);
dataDF$test=sub("-false-true","-stich",dataDF$test);
dataDF$test=sub("Ivar","-ivar",dataDF$test);

#**********************************************************
# Sec-02 Sub-03:
#  - Refactor data (output order)
#**********************************************************

dataDF$test =
   factor(
      x = dataDF$test, 
      levels =
         c(
            "artic-NA",
            "ivar-NA",
            "ivarTrim-NA",
            "Lilo-NA",
            "Lilo-ivar-NA",
            "Lilo-scaffold",
            "Lilo-ivar-scaffold",
            "Lilo-stich",
            "Lilo-ivar-stich",
            "buildCon-stich",
            "buildCon-ivar-stich"
      ) # Levels
); # Orgainzie my data (factor)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03: 
#  - Graph Lilo snp (% reference)
#  o sec-03 sub-01:
#    - Make the graph
#  o sec-03 sub-02:
#    - Add in percent error markers
#  o sec-03 sub-03:
#    - Add in the points
#  o sec-03 sub-04:
#    - Add in the colors and shapes
#  o sec-03 sub-05:
#    - Add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-03 Sub-01:
#  - Make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "Lilo" |
              dataDF$program == "LiloIvar") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-03 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.002 * 29803, lty=3);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.002 * 29803,
      label = "0.2% difference"
); # add an 0.2% line marker

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29803,lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29803,lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-03 Sub-03:
#  - Add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 0.3)
); # Plot points

#**********************************************************
# Sec-03 Sub-04:
#  - Add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "Program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

#**********************************************************
# Sec-03 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Reference % difference from sequence");
graphObj = graphObj + ylab("Number of false snps");

graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/lilo-snps-ref", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-04: 
#  - Graph Lilo snp (by depth)
#  o sec-04 sub-01:
#    - Make the graph
#  o sec-04 sub-02:
#    - Add in percent error markers
#  o sec-04 sub-03:
#    - Add in the points
#  o sec-04 sub-04:
#    - Add in the colors and shapes
#  o sec-04 sub-05:
#    - Add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-04 Sub-01:
#  - Make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "Lilo" |
              dataDF$program == "LiloIvar") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth),
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-04 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.002 * 29803, lty=3);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.002 * 29803,
      label = "0.2% difference"
); # add an 0.2% line marker

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29803, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29803, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-04 Sub-03:
#  - Add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 10)
); # Plot points

#**********************************************************
# Sec-04 Sub-04:
#  - Add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "Program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

#**********************************************************
# Sec-04 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Read depth");
graphObj = graphObj + ylab("Number of false snps");
graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/lilo-snps-depth", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-05 
#  - graph lilo snp / indels
#  o sec-05 sub-01:
#    - make the graph
#  o sec-05 sub-02:
#    - add in percent error markers
#  o sec-05 sub-03:
#    - add in the points
#  o sec-05 sub-04:
#    - add in the colors and shapes
#  o sec-05 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-05 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "Lilo" |
              dataDF$program == "LiloIvar") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(ins) + as.numeric(del),
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-04 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj +
   geom_hline(yintercept = 0.002 * 29803, lty = 3);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.002 * 29803, lty = 3);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.002 * 29803,
      y = 0.002 * 29803,
      label = "0.2% difference"
); # add an 0.2% line marker

graphObj =
   graphObj +
   geom_hline(yintercept = 0.001 * 29803, lty = 2);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.001 * 29803, lty = 2);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.001 * 29803,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj +
   geom_hline(yintercept=0.0005 * 29803, lty = 5);

graphObj =
   graphObj +
   geom_vline(xintercept=0.0005 * 29803, lty = 5);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.0005 * 29803 + 2,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-04 Sub-03:
#  - Add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 10)
); # Plot points

#**********************************************************
# Sec-05 Sub-03:
#  - add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 0.3)
); # plot points

#**********************************************************
# Sec-05 Sub-04:
#  - add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "program"
); # add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = 1,
       option = "D",
      name = "program"
); # add in color scheme

#**********************************************************
# Sec-05 Sub-05:
#  - add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj + xlab("false insertions + false deletions");
graphObj = graphObj + ylab("number of false snps");

graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/lilo-snps-indels", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-06: 
#  - Graph buildCon snp (% reference)
#  o sec-06 sub-01:
#    - Make the graph
#  o sec-06 sub-02:
#    - Add in percent error markers
#  o sec-06 sub-03:
#    - Add in the points
#  o sec-06 sub-04:
#    - Add in the colors and shapes
#  o sec-06 sub-05:
#    - Add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-06 Sub-01:
#  - Make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$test == "buildCon-stich" |
              dataDF$test == "buildCon-ivar-stich") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-06 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29803, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29803, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-06 Sub-03:
#  - Add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 0.3)
); # Plot points

#**********************************************************
# Sec-06 Sub-04:
#  - Add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "Program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = -1,
      option = "D",
      name = "Program"
); # Add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = -1,
      option = "D",
      name = "Program"
); # Add in color scheme

#**********************************************************
# Sec-06 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Reference % difference from sequence");
graphObj = graphObj + ylab("Number of false snps");

graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/buildCon-snps-ref", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-07: 
#  - Graph buildCon snp (by depth)
#  o sec-07 sub-01:
#    - Make the graph
#  o sec-07 sub-02:
#    - Add in percent error markers
#  o sec-07 sub-03:
#    - Add in the points
#  o sec-07 sub-04:
#    - Add in the colors and shapes
#  o sec-07 sub-05:
#    - Add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-07 Sub-01:
#  - Make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "buildCon" |
              dataDF$program == "buildConIvar") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth),
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-07 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj + geom_hline(yintercept = 0.001 * 29803,lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj + geom_hline(yintercept=0.0005 * 29803,lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-07 Sub-03:
#  - Add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 10)
); # Plot points

#**********************************************************
# Sec-07 Sub-04:
#  - Add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "Program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = -1,
      option = "D",
      name = "Program"
); # Add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = -1,
      option = "D",
      name = "Program"
); # Add in color scheme

#**********************************************************
# Sec-07 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Read depth");
graphObj = graphObj + ylab("Number of false snps");
graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/buildCon-snps-depth", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-08: 
#  - graph buildCon snp / indels
#  o sec-08 sub-01:
#    - make the graph
#  o sec-08 sub-02:
#    - add in percent error markers
#  o sec-08 sub-03:
#    - add in the points
#  o sec-08 sub-04:
#    - add in the colors and shapes
#  o sec-08 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-08 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "buildCon" |
              dataDF$program == "buildConIvar") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(ins) + as.numeric(del),
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-04 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj +
   geom_hline(yintercept = 0.002 * 29803, lty = 3);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.002 * 29803, lty = 3);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.002 * 29803,
      y = 0.002 * 29803,
      label = "0.2% difference"
); # add an 0.2% line marker

graphObj =
   graphObj +
   geom_hline(yintercept = 0.001 * 29803, lty = 2);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.001 * 29803, lty = 2);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.001 * 29803,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj +
   geom_hline(yintercept=0.0005 * 29803, lty = 5);

graphObj =
   graphObj +
   geom_vline(xintercept=0.0005 * 29803, lty = 5);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.0005 * 29803 + 2,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker


#**********************************************************
# Sec-08 Sub-03:
#  - add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 0.3)
); # plot points

#**********************************************************
# Sec-08 Sub-04:
#  - add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = -1,
       option = "D",
      name = "program"
); # add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = -1,
       option = "D",
      name = "program"
); # add in color scheme

#**********************************************************
# Sec-08 Sub-05:
#  - add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("false insertions + false deletions");
graphObj = graphObj + ylab("number of false snps");

graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/buildCon-snps-indels", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-09: 
#  - Graph ivar snp (% reference)
#  o sec-09 sub-01:
#    - Make the graph
#  o sec-09 sub-02:
#    - Add in percent error markers
#  o sec-09 sub-03:
#    - Add in the points
#  o sec-09 sub-04:
#    - Add in the colors and shapes
#  o sec-09 sub-05:
#    - Add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-09 Sub-01:
#  - Make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "ivarTrim") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-09 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29803, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29803, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-09 Sub-03:
#  - Add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 0.3)
); # Plot points

#**********************************************************
# Sec-09 Sub-04:
#  - Add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "Program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = -1,
      option = "D",
      name = "Program"
); # Add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = -1,
      option = "D",
      name = "Program"
); # Add in color scheme

#**********************************************************
# Sec-09 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Reference % difference from sequence");
graphObj = graphObj + ylab("Number of false snps");

graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/ivar-snps-ref", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-10: 
#  - Graph ivar snp (by depth)
#  o sec-10 sub-01:
#    - Make the graph
#  o sec-10 sub-02:
#    - Add in percent error markers
#  o sec-10 sub-03:
#    - Add in the points
#  o sec-10 sub-04:
#    - Add in the colors and shapes
#  o sec-10 sub-05:
#    - Add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-10 Sub-01:
#  - Make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "ivarTrim") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth),
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-10 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj + geom_hline(yintercept = 0.001 * 29803,lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj + geom_hline(yintercept=0.0005 * 29803,lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-10 Sub-03:
#  - Add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 10)
); # Plot points

#**********************************************************
# Sec-10 Sub-04:
#  - Add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "Program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = -1,
      option = "D",
      name = "Program"
); # Add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = -1,
      option = "D",
      name = "Program"
); # Add in color scheme

#**********************************************************
# Sec-10 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Read depth");
graphObj = graphObj + ylab("Number of false snps");
graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/ivar-snps-depth", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-11: 
#  - graph ivar snp / indels
#  o sec-11 sub-01:
#    - make the graph
#  o sec-11 sub-02:
#    - add in percent error markers
#  o sec-11 sub-03:
#    - add in the points
#  o sec-11 sub-04:
#    - add in the colors and shapes
#  o sec-11 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-11 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "ivarTrim") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(ins) + as.numeric(del),
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-04 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj +
   geom_hline(yintercept = 0.002 * 29803, lty = 3);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.002 * 29803, lty = 3);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.002 * 29803,
      y = 0.002 * 29803,
      label = "0.2% difference"
); # add an 0.2% line marker

graphObj =
   graphObj +
   geom_hline(yintercept = 0.001 * 29803, lty = 2);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.001 * 29803, lty = 2);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.001 * 29803,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj +
   geom_hline(yintercept=0.0005 * 29803, lty = 5);

graphObj =
   graphObj +
   geom_vline(xintercept=0.0005 * 29803, lty = 5);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.0005 * 29803 + 2,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-11 Sub-03:
#  - add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 10)
); # plot points

#**********************************************************
# Sec-11 Sub-04:
#  - add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = -1,
       option = "D",
      name = "program"
); # add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = -1,
       option = "D",
      name = "program"
); # add in color scheme

#**********************************************************
# Sec-11 Sub-05:
#  - add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("false insertions + false deletions");
graphObj = graphObj + ylab("number of false snps");

graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/ivar-snps-indels", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-12: 
#  - graph each pipelines snps by reference
#  o sec-12 sub-01:
#    - make the graph
#  o sec-12 sub-02:
#    - add in percent error markers
#  o sec-12 sub-03:
#    - add in the points
#  o sec-12 sub-04:
#    - add in the colors and shapes
#  o sec-12 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-12 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-ivar-stich" |
              dataDF$test == "Lilo-stich" |
              dataDF$program == "buildConIvar") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-12 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29803, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29803, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-12 Sub-03:
#  - Add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 0.3)
); # Plot points

#**********************************************************
# Sec-12 Sub-04:
#  - Add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "Program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

#**********************************************************
# Sec-12 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Reference % difference from sequence");
graphObj = graphObj + ylab("Number of false snps");

graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/all-snps-ref", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-13: 
#  - graph each pipelines snps by depth
#  o sec-13 sub-01:
#    - make the graph
#  o sec-13 sub-02:
#    - add in percent error markers
#  o sec-13 sub-03:
#    - add in the points
#  o sec-13 sub-04:
#    - add in the colors and shapes
#  o sec-13 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-13 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-ivar-stich" |
              dataDF$test == "Lilo-stich" |
              dataDF$program == "buildConIvar") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth), # 100 = 1%
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-13 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29803, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29803, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-13 Sub-03:
#  - Add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 10)
); # Plot points

#**********************************************************
# Sec-13 Sub-04:
#  - Add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "Program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

#**********************************************************
# Sec-13 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Read depth");
graphObj = graphObj + ylab("Number of false snps");

graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/all-snps-depth", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-14: 
#  - graph each pipelines snps by depth
#  o sec-14 sub-01:
#    - make the graph
#  o sec-14 sub-02:
#    - add in percent error markers
#  o sec-14 sub-03:
#    - add in the points
#  o sec-14 sub-04:
#    - add in the colors and shapes
#  o sec-14 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-14 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-ivar-stich" |
              dataDF$test == "Lilo-stich" |
              dataDF$program == "buildConIvar") &
            dataDF$lenRef > 20000 &
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(ins) + as.numeric(del),
         y = as.numeric(snp),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-14 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj +
   geom_hline(yintercept = 0.002 * 29803, lty = 3);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.002 * 29803, lty = 3);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.002 * 29803,
      y = 0.002 * 29803,
      label = "0.2% difference"
); # add an 0.2% line marker

graphObj =
   graphObj +
   geom_hline(yintercept = 0.001 * 29803, lty = 2);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.001 * 29803, lty = 2);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.001 * 29803,
      y = 0.001 * 29803,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj +
   geom_hline(yintercept=0.0005 * 29803, lty = 5);

graphObj =
   graphObj +
   geom_vline(xintercept=0.0005 * 29803, lty = 5);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.0005 * 29803 + 2,
      y = 0.0005 * 29803,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-14 Sub-03:
#  - Add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 0.3)
); # Plot points

#**********************************************************
# Sec-14 Sub-04:
#  - Add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "Program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

#**********************************************************
# Sec-14 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Number of false indels");
graphObj = graphObj + ylab("Number of false snps");

graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/all-snps-indels", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-15: 
#  - graph each pipelines snps by depth
#  o sec-15 sub-01:
#    - make the graph
#  o sec-15 sub-02:
#    - add in percent error markers
#  o sec-15 sub-03:
#    - add in the points
#  o sec-15 sub-04:
#    - add in the colors and shapes
#  o sec-15 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-15 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-ivar-stich" |
              dataDF$test == "Lilo-stich" |
              dataDF$program == "buildConIvar") &
            dataDF$lenRef > 20000 &
            dataDF$Ns < 29000 & # Removes a few artic runs
            dataDF$rep > lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth),
         y = as.numeric(Ns),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-15 Sub-02:
#  - Add in percent error markers
#**********************************************************

#**********************************************************
# Sec-15 Sub-03:
#  - Add in the points
#**********************************************************

graphObj =
   graphObj +
   geom_point(
      size = 3,
      alpha = 0.5,
      position = position_jitter(height = 0, width = 10)
); # Plot points

#**********************************************************
# Sec-15 Sub-04:
#  - Add in the colors and shapes
#**********************************************************

graphObj=graphObj + scale_shape_discrete(name = "Program");

graphObj =
   graphObj +
   scale_fill_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

graphObj =
   graphObj +
   scale_color_viridis(
      discrete = TRUE,
      direction = 1,
      option = "D",
      name = "Program"
); # Add in color scheme

#**********************************************************
# Sec-15 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Depth");
graphObj = graphObj + ylab("Number of masked bases");

graphObj = graphObj + theme_pubr();

saveggplot("10-new-bench/all-masked", asSvgBl);
