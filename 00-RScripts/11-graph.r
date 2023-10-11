# This should be called as "Rscript 00-Rscripts/11-graph.r"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 11-graph.r TOC: Table Of Contents
#  o sec-00: 
#    - Functions and libraries
#  o sec-01: 
#    - Variable declerations
#  o sec-02: 
#    - Prepare data
#  o sec-03: 
#    - Graph Lilo snp (% reference)
#  o sec-04: 
#    - Graph Lilo snp (by depth)
#  o sec-05: 
#    - graph lilo snp / indels
#  o sec-06: 
#    - Graph buildCon snp (% reference)
#  o sec-07: 
#    - Graph buildCon snp (by depth)
#  o sec-08: 
#    - graph buildCon snp / indels
#  o sec-09: 
#    - Graph ivar snp (% reference)
#  o sec-10: 
#    - Graph ivar snp (by depth)
#  o sec-11: 
#    - graph ivar snp / indels
#  o sec-12: 
#    - graph each pipelines snps by reference
#  o sec-13: 
#    - graph each pipelines snps by depth
#  o sec-14: 
#    - graph each pipelines snps by depth
#  o sec-15: 
#    - graph each pipelines number of masked bases
#  o sec-16: 
#    - graph each pipelines (over 30x depth) snps & indels 
#  o sec-17: 
#    - graph each pipelines snps by ref (over 30x depth)
#  o sec-18: 
#    - graph each pipelines snps by depth (over 30x)
#  o sec-19: 
#    - graph each pipelines (50x depth) snps and indels 
#  o sec-20: 
#    - graph each pipelines snps by ref (50x read depth)
#  o sec-21: 
#    - graph each pipelines (100x depth) snps and indels 
#  o sec-22: 
#    - graph each pipelines snps by ref (100x read depth)
#  o sec-23: 
#    - graph each pipelines (300x depth) snps and indels 
#  o sec-24: 
#    - graph each pipelines snps by ref (300x read depth)
#  o sec-25: 
#    - graph each pipelines (500x depth) snps and indels 
#  o sec-26: 
#    - graph each pipelines snps by ref (500x read depth)
#  o sec-27: 
#    - graph each pipelines (1000x depth) snps and indels 
#  o sec-28: 
#    - graph each pipelines snps by ref (1000x read depth)
#  o sec-29: 
#    - graph each consensus length by depth (>30x, all)
#  o sec-30: 
#    - graph each consensus length by ref (depth >30x, all)
#  o sec-31: 
#    - graph the time usage for each method
#  o sec-32: 
#    - graph the time usage for each method
#  o sec-33: 
#    - graph lilo by consensus length and reference
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-00: 
#  - Functions and libraries
#  o sec-00 sub-01:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

fileStr="11-quast-bench/11-quast-bench-cleanup.tsv";
dataDF = read.table(fileStr, sep = "\t", header = TRUE);
asSvgBl = 1; # 1: save files as svg; 0: save as png

graphObj = NULL;
lastCompleteRepI = 1;

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02: 
#  - Prepare data
#  o sec-02 sub-01:
#    - Clean up invalid entries (errored out)
#  o sec-02 sub-02:
#    - Convert trues to actual program names
#  o sec-02 sub-03:
#    - Set up test datatype (program-ivar-scaffoldMethod)
#  o sec-02 sub-04:
#    - Re-factor data (output order)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-02 Sub-01:
#  - Clean up invalid entries (errored out)
#**********************************************************

dataDF = dataDF[dataDF$snp != "ERR",];

#**********************************************************
# Sec-02 Sub-02:
#  - Convert trues to actual program names
#**********************************************************

dataDF$usedMedaka =gsub("TRUE","medaka",dataDF$usedMedaka);
dataDF$usedMedaka = gsub("FALSE", "", dataDF$usedMedaka);

dataDF$usedMajcon =gsub("TRUE","majcon",dataDF$usedMajcon);
dataDF$usedMajcon =gsub("FALSE","",dataDF$usedMajcon);

dataDF$usedRacon =gsub("TRUE","racon",dataDF$usedRacon);
dataDF$usedRacon =gsub("FALSE","",dataDF$usedRacon);

dataDF$usedIvar =gsub("TRUE","ivar",dataDF$usedIvar);
dataDF$usedRacon =gsub("FALSE","",dataDF$usedRacon);

dataDF$ivarPolish =
   gsub("TRUE", "ivarPolish", dataDF$ivarPolish);
dataDF$ivarPolish = gsub("FALSE","",dataDF$ivarPolish);

dataDF$medakaPolish =
   gsub("TRUE", "medakaPolish", dataDF$medakaPolish);
dataDF$medakaPolish = gsub("FALSE","",dataDF$medakaPolish);

dataDF$conDef = gsub("TRUE", "conDefault", dataDF$conDef);
dataDF$conDef = gsub("FALSE", "", dataDF$conDef);

dataDF$conStich = gsub("TRUE", "stich", dataDF$conStich);
dataDF$conStich = gsub("FALSE", "", dataDF$conStich);

dataDF$conScaffoldBuilder =
   gsub("TRUE", "scaffold", dataDF$conScaffoldBuilder);
dataDF$conScaffoldBuilder =
   gsub("FALSE", "", dataDF$conScaffoldBuilder);

#**********************************************************
# Sec-02 Sub-02:
#  - Set up the test names
#**********************************************************

# Add in if used stich/scaffold builder
dataDF$test =
   paste(
      dataDF$program,
      dataDF$usedMedaka,
      dataDF$usedMajcon,
      dataDF$useIvar,
      dataDF$useRacon,
      dataDF$medakaPolish,
      dataDF$ivarPolish,
      dataDF$conScaffoldBuilder,
      dataDF$conStich,
      sep = "-"
); # paste test together

# Remove any exta dashes
dataDF$test = gsub("--*", "-", dataDF$test);
dataDF$test = sub("-$", "", dataDF$test);

#**********************************************************
# Sec-02 Sub-03:
#  - Refactor data (output order)
#**********************************************************

dataDF$test =
   factor(
      x = dataDF$test, 
      levels =
         c(
            "artic-medaka",
            "ivar",
            "ivarTrim",
            "Lilo-medaka",
            "Lilo-medaka-ivarPolish",
            "Lilo-medaka-scaffold",
            "Lilo-medaka-ivarPolish-scaffold",
            "Lilo-medaka-stich",
            "Lilo-medaka-ivarPolish-stich",
            "buildCon-majcon-stich",
            "buildCon-majcon-ivarPolish-stich"
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
            dataDF$program == "Lilo" &
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-03 Sub-02:
#  - Add in percent error markers
#**********************************************************

#graphObj =
#   graphObj+geom_hline(yintercept = 0.002 * 29903, lty=3);
#
#graphObj =
#   graphObj +
#   annotate(
#      "label", # Adding an text with box annotation
#      x = 1,
#      y = 0.002 * 29903,
#      label = "0.2% difference"
#); # add an 0.2% line marker
#
#graphObj =
#   graphObj+geom_hline(yintercept = 0.001 * 29903,lty=2);
#graphObj =
#   graphObj +
#   annotate(
#      "label", # Adding an text with box annotation
#      x = 1,
#      y = 0.001 * 29903,
#      label = "0.1% difference"
#); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903,lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
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

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/lilo-snps-ref", asSvgBl);

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
            dataDF$program == "Lilo" &
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth),
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-04 Sub-02:
#  - Add in percent error markers
#**********************************************************

#graphObj =
#   graphObj+geom_hline(yintercept = 0.002 * 29903, lty=3);
#graphObj =
#   graphObj +
#   annotate(
#      "label", # Adding an text with box annotation
#      x = 750,
#      y = 0.002 * 29903,
#      label = "0.2% difference"
#); # add an 0.2% line marker
#
#graphObj =
#   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
#graphObj =
#   graphObj +
#   annotate(
#      "label", # Adding an text with box annotation
#      x = 750,
#      y = 0.001 * 29903,
#      label = "0.1% difference"
#); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.0005 * 29903,
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

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/lilo-snps-depth", asSvgBl);

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
            dataDF$program == "Lilo" &
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(indel) * (refLen / 1000),
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-05 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj +
   geom_vline(xintercept = 0.002 * 29903, lty = 3);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.002 * 29903,
      y = 0.0005 * 29903 + 1,
      label = "0.2% difference"
); # add an 0.2% line marker

graphObj =
   graphObj +
   geom_vline(xintercept = 0.001 * 29903, lty = 2);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.001 * 29903,
      y = 0.0005 * 29903 + 1,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj +
   geom_hline(yintercept=0.0005 * 29903, lty = 5);

graphObj =
   graphObj +
   geom_vline(xintercept=0.0005 * 29903, lty = 5);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.0005 * 29903 + 2,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-05 Sub-03:
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

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/lilo-snps-indels", asSvgBl);

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
            dataDF$program == "buildCon" &
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
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
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
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

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();
saveggplot("11-quast-bench/buildCon-snps-ref", asSvgBl);

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
            dataDF$program == "buildCon" &
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth),
         y = as.numeric(snp) * (refLen / 1000),
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
   graphObj + geom_hline(yintercept = 0.001 * 29903,lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj + geom_hline(yintercept=0.0005 * 29903,lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.0005 * 29903,
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

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/buildCon-snps-depth", asSvgBl);

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
            dataDF$program == "buildCon" &
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(indel) * (refLen / 1000),
         y = as.numeric(snp) * (refLen / 1000),
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
   geom_hline(yintercept = 0.002 * 29903, lty = 3);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.002 * 29903, lty = 3);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.002 * 29903,
      y = 0.002 * 29903,
      label = "0.2% difference"
); # add an 0.2% line marker

graphObj =
   graphObj +
   geom_hline(yintercept = 0.001 * 29903, lty = 2);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.001 * 29903, lty = 2);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.001 * 29903,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj +
   geom_hline(yintercept=0.0005 * 29903, lty = 5);

graphObj =
   graphObj +
   geom_vline(xintercept=0.0005 * 29903, lty = 5);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.0005 * 29903 + 2,
      y = 0.0005 * 29903,
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

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/buildCon-snps-indels", asSvgBl);

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
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
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
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
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
graphObj = graphObj + ylim(0,NA); # Start at 0

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/ivar-snps-ref", asSvgBl);

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
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth),
         y = as.numeric(snp) * (refLen / 1000),
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
   graphObj + geom_hline(yintercept=0.0005 * 29903,lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.0005 * 29903,
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
graphObj = graphObj + ylim(0,NA); # Start at 0

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/ivar-snps-depth", asSvgBl);

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
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(indel) * (refLen / 1000),
         y = as.numeric(snp) * (refLen / 1000),
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
   geom_hline(yintercept=0.0005 * 29903, lty = 5);

graphObj =
   graphObj +
   geom_vline(xintercept=0.0005 * 29903, lty = 5);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.0005 * 29903 + 2,
      y = 0.0005 * 29903,
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
      position = position_jitter(height = 0, width = 0.3)
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

graphObj = graphObj + xlim(0,NA); # Start at 0

graphObj = graphObj + ylab("number of false snps");
graphObj = graphObj + ylim(0,NA); # Start at 0

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/ivar-snps-indels", asSvgBl);

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
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
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
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
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

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/all-snps-ref", asSvgBl);

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
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth), # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
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
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.0005 * 29903,
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

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/all-snps-depth", asSvgBl);

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
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(indel) * (refLen / 1000),
         y = as.numeric(snp) * (refLen / 1000),
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
   geom_hline(yintercept = 0.002 * 29903, lty = 3);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.002 * 29903, lty = 3);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.002 * 29903,
      y = 0.002 * 29903,
      label = "0.2% difference"
); # add an 0.2% line marker

graphObj =
   graphObj +
   geom_hline(yintercept = 0.001 * 29903, lty = 2);

graphObj =
   graphObj +
   geom_vline(xintercept = 0.001 * 29903, lty = 2);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.001 * 29903,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj +
   geom_hline(yintercept=0.0005 * 29903, lty = 5);

graphObj =
   graphObj +
   geom_vline(xintercept=0.0005 * 29903, lty = 5);

graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 0.0005 * 29903 + 2,
      y = 0.0005 * 29903,
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

graphObj = graphObj + ylim(0,NA);

graphObj =
   graphObj +
   xlab("Number of false indels");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/all-snps-indels", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-15: 
#  - graph each pipelines number of masked bases
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
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            #dataDF$Ns < 29000 & # Removes a few artic runs
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth),
         y =
            as.numeric(Ns) * (as.numeric(refLen) / 1000) +
            ifelse(
               as.numeric(refLen) - as.numeric(conLen) > 0,
               as.numeric(refLen) - as.numeric(conLen),
               0
             ),
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

graphObj = graphObj + ylim(0,NA);

graphObj =
   graphObj +
   xlab("Depth");
graphObj =
   graphObj +
   ylab("Number of masked bases + missing bases");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/all-masked", asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-16: 
#  - graph each pipelines (over 30x depth) snps and indels 
#  o sec-16 sub-01:
#    - make the graph
#  o sec-16 sub-02:
#    - add in percent error markers
#  o sec-16 sub-03:
#    - add in the points
#  o sec-16 sub-04:
#    - add in the colors and shapes
#  o sec-16 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-16 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth > 30 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-16 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-16 Sub-03:
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
# Sec-16 Sub-04:
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
# Sec-16 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Reference % difference from sequence");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/all-over30x-snps-ref",asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-17: 
#  - graph each pipelines snps by ref (over 30x depth)
#  o sec-17 sub-01:
#    - make the graph
#  o sec-17 sub-02:
#    - add in percent error markers
#  o sec-17 sub-03:
#    - add in the points
#  o sec-17 sub-04:
#    - add in the colors and shapes
#  o sec-17 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-17 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth > 30 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(indel) * (refLen / 1000),
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-17 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-17 Sub-03:
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
# Sec-17 Sub-04:
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
# Sec-17 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Number of false indels");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
    "11-quast-bench/all-over30x-snps-indels",
    asSvgBl
);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-18: 
#  - graph each pipelines snps by depth (over 30x)
#  o sec-18 sub-01:
#    - make the graph
#  o sec-18 sub-02:
#    - add in percent error markers
#  o sec-18 sub-03:
#    - add in the points
#  o sec-18 sub-04:
#    - add in the colors and shapes
#  o sec-18 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-18 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth > 30 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth), # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-18 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-18 Sub-03:
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
# Sec-18 Sub-04:
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
# Sec-18 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Read depth");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
   "11-quast-bench/all-over30x-snps-depth",
   asSvgBl
);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-19: 
#  - graph each pipelines (50x depth) snps and indels 
#  o sec-19 sub-01:
#    - make the graph
#  o sec-19 sub-02:
#    - add in percent error markers
#  o sec-19 sub-03:
#    - add in the points
#  o sec-19 sub-04:
#    - add in the colors and shapes
#  o sec-19 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-19 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth == 50 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-19 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-19 Sub-03:
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
# Sec-19 Sub-04:
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
# Sec-19 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Reference % difference from sequence");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/all-50x-snps-ref",asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-20: 
#  - graph each pipelines snps by reference(50x read depth)
#  o sec-20 sub-01:
#    - make the graph
#  o sec-20 sub-02:
#    - add in percent error markers
#  o sec-20 sub-03:
#    - add in the points
#  o sec-20 sub-04:
#    - add in the colors and shapes
#  o sec-20 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-20 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth == 50 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(indel) * (refLen / 1000),
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-20 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-20 Sub-03:
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
# Sec-20 Sub-04:
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
# Sec-20 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Number of false indels");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
    "11-quast-bench/all-50x-snps-indels",
    asSvgBl
);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-21: 
#  - graph each pipelines (100x depth) snps and indels 
#  o sec-21 sub-01:
#    - make the graph
#  o sec-21 sub-02:
#    - add in percent error markers
#  o sec-21 sub-03:
#    - add in the points
#  o sec-21 sub-04:
#    - add in the colors and shapes
#  o sec-21 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-21 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth == 100 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-21 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-21 Sub-03:
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
# Sec-21 Sub-04:
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
# Sec-21 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Reference % difference from sequence");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/all-100x-snps-ref",asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-22: 
#  - graph each pipelines snps by ref (100x read depth)
#  o sec-22 sub-01:
#    - make the graph
#  o sec-22 sub-02:
#    - add in percent error markers
#  o sec-22 sub-03:
#    - add in the points
#  o sec-22 sub-04:
#    - add in the colors and shapes
#  o sec-22 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-22 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth == 100 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(indel) * (refLen / 1000),
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-22 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-22 Sub-03:
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
# Sec-22 Sub-04:
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
# Sec-22 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Number of false indels");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
    "11-quast-bench/all-100x-snps-indels",
    asSvgBl
);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-23: 
#  - graph each pipelines (300x depth) snps and indels 
#  o sec-23 sub-01:
#    - make the graph
#  o sec-23 sub-02:
#    - add in percent error markers
#  o sec-23 sub-03:
#    - add in the points
#  o sec-23 sub-04:
#    - add in the colors and shapes
#  o sec-23 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-23 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth == 300 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-23 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-23 Sub-03:
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
# Sec-23 Sub-04:
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
# Sec-23 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Reference % difference from sequence");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/all-300x-snps-ref",asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-24: 
#  - graph each pipelines snps by ref (300x read depth)
#  o sec-24 sub-01:
#    - make the graph
#  o sec-24 sub-02:
#    - add in percent error markers
#  o sec-24 sub-03:
#    - add in the points
#  o sec-24 sub-04:
#    - add in the colors and shapes
#  o sec-24 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-24 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth == 300 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(indel) * (refLen / 1000),
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-24 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-24 Sub-03:
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
# Sec-24 Sub-04:
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
# Sec-24 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Number of false indels");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
    "11-quast-bench/all-300x-snps-indels",
    asSvgBl
);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-25: 
#  - graph each pipelines (500x depth) snps and indels 
#  o sec-25 sub-01:
#    - make the graph
#  o sec-25 sub-02:
#    - add in percent error markers
#  o sec-25 sub-03:
#    - add in the points
#  o sec-25 sub-04:
#    - add in the colors and shapes
#  o sec-25 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-25 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth == 500 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-25 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-25 Sub-03:
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
# Sec-25 Sub-04:
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
# Sec-25 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Reference % difference from sequence");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/all-500x-snps-ref",asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-26: 
#  - graph each pipelines snps by ref (500x read depth)
#  o sec-26 sub-01:
#    - make the graph
#  o sec-26 sub-02:
#    - add in percent error markers
#  o sec-26 sub-03:
#    - add in the points
#  o sec-26 sub-04:
#    - add in the colors and shapes
#  o sec-26 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-26 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth == 500 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(indel) * (refLen / 1000),
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-26 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-26 Sub-03:
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
# Sec-26 Sub-04:
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
# Sec-26 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Number of false indels");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
    "11-quast-bench/all-500x-snps-indels",
    asSvgBl
);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-27: 
#  - graph each pipelines (1000x depth) snps and indels 
#  o sec-27 sub-01:
#    - make the graph
#  o sec-27 sub-02:
#    - add in percent error markers
#  o sec-27 sub-03:
#    - add in the points
#  o sec-27 sub-04:
#    - add in the colors and shapes
#  o sec-27 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-27 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth == 1000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-27 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-27 Sub-03:
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
# Sec-27 Sub-04:
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
# Sec-27 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Reference % difference from sequence");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot("11-quast-bench/all-1000x-snps-ref",asSvgBl);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-28: 
#  - graph each pipelines snps by ref (1000x read depth)
#  o sec-28 sub-01:
#    - make the graph
#  o sec-28 sub-02:
#    - add in percent error markers
#  o sec-28 sub-03:
#    - add in the points
#  o sec-28 sub-04:
#    - add in the colors and shapes
#  o sec-28 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-28 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth == 1000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(indel) * (refLen / 1000),
         y = as.numeric(snp) * (refLen / 1000),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-28 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj =
   graphObj+geom_hline(yintercept = 0.001 * 29903, lty=2);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.001 * 29903,
      label = "0.1% difference"
); # add an 0.1% line marker

graphObj =
   graphObj+geom_hline(yintercept=0.0005 * 29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 1,
      y = 0.0005 * 29903,
      label = "0.05% difference"
); # add an 0.05% line marker

#**********************************************************
# Sec-28 Sub-03:
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
# Sec-28 Sub-04:
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
# Sec-28 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Number of false indels");
graphObj = graphObj + ylab("Number of false snps");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
    "11-quast-bench/all-1000x-snps-indels",
    asSvgBl
);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-29: 
#  - graph each consensus length by depth (>30x, all)
#  o sec-29 sub-01:
#    - make the graph
#  o sec-29 sub-02:
#    - add in percent error markers
#  o sec-29 sub-03:
#    - add in the points
#  o sec-29 sub-04:
#    - add in the colors and shapes
#  o sec-29 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-29 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth > 30 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth), # 100 = 1%
         y = as.numeric(conLen),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-29 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj = graphObj + geom_hline(yintercept=29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 750,
      y = 29903,
      label = "Reference length"
); # add an 0.05% line marker

#**********************************************************
# Sec-29 Sub-03:
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
# Sec-29 Sub-04:
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
# Sec-29 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj =
   graphObj +
   xlab("Read depth");
graphObj = graphObj + ylab("Length of consensus");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
   "11-quast-bench/all-over30x-length-depth",
   asSvgBl
);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-30: 
#  - graph each consensus length by reference (>30x, all)
#  o sec-30 sub-01:
#    - make the graph
#  o sec-30 sub-02:
#    - add in percent error markers
#  o sec-30 sub-03:
#    - add in the points
#  o sec-30 sub-04:
#    - add in the colors and shapes
#  o sec-30 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-30 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth > 30 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(conLen),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-30 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj = graphObj + geom_hline(yintercept=29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 7.5,
      y = 29903,
      label = "Reference length"
); # add an 0.05% line marker

#**********************************************************
# Sec-30 Sub-03:
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
# Sec-30 Sub-04:
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
# Sec-30 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Reference % difference");
graphObj = graphObj + ylab("Length of consensus");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
   "11-quast-bench/all-over30x-length-ref",
   asSvgBl
);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-31: 
#  - graph the time usage for each method
#  o sec-31 sub-01:
#    - make the graph
#  o sec-31 sub-02:
#    - add in the points
#  o sec-31 sub-03:
#    - add in the colors and shapes
#  o sec-31 sub-04:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-31 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth > 30 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth),
         y = as.numeric(elpTime) / 60,
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-31 Sub-02:
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
# Sec-31 Sub-03:
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
# Sec-31 Sub-04:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Read depth");
graphObj = graphObj + ylab("Time in minutes");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
   "11-quast-bench/all-over30x-time",
   asSvgBl
);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-32: 
#  - graph the time usage for each method
#  o sec-32 sub-01:
#    - make the graph
#  o sec-32 sub-02:
#    - add in the points
#  o sec-32 sub-03:
#    - add in the colors and shapes
#  o sec-32 sub-04:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-32 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            ( dataDF$program == "ivar" |
              dataDF$program == "artic" |
              dataDF$test == "Lilo-medaka-ivarPolish-stich" |
              dataDF$test == "Lilo-medaka-ivarPolish-scaffold" |
              dataDF$test == "buildCon-majcon-ivarPolish-stich"
             ) &
            dataDF$refLen > 20000 &
            dataDF$depth > 30 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(depth),
         y = as.numeric(memKb) / 1000,
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-32 Sub-02:
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
# Sec-32 Sub-03:
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
# Sec-32 Sub-04:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Read depth");
graphObj = graphObj + ylab("Memory usage in megabytes");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
   "11-quast-bench/all-over30x-memory",
   asSvgBl
);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-33: 
#  - graph lilo by consensus length and reference
#  o sec-33 sub-01:
#    - make the graph
#  o sec-33 sub-02:
#    - add in percent error markers
#  o sec-33 sub-03:
#    - add in the points
#  o sec-33 sub-04:
#    - add in the colors and shapes
#  o sec-33 sub-05:
#    - add in labels, axis formating, and save graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-33 Sub-01:
#  - make the graph
#**********************************************************

graphObj = 
   ggplot(
      data =
         dataDF[
            dataDF$program == "Lilo" &
            dataDF$refLen > 20000 &
            dataDF$rep >= lastCompleteRepI
            ,
      ],
      aes(
         x = as.numeric(percMutate) / 100, # 100 = 1%
         y = as.numeric(conLen),
         color = as.character(test),
         fill = as.character(test),
         shape = as.character(test)
     ) # aes
); # make graph

#**********************************************************
# Sec-33 Sub-02:
#  - Add in percent error markers
#**********************************************************

graphObj = graphObj + geom_hline(yintercept=29903, lty=5);
graphObj =
   graphObj +
   annotate(
      "label", # Adding an text with box annotation
      x = 7.5,
      y = 29903,
      label = "Reference length"
); # add an 0.05% line marker

#**********************************************************
# Sec-33 Sub-03:
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
# Sec-33 Sub-04:
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
# Sec-33 Sub-05:
#  - Add in labels, axis formating, and save graph
#**********************************************************

graphObj = graphObj + xlab("Reference % difference");
graphObj = graphObj + ylab("Length of consensus");

# Make sure legend is readable
graphObj =
   graphObj +
   guides(
      color = guide_legend(nrow = 3),
      fill = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 3),
      linetype = guide_legend(nrow = 3),
);

graphObj = graphObj + theme_pubr();

saveggplot(
   "11-quast-bench/lilo-length-ref",
   asSvgBl
);

