

#tidying up the call function for getting MGtoUTI blast data 
MGhomology <- function(){
  MGtoUTI <- read.table('/Users/sylvester/Dropbox/NUS Mod Stuff/FYP/Programming @ GIS/SLComp1-24_MAE Visualisation/MG1655-against-UTI89.blastn', 
                        header = FALSE,  
                        sep = "\t", 
                        as.is = TRUE, 
                        comment.char = "#", 
                        quote = "", 
                        col.names = c("queryseqid", "subjectseqid", "pIdentity", "length", "mismatch", "gapopen", "querystart", "queryend",
                                      "subjectstart", "subjectend", "eval", "bitscore"))
  MGtoUTI <-  MGtoUTI[order(MGtoUTI$querystart),] #sort in ascending order
  return(MGtoUTI)
}

# function to filter the MGtoUTI homology data, returns a filtered dataframe 
filtercoords <- function(start, stop){
  output <- filter(MGtoUTI, (MGtoUTI$querystart > start & MGtoUTI$queryend < stop))
  #output <- filter(output, output$length > 1000 & output$pIdentity>80)
  return(output)
}

# function for removing overlapping homologous sequences, returns a dataframe without overlaps
remove_overlap <- function(interestdf){
  interestdf <-  interestdf[order(interestdf$querystart),] #sort in ascending order
  for(i in c(1:nrow(interestdf))){
    if(i == nrow(interestdf)){
      break
    }
    else if(interestdf[i,8] > interestdf[i+1,7]){
      interestdf[i,1] = "cull"
    }
  }
  interestdf <- filter(interestdf, (interestdf$queryseqid != "cull"))
  return(interestdf)
}

#function for adding the mg1655pos and uti89pos, then adding a mapdiff that accounts for cross mappers, returns a dataframe with mg1655pos, uti89pos and mapdiff
mappingf <- function(interestdf){
  mg1655pos <- (((interestdf$queryend - interestdf$querystart)/2) + interestdf$querystart)/4641652
  uti89pos <- (((interestdf$subjectend - interestdf$subjectstart)/2) + interestdf$subjectstart)/5065742
  interestdf <- cbind(interestdf, mg1655pos)
  interestdf <- cbind(interestdf, uti89pos)
  mapdiff <- abs(interestdf$uti89pos - interestdf$mg1655pos)
  interestdf <- cbind(interestdf, mapdiff)
  # calculate crossover mapping mapdiff properly
  #filtering for distances greater than an arbitrary threshold
  interestdf1 <- filter(interestdf, (interestdf$mapdiff <= 0.5)) #filtering for non-crossed mapping
  interestdf <- filter(interestdf, (interestdf$mapdiff > 0.5)) #filtering for crossed mapping 
  interestdf2 <- filter(interestdf, (interestdf$mg1655pos <= 0.5)) #filtering for calculating the distance amongst crossed mappers
  interestdf3 <- filter(interestdf, (interestdf$uti89pos <= 0.5)) #filtering for calculating the distance amongst crossed mappers
  temp2 <- c(interestdf2$mg1655pos + (1-interestdf2$uti89pos))
  temp3 <- c(interestdf3$uti89pos + (1-interestdf3$mg1655pos))
  interestdf2$mapdiff <- temp2
  interestdf3$mapdiff <- temp3
  interestdf_cross <- rbind(interestdf2, interestdf3)
  total <- rbind(interestdf1, interestdf_cross)
  total <- filter(total, total$pIdentity >= 80 & total$length >= 1000)
  return(total)
}

mappinguh <- function(interestdf){
  mg1655pos <- (((interestdf$queryend - interestdf$querystart)/2) + interestdf$querystart)/4641652
  uti89pos <- (((interestdf$subjectend - interestdf$subjectstart)/2) + interestdf$subjectstart)/5065742
  interestdf <- cbind(interestdf, mg1655pos)
  interestdf <- cbind(interestdf, uti89pos)
  mapdiff <- abs(interestdf$uti89pos - interestdf$mg1655pos)
  interestdf <- cbind(interestdf, mapdiff)
  # calculate crossover mapping mapdiff properly
  #filtering for distances greater than an arbitrary threshold
  interestdf1 <- filter(interestdf, (interestdf$mapdiff <= 0.5)) #filtering for non-crossed mapping
  interestdf <- filter(interestdf, (interestdf$mapdiff > 0.5)) #filtering for crossed mapping 
  interestdf2 <- filter(interestdf, (interestdf$mg1655pos <= 0.5)) #filtering for calculating the distance amongst crossed mappers
  interestdf3 <- filter(interestdf, (interestdf$uti89pos <= 0.5)) #filtering for calculating the distance amongst crossed mappers
  temp2 <- c(interestdf2$mg1655pos + (1-interestdf2$uti89pos))
  temp3 <- c(interestdf3$uti89pos + (1-interestdf3$mg1655pos))
  interestdf2$mapdiff <- temp2
  interestdf3$mapdiff <- temp3
  interestdf_cross <- rbind(interestdf2, interestdf3)
  total <- rbind(interestdf1, interestdf_cross)
  return(total)
}


#function for plotting mapping
map_plot <- function(interestdf){
  mapout <- ggplot(data = interestdf) +
    geom_segment(aes(x = mg1655pos, y=0, xend = uti89pos, yend = 1)) +
    ggtitle("Homology Mapping for Hybridized Regions") +
    scale_x_continuous(name = "Relative MG1655 Genomic Coordinate", breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), limits = c(0,1), sec.axis = sec_axis(trans = ~., name = "Relative UTI89 Genomic Coordinate")) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          text = element_text(size = 12),
          axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) +
    ylim(0,1)
  show(mapout)
  return(mapout)
}

mgcds <- read.table('/Users/sylvester/Desktop/GenomeVizTest/mgcds.txt', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("location", "strand", "length", "PID", "gene", "synonym", "code", "COG",
                                    "product"))

uticds <- read.table('/Users/sylvester/Desktop/GenomeVizTest/uticds.txt', 
                    header = FALSE,  
                    sep = "\t", 
                    as.is = TRUE, 
                    comment.char = "#", 
                    quote = "", 
                    col.names = c("location", "strand", "length", "PID", "gene", "synonym", "code", "COG",
                                  "product"))

# syntenic mapping only for hybridized sites  

MGtoUTI <- MGhomology()
MGtoUTI <- MGtoUTI[order(MGtoUTI$querystart),]


## violin plot for all the hybridized region, this needs working on because there are mapdiffs > 0.5. 

UH1cf <- filter(MGtoUTI, (MGtoUTI$querystart > 392 & MGtoUTI$queryend <= 566332))
UH2cf <- filter(MGtoUTI, (MGtoUTI$querystart >= 582018 & MGtoUTI$queryend <= 1173242))
UH3cf <- filter(MGtoUTI, (MGtoUTI$querystart >= 1319462 & MGtoUTI$queryend <= 1412029))
UH4cf <- filter(MGtoUTI, (MGtoUTI$querystart >= 1435239 & MGtoUTI$queryend <= 1631056))
UH5cf <- filter(MGtoUTI, (MGtoUTI$querystart >= 1652812 & MGtoUTI$queryend <= 4500714))
UH6cf <- filter(MGtoUTI, (MGtoUTI$querystart >= 4561492 & MGtoUTI$queryend <= 4641440))
UH7cf <- filter(MGtoUTI, (MGtoUTI$querystart >= 4641441))
all <- rbind(UH1cf, UH2cf, UH3cf, UH4cf, UH5cf, UH6cf, UH7cf)
all <- filter(all, (all$length > 1000 & all$pIdentity >= 80))
all <- mappingf(all)
all$mapdiff <- round(all$mapdiff, digits = 2) #rounding off to 2 significant figures 

box <- ggplot(all, aes(x = "", y = mapdiff)) +
  geom_violin(trim = FALSE, fill = "lightsteelblue3") +
  geom_boxplot(width = 0.07) +
  geom_jitter(position = position_jitter(0.1))+
  theme_bw() +
  ggtitle("Mapping Distance of Homologous Sites") +
  ylab("Mapping Distance (Relative Units)") +
  xlab("MG1655 : UTI89 Homologs") +
  ylim(0, 0.5) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 10, hjust=1))
show(box)

## redoing the syntenic mapping
UH2heatmap <- function(){
  hybstart <- 566333
  hybend <- 582018
  gap <- abs(hybend-hybstart)
  df <- filter(MGtoUTI, MGtoUTI$length > 500, MGtoUTI$pIdentity > 80)
  mguuh2 <- ggplot(df) +
    geom_rect(aes(xmin = querystart, xmax = queryend, ymin = 0, ymax = 1, fill = pIdentity), alpha=0.8) +
    theme_bw() +
    scale_x_continuous("MG1655 coordinates", breaks = c(hybstart - gap, hybstart, hybend, hybend + gap), limits = c(hybstart-gap, hybend + gap)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), color = c("black", "maroon3", "maroon3", "black")))+
    ylim(0,1) +
    scale_fill_gradient(low ="black", high="red") +
    ggtitle("Homology Hit for UH2") +
    xlab("MG1655 coordinates")
  return(mguuh2)
}
UH3heatmap <- function(){
  hybstart <- 1173243
  hybend <- 1435238
  gap <- abs(hybend-hybstart)
  df <- filter(MGtoUTI, MGtoUTI$length > 1000, MGtoUTI$pIdentity > 80)
  mguuh3 <- ggplot(df) +
    geom_rect(aes(xmin = querystart, xmax = queryend, ymin = 0, ymax = 1, fill = pIdentity), alpha=0.8) +
    theme_bw() +
    scale_x_continuous("MG1655 coordinates", breaks = c(hybstart - gap, hybstart, hybend, hybend + gap), limits = c(hybstart-gap, hybend + gap)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), color = c("black", "maroon3", "maroon3", "black")))+
    ylim(0,1) +
    scale_fill_gradient(low ="black", high="red") +
    ggtitle("Homology Hit for UH3") +
    xlab("MG1655 coordinates")
  return(mguuh3)
}
UH4heatmap <- function(){
  hybstart <- 1412030
  hybend <- 1435238
  gap <- abs(hybend-hybstart)
  df <- filter(MGtoUTI, MGtoUTI$length > 1000, MGtoUTI$pIdentity > 80)
  mguuh2 <- ggplot(df) +
    geom_rect(aes(xmin = querystart, xmax = queryend, ymin = 0, ymax = 1, fill = pIdentity), alpha=0.8) +
    theme_bw() +
    scale_x_continuous("MG1655 coordinates", breaks = c(hybstart - gap, hybstart, hybend, hybend + gap), limits = c(hybstart-gap, hybend + gap)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), color = c("black", "maroon3", "maroon3", "black")))+
    ylim(0,1) +
    scale_fill_gradient(low ="black", high="red") +
    ggtitle("Homology Hit for UH4") +
    xlab("MG1655 coordinates")
  return(mguuh2)
}
UH5heatmap <- function(){
  hybstart <- 1631057
  hybend <- 1652811
  gap <- abs(hybend-hybstart)
  df <- filter(MGtoUTI, MGtoUTI$length > 1000, MGtoUTI$pIdentity > 80)
  mguuh2 <- ggplot(df) +
    geom_rect(aes(xmin = querystart, xmax = queryend, ymin = 0, ymax = 1, fill = pIdentity), alpha=0.8) +
    theme_bw() +
    scale_x_continuous("MG1655 coordinates", breaks = c(hybstart - gap, hybstart, hybend, hybend + gap), limits = c(hybstart-gap, hybend + gap)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), color = c("black", "maroon3", "maroon3", "black")))+
    ylim(0,1) +
    scale_fill_gradient(low ="black", high="red") +
    ggtitle("Homology Hit for UH5") +
    xlab("MG1655 coordinates")
  return(mguuh2)
}
UH6heatmap <- function(){
  hybstart <- 4500715
  hybend <- 4561491
  gap <- abs(hybend-hybstart)
  df <- filter(MGtoUTI, MGtoUTI$length > 1000, MGtoUTI$pIdentity > 80)
  mguuh2 <- ggplot(df) +
    geom_rect(aes(xmin = querystart, xmax = queryend, ymin = 0, ymax = 1, fill = pIdentity), alpha=0.8) +
    theme_bw() +
    scale_x_continuous("MG1655 coordinates", breaks = c(hybstart - 80161, hybstart, hybend, 4641652), limits = c(hybstart-80161, 4641652)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), color = c("black", "maroon3", "maroon3", "black")))+
    ylim(0,1) +
    scale_fill_gradient(low ="black", high="red") +
    ggtitle("Homology Hit for UH6") +
    xlab("MG1655 coordinates")
  return(mguuh2)
}

UH2mapgraph <- function(){
  uh2f1 <- filtercoords(534963,566332)
  uh2f1 <- remove_overlap(uh2f1)
  uh2f1 <- mappingf(uh2f1)
  uh2 <- filtercoords(566333, 582018)
  #uh2 <- remove_overlap(uh2)
  uh2 <- mappinguh(uh2)
  uh2f2 <- filtercoords(582019, 613388)
  uh2f2 <- remove_overlap(uh2f2)
  uh2f2 <- mappingf(uh2f2)
  uh2flanks <- rbind(uh2f1, uh2f2)
  uh2 <- filter(uh2, uh2$length > 150)
  
  UH2 <- ggplot() +
    geom_segment(data = uh2flanks, aes(x = mg1655pos, y=0, xend = uti89pos, yend=1, color = length)) + #for charting within the flanking sites
    geom_segment(data = uh2, aes(x = mg1655pos, y=0, xend = uti89pos, yend=1)) + #for charting within the homologous site
    scale_color_gradient(low ="yellow", high="red") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Homology Mapping for UH2") +
    scale_x_continuous(name = "Relative MG1655 Genomic Coordinate", breaks = c(0,0.25,0.5,0.75, 1), limits = c(0,1), sec.axis = sec_axis(trans = ~., name = "Relative UTI89 Genomic Coordinate")) +
    ylim(0,1)
  return(UH2)
}
UH3mapgraph <- function(){
  uh3f1 <- filtercoords(1173243,1435238)
  uh3f1 <- remove_overlap(uh3f1)
  uh3f1 <- mappingf(uh3f1)
  uh3 <- filtercoords(1173243,1435238)
  #uh3 <- remove_overlap(uh3)
  uh3 <- mappinguh(uh3)
  uh3f2 <- filtercoords(1173243,1435238)
  uh3f2 <- remove_overlap(uh3f2)
  uh3f2 <- mappingf(uh3f2)
  uh3flanks <- rbind(uh3f1, uh3f2)
  uh3 <- filter(uh3, uh3$length > 150)
  UH3 <- ggplot() +
    geom_segment(data = uh3flanks, aes(x = mg1655pos, y=0, xend = uti89pos, yend=1, color = length)) + #for charting within the flanking sites
    geom_segment(data = uh3, aes(x = mg1655pos, y=0, xend = uti89pos, yend=1)) + #for charting within the homologous site
    scale_color_gradient(low ="yellow", high="red") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Homology Mapping for UH3") +
    scale_x_continuous(name = "Relative MG1655 Genomic Coordinate", breaks = c(0,0.25,0.5,0.75, 1), limits = c(0,1), sec.axis = sec_axis(trans = ~., name = "Relative UTI89 Genomic Coordinate")) +
    ylim(0,1)
  return(UH3)
}
UH4mapgraph <- function(){
  uh4f1 <- filtercoords(1365614,1412029)
  uh4f1 <- remove_overlap(uh4f1)
  uh4f1 <- mappingf(uh4f1)
  uh4 <- filtercoords(1412030, 1435238)
  #uh4 <- remove_overlap(uh4)
  uh4 <- mappinguh(uh4)
  uh4f2 <- filtercoords(1435239, 1481654)
  uh4f2 <- remove_overlap(uh4f2)
  uh4f2 <- mappingf(uh4f2)
  uh4flanks <- rbind(uh4f1, uh4f2)
  uh4 <- filter(uh4, uh4$length > 150)
  
  UH4 <- ggplot() +
    geom_segment(data = uh4flanks, aes(x = mg1655pos, y=0, xend = uti89pos, yend=1, color = length)) + #for charting within the flanking sites
    geom_segment(data = uh4, aes(x = mg1655pos, y=0, xend = uti89pos, yend=1)) + #for charting within the homologous site
    scale_color_gradient(low ="yellow", high="red") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Homology Mapping for UH4") +
    scale_x_continuous(name = "Relative MG1655 Genomic Coordinate", breaks = c(0,0.25,0.5,0.75, 1), limits = c(0,1), sec.axis = sec_axis(trans = ~., name = "Relative UTI89 Genomic Coordinate")) +
    ylim(0,1)
  return(UH4)
}
UH5mapgraph <- function(){
  uh5f1 <- filtercoords(1587549,1631056)
  uh5f1 <- remove_overlap(uh5f1)
  uh5f1 <- mappingf(uh5f1)
  uh5 <- filtercoords(1631057,1652811)
  #uh5 <- remove_overlap(uh5)
  uh5 <- mappinguh(uh5)
  uh5f2 <- filtercoords(1652812, 1696319)
  uh5f2 <- remove_overlap(uh5f2)
  uh5f2 <- mappingf(uh5f2)
  uh5flanks <- rbind(uh5f1, uh5f2)
  uh5flanks <- filter(uh5flanks, uh5flanks$pIdentity > 81)
  uh5 <- filter(uh5, uh5$length > 150)
  
  UH5 <- ggplot() +
    geom_segment(data = uh5flanks, aes(x = mg1655pos, y=0, xend = uti89pos, yend=1, color = length)) + #for charting within the flanking sites
    geom_segment(data = uh5, aes(x = mg1655pos, y=0, xend = uti89pos, yend=1)) + #for charting within the homologous site
    scale_color_gradient(low ="yellow", high="red") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Homology Mapping for UH5") +
    scale_x_continuous(name = "Relative MG1655 Genomic Coordinate", breaks = c(0,0.25,0.5,0.75, 1), limits = c(0,1), sec.axis = sec_axis(trans = ~., name = "Relative UTI89 Genomic Coordinate")) +
    ylim(0,1)
  return(UH5)
}
UH6mapgraph <- function(){
  uh6f1 <- filtercoords(4420553,4500714)
  uh6f1 <- remove_overlap(uh6f1)
  uh6f1 <- mappingf(uh6f1)
  uh6 <- filtercoords(4500715,4561491)
  #uh6 <- remove_overlap(uh6)
  uh6 <- mappinguh(uh6)
  uh6f2 <- filtercoords(4561492, 4641652)
  uh6f2 <- remove_overlap(uh6f2)
  uh6f2 <- mappingf(uh6f2)
  uh6flanks <- rbind(uh6f1, uh6f2)
  uh6flanks <- filter(uh6flanks, uh6flanks$length> 1000)
  uh6 <- filter(uh6, uh6$length > 150)
  UH6 <- ggplot() +
    geom_segment(data = uh6flanks, aes(x = mg1655pos, y=0, xend = uti89pos, yend=1, color = length)) + #for charting within the flanking sites
    geom_segment(data = uh6, aes(x = mg1655pos, y=0, xend = uti89pos, yend=1)) + #for charting within the homologous site
    scale_color_gradient(low ="yellow", high="red") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Homology Mapping for UH6") +
    scale_x_continuous(name = "Relative MG1655 Genomic Coordinate", breaks = c(0,0.25,0.5,0.75, 1), limits = c(0,1), sec.axis = sec_axis(trans = ~., name = "Relative UTI89 Genomic Coordinate")) +
    ylim(0,1)
  return(UH6)
}

UH2heat <- UH2heatmap()
UH3heat <- UH3heatmap()
UH4heat <- UH4heatmap()
UH5heat <- UH5heatmap()
UH6heat <- UH6heatmap()

UH2map <- UH2mapgraph()
UH3map <- UH3mapgraph()
UH4map <- UH4mapgraph()
UH5map <- UH5mapgraph()
UH6map <- UH6mapgraph()

outheat <- ggarrange(UH2heat, UH3heat,UH4heat, 
                 UH5heat, UH6heat,
                 nrow = 5, ncol = 1,
                 common.legend = TRUE,
                 legend = "right")
outmap <- ggarrange(UH2map, UH3map, UH4map,
                    UH5map, UH6map,
                    nrow = 5, ncol = 1,
                    common.legend = TRUE,
                    legend = "right")
out <- ggarrange(outheat, outmap,
                 nrow = 1, ncol = 2)
show(out)