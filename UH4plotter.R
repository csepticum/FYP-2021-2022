WEB5984 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH4gcov/WEB5984.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5985 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH4gcov/WEB5985.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5986 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH4gcov/WEB5986.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5990 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH4gcov/WEB5990.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5991 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH4gcov/WEB5991.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))



#given the gcov alignment, takes out the coordinates i want and return as a dataframe ready for ggplotting
makeUHdf <- function(df, star, sto){
  star <- star - 20000
  sto <- sto + 20000
  UHout <- df[star:sto,]
  return(UHout)
}

#plotting single graphs 
utiplotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "tomato3") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), 1412030, 1435238, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("UTI89 reads aligned against MG1655 genome (UH4 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

mgplotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "skyblue3") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), 1412030, 1435238, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("MG1655 reads aligned against MG1655 genome (UH4 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

plotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "grey38") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), 1412030, 1435238, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("WEB5984 reads aligned against MG1655 genome (UH4 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

uh4a <- makeUHdf(WEB5984, 1412030, 1435238)
uh4b <- makeUHdf(WEB5985, 1412030, 1435238)
uh4c <- makeUHdf(WEB5986, 1412030, 1435238)
uh4d <- makeUHdf(WEB5990, 1412030, 1435238)
uh4e <- makeUHdf(WEB5991, 1412030, 1435238)

#plotting multiple graphs simultaneously 
out4 <- ggplot() +
  geom_line(data = uh4a, aes(x = x0, y = y0, color = "WEB5984"), alpha = 0.8) +
  geom_line(data = uh4b, aes(x = x0, y = y0, color = "WEB5985"), alpha = 0.8) +
  geom_line(data = uh4c, aes(x = x0, y = y0, color = "WEB5986"), alpha = 0.8) +
  geom_line(data = uh4d, aes(x = x0, y = y0, color = "WEB5990"), alpha = 0.8) +
  geom_line(data = uh4e, aes(x = x0, y = y0, color = "WEB5991"), alpha = 0.8) +
  scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(uh4a$x0),1412030, 1435238, max(uh4a$x0))) +
  ylab("Sequencing Depth") +
  ggtitle("UH4 Clones Sequencing Depth") +
  theme_light() +
  theme(
    axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
show(out4)



uti4 <- utiplotter(makeUHdf(uti, 1412030, 1435238))
mg4 <- mgplotter(makeUHdf(mg, 1412030,1435238))
uh4a <- plotter(makeUHdf(WEB5984, 1412030, 1435238))


out <- ggarrange(uti4, uti5,
                 mg4, mg5,
                 uh4a, uh5a,
                 nrow = 3, ncol = 2)
show(out)
