WEB6004 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/WEB6004.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5998 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/WEB5998.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB6006 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/WEB6006.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB6007 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/WEB6007.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5997 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/WEB5997.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5992 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/WEB5992.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5999 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/WEB5999.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5995 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/WEB5995.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB6004 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/WEB6004.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5993 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/WEB5993.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5994 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/WEB5994.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

t <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/LOY7-23-10.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

e <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/LOY7-23-18.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

n <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH5gcov/LOY7-23-19.gcov', 
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
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), 1631057, 1652811, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("UTI89 reads aligned against MG1655 genome (UH5 region)") +
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
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), 1631057, 1652811, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("MG1655 reads aligned against MG1655 genome (UH5 region)") +
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
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), 1631057, 1652811, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("WEB6004 reads aligned against MG1655 genome (UH5 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

uh5a <- makeUHdf(WEB5998, 1631057, 1652811)
uh5b <- makeUHdf(WEB5997, 1631057, 1652811)
uh5c <- makeUHdf(WEB5992, 1631057, 1652811)
uh5d <- makeUHdf(WEB6006, 1631057, 1652811)
uh5e <- makeUHdf(WEB6007, 1631057, 1652811)
uh5f <- makeUHdf(WEB5999, 1631057, 1652811)
uh5g <- makeUHdf(WEB5995, 1631057, 1652811)
uh5h <- makeUHdf(WEB5994, 1631057, 1652811)
uh5i <- makeUHdf(WEB6004, 1631057, 1652811)
uh5j <- makeUHdf(WEB5993, 1631057, 1652811)

#plotting multiple graphs simultaneously 
out5 <- ggplot() +
  geom_line(data = uh5a, aes(x = x0, y = y0, color = "WEB5998"), alpha = 0.8) +
  geom_line(data = uh5b, aes(x = x0, y = y0, color = "WEB5997"), alpha = 0.8) +
  geom_line(data = uh5c, aes(x = x0, y = y0, color = "WEB5992"), alpha = 0.8) +
  geom_line(data = uh5d, aes(x = x0, y = y0, color = "WEB6006"), alpha = 0.8) +
  geom_line(data = uh5e, aes(x = x0, y = y0, color = "WEB6007"), alpha = 0.8) +
  geom_line(data = uh5f, aes(x = x0, y = y0, color = "WEB5999"), alpha = 0.8) +
  geom_line(data = uh5g, aes(x = x0, y = y0, color = "WEB5995"), alpha = 0.8) +
  geom_line(data = uh5h, aes(x = x0, y = y0, color = "WEB5994"), alpha = 0.8) +
  geom_line(data = uh5i, aes(x = x0, y = y0, color = "WEB6004"), alpha = 0.8) +
  geom_line(data = uh5j, aes(x = x0, y = y0, color = "WEB5993"), alpha = 0.8) +
  scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(uh5a$x0),1631057, 1652811, max(uh5a$x0))) +
  ylab("Sequencing Depth") +
  ggtitle("UH5 Clones Sequencing Depth") +
  theme_light() +
  theme(
    axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
show(out5)

uti5 <- utiplotter(makeUHdf(uti, 1631057, 1652811))
mg5 <- mgplotter(makeUHdf(mg, 1631057, 1652811))
uh5a <- plotter(makeUHdf(n, 1631057, 1652811))

out <- ggarrange(uti5, mg5, uh5a,
                 nrow = 3, ncol = 1)

out5 <- ggarrange(uti2 +rremove("xlab"), uti4 +rremove("ylab") + rremove("xlab"), uti5 +rremove("ylab") +rremove("xlab"),
                  mg2 + rremove("xlab"), mg4 +rremove("ylab") + rremove("xlab"), mg5 +rremove("ylab") + rremove("xlab"),
                  loy7, uh4a + rremove("ylab"), uh5a + rremove("ylab"),
                 nrow = 3, ncol = 3)
show(out5)
