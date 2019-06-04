#Essa funcao produz o primeiro grafico da aba plotoutput 
#que e o grafico dos dados em funcao das coordenadas na area experimental

chart1<-function(dados,resp,trat,coordX,coordY){
  coords<-cbind(dados[,coordX],dados[,coordY])
  geodados <- as.geodata(data.frame(coords,dados[,resp]),coords.col = 1:2,data.col = 3)
  #plot(geodados)
  layout(matrix(c(1,2),nrow=2),widths= c(1,1), heights=c(5.5,1))
  #layout.show(n = 2)
  #
  #Primeiro gráfico do geoR
  coords.lims <- set.coords.lims(coords = coords)
  plot(coords, xlab = "X Coord", ylab = "Y Coord ", type = "n", xlim = coords.lims[, 1], 
       ylim = coords.lims[, 2])
  data.lab <- "Data"
  qt.col <- c("blue", "green", "yellow2", "red")
  data.breaks <- unique(quantile(dados[,resp]))
  data.cut <- cut(dados[,resp], breaks = data.breaks, include.l = TRUE, 
                  labels = FALSE)
  points(coords, pch = 20, col = qt.col[data.cut])
  par(mar=c(0,0,0,0))
  plot(1,1,pch=NA, axes=F)
  legend(x="center",title = "Quartile",pch = rep(20,4),
         col = qt.col,legend = c(expression(Q[1]),expression(Q[2]),
           expression(Q[3]),expression(Q[4])),ncol=4)
  
}

chart2<-function(dados,resp,trat,coordX,coordY){
  coords<-cbind(dados[,coordX],dados[,coordY])
  coords.lims <- set.coords.lims(coords = coords)
  data<-dados[,resp]
  plot(data, coords[, 2], ylab = "Y Coord", xlab = "Data", 
       cex = 1, ylim = coords.lims[, 2])
}

chart3<-function(dados,resp,trat,coordX,coordY){
  coords<-cbind(dados[,coordX],dados[,coordY])
  coords.lims <- set.coords.lims(coords = coords)
  data<-dados[,resp]
  plot(coords[, 1], data, xlab = "X Coord", ylab = "Data",
       cex = 1, xlim = coords.lims[, 1])
}


# resp<-"H"
# trat<-"Tratamento"  
# Coord_X<-"Coord_X"
# Coord_Y<-"Coord_Y"
# chart1(dados,resp,trat,coordX=Coord_X,coordY=Coord_Y)