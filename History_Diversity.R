#Diferentes analisis de diversidad y rarefracción usando funciones del paquete Vegan


#Para Cargar el paquete
library(vegan) #se requiere esta libreria para bajar la base de datos BCI
data(BCI)
TablaH <- BCI[1:23,] #Solo escojo unos renglones (En este caso sitios)
nrow(TablaH)
UU <- rep("Veg",nrow(TablaH))
UUU <- seq(1:nrow(TablaH))
FactoresH <- paste(UU,UUU, sep=".")
FactoresH
#head(TablaH,3)
apply(TablaH,1,sum)
nrow(TablaH)


#Este seria el indicador de los ocho renglones de la TablaH


########################################
######## Diversidad Varios analisis
#Aqui analizamos varios parametros de diversidad: Shannon, Inversa de Simpson, equidad de Pielou y el numero total de especies
#Factor seria las variables alfanumericas que escogeriamos para hacer el analisis (e.g. vegetacion, sitios de muestreo)

#Variable "Tabla" debe tener el mismo numero de renglones que la variable "factor"
#Variable "factor" debe tener distinto nombre

DiversidadCC <- function(Tabla,factor){
		require(vegan)
		SpecNum <- specnumber(Tabla) ## #rowSums(BCI > 0)#  Species 		richness
		ShannonD <- diversity(Tabla)#Shannon entropy
  	Pielou <- ShannonD/log(SpecNum)#Pielou's evenness
		Simp <- diversity(Tabla, "simpson")#        Indice de dominacia de Simpson
		TablaF <- data.frame(factor,SpecNum, ShannonD, Simp, Pielou)
		print("Indicadors de Diversidad")
		print(TablaF)
}
DiversidadCC(TablaH,FactoresH)

#############################################
##############Para Rarefraccion
#Variable "Tabla" debe tener el mismo numero de renglones que la variable "factor"
#Variable "factor" debe tener distinto nombre


RarefraccionCC <- function(Tabla,factor){
	require(vegan)
	Tabla1 <- data.frame(Tabla, row.names = factor)
	raremax <- min(rowSums(Tabla1))
	col1 <- seq(1:nrow(Tabla1)) #Para poner color a las lineas
	lty1 <- c("solid","dashed","longdash","dotdash")
	rarecurve(Tabla1, sample = raremax, col = "black", lty = lty1, cex = 0.6)
  #Para calcular el numero de especies de acuerdo a rarefraccion
  UUU <- rarefy(Tabla1, raremax)
  print(UUU)
}


RarefraccionCC(TablaH,FactoresH)

#############################################
####Para Calcular Renyi#######
#Variable "Tabla" debe tener el mismo numero de renglones que la variable "factor"
#Variable "factor" debe tener distinto nombre

RenyiCC <- function(Tabla, factor){
		require(vegan) #Paquete para la funcion "renyi"
		require(ggplot2)#Paquete para hacer la funcion "qplot"
		require(reshape)#Paquete para la funcion "melt"
		Tabla <- data.frame(Tabla, row.names = factor)
		mod <- renyi(Tabla)
		vec <- seq(1:11)
		mod1 <- data.frame(vec,t(mod))
		mod2 <- melt(mod1, id = c("vec"))
		mod2
		#mod2$variable <- as.numeric(mod2$variable)
	orange <- qplot(vec, value, data = mod2, colour = variable, geom = 	"line") + theme_bw()
	orange + scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11), 		labels = c("0","0.25","0.5","1","2","4","8","16","32","64","Inf"))
}

RenyiCC(TablaH, FactoresH)

library(vegan)

data(BCI)
dim(BCI)
BCI1 <- BCI[1:7,]
vec1 <- paste(rep("site",7),seq(1:7), sep = ".")

DiversidadCC(BCI1,vec1)
RenyiCC(BCI1, vec1)

########################################################

# Calcular Diversidad Beta para presencia y ausencia o abundancia
# Tabla, debera de tener los valores de cada sitio
# variable factor debe tener distinto nombre
# el valor "n" es el tipo de metodo de analisis de la diversidad beta  
# el valor de "n" debe de ir del 1 al 24 de acuerdo al articulo de
# Koleff et al., 2003 o buscar el numero en: betadiver(help = TRUE)

#Alex quiza hacer una opcion para que si el usuario quiere abundancia
# o los datos en presencia y ausencia?

#Para presencia y ausencia
DivBetaPA <- function(Tabla, Factor, n1){
  require(vegan) #Paquete para la funcion "betadiver"
  #Diversidad Beta de acuerdo a Koleff et al., con presencia y ausencia
  Tabla <- data.frame(Tabla, row.names = Factor)
  Tabla
  Tabla <- decostand(Tabla, "pa")
  mod <- betadiver(Tabla, method = n1, binary = F)
  mod
  hca <- hclust(mod, method = "ward.D")
  plot(as.dendrogram(hca), horiz = T, main = "Diversidad Beta con Presencia-Ausencia")
  print(mod)
}
ls()


#Para abundancia
DivBetaAbun <- function(Tabla, Factor, n1){
  require(vegan) #Paquete para la funcion "vegdist"
  Tabla <- data.frame(Tabla, row.names = Factor)
  Tabla
  
  #Los distintos tipos de analisis
  if (n1 == 1) {
    n2 <- "manhattan"
  } else if (n1 == 2) {
    n2 <- "euclidean"
  } else if (n1 == 3) {
    n2 <- "canberra"
  } else if (n1 == 4) {
    n2 <- "bray"
  } else if (n1 == 5) {
    n2 <- "kulczynski"
  } else if (n1 == 6) {
    n2 <- "jaccard"
  } else if (n1 == 7) {
    n2 <- "gower"
  } else if (n1 == 8) {
    n2 <- "altGower"
  } else if (n1 == 9) {
    n2 <- "morosita"
  } else if (n1 == 10) {
    n2 <- "horn"
  } else if (n1 == 11) {
    n2 <- "mountford"
  }else if (n1 == 12) {
    n2 <- "raup"
  } else if (n1 == 13) {
    n2 <- "binomial"
  } else if (n1 == 14) {
    n2 <- "chao"
  } else if (n1 == 15) {
    n2 <- "cao"
  } else if (n1 == 16) {
    n2 <- "mahalanobis"
  }
  
  mod <- vegdist(Tabla, method = n2)
  mod
  hca <- hclust(mod, method = "ward.D")
  plot(as.dendrogram(hca), horiz = T, main = "Diversidad Beta con Abundancia")
  print(mod)
}

ls()
par(mfrow = c(1,1))
DivBetaPA(TablaH, FactoresH, 9)
DivBetaAbun(TablaH, FactoresH, 4)




########################################
#Chao data
#Many species will always remain unseen or undetected in a 
#collection of sample plots. 
#The function uses some popular ways of estimating 
#the number of these unseen species and adding them 
#to the observed species richness (Palmer 1990, Colwell & Coddington 1994).
#########################################
ChaoF <- function(Tabla, Factor){
  require(vegan)
  require(ggplot2)
  require(reshape)
  Tabla <- data.frame(Tabla, row.names = Factor)
  Tabla
  UNO <- estimateR(Tabla)[1:2,]
  UNO1 <- melt(UNO)
  names(UNO1)[1] <- c("Var")
  names(UNO1)[2] <- c("Veg")
  LL <- ggplot(UNO1, aes(x = Veg, y = value, fill = Var)) 
  LL <- LL + geom_bar(position = "stack", stat = "identity") 
  LL <- LL + facet_wrap( ~ Var) + theme_bw() + coord_flip()
  print(LL)
  print(UNO)
}
ls()
ChaoF(TablaH, FactoresH)


##############################################
#Curvas de Whitaker

#legend("bottomright","foo")

Whitaker1 <- function(Tabla, Factor){
  require(vegan)
  Tabla <- data.frame(Tabla, row.names = Factor)
  Tabla
  mod <- radfit(Tabla)
  mod
  plot(mod, pch = 19)
  }
ls()
Whitaker1(TablaH, FactoresH)

