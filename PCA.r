#---------------------------------------------
# Lectura de la base de datos "FQmarino"
datos <-read.csv2("FQmarino.csv",row.names=1)

#---------------------------
# LIBRERÍAS REQUERIDAS
library(lattice)
library(ellipse)
require(SciViews)
library(FactoMineR)

# 1. Analisis numérico
acp <-princomp(datos,cor=T)

# Información sobre el análisis
help(princomp)

# Resultados del PCA
summary(acp)

# Insumos
names(acp)
# Pesos de las variables (autovectores)
acp$loadings
#
acp$scores

# 2. Figura de bahías
with(acp,plot(acp$scores,pch=8,col=4,xlab="CP1", ylab="CP2"),xlim(c(-4,4)))

# Colocar los rotulos de las bahías (observaciones) 
with(acp,text(acp$scores,row.names(acp$scores),pos=3))
abline(v=0,lty=2, col=2, lwd=2)
abline(h=0,lty=2, col=2, lwd=2)


# 3. Biplot para bahías con sus variables ambientales
biplot(acp,cex=0.9,xlab="CP1", ylab="CP2")
abline(v=0,lty=2, col=4)
abline(h=0,lty=2, col=4)


#-------------------------------------------
# MODELO LINEAL A PARTIR DEL PCA
datos.PCA<-datos[, c("pH", "Cond", "Turbidez", "Temp", "Salinidad",
                     "CapaFotica", "Oxigeno")]
#               
pca<-PCA(datos.PCA , scale.unit=TRUE, ncp=5, graph = FALSE)
#
plot.PCA(pca, axes=c(1, 2), choix="ind")
#                                                                     
plot.PCA(pca, axes=c(1, 2), choix="var", col.var="#ff0000", new.plot=T,
         col.quanti.sup="blue", label=c("var", "quanti.sup"), lim.cos2.var=0)
#
pca$eig
#
pca$var

#
dimdesc(pca)
#
remove(datos.PCA)

#------------------------
# PCA de datos ambientales con escalamiento de las variables(estandarización)
# Se usa la función "rda" del análisis de redundancia
pca1 <- rda(datos, scale=TRUE)

summary(pca1)
# Varianza de los dos primeros ejes= 48%, se requieren 4 ejes para 70%

# Comparación gráfica de los dos Scaling
windows(title="Biplot de datos ambientales", 12, 6) par(mfrow=c(1,2))
par(mfrow=c(1,2))
biplot(pca1, scaling=1, main="PCA - scaling 1")
biplot(pca1, main="PCA - scaling 2")  
# Scaling 1 =  Se observa gradiente izquierda-derecha en términos de contaminación
# ej. grupo 1= mayor altitud, grupo 2= mayor nivel de oxigeno,..., grupo4= mayor cont.
# Scalig 2 = muestra mayor asociación de variables con los grupos

#---------------
# Dos métodos para seleccionar el número de ejes del PCA
# Autovalores del pca (av)
(av <- pca1$CA$eig)

# Criterio Kaiser-Guttman para seleccionar el número de ejes
# Calcula el promedio de los autovalores y escoge los mayores a la media
av[av > mean(av)]
# Define tres ejes

# Criterio de Vara quebrada (Broken stick model) para seleccionar núm. ejes
# Divide en una rama de longitud unitaria, al mismo número de piezas (p) 
# en cada eje del PCA 
# Se escogen los ejes con autovalores > a su pieza correspondiente
n <- length(av)
vq <- data.frame(j=seq(1:n), p=0)
vq$p[1] <- 1/n
for (i in 2:n) {
	vq$p[i] = vq$p[i-1] + (1/(n + 1 - i))
}
vq$p <- 100*vq$p/n
vq

# Comparación gráfica de los dos modelos
# Plot eigdatosalues and % of variance for each axis
windows(title="Autovalores del PCA")
par(mfrow=c(2,1))
barplot(ev, main="Valores Propios", col="bisque", las=2)
abline(h=mean(av), col="red")	
# autovalores promedio
legend("topright", "Promedio de autovalores", lwd=1, col=2, bty="n")
barplot(t(cbind(100*av/sum(av),vq$p[n:1])), beside=TRUE, 
	main="% Varianza", col=c("bisque",2), las=2)
legend("topright", c("% Autovalores", "Modelo Vara Quebrada"), 
	pch=15, col=c("bisque",2), bty="n")

# Criterio 1= selecciona a 2 ejes.
# Criterio 2= selecciona a 3 ejes (autovalores > Particiones)	

