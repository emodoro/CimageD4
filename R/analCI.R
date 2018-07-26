#' @title summary from those data which fulfill the condition of interest
#' @description This functions provide a summary, for each experiment, and all of them into a only one .csv, of data which fulfill a condition of interest.
#' @author Enrique Perez_Riesgo
#' @param grupos
#' @return .csv
#' @export recopilationROI

recopilationROI <- function(column = "oscilation.index", variables = "oscilation.index", threshold = 1.05, category = category, centr.par = "median", disp.par = "mad"){
  directory <- getwd()
  datosfich <- file.path(directory, "resultados")
  ficheros <- dir(datosfich)
  ficheros <- ficheros[-grep("datosO", ficheros)]
  ficheros.datos <- ficheros[grep("datos", ficheros)]
  datos <- data.frame(matrix(0, nrow = length(ficheros.datos), ncol = 4))
  colnames(datos) <- c("Experiment", "category", "No", "Yes")
  media <- as.vector(matrix(0, nrow = length(ficheros.datos), ncol = 1))
  for(i in 1:length(ficheros.datos)){
    datos.tabla <- read.csv2(file.path(datosfich, ficheros.datos[i]))
    #Nombre del experimento
    exp.name <- sub("datos", x = ficheros.datos[i], replacement = "")
    exp.name <- sub(".csv", x = exp.name, replacement = "")
    #Seleccion variable de interes where asses the condition (column)
    variable <- datos.tabla[, variables]
    tabla <- c(sum(variable < threshold), sum(variable >= threshold))
    #column is the column, related to variable, where evaluate the mean. Sometimes could be the same
    if(centr.par == "median"){
      media[i] <- median(datos.tabla[variable >= threshold, column])
    }
    if(centr.par == "mean"){
      media[i] <- mean(datos.tabla[variable >= threshold, column])
    }
    datos[i, ] <- c(exp.name, as.character(category[i]), tabla)
  }

  if(centr.par == "median"){
    datos <- cbind(datos, Median = media)
    datos <- rbind(datos, n= c(NA, NA, sum(as.numeric(datos$No)), sum(as.numeric(datos$Yes)), NA), Median = c(NA, NA, NA, NA, median(datos$Median)), Sd = c(NA, NA, NA, NA, mad(datos$Median)))
    rownames(datos) <- c(1: (nrow(datos) - 3), "n", "Median", "mad")
  }
  if(centr.par == "mean"){
    datos <- cbind(datos, Mean = media)
    datos <- rbind(datos, n= c(NA, NA, sum(as.numeric(datos$No)), sum(as.numeric(datos$Yes)), NA), Mean = c(NA, NA, NA, NA, mean(datos$Mean)), Sd = c(NA, NA, NA, NA, sd(datos$Mean)))
    rownames(datos) <- c(1: (nrow(datos) - 3), "n", "Mean", "Sd")
  }
  write.csv2(datos, file = file.path(directory, paste("resumen", variables, centr.par, ".csv", sep = "")))
}

#' @title wave length
#' @description This functions provide the number of peaks, oscilations.
#' @author Enrique Perez_Riesgo
#' @param grupos
#' @return plots
#' @export wave.length

wave.length <- function(interval = NULL, data){

  if(is.null(interval)){ #En caso de ser NULL el intervalo, por defecto se toma el intervalo [0, 1] para estimar el índice de oscilaciones
    interval <- c(0, 1)
  }
  data <- data[data[, 1] <= interval[2] &  data[, 1] >= interval[1], ]
  sigma <- apply(data[, -1], MARGIN = 2, function(X){summary(lm(X ~ data[, 1]))$sigma})
  resid <- apply(data[, -1], MARGIN = 2, function(X){residuals(lm(X ~ data[, 1]))})
  DW <- apply(data[, -1], MARGIN = 2, function(X){lmtest::dwtest(data[, 1] ~ X)$p.value})
  autoc <- DW <= 0.05

  resid.lag <- rbind(c(rep(NA, ncol(resid))), resid[- nrow(resid), ])
  resid.sign <- sign(resid.lag) + sign(resid)
  resid.sign2 <- apply(resid.sign[-1, ], 2, function(x){
    results <- NULL
    for(i in 1:(length(x))){results <- c(results, ifelse(x[i] == 0 & x[(i + 1)] == 0,x[(i - 1)] , x[(i)]))}
    return(results)
  })
  resid.sign2[nrow(resid.sign2), ] <- resid.sign[nrow(resid.sign), ]
  wave.length <- apply(resid.sign2, 2, function(x){sum(x != 0)/sum(x == 0)})*(data[2, 1] - data[1, 1])
  amplitud <-  apply(resid, 2, function(x){abs(max(x) - min(x))})


  areas <- apply(abs(resid) , 2, function(x){trapz(x = data[, 1], y = x)})

  return(data.frame(WL = wave.length, Amplitud = amplitud, OA = areas, Dispersion = sigma))
}



#' @title Oscilations number.
#' @description This functions provide the number of peaks, oscilations.
#' @author Enrique Perez_Riesgo
#' @param grupos
#' @return plots
#' @export OI

OI <- function(interval = NULL, data){

  if(is.null(interval)){ #En caso de ser NULL el intervalo, por defecto se toma el intervalo [0, 1] para estimar el índice de oscilaciones
    interval <- c(0, 1)
  }
  data <- data[data[, 1] <= interval[2] &  data[, 1] >= interval[1], ]
  data.lag <- data[- nrow(data), ]
  diferencia <- data[-1, ]- data.lag
  diferencia2 <- diferencia ^ 2
  d2 <- diferencia2[, -1] + diferencia2[, 1]
  distancia <- sqrt(d2)
  IO <- colSums(distancia)/(data[nrow(data), 1]- data[1, 1])
  return(IO)
  #prueba <- NULL
  #for(i in 2:nrow(as.matrix(dist(data[, c(1, 3)])))){prueba <- c(prueba, as.matrix(dist(data[, c(1, 3)]))[i, (i-1)])}
  #sum(prueba)
}




#' @title Multivariate outliers.
#' @description An oultier detection is carried out by applying Mahalanobis distances, where means vetor and correlation matrix is computed with those data whose mahalanobis distance is lower than median of mahalanobis distances.
#' @author Enrique Perez_Riesgo
#' @param grupos
#' @return plots
#' @export mahOutlier

mahOutlier <- function(X){
  #rango
  QR.desc <- qr(X)
  rangoX <- QR.desc$rank
  # Media de todos los datos, el centro
  centro <- colMeans(X)
  # Distancia de Mahalanobis al centro
  #if(rangoX < ncol(X)){
    dm <- mahalanobis(X, centro, MASS::ginv(cov(X)), inverted = TRUE)
  #}else{
    #dm <- mahalanobis(X, centro, cov(X))
  #}
  # Selección del 50% de los datos con menor dm
  X.menordm <- X[dm <= quantile(dm)[4],]
  # Estimadores reducidos
  centroR <- colMeans(X.menordm)
  covarianzaR <- cov(X.menordm)
  # distancias de mahalanobis al centro reducido
  p <- ncol(X)
  #if(rangoX < ncol(X)){
    dmR <- mahalanobis(X, centroR, MASS::ginv(covarianzaR), inverted = TRUE)
  #}else{
   # dmR <- mahalanobis(X, centroR, covarianzaR)
  #}

  outliers <- as.numeric(which(dmR > p + 3 * sqrt(2 * p)))
  return(outliers)
}




#' @title Calcium Image data analysis.
#' @description this function allows arring out a comprensive calcium image data
#'   analysis. Notably, the folder which stores the .txt also has to store a
#'   .csv file where each row correspond to a stimuli, shuch that the first
#'   colomn store the name of each stimuli, the second column the stimuli start
#'   time, in seconds, and the third column the stimuli-end-time. Aditional rows
#'   can be typing: cut: This row meaning the interval of experiment which you
#'   want to remove.
#'   IO: Interval where asses oscilation index.
#'   Nisoldipina: The interval of experiment where Nisoldipina,
#'   an antagonist of VOCCs, is emplyed in order to remove the masking efect of
#'   Calcium entry through VOCCs when you are concern with asses SOCE
#' @author Enrique Perez_Riesgo
#' @param grupos Allows stablishing a number of groups deshired
#' @param
#' @return plots
#' @export analCI

#Analisis de imagen

analCI <- function(grupos = NULL, agrupacion = "silueta", modo = "Kmedioids", outlier = TRUE,directory = NULL, skip = 5, data.scale = TRUE, legend.ROIs = TRUE, interval = NULL, Units = "ms", Smooth. = TRUE, y.int =c(0, 1.5), min.threshold = 0) {
  #directories
  if(is.null(directory)){
    directory <- getwd()
  }
  archivos <- dir(directory)

  if(outlier == TRUE){
    if(length(grep(pattern = "resultadosOUT", archivos)) == 0){
      dir.create(file.path(directory, "resultadosOUT"))
    }
    results.dir <- file.path(directory, "resultadosOUT")
  }
  if(outlier == FALSE){
    if(length(grep(pattern = "resultados", archivos)) == 0){
      dir.create(file.path(directory, "resultados"))
    }
    results.dir <- file.path(directory, "resultados")
  }


  if(length(grep("resultados", archivos)) != 0){archivos <- archivos[-(grep("resultados", archivos))]}
  if(length(grep(".Rmd", archivos)) != 0){archivos <- archivos[-(grep(".Rmd", archivos))]}
  if(length(grep(".csv", archivos)) != 0){archivos <- archivos[-(grep(".csv", archivos))]}
  if(length(grep(".xls", archivos)) != 0){archivos <- archivos[-(grep(".xls", archivos))]}
  if(length(grep(".txt", archivos)) != 0){archivos <- archivos[-(grep(".txt", archivos))]}
  grupo.numero <- 0
  #read .txt
  for(z in archivos){
    tequiste <- dir(file.path(directory, z))
    tequiste <-tequiste[grep(".txt", tequiste)]
    datos <- read.table(file.path(file.path(directory, z), tequiste), header = FALSE, skip = skip)
    colnames(datos) <- c("Time", paste("ROI", 1:(dim(datos)[2]-1)))
    #Unidades
    Unidades <- c("ms", "s")
    unidades <- c(6*10^4, 6*10)
    datos$Time <- datos$Time/ unidades[grep(paste("^", Units, sep = ""), Unidades)]
    #remove those ROIs whose response is bad
    if(length(grep("remove", dir(file.path(directory, z)))) != 0){
      remove.exp <- read.csv2(file.path(file.path(directory, z), "remove.csv"), header = TRUE)
      datos <- datos[, -(as.numeric(remove.exp$remove)+1)]
    }

    #Estímulos
    estimulos <- read.csv2(file.path(file.path(directory, z), "estimulos.csv"), header = TRUE)

    #Cortar el registro hasta donde interese, denominado como cut en estimulos
    if(length(grep("cut", estimulos[,1])) != 0){
      cut <- estimulos[grep("cut", estimulos[, 1]), ]
      estimulos <- estimulos[- grep("cut", estimulos[,1]), ]
      datos <- datos[datos[, 1] <= as.numeric(cut[2]), ]

    }
    if(length(grep("IO", estimulos[,1])) != 0){
      interval <- as.numeric(estimulos[grep("IO", estimulos[, 1]), -1])
      estimulos <- estimulos[- grep("IO", estimulos[,1]), ]
    }
    Nisoldipina <- NULL
    if(length(grep("Nisoldipin", estimulos[,1])) != 0){
      Nisoldipina <- as.numeric(estimulos[grep("Nisoldipina", estimulos[, 1]), -1])
      Nisoldipina <- c(datos[sum(datos$Time <= Nisoldipina[1]), 1], datos[sum(datos$Time < Nisoldipina[2]) + 1, 1])
      estimulos <- estimulos[- grep("Nisoldipina", estimulos[,1]), ]
    }

    #elimnar segun min.threshold
    encima <- apply(datos[, -1] <= min.threshold, MARGIN = 2, sum)
    datos <- datos[, c(TRUE, encima == 0)]


    #pdf datos sin suavizar
    pdf(paste(results.dir,"/Graficos", z, ".pdf", sep = ""))
    plot(datos$Time, datos[,2], type = "l", col = 2, xlab = "tiempo", ylab = "Ratio F340/380", ylim = c(ifelse(is.null(y.int[1]), -0.1, y.int[1]), ifelse(is.null(y.int[2]), max(datos[,-1])*1.25, y.int[2]*1.25)) , main = z, axes = FALSE)
    axis(side = 2, at = seq(0, round(max(datos[,-1])+max(datos[,-1])*0.25, 0), by = 0.1))
    for(i in 3:dim(datos)[2]){
      lines(datos$Time, datos[,i], type = "l", col = i, lwd = 2)
    }
    if(legend.ROIs == TRUE){
      legend("topleft", legend = colnames(datos)[-1], col = 2:dim(datos)[2], cex = 0.5, lty = 1, ncol = 2)
    }

    color <- estimulos[,1]
    for(i in 1:dim(estimulos)[1]){
      lines(estimulos[i,2:3], c(0,0), lty = 1, col = as.numeric(color)[i], lwd = 10)
      text(mean(as.numeric(estimulos[i,2:3])), c(-0.05, -0.05), labels = estimulos[i,1])
    }
    lines(c(max(datos[,1])-1.5, max(datos[,1])-0.5), c(max(datos[,-1])+0.15*max(datos[,-1]), max(datos[,-1])+max(datos[,-1])*0.15), lty = 1, col = "black", lwd = 10)
    text(mean(c(max(datos[,1])-1.5, max(datos[,1])-0.5)),c(max(datos[,-1])+0.20*max(datos[,-1]), max(datos[,-1])+max(datos[,-1])*0.20), labels = "1min")
    dev.off()

    #Datos Suavizados (datos) y sin suavizar (datosraw)
    datosraw <- datos
    if(Smooth. == TRUE){
      datos <- data.frame(apply(datos, 2, function(x){smooth(x)}))
    }

    #ejes estimulos
    require(pracma)
    y.estimulosS <- data.frame(matrix(0,nrow=dim(estimulos)[1], ncol=dim(datos)[2]-1)) #matriz de longitud estiulosXROIS
    y.estimulosE <- data.frame(matrix(0,nrow=dim(estimulos)[1], ncol=dim(datos)[2]-1)) #matriz de longitud estiulosXROIS
    if(estimulos[dim(estimulos)[1], dim(estimulos)[2]] > datos$Time[length(datos$Time)]){
      estimulos[dim(estimulos)[1], dim(estimulos)[2]] <- datos$Time[length(datos$Time)]
    } #En caso de que el estimulo final termine mas tarde del tiempo de los datos, se fija que termina cuando acaba el tiempo de los datos
    if(estimulos[1, 2] == 0){
      estimulos[1, 2] <- datos$Time[2]
    } #En caso de que el estimulo inicial empiece a tiempo cero, se corrige para que sea el instante posterior y no de posteriormente errores
    for(i in 1:dim(estimulos)[1]){
      posicionS <- datos[datos$Time < estimulos[i,2],] #TRUE todos los valores cuyo tiempo sea interior al del estimulo
      y.estimulosS[i,] <- as.numeric(posicionS[dim(posicionS)[1],-1]) #Se fija el valor de respuesta correspondiente al ultimo TRUE de posicionS, y sera el valor de respuesta cuando comienza el estimulo (basal por estimulo)
      rownames(y.estimulosS)[i] <- as.character(posicionS[dim(posicionS)[1],1]) #Se guarda a que tiempo comienza el estimulo segun los datos temporales registrados
      y.estimulosE[i,] <- as.numeric(datos[datos$Time > estimulos[i,3],-1][1,]) #Aqui el primer true corresponde con el fin del estimulo
      if(sum(is.na(y.estimulosE[i,])) != 0){
        y.estimulosE[i,] <- as.numeric(datos[datos$Time >= estimulos[i,3],-1][1,])
        rownames(y.estimulosE)[i] <- as.character(datos[datos$Time >= estimulos[i,3],1][1]) #En caso de que no haya tiempos mayores se coge el igual o mayor en vez de mayor
      }else {
        rownames(y.estimulosE)[i] <- as.character(datos[datos$Time > estimulos[i,3],1][1])
      }
    }
    #Alturas
    alturas <- data.frame(matrix(0,ncol=(dim(estimulos)[1]+2), nrow = dim(datos)[2]-1)) #Matriz de dimension ROIS X Estimulos+1
    y.estimulos <- rbind(y.estimulosS, y.estimulosE) #Se une la matriz de registro para principio del estimulo con la del final del estimulo para cada ROI
    tiempos <- as.numeric(rownames(y.estimulos)) #Los tiempos de principio y fin de estimulo segun tiempo regustrado
    alturas[,(dim(alturas)[2]-1)] <- as.numeric(datos[1,-1]) #datos basales inicio
    alturas[,dim(alturas)[2]] <- as.numeric(datos[nrow(datos),-1]) #datos basales final
    colnames(alturas)[(dim(alturas)[2]-1)] <- "BASAL.Principio"
    colnames(alturas)[dim(alturas)[2]] <- "BASAL.Final"
    for(i in 1:dim(estimulos)[1]){
      #Selecciona el intervalo del estimulo y busca el max
      altura.total <- apply(datos[sum(datos$Time < estimulos[i,2]):(sum(datos$Time < estimulos[i,3])+1),-1], 2, max)
      if(!is.null(Nisoldipina)){ #Primero se evalua si Nisoldipina es NULL, y en caso negativo, si el estimulo i esta dentro del intervalo de accion de la nisoldipina
        print(Nisoldipina)
        if(Nisoldipina[1] <= datos[sum(datos$Time < estimulos[i,2]), 1]  & Nisoldipina[2] >= datos[sum(datos$Time < (estimulos[i,3])+1), 1])
        {
          altura.total <- as.numeric(datos[(sum(datos$Time < estimulos[i,3])+1),-1])
        }
      }
      #Segundas fases
      if(length(grep("2f", estimulos[i, 1])) != 0)
      {
        altura.total <- as.numeric(datos[(sum(datos$Time < estimulos[i,3])+1),-1])
      }
      #Selecciona el intervalo del estimulo y busca el min
      altura.total.min <- apply(datos[sum(datos$Time < estimulos[i,2]):(sum(datos$Time < estimulos[i,3])+1),-1], 2, min)
      #Aquí tiene en cuenta si el estímulo hace bajar o subir la señal
      alturas[,i] <- apply(rbind(altura.total, altura.total.min, as.numeric(y.estimulosS[i,])), 2, FUN = function(x){ifelse(abs(x[1] - x[3]) >= abs(x[2] - x[3]), x[1] - x[3], x[2] - x[3])}) #corrige por el valor de respuesta (min) previo al estimulo
      colnames(alturas)[i] <- paste("ALTURA", i, sep="")
    }
    rownames(alturas) <- colnames(datos)[-1] #Nombra segun los ROIs
    #variables alturas
    datos.basal <- datos
    #variables área
    areas <- data.frame(matrix(0,ncol=dim(estimulos)[1], nrow = dim(datos)[2]-1))
    for(i in 1:dim(estimulos)[1]){
      area.total <- apply(datos.basal[sum(datos.basal$Time < estimulos[i,2]):(sum(datos.basal$Time < estimulos[i,3])+1),-1], 2, function(x){trapz(x = datos.basal[sum(datos.basal$Time < estimulos[i,2]):(sum(datos.basal$Time < estimulos[i,3])+1),1], y = x)})
      area.restar <- apply(y.estimulos[c(i,(i+dim(estimulos)[1])),], 2, function(x){trapz(x = tiempos[c(i,(i+dim(estimulos)[1]))], y = x)})
      areas[,i] <- area.total- area.restar
      colnames(areas)[i] <- paste("AREA", i, sep="")
    }
    rownames(areas) <- colnames(datos)[-1]

    #Oscilations Index
    oscilation.index <- OI(interval = interval, data = datosraw)
    longitud.onda <- wave.length(interval = interval, data = datosraw)


    #table
    write.csv2(cbind(areas, alturas, longitud.onda, oscilation.index), paste(results.dir,"/datos", z, ".csv", sep = ""))


    #Decidir si hay o no señal
    dispersion <- longitud.onda$Dispersion
    datos.responden <- cbind(alturas[, -c((ncol(alturas) - 1): ncol(alturas))], dispersion)
    colnames(datos.responden)[1:length(estimulos[, 1])] <- as.character(estimulos[, 1])
    phase2 <- grep("2f", estimulos[, 1])
    decision <- t(apply(datos.responden, 1, function(x){
      signal <- x[-length(x)] >= as.numeric(x[length(x)])*3.29
      signal.lag <- x[-length(x)][phase2-1] >= as.numeric(x[length(x)])*3.29 #deteccion de la primera fase
      signal.2 <- x[-length(x)][phase2] >= as.numeric(x[length(x)])*1.645 #decision de si hay segunda fase
      signal[phase2] <- signal.lag * signal.2 #para que haya segunda fase ha de cumplirse que haya primera
      return(signal)
    }))


    #tabla total
    write.csv2(cbind(areas, alturas, longitud.onda, oscilation.index, decision), paste(results.dir,"/datos", z, ".csv", sep = ""))

    #Outliers Se tiene en cuenta que el numero de observaciones no sea menor que el numero de variables. En ese caso, no se obtienen los outliers
    datosO <- datos #En el gráfico del final sin outliers se usará datosO en lugar de datos, por si a caso se han eliminado outliers cumpliendose que p >= n
    if(outlier == TRUE && nrow(areas) > ncol(cbind(areas, alturas))){
      X = cbind(areas, alturas)
      outliers <- mahOutlier(X)
      if(length(outliers) != 0){ #Si no hay outliers no se quitan
        areas <- areas[-outliers, ]
        alturas <- alturas[-outliers, ]
        longitud.onda <- longitud.onda[-outliers, ]
        oscilation.index <- oscilation.index[-outliers]
        datosO <- datos[,-(outliers+1)]
        decision <- decision[-outliers, ]
      }
    }


    #distancias
    distancias <- dist(scale(cbind(areas, alturas, oscilation.index)), method = "euclidean")
    grupo.numero <- grupo.numero + 1
    if(is.null(grupos)){
      grupos = rep(3, length(archivos))
    }
    if(agrupacion == "silueta"){
      grupos = rep(0, length(archivos))
      tope <- ifelse(ncol(datosO[,-1]) >= 6, 5, ncol(datosO)-2)
      siluetas <- as.numeric(vector(length = tope-1))
      for (i in 2:tope) {
        silueta.media <- cluster::silhouette(cluster::pam(scale(cbind(areas, alturas)), k = i), grupos[grupo.numero], dist = distancias)
        siluetas[(i-1)] <- mean(silueta.media[,3])
      }
      grupos[grupo.numero] <- which(siluetas == max(siluetas))+1
    }
    pdf(paste(results.dir,"/Cluster", z, ".pdf", sep = ""))
    plot(hclust(distancias, method = "ward.D2"), main = z)
    dev.off()
    kmedioides <- cluster::pam(scale(cbind(areas, alturas, oscilation.index)), grupos[grupo.numero])
    grupos2 <- cutree(hclust(distancias), grupos[grupo.numero])

    #tabla total
    write.csv2(cbind(areas, alturas, longitud.onda, oscilation.index, grupos.Kmedioids = kmedioides$clustering, grupos.Cluster = grupos2, decision), paste(results.dir,"/datosOut", z, ".csv", sep = ""))

    #Grupos segun seleccion en argumento
    if(modo == "Kmedioids"){
      grupitos = as.numeric(table(kmedioides$clustering))
      asignacion <- kmedioides$clustering
    }
    if(modo == "cluster"){
      grupitos = as.numeric(table(grupos2))
      asignacion <- grupos2
    }
    require(ggplot2)
    require(ggfortify)
    require(cluster)
    PCA <- prcomp(cbind(areas, alturas, oscilation.index), scale. = data.scale)
    pdf(paste(results.dir,"/PCA", z, ".pdf", sep = ""))
    #autoplot(PCA, shape = FALSE, label.size = 3)
    plot(PCA$x[,1], PCA$x[,2], xlab = paste("PC1", round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100,2),"%"), ylab = paste("PC2", round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100,2),"%"), main = z, col = 0)
    text(PCA$x[,1], PCA$x[,2], labels = colnames(datosO)[-1], col = as.numeric(asignacion))
    dev.off()
    tabla.medias <- apply(cbind(areas, alturas, longitud.onda, oscilation.index), MARGIN = 2, FUN = tapply, INDEX=asignacion, mean)
    tabla.desviaciones <- apply(cbind(areas, alturas, longitud.onda, oscilation.index), MARGIN = 2, FUN = tapply, INDEX=asignacion, sd)/sqrt(grupitos)

    #ALTURAS
    pdf(paste(results.dir,"/Barras.Altura", z, ".pdf", sep = ""))
    media.altura <- tabla.medias[,grep("ALTURA", colnames(tabla.medias))]
    media.altura[is.na(media.altura)] <- 0
    desviacion.altura <- tabla.desviaciones[,grep("ALTURA", colnames(tabla.desviaciones))]
    desviacion.altura[is.na(desviacion.altura)] <- 0

    #diagramas de barras
    barras <- barplot(media.altura, beside = TRUE, las = 2, cex.names = 1, ylim = c(min(pretty(media.altura - desviacion.altura)), max(media.altura + desviacion.altura)*1.2), col = 2:(length(grupitos)+1), main = z, names.arg = estimulos[,1])
    arrows(barras, media.altura + desviacion.altura, barras, media.altura - desviacion.altura, angle = 90, code = 3)
    legend("topright", legend = paste("n = ", grupitos), fill = 2:(length(grupitos)+1))
    box(bty = "l")
    dev.off()

    #OI
    pdf(paste(results.dir,"/Barras.OI", z, ".pdf", sep = ""))
    media.OI <- tabla.medias[,grep("oscilation.index", colnames(tabla.medias))]
    media.OI[is.na(media.OI)] <- 0
    desviacion.OI <- tabla.desviaciones[,grep("oscilation.index", colnames(tabla.desviaciones))]
    desviacion.OI[is.na(desviacion.OI)] <- 0

    #diagramas de barras
    barras <- barplot(media.OI, beside = TRUE, las = 2, cex.names = 1, ylim = c(min(pretty(media.OI - desviacion.OI)), max(pretty(media.OI + desviacion.OI*1.2))), col = 2:(length(grupitos)+1), main = z, xpd = F)
    arrows(barras, media.OI + desviacion.OI, barras, media.OI - desviacion.OI, angle = 90, code = 3)
    legend("topright", legend = paste("n = ", grupitos), fill = 2:(length(grupitos)+1))
    box(bty = "l")
    dev.off()

    #OA
    pdf(paste(results.dir,"/Barras.OA", z, ".pdf", sep = ""))
    media.OA <- tabla.medias[,grep("OA", colnames(tabla.medias))]
    media.OA[is.na(media.OA)] <- 0
    desviacion.OA <- tabla.desviaciones[,grep("oscilation.index", colnames(tabla.desviaciones))]
    desviacion.OA[is.na(desviacion.OA)] <- 0

    #diagramas de barras
    barras <- barplot(media.OA, beside = TRUE, las = 2, cex.names = 1, ylim = c(min(pretty(media.OA - desviacion.OA)), max(pretty(media.OA + desviacion.OA*1.2))), col = 2:(length(grupitos)+1), main = z, xpd = F)
    arrows(barras, media.OA + desviacion.OA, barras, media.OA - desviacion.OA, angle = 90, code = 3)
    legend("topright", legend = paste("n = ", grupitos), fill = 2:(length(grupitos)+1))
    box(bty = "l")
    dev.off()

    #Dispersion
    pdf(paste(results.dir,"/Barras.sigma", z, ".pdf", sep = ""))
    media.Dispersion <- tabla.medias[,grep("Dispersion", colnames(tabla.medias))]
    media.Dispersion[is.na(media.Dispersion)] <- 0
    desviacion.Dispersion <- tabla.desviaciones[,grep("Dispersion", colnames(tabla.desviaciones))]
    desviacion.Dispersion[is.na(desviacion.Dispersion)] <- 0

    #diagramas de barras
    barras <- barplot(media.Dispersion, beside = TRUE, las = 2, cex.names = 1, ylim = c(min(pretty(media.Dispersion - desviacion.Dispersion)), max(pretty(media.Dispersion + desviacion.Dispersion*1.2))), col = 2:(length(grupitos)+1), main = z, xpd = F)
    arrows(barras, media.Dispersion + desviacion.Dispersion, barras, media.Dispersion - desviacion.Dispersion, angle = 90, code = 3)
    legend("topright", legend = paste("n = ", grupitos), fill = 2:(length(grupitos)+1))
    box(bty = "l")
    dev.off()


    #tablas medias desviaciones
    descriptiva <- matrix(0, ncol = (2*length(grupitos)+2), nrow = (length(estimulos[,1])+4))
    rownames(descriptiva) <- c(as.character(estimulos[,1]), "n", "OI", "OA", "sigma")
    colnames(descriptiva) <- c(paste(rep(c("media", "desviación"), length(grupitos)), rep(1:length(grupitos), each = 2)), "Media Global", "Desviación Global")
    descriptiva <- data.frame(descriptiva)
    for(i in 1:length(grupitos)){
      descriptiva[,(2*(i-1)+1)] <- c(signif(t(media.altura)[,i],2),grupitos[i], t(media.OI)[,i], t(media.OA)[,i], t(media.Dispersion)[,i])
      descriptiva[,(2*(i))] <- c(signif(t(desviacion.altura)[,i],2)," ", t(desviacion.OI)[,i], t(desviacion.OA)[,i] , t(desviacion.Dispersion)[,i])
    }
    if(length(estimulos[,1]) > 1){
      descriptiva[,(dim(descriptiva)[2]-1)] <- c(signif(apply(alturas[,grep("ALTURA", colnames(alturas))], MARGIN = 2, mean), 2), sum(grupitos), mean(oscilation.index), mean(longitud.onda$OA), mean(longitud.onda$Dispersion))
      descriptiva[,(dim(descriptiva)[2])] <- c(signif(apply(alturas[,grep("ALTURA", colnames(alturas))], MARGIN = 2, sd),2), "", sd(oscilation.index), sd(longitud.onda$OA), sd(longitud.onda$Dispersion))
    }else{
      descriptiva[,(dim(descriptiva)[2]-1)] <- c(signif(mean(alturas[,grep("ALTURA", colnames(alturas))]),2), sum(grupitos), mean(oscilation.index), mean(longitud.onda$OA), mean(longitud.onda$Dispersion))
      descriptiva[,(dim(descriptiva)[2])] <- c(signif(sd(alturas[,grep("ALTURA", colnames(alturas))]),2), "", sd(oscilation.index), sd(longitud.onda$OA), sd(longitud.onda$Dispersion))
    }

    write.csv2(descriptiva, paste(results.dir,"/descriptiva", z, ".csv", sep = ""))

    #pdf raw data agrupados y sin outliers
    pdf(paste(results.dir,"/Graficos_Grupos", z, ".pdf", sep = ""))
    color <- as.numeric(asignacion) + 1
    plot(datos$Time, datosO[,2], type = "l", col = color[1], xlab = "tiempo", ylab = "Ratio F340/380", ylim = c(ifelse(is.null(y.int[1]), -0.1, y.int[1]), ifelse(is.null(y.int[2]), max(datosO[,-1])*1.25, y.int[2]*1.25)), main = z, axes = FALSE)
    axis(side = 2, at = seq(0, round(max(datosO[,-1])+max(datosO[,-1])*0.25, 1), by = 0.1))
    for(i in 3:dim(datosO)[2]){
      lines(datosO$Time, datosO[,i], type = "l", col = color[(i-1)], lwd = 2)
    }
    legend("topleft", legend = paste("n = ", grupitos), fill = unique(color))

    color <- estimulos[,1]
    for(i in 1:dim(estimulos)[1]){
      lines(estimulos[i,2:3], c(0,0), lty = 1, col = as.numeric(color)[i], lwd = 10)
      text(mean(as.numeric(estimulos[i,2:3])), c(-0.05, -0.05), labels = estimulos[i,1])
    }
    lines(c(max(datosO[,1])-1.5, max(datosO[,1])-0.5), c(max(datosO[,-1])+0.15*max(datosO[,-1]), max(datosO[,-1])+max(datosO[,-1])*0.15), lty = 1, col = "black", lwd = 10)
    text(mean(c(max(datosO[,1])-1.5, max(datosO[,1])-0.5)),c(max(datosO[,-1])+0.20*max(datosO[,-1]), max(datosO[,-1])+max(datosO[,-1])*0.20), labels = "1min")
    dev.off()
    print(paste("Exp", z, "(", grep(z, archivos), "of", length(archivos),")"))
  }
}
