library(shiny)

ui <- shinyUI(fluidPage(

  titlePanel(title = "Calcium Image D4 SuperSoftware"),
  sidebarLayout(position = "left",
    sidebarPanel(h2("Data Inputs"),
                 fileInput("txt", h3("Introducir .txt"), accept = c(
                   ".txt")),
                 fileInput("estimulos", h3("Introducir estimulos.csv"), accept = c(".csv")),
                 textInput("skip", "Lineas que quitar", 5),
                 radioButtons("legends", "Show Legends?", list("TRUE", "FALSE")),
                 radioButtons("Units", "Units", list("ms", "s")),
                 radioButtons("remove", "Remove any ROI?", list("TRUE", "FALSE")),
                 fileInput("removefile", h3("Introducir remove.csv"), accept = c(".csv")),
                 radioButtons("cellone", h3("One Cell?"), list("TRUE", "FALSE")),
                 numericInput("cell", label = "cell number", value = 1),
                 radioButtons("Outliers", "Remove Multivariate Outliers?", list("TRUE", "FALSE")),
                 radioButtons("Smooth", "Signal Smooth?", list("TRUE", "FALSE")),
                 radioButtons("IntervalSel", "Select Interval?", list("TRUE", "FALSE")),
                 sliderInput("interval", "Intervalo para Indice de Oscilaciones", min = 0, max = 20, value = 1)),
    mainPanel(h3("Results"),
                 plotOutput("cabecera"),
                 tableOutput("tabla"))
    )
  )

)

server <- shinyServer(
  function(input, output){
    output$cabecera <- renderPlot({
      datos <- input$txt$datapath
      datos <- read.table(datos, header = FALSE, skip = input$skip)
      colnames(datos) <- c("Time", paste("ROI", 1:(dim(datos)[2]-1)))
      Unidades <- c("ms", "s")
      unidades <- c(6*10^4, 6*10)
      datos$Time <- datos$Time/unidades[grep(paste("^", input$Units, sep = ""), Unidades)]

      #remove those ROIs whose response is bad
      remove <- input$removefile$datapath
      if(input$remove == TRUE & !is.null(remove)){
        remove.exp <- read.csv2(remove, header = TRUE)
        datos <- datos[, -(as.numeric(remove.exp$remove)+1)]
      }
      #Remove Outliers
      mahOutlier <- function(X){
        QR.desc <- qr(X)
        rangoX <- QR.desc$rank
        centro <- colMeans(X)
        dm <- mahalanobis(X, centro, MASS::ginv(cov(X)), inverted = TRUE)
        X.menordm <- X[dm <= quantile(dm)[4],]
        centroR <- colMeans(X.menordm)
        covarianzaR <- cov(X.menordm)
        p <- ncol(X)
        dmR <- mahalanobis(X, centroR, MASS::ginv(covarianzaR), inverted = TRUE)
        outliers <- as.numeric(which(dmR > p + 3 * sqrt(2 * p)))
        return(outliers)
      }

      #Estimulos
      estimulos <- input$estimulos$datapath
      estimulos <- read.csv2(estimulos, header = TRUE)
      #Cortar el registro hasta donde interese, denominado como cut en estimulos
      if(length(grep("cut", estimulos[,1])) != 0){
        cut <- estimulos[grep("cut", estimulos[, 1]), ]
        estimulos <- estimulos[- grep("cut", estimulos[,1]), ]
        datos <- datos[datos[, 1] <= as.numeric(cut[2]), ]
      }
      Nisoldipina <- NULL
      if(length(grep("Nisoldipin", estimulos[,1])) != 0){
        Nisoldipina <- as.numeric(estimulos[grep("Nisoldipina", estimulos[, 1]), -1])
        Nisoldipina <- c(datos[sum(datos$Time <= Nisoldipina[1]), 1], datos[sum(datos$Time < Nisoldipina[2]) + 1, 1])
        estimulos <- estimulos[- grep("Nisoldipina", estimulos[,1]), ]
      }
      interval <- c(0, input$interval)
      if(length(grep("IO", estimulos[,1])) != 0 &  input$IntervalSel == FALSE){
        interval <- as.numeric(estimulos[grep("IO", estimulos[, 1]), -1])
        estimulos <- estimulos[- grep("IO", estimulos[,1]), ]
      }else{
        if(length(grep("IO", estimulos[,1])) != 0){
          estimulos <- estimulos[- grep("IO", estimulos[,1]), ]
        }
      }

      #Datos Suavizados (datos) y sin suavizar (datosraw)
      datosraw <- datos
      if(input$Smooth == TRUE){
        datos <- data.frame(apply(datos, 2, function(x){smooth(x)}))
      }
      output$tabla <- renderTable(input$interval)
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
          if(Nisoldipina[1] <= datos[sum(datos$Time < estimulos[i,2]), 1]  & Nisoldipina[2] >= datos[sum(datos$Time < (estimulos[i,3])+1), 1]){
            altura.total <- as.numeric(datos[(sum(datos$Time < estimulos[i,3])+1),-1])
            print(altura.total)
          }
        }
        if(length(grep("2f", estimulos[i, 1])) != 0)
        {
          altura.total <- as.numeric(datos[(sum(datos$Time < estimulos[i,3])+1),-1])
          print(altura.total)
        }

        #Selecciona el intervalo del estimulo y busca el min
        altura.total.min <- apply(datos[sum(datos$Time < estimulos[i,2]):(sum(datos$Time < estimulos[i,3])+1),-1], 2, min)
        #Aqui tiene en cuenta si el estimulo hace bajar o subir la senal
        alturas[,i] <- apply(rbind(altura.total, altura.total.min, as.numeric(y.estimulosS[i,])), 2, FUN = function(x){ifelse(abs(x[1] - x[3]) >= abs(x[2] - x[3]), x[1] - x[3], x[2] - x[3])}) #corrige por el valor de respuesta (min) previo al estimulo
        colnames(alturas)[i] <- paste("ALTURA", i, sep="")
      }
      rownames(alturas) <- colnames(datos)[-1] #Nombra segun los ROIs
      #variables alturas
      datos.basal <- datos
      #variables area
      areas <- data.frame(matrix(0,ncol=dim(estimulos)[1], nrow = dim(datos)[2]-1))
      for(i in 1:dim(estimulos)[1]){
        area.total <- apply(datos.basal[sum(datos.basal$Time < estimulos[i,2]):(sum(datos.basal$Time < estimulos[i,3])+1),-1], 2, function(x){trapz(x = datos.basal[sum(datos.basal$Time < estimulos[i,2]):(sum(datos.basal$Time < estimulos[i,3])+1),1], y = x)})
        area.restar <- apply(y.estimulos[c(i,(i+dim(estimulos)[1])),], 2, function(x){trapz(x = tiempos[c(i,(i+dim(estimulos)[1]))], y = x)})
        areas[,i] <- area.total- area.restar
        colnames(areas)[i] <- paste("AREA", i, sep="")
      }
      rownames(areas) <- colnames(datos)[-1]

      #Oscilations Index
      OI <- function(interval = NULL, data){

        if(is.null(interval)){
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
      }
      oscilation.index <- OI(interval = interval, data = datosraw)


      #Outliers Se tiene en cuenta que el numero de observaciones no sea menor que el numero de variables. En ese caso, no se obtienen los outliers
      datosO <- datos #En el grafico del final sin outliers se usara datosO en lugar de datos, por si a caso se han eliminado outliers cumpliendose que p >= n
      if(input$Outliers == TRUE && nrow(areas) > ncol(cbind(areas, alturas))){
        X = cbind(areas, alturas, oscilation.index)
        outliers <- mahOutlier(X)
        if(length(outliers) != 0){ #Si no hay outliers no se quitan
          areas <- areas[-outliers, ]
          alturas <- alturas[-outliers, ]
          oscilation.index <- oscilation.index[-outliers]
          datos <- datos[,-(outliers+1)]
        }
      }



      #Plot
      if(input$cellone == TRUE){
        plot(datos$Time, datos[,(input$cell + 1)], type = "l", col = 1, xlab = " ", ylab = "Ratio F340/380", ylim = c(-0.1,max(datos[,(input$cell + 1)])+max(datos[,(input$cell + 1)])*0.25), axes = FALSE)
        axis(side = 2, at = seq(0, round(max(datos[,(input$cell + 1)])+max(datos[,(input$cell + 1)])*0.25, 1), by = 0.1))
        color <- estimulos[,1]
        for(i in 1:dim(estimulos)[1]){
          lines(estimulos[i,2:3], c(0,0), lty = 1, col = as.numeric(color)[i], lwd = 10)
          text(mean(as.numeric(estimulos[i,2:3])), c(-0.05, -0.05), labels = estimulos[i,1])
          if(input$IntervalSel == TRUE){
            abline(v = interval[2], col = "red")
          }
        }
      }else{
        plot(datos$Time, datos[,2], type = "l", col = 2, xlab = " ", ylab = "Ratio F340/380", ylim = c(-0.1,max(datos[,-1])+max(datos[,-1])*0.25), axes = FALSE)
        axis(side = 2, at = seq(0, round(max(datos[,-1])+max(datos[,-1])*0.25, 1), by = 0.1))
        for(i in 3:dim(datos)[2]){
          lines(datos$Time, datos[,i], type = "l", col = i, lwd = 2)
        }
        if(input$legends == TRUE){
          legend("topleft", legend = colnames(datos)[-1], col = 2:dim(datos)[2], cex = 0.5, lty = 1, ncol = 2)
        }
        color <- estimulos[,1]
        for(i in 1:dim(estimulos)[1]){
          lines(estimulos[i,2:3], c(0,0), lty = 1, col = as.numeric(color)[i], lwd = 10)
          text(mean(as.numeric(estimulos[i,2:3])), c(-0.05, -0.05), labels = estimulos[i,1])
        }
        lines(c(max(datos[,1])-1.5, max(datos[,1])-0.5), c(max(datos[,-1])+0.15*max(datos[,-1]), max(datos[,-1])+max(datos[,-1])*0.15), lty = 1, col = "black", lwd = 10)
        text(mean(c(max(datos[,1])-1.5, max(datos[,1])-0.5)),c(max(datos[,-1])+0.20*max(datos[,-1]), max(datos[,-1])+max(datos[,-1])*0.20), labels = "1min")
        if(input$IntervalSel == TRUE){
          abline(v = interval[2], col = "red")
        }
      }


      output$tabla <- renderTable({cbind(alturas, oscilation.index)}, rownames = TRUE)
      })

    output$legends <- renderText(as.character(input$legends))


  }
)

shinyApp(ui = ui, server = server)
