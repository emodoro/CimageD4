library(shiny)

ui <- shinyUI(fluidPage(

  titlePanel(title = "Calcium Image interface"),
  sidebarLayout(position = "left",
    sidebarPanel(h2("Data Inputs"),
                 textInput("directory", "PATH de los Datos", "G:/ImagenCalcio/jelen/ANALISIS DE DATOS IMAGEN DE CALCIO/ME6/ME6_01"),
                 textInput("skip", "Lineas que quitar", 5),
                 radioButtons("legends", "Show Legends?", list("TRUE", "FALSE")),
                 radioButtons("remove", "Remove any ROI?", list("TRUE", "FALSE"))),
    mainPanel(h3("Results"),
                 plotOutput("cabecera"),
                 textOutput("legends"))
    )
  )

)

server <- shinyServer(
  function(input, output){

    output$cabecera <- renderPlot({
      tequiste <- dir(file.path(input$directory))
      tequiste <-tequiste[grep(".txt", tequiste)]
      datos <- read.table(file.path(input$directory, tequiste), header = FALSE, skip = input$skip)
      colnames(datos) <- c("Time", paste("ROI", 1:(dim(datos)[2]-1)))
      datos$Time <- datos$Time/60000
      #remove those ROIs whose response is bad
      if(input$remove == TRUE & length(grep("remove", dir(file.path(input$directory)))) != 0){
        remove.exp <- read.csv2(file.path(file.path(input$directory), "remove.csv"), header = TRUE)
        datos <- datos[, -(as.numeric(remove.exp$remove)+1)]
      }
      plot(datos$Time, datos[,2], type = "l", col = 2, xlab = "tiempo", ylab = "Ratio F340/380", ylim = c(-0.1,max(datos[,-1])+max(datos[,-1])*0.25), axes = FALSE)
      axis(side = 2, at = seq(0, round(max(datos[,-1])+max(datos[,-1])*0.25, 0), by = 0.1))
      for(i in 3:dim(datos)[2]){
        lines(datos$Time, datos[,i], type = "l", col = i, lwd = 2)
      }
      if(input$legends == TRUE){
        legend("topleft", legend = colnames(datos)[-1], col = 2:dim(datos)[2], cex = 0.5, lty = 1, ncol = 2)
      }
      #Estimulos
      estimulos <- read.csv2(file.path(input$directory, "estimulos.csv"), header = TRUE)
      #Cortar el registro hasta donde interese, denominado como cut en estimulos
      if(length(grep("cut", estimulos[,1])) != 0){
        cut <- estimulos[grep("cut", estimulos[, 1]), ]
        estimulos <- estimulos[- grep("cut", estimulos[,1]), ]
        datos <- datos[datos[, 1] <= as.numeric(cut[2]), ]
      }
      color <- estimulos[,1]
      for(i in 1:dim(estimulos)[1]){
        lines(estimulos[i,2:3], c(0,0), lty = 1, col = as.numeric(color)[i], lwd = 10)
        text(mean(as.numeric(estimulos[i,2:3])), c(-0.05, -0.05), labels = estimulos[i,1])
      }
      lines(c(max(datos[,1])-1.5, max(datos[,1])-0.5), c(max(datos[,-1])+0.15*max(datos[,-1]), max(datos[,-1])+max(datos[,-1])*0.15), lty = 1, col = "black", lwd = 10)
      text(mean(c(max(datos[,1])-1.5, max(datos[,1])-0.5)),c(max(datos[,-1])+0.20*max(datos[,-1]), max(datos[,-1])+max(datos[,-1])*0.20), labels = "1min")



      })

    output$legends <- renderText(as.character(input$legends))

  }
)

shinyApp(ui = ui, server = server)
