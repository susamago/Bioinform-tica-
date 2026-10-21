###############################################################################
# TALLER R: expresión diferencial en GSE5583
# Referencia: TallerR.R, adaptado del código de Ahmed Moustafa.
# Datos: microarray de ratón, 3 muestras WT y 3 muestras HDAC1 KO.
# Abre este archivo en RStudio y pulsa Source para ejecutar todo el taller.
# Las respuestas están en los comentarios de cada apartado.
###############################################################################

# PREPARACIÓN: localizar esta carpeta sin cambiar el directorio de trabajo.
ruta_script <- NULL
for (marco in sys.frames()) {
  if (!is.null(marco$ofile)) ruta_script <- marco$ofile
}
if (is.null(ruta_script)) {
  argumento <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argumento)) ruta_script <- sub("^--file=", "", argumento[1])
}
carpeta <- if (!is.null(ruta_script)) {
  dirname(normalizePath(ruta_script, mustWork = TRUE))
} else {
  getwd()
}
carpeta_datos <- file.path(carpeta, "datos")
carpeta_resultados <- file.path(carpeta, "resultados")
carpeta_graficos <- file.path(carpeta_resultados, "graficos")
dir.create(carpeta_datos, showWarnings = FALSE, recursive = TRUE)
dir.create(carpeta_graficos, showWarnings = FALSE, recursive = TRUE)

# Paquetes del heatmap avanzado. Se instalan solo si faltan, desde CRAN,
# en una biblioteca de esta carpeta. No se pide actualizar otros paquetes.
biblioteca <- file.path(carpeta, "R_lib")
dir.create(biblioteca, showWarnings = FALSE)
.libPaths(c(biblioteca, .libPaths()))
for (paquete in c("gplots", "RColorBrewer")) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    install.packages(paquete, lib = biblioteca,
                     repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace(paquete, quietly = TRUE)) {
    stop("No se ha podido cargar el paquete: ", paquete)
  }
}

# Esta función guarda cada gráfico en PNG y, en una sesión interactiva,
# también lo muestra en Plots. Solo cierra el dispositivo del archivo PNG.
guardar_grafico <- function(nombre, dibujar, ancho = 1800, alto = 1400) {
  png(file.path(carpeta_graficos, paste0(nombre, ".png")),
      width = ancho, height = alto, res = 160,
      type = if (capabilities("aqua")) "quartz" else getOption("bitmapType"))
  tryCatch(dibujar(), finally = dev.off())
  if (interactive()) dibujar()
  invisible(NULL)
}

# 1. CARGAR LOS DATOS Y EXAMINARLOS
# Se incluye la tabla del enlace original del taller. No hace falta RCurl.
archivo <- file.path(carpeta_datos, "GSE5583.txt")
url_datos <- paste0(
  "https://gist.githubusercontent.com/ahmedmoustafa/5520526/raw/",
  "0bd8677a0b7f713d4dc7b9a72775d1ec474c77a8/GSE5583.txt"
)
if (!file.exists(archivo)) download.file(url_datos, archivo, mode = "wb")
data <- as.matrix(read.table(archivo, header = TRUE, row.names = 1))
stopifnot(is.numeric(data), ncol(data) == 6, all(is.finite(data)),
          all(data > 0), !anyDuplicated(rownames(data)))
print(dim(data))
print(head(data))
print(tail(data))
# Respuesta: 12488 filas y 6 columnas. Cada fila es una sonda/identificador
# de Affymetrix; distintas sondas pueden corresponder al mismo gen.
# Las tres primeras columnas son WT y las tres últimas, HDAC1 KO.

guardar_grafico("01_histograma_original", function() {
  hist(data, breaks = 60, col = "gray75", border = "white",
       main = "GSE5583: intensidades originales",
       xlab = "Intensidad", ylab = "Frecuencia")
})

# 2. TRANSFORMACIÓN LOGARÍTMICA
data2 <- log2(data)
guardar_grafico("02_histograma_log2", function() {
  hist(data2, breaks = 60, col = "gray75", border = "white",
       main = "GSE5583: intensidades transformadas",
       xlab = "log2(intensidad)", ylab = "Frecuencia")
})
# ¿Qué ocurre y para qué sirve?
# log2 comprime las intensidades altas y reduce la asimetría hacia la derecha.
# Convierte cocientes en diferencias: duplicar la intensidad suma 1 en log2.
# Puede ayudar a estabilizar la variabilidad, pero no garantiza normalidad.
# Un histograma global no comprueba la normalidad de las réplicas de cada gen.
# Todos los valores son positivos, por lo que no necesitamos sumar 1.
# Si log(X) fuese normal, X sería lognormal; no se llama lognormal a log(X).

# 3. BOXPLOT Y SUS PARÁMETROS
colores_grupo <- c(rep("#0072B2", 3), rep("#E69F00", 3))
etiquetas <- sub("\\.", "\n", colnames(data))
guardar_grafico("03_boxplot", function() {
  par(mar = c(7, 4, 4, 2))
  boxplot(data2, col = colores_grupo, names = etiquetas,
          main = "GSE5583: distribución por muestra", las = 2,
          ylab = "log2(intensidad)")
  legend("topleft", c("WT", "HDAC1 KO"),
         fill = c("#0072B2", "#E69F00"), bty = "n")
})
# col asigna los colores; main, el título; las = 2 coloca las etiquetas
# perpendiculares a los ejes. La caja va del primer al tercer cuartil;
# la línea interior es la mediana. Los bigotes llegan al dato más extremo
# que no supera 1.5 veces el rango intercuartílico desde la caja.
# Los puntos restantes se dibujan como valores atípicos.

# 4. AGRUPAMIENTO JERÁRQUICO DE LAS MUESTRAS
hc <- hclust(as.dist(1 - cor(data2)))
guardar_grafico("04_clustering_muestras", function() {
  plot(hc, labels = etiquetas, main = "GSE5583: agrupamiento de muestras",
       xlab = "", sub = "Distancia: 1 - correlación de Pearson; enlace completo",
       ylab = "Distancia")
  rect.hclust(hc, k = 2, border = c("#0072B2", "#E69F00"))
})
grupos_cluster <- cutree(hc, k = 2)
print(grupos_cluster)
# ¿Es correcta la separación? Sí: al cortar el árbol en dos grupos,
# las tres muestras WT quedan juntas y las tres KO forman el otro grupo.
# cor calcula correlaciones entre columnas; as.dist crea una distancia;
# hclust construye el árbol (enlace completo por defecto).

# 5. MATRICES WT/KO, MEDIAS Y GRÁFICO DE DISPERSIÓN
wt <- data2[, 1:3]
ko <- data2[, 4:6]
print(class(wt))
# Son matrices numéricas de 12488 filas por 3 columnas.
wt.mean <- apply(wt, 1, mean)
ko.mean <- apply(ko, 1, mean)
print(head(wt.mean))
print(head(ko.mean))
# En apply, 1 indica que se aplica mean por filas, usando las tres réplicas.
limit <- max(wt.mean, ko.mean)
print(limit)
print(c(WT = max(wt.mean), KO = max(ko.mean)))
print(c(WT = names(which.max(wt.mean)), KO = names(which.max(ko.mean))))
# Media más alta: KO, 15.14874, sonda 100213_f_at.
# Máxima media WT: 15.09946, sonda 94766_at. Son medias en escala log2.
guardar_grafico("05_dispersion_medias", function() {
  limites <- range(wt.mean, ko.mean)
  plot(wt.mean, ko.mean, xlab = "Media WT (log2)", ylab = "Media KO (log2)",
       main = "GSE5583: comparación de medias", xlim = limites, ylim = limites,
       pch = 16, cex = 0.4, col = adjustcolor("black", alpha.f = 0.3))
  grid()
  abline(0, 1, col = "red", lwd = 2)
})
# grid() añade la cuadrícula. abline(0, 1) dibuja y = x.
# Por encima de la diagonal la media es mayor en KO; por debajo, en WT.
# Los límites incluyen también las medias negativas de log2, si las hay.

# 6. DIFERENCIA DE MEDIAS
diff.mean <- wt.mean - ko.mean
guardar_grafico("06_diferencias_medias", function() {
  hist(diff.mean, breaks = 60, col = "gray75", border = "white",
       main = "GSE5583: diferencias de medias",
       xlab = "Media log2(WT) - media log2(KO)", ylab = "Frecuencia")
})
# diff.mean es el log2 del cociente de las medias geométricas WT/KO.
# Un valor positivo indica mayor expresión en WT; uno negativo, en KO.

# 7. TEST POR SONDA
# El enunciado dice "SIN TRANSFORMAR", pero su código define wt y ko a
# partir de data2. Aquí seguimos ese código: el test usa valores log2.
# La transformación se usa tanto para las medias como para el test.
# Cada sonda tiene 3 valores WT y 3 valores KO: seis observaciones en total.
# t.test usa por defecto Welch bilateral para grupos independientes;
# no exige varianzas iguales. Con n = 3 por grupo es difícil comprobar
# la normalidad de las réplicas y los resultados requieren cautela.
pvalue <- numeric(nrow(data))
tstat <- numeric(nrow(data))
for (i in seq_len(nrow(data))) {
  x <- wt[i, ]
  y <- ko[i, ]
  test <- t.test(x, y)
  pvalue[i] <- test$p.value
  tstat[i] <- unname(test$statistic)
}
names(pvalue) <- names(tstat) <- rownames(data)
stopifnot(length(pvalue) == nrow(data), all(is.finite(pvalue)))
print(head(pvalue))
print(length(pvalue))
# Se obtienen 12488 p-valores y 12488 estadísticos, en vectores numéricos.
# No se aplica la transformación inversa a los p-valores.

# 8. HISTOGRAMAS DE LOS P-VALORES
guardar_grafico("07_histograma_pvalores", function() {
  hist(pvalue, breaks = 50, col = "gray75", border = "white",
       main = "GSE5583: p-valores", xlab = "p-valor", ylab = "Frecuencia")
})
guardar_grafico("08_histograma_menos_log10p", function() {
  hist(-log10(pvalue), breaks = 50, col = "gray75", border = "white",
       main = "GSE5583: p-valores transformados",
       xlab = "-log10(p-valor)", ylab = "Frecuencia")
})
# -log10(p) convierte p-valores pequeños en valores grandes y facilita
# distinguirlos. p = 0.01 se representa como 2; p = 0.001, como 3.
# Es una transformación para representar p, no una corrección estadística.

# 9. VOLCANO PLOT Y FILTROS
diff.mean_cutoff <- 2
pvalue_cutoff <- 0.01
filter_by_diff.mean <- abs(diff.mean) >= diff.mean_cutoff
filter_by_pvalue <- pvalue <= pvalue_cutoff
filter_combined <- filter_by_diff.mean & filter_by_pvalue
filtered <- data[filter_combined, , drop = FALSE]
print(dim(data[filter_by_diff.mean, , drop = FALSE])) # 210 x 6
print(dim(data[filter_by_pvalue, , drop = FALSE]))    # 429 x 6
print(dim(filtered))                              # 42 x 6
print(head(filtered))
# Cumplen ambos criterios 42 sondas. & exige que se cumplan los dos filtros.
# |diferencia log2| >= 2 equivale a un cambio de al menos cuatro veces.
# Son p-valores sin ajustar, como pide el taller. Como información adicional:
pvalue_BH <- p.adjust(pvalue, method = "BH")
print(sum(filter_by_diff.mean & pvalue_BH <= pvalue_cutoff)) # 0
# Ninguna cumple simultáneamente el umbral de cambio y p ajustado BH <= 0.01.
# Por tanto, no se atribuye control de FDR a la lista de 42 del ejercicio.

volcano_base <- function(titulo) {
  plot(diff.mean, -log10(pvalue), pch = 16, cex = 0.5,
       col = "gray65", main = titulo,
       xlab = "Media log2(WT) - media log2(KO)", ylab = "-log10(p-valor)")
  abline(v = c(-diff.mean_cutoff, diff.mean_cutoff),
         col = "#0072B2", lwd = 1.5, lty = 2)
  abline(h = -log10(pvalue_cutoff), col = "#009E73", lwd = 1.5, lty = 2)
}
guardar_grafico("09_volcano_umbrales", function() {
  volcano_base("GSE5583: umbrales de selección")
})
guardar_grafico("10_volcano_seleccionados", function() {
  volcano_base("GSE5583: 42 sondas seleccionadas")
  points(diff.mean[filter_combined], -log10(pvalue[filter_combined]),
         col = "red", pch = 16, cex = 0.7)
})
sobreexpresados_KO <- filter_combined & diff.mean < 0
reprimidos_KO <- filter_combined & diff.mean > 0
print(c(mayor_en_KO = sum(sobreexpresados_KO),
        menor_en_KO = sum(reprimidos_KO))) # 26 y 16
guardar_grafico("11_volcano_direccion", function() {
  volcano_base("GSE5583: dirección del cambio respecto a WT")
  points(diff.mean[sobreexpresados_KO], -log10(pvalue[sobreexpresados_KO]),
         col = "red", pch = 16, cex = 0.7)
  points(diff.mean[reprimidos_KO], -log10(pvalue[reprimidos_KO]),
         col = "blue", pch = 16, cex = 0.7)
  legend("topright", c("Mayor en KO: 26", "Menor en KO: 16"),
         col = c("red", "blue"), pch = 16, bty = "n")
})
# ¿Por qué parecen al revés? Hemos restado WT - KO, no KO - WT.
# Los 26 sobreexpresados en KO están a la izquierda (diferencia negativa).
# Los 16 reprimidos en KO están a la derecha (diferencia positiva).

# 10. HEATMAP BASE Y CAMBIO DE COLORES
# Conservamos la definición del taller: filtered contiene intensidades
# originales. Los dendrogramas siguientes usan esos valores; los colores
# del heatmap se estandarizan por fila. El test anterior sigue usando log2.
rowv <- as.dendrogram(hclust(as.dist(1 - cor(t(filtered)))))
colv <- as.dendrogram(hclust(as.dist(1 - cor(filtered))))
guardar_grafico("12_heatmap_correlacion", function() {
  heatmap(filtered, Rowv = rowv, Colv = colv, scale = "row",
          cexCol = 0.7, labRow = NA, margins = c(9, 4),
          main = "GSE5583: heatmap con distancia por correlación")
})
guardar_grafico("13_heatmap_por_defecto", function() {
  heatmap(filtered, cexRow = 0.45, cexCol = 0.7, margins = c(9, 8),
          main = "GSE5583: heatmap con agrupamiento por defecto")
}, alto = 1800)
guardar_grafico("14_heatmap_colores_hcl", function() {
  heatmap(filtered, Rowv = rowv, Colv = colv, scale = "row",
          col = hcl.colors(100, "Blue-Red 3"),
          cexCol = 0.7, labRow = NA, margins = c(9, 4),
          main = "GSE5583: cambio de paleta con hcl.colors")
})
# filtered: matriz de expresión; Rowv: dendrograma de sondas;
# Colv: dendrograma de muestras; cexCol: tamaño relativo de sus etiquetas;
# labRow = NA: oculta los nombres de las filas; col: paleta de colores.
# scale = "row": resta la media y divide por la desviación típica de cada
# fila. Así los colores muestran expresión relativa dentro de cada sonda.
# heatmap(filtered) recalcula los árboles con distancia euclídea y sus
# opciones por defecto; por eso no tiene por qué conservar el mismo orden.

# 11. HEATMAP AVANZADO CON gplots Y RColorBrewer
paleta <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(256)
dibujar_heatmap <- function(colores = paleta) {
  gplots::heatmap.2(
    filtered, Rowv = rowv, Colv = colv, scale = "row",
    col = colores, trace = "none", density.info = "none",
    ColSideColors = colores_grupo, labCol = etiquetas,
    cexCol = 0.8, cexRow = 0.5, margins = c(7, 9),
    key.title = "Expresión relativa", key.xlab = "Z por sonda",
    keysize = 1.3, symkey = TRUE,
    main = "GSE5583: 42 sondas seleccionadas"
  )
}
guardar_grafico("15_heatmap_gplots", dibujar_heatmap, alto = 1900)
# Azul: expresión por debajo de la media de esa sonda; rojo: por encima.
# La barra superior identifica WT (azul) y KO (naranja), con independencia
# del orden en que el agrupamiento coloque las columnas.

# Variante de colores propuesta por el taller (rojo/verde).
guardar_grafico("16_heatmap_rojo_verde", function() {
  dibujar_heatmap(gplots::redgreen(75))
}, alto = 1900)
# Para la figura principal se usa azul/rojo para facilitar su lectura.

# Guardar el heatmap en PDF sin cerrar el panel Plots de RStudio.
local({
  pdf(file.path(carpeta_resultados, "GSE5583_DE_Heatmap.pdf"),
      width = 10, height = 11, useDingbats = FALSE)
  on.exit(dev.off())
  dibujar_heatmap()
})
# Dejar como último gráfico de Plots el heatmap principal.
if (interactive()) dibujar_heatmap()

# 12. EXPORTAR TABLAS Y RESUMEN
# ID se guarda como columna explícita para que el encabezado quede alineado.
tabla_filtrada <- data.frame(ID = rownames(filtered), filtered,
                             row.names = NULL, check.names = FALSE)
write.table(tabla_filtrada,
            file.path(carpeta_resultados, "GSE5583_DE.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
resultados <- data.frame(
  ID = rownames(data), Media_WT_log2 = wt.mean, Media_KO_log2 = ko.mean,
  Diferencia_WT_menos_KO = diff.mean, Cambio_KO_sobre_WT = 2^(-diff.mean),
  Estadistico_t = tstat, p_value = pvalue, p_ajustado_BH = pvalue_BH,
  Seleccion_taller = filter_combined,
  Direccion_KO = ifelse(diff.mean < 0, "Mayor", ifelse(diff.mean > 0, "Menor", "Igual")),
  row.names = NULL
)
write.table(resultados, file.path(carpeta_resultados, "GSE5583_resultados_completos.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(resultados[filter_combined, ],
            file.path(carpeta_resultados, "GSE5583_seleccion_con_estadisticas.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

resumen <- c(
  "TALLER R - RESULTADOS GSE5583",
  sprintf("Sondas: %d. Muestras: %d (3 WT y 3 HDAC1 KO).", nrow(data), ncol(data)),
  sprintf("Máxima media WT en log2: %.5f (%s).", max(wt.mean), names(which.max(wt.mean))),
  sprintf("Máxima media KO en log2: %.5f (%s).", max(ko.mean), names(which.max(ko.mean))),
  sprintf("Sondas con |diferencia log2| >= 2: %d.", sum(filter_by_diff.mean)),
  sprintf("Sondas con p <= 0.01 sin ajustar: %d.", sum(filter_by_pvalue)),
  sprintf("Sondas que cumplen ambos criterios: %d.", nrow(filtered)),
  sprintf("Mayor expresión en KO respecto a WT: %d.", sum(sobreexpresados_KO)),
  sprintf("Menor expresión en KO respecto a WT: %d.", sum(reprimidos_KO)),
  sprintf("Con diferencia >= 2 en valor absoluto y p BH <= 0.01: %d.",
          sum(filter_by_diff.mean & pvalue_BH <= pvalue_cutoff)),
  "",
  "El test es Welch bilateral sobre log2, siguiendo el código del taller.",
  "El comentario original 'SIN TRANSFORMAR' contradice su definición de wt/ko.",
  "La selección de 42 emplea p-valores sin ajustar; no es una lista con control de FDR.",
  "La unidad analizada es la sonda del microarray, no necesariamente un gen único.",
  "El heatmap sigue filtered = data[filter_combined, ] del taller: intensidades",
  "originales, distancia 1-cor para los árboles y colores estandarizados por fila.",
  "",
  "Fuentes:",
  "TallerR.R proporcionado como referencia.",
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583",
  "https://gist.github.com/ahmedmoustafa/5520526",
  "https://stat.ethz.ch/R-manual/R-devel/library/stats/html/t.test.html",
  "https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html",
  "https://cran.r-project.org/web/packages/gplots/refman/gplots.html",
  paste("Tabla utilizada:", url_datos),
  paste("MD5 de la tabla:", unname(tools::md5sum(archivo)))
)
writeLines(resumen, file.path(carpeta_resultados, "Resumen_resultados.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(carpeta_resultados, "Sesion_R.txt"))
cat(paste(resumen[1:11], collapse = "\n"), "\n")
message("Taller terminado. Resultados guardados en: ", carpeta_resultados)
