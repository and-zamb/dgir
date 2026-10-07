# ---- 0. PAQUETES Y CONEXIÓN -------------------------------------------------
paquetes <- c("reticulate", "leaflet", "sf", "jsonlite")
faltan <- paquetes[!paquetes %in% rownames(installed.packages())]
if (length(faltan) > 0) install.packages(faltan)
library(reticulate)
library(leaflet)
use_virtualenv("ee-env", required = TRUE)
ee <- import("ee")
# 1. Iniciar sesión (solo esta vez)
ee$Authenticate(auth_mode = "localhost")
# ee$Authenticate()                                # solo la primera vez
ee$Initialize(project = "curso-big-data-501201")
# =============================================================================
# ---- 1. PARÁMETROS (lo único que normalmente hay que editar) ----------------
# =============================================================================

# Área de estudio ------------------------------------------------------------
# nivel = "provincia" -> nombres de provincia (ADM1 de GAUL)
# nivel = "canton"    -> nombres de cantón; si el nombre se repite en otra
#                        provincia, usa "Provincia/Canton" (ej. "Guayas/Daule")
# nivel = "shapefile" -> usa ruta_shp (y campo_nombre si quieres una región
#                        por cada valor de ese campo)
nivel<-"provincia"
regiones<-c("Guayas")
ruta_shp<-NULL #ej."C:/Users/andrea.zambrano/Documents/limites/cuenca_guayas.shp"
campo_nombre <- NULL   # ej. "DPA_DESCAN"; NULL = todo el shp como una región
# Periodos ---------------------------------------------------------------------
# Desde la pérdida de Sentinel-1B (dic. 2021) la revisita en Ecuador suele ser
# de 12 días: usa ventanas de al menos 12-15 días para asegurar imágenes.
# Las fechas de abajo son solo de ejemplo (época seca vs. época lluviosa).
before_start <- "2024-11-01"
before_end   <- "2024-12-15"
after_start  <- "2025-03-01"
after_end    <- "2025-03-31"
#Parametros SAR
polarization         <- "VV"
pass_direction       <- "AUTO"   # "AUTO", "ASCENDING" o "DESCENDING"
difference_threshold <- 1.25     # umbral sobre el cociente después/antes
umbral_agua_db       <- NULL     # opcional para VV: exigir además que el píxel
# "después" sea agua abierta, ej. -15 (dB).
# Reduce falsos positivos. NULL = desactivado.
smoothing_radius     <- 50       # filtro speckle (m)

