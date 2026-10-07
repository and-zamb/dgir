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
# Salidas ----------------------------------------------------------------------
exportar        <- TRUE                 # raster de inundación y población a Drive
exportar_vector <- FALSE                # polígonos SHP (lento en provincias grandes)
escala_vector   <- 30L                  # m; 10 es muy pesado para áreas grandes
carpeta_drive   <- "Inundaciones_SNGR"
crs_salida      <- "EPSG:32717"         # WGS84 / UTM 17S
archivo_csv     <- "resultados_inundacion.csv"
# =============================================================================
# ---- 2. FUNCIONES AUXILIARES (no editar) ------------------------------------
# =============================================================================

gaul1 <- ee$FeatureCollection("FAO/GAUL/2015/level1")$
  filter(ee$Filter$eq("ADM0_NAME", "Ecuador"))
gaul2 <- ee$FeatureCollection("FAO/GAUL/2015/level2")$
  filter(ee$Filter$eq("ADM0_NAME", "Ecuador"))

listar_regiones <- function(nivel = "provincia") {
  if (nivel == "provincia") {
    sort(unlist(gaul1$aggregate_array("ADM1_NAME")$distinct()$getInfo()))
  } else {
    prov <- unlist(gaul2$aggregate_array("ADM1_NAME")$getInfo())
    cant <- unlist(gaul2$aggregate_array("ADM2_NAME")$getInfo())
    sort(paste(prov, cant, sep = "/"))
  }
}

region_gaul <- function(nombre, nivel) {
  if (nivel == "provincia") {
    fc <- gaul1$filter(ee$Filter$eq("ADM1_NAME", nombre))
  } else {
    partes <- strsplit(nombre, "/")[[1]]
    if (length(partes) == 2) {
      fc <- gaul2$
        filter(ee$Filter$eq("ADM1_NAME", partes[1]))$
        filter(ee$Filter$eq("ADM2_NAME", partes[2]))
    } else {
      fc <- gaul2$filter(ee$Filter$eq("ADM2_NAME", nombre))
    }
  }
  if (fc$size()$getInfo() == 0) {
    stop("No se encontró '", nombre, "'. Revisa los nombres con listar_regiones(\"",
         nivel, "\").")
  }
  fc$geometry()
}

# Convierte un objeto sf (shapefile leído en R) en geometría de Earth Engine
sf_a_ee <- function(capa) {
  capa <- sf::st_zm(capa)
  capa <- sf::st_transform(capa, 4326)
  capa <- sf::st_make_valid(sf::st_union(capa))
  capa <- sf::st_simplify(capa, dTolerance = 30)   # aligera la geometría (m)
  tmp  <- tempfile(fileext = ".geojson")
  sf::st_write(sf::st_sf(geometry = capa), tmp, quiet = TRUE, delete_dsn = TRUE)
  gj <- jsonlite::read_json(tmp, simplifyVector = FALSE)
  ee$Geometry(gj$features[[1]]$geometry, NULL, FALSE)
}

construir_regiones <- function() {
  if (nivel == "shapefile") {
    if (is.null(ruta_shp)) stop("Define ruta_shp.")
    shp <- sf::st_read(ruta_shp, quiet = TRUE)
    if (is.null(campo_nombre)) {
      out <- list(sf_a_ee(shp))
      names(out) <- tools::file_path_sans_ext(basename(ruta_shp))
    } else {
      nombres <- unique(as.character(shp[[campo_nombre]]))
      if (length(regiones) > 0 && !identical(regiones, "")) {
        sel <- nombres[nombres %in% regiones]
        if (length(sel) > 0) nombres <- sel
      }
      out <- lapply(nombres, function(n) sf_a_ee(shp[shp[[campo_nombre]] == n, ]))
      names(out) <- nombres
    }
    out
  } else {
    out <- lapply(regiones, region_gaul, nivel = nivel)
    names(out) <- regiones
    out
  }
}

slug <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x))
}

dates <- function(imgcol) {
  rango <- imgcol$reduceColumns(ee$Reducer$minMax(), list("system:time_start"))
  ee$String("from ")$
    cat(ee$Date(rango$get("min"))$format("YYYY-MM-dd"))$
    cat(" to ")$
    cat(ee$Date(rango$get("max"))$format("YYYY-MM-dd"))
}

elegir_orbita <- function(base) {
  if (pass_direction != "AUTO") return(pass_direction)
  conteos <- sapply(c("ASCENDING", "DESCENDING"), function(d) {
    col <- base$filter(ee$Filter$eq("orbitProperties_pass", d))
    nb  <- col$filterDate(before_start, before_end)$size()$getInfo()
    na  <- col$filterDate(after_start, after_end)$size()$getInfo()
    min(nb, na)
  })
  if (all(conteos == 0)) return(NA_character_)
  names(conteos)[which.max(conteos)]
}

# Añade una imagen de Earth Engine como capa de leaflet
capa_ee <- function(mapa, imagen, vis, nombre) {
  mid <- imagen$getMapId(vis)
  addTiles(mapa, urlTemplate = mid$tile_fetcher$url_format, group = nombre)
}


# =============================================================================
# ---- 3. FUNCIÓN PRINCIPAL ----------------------------------------------------
# =============================================================================

mapear_inundacion <- function(nombre, aoi) {
  
  message("\n=== ", nombre, " ===")
  
  # -- Selección de datos Sentinel-1 GRD --
  base <- ee$ImageCollection("COPERNICUS/S1_GRD")$
    filter(ee$Filter$eq("instrumentMode", "IW"))$
    filter(ee$Filter$listContains("transmitterReceiverPolarisation", polarization))$
    filter(ee$Filter$eq("resolution_meters", 10L))$
    filterBounds(aoi)$
    select(polarization)
  
  orbita <- elegir_orbita(base)
  if (is.na(orbita)) {
    warning(nombre, ": no hay imágenes antes y después en ninguna órbita. Amplía las fechas.")
    return(NULL)
  }
  
  collection        <- base$filter(ee$Filter$eq("orbitProperties_pass", orbita))
  before_collection <- collection$filterDate(before_start, before_end)
  after_collection  <- collection$filterDate(after_start, after_end)
  
  n_before <- before_collection$size()$getInfo()
  n_after  <- after_collection$size()$getInfo()
  message("Órbita: ", orbita, " | imágenes antes: ", n_before, " | después: ", n_after)
  if (n_before == 0 || n_after == 0) {
    warning(nombre, ": colección vacía con la órbita ", orbita, ". Amplía las fechas.")
    return(NULL)
  }
  
  # -- Mosaico, recorte y filtro speckle --
  before <- before_collection$mosaic()$clip(aoi)
  after  <- after_collection$mosaic()$clip(aoi)
  
  before_filtered <- before$focal_mean(radius = smoothing_radius,
                                       kernelType = "circle", units = "meters")
  after_filtered  <- after$focal_mean(radius = smoothing_radius,
                                      kernelType = "circle", units = "meters")
  
  # -- Detección de cambios --
  difference        <- after_filtered$divide(before_filtered)
  difference_binary <- difference$gt(difference_threshold)
  
  if (!is.null(umbral_agua_db)) {
    difference_binary <- difference_binary$And(after_filtered$lt(umbral_agua_db))
  }
  
  # -- Refinamiento --
  # Agua permanente (> 10 meses/año)
  swater      <- ee$Image("JRC/GSW1_4/GlobalSurfaceWater")$select("seasonality")
  swater_mask <- swater$gte(10)$updateMask(swater$gte(10))
  flooded_mask <- difference_binary$where(swater_mask, 0)
  flooded      <- flooded_mask$updateMask(flooded_mask)
  
  # Conectividad: elimina grupos de 8 píxeles o menos
  connections <- flooded$connectedPixelCount()
  flooded     <- flooded$updateMask(connections$gte(8L))
  
  # Pendiente > 5 %
  DEM   <- ee$Image("WWF/HydroSHEDS/03VFDEM")
  slope <- ee$Terrain$slope(DEM)
  flooded <- flooded$updateMask(slope$lt(5))
  
  # -- Área inundada --
  flood_stats <- flooded$select(polarization)$
    multiply(ee$Image$pixelArea())$
    reduceRegion(reducer = ee$Reducer$sum(), geometry = aoi, scale = 10L,
                 maxPixels = 1e13, bestEffort = TRUE)
  flood_area_ha <- flood_stats$getNumber(polarization)$divide(10000)$round()
  
  # -- Población expuesta (GHSL 2020, 100 m) --
  population_count <- ee$Image("JRC/GHSL/P2023A/GHS_POP/2020")$
    select("population_count")$clip(aoi)
  GHSLprojection <- population_count$projection()
  flooded_res1   <- flooded$reproject(crs = GHSLprojection)
  population_exposed <- population_count$
    updateMask(flooded_res1)$
    updateMask(population_count)
  pop_stats <- population_exposed$reduceRegion(
    reducer = ee$Reducer$sum(), geometry = aoi, scale = 100L,
    maxPixels = 1e13, bestEffort = TRUE)
  number_pp_exposed <- pop_stats$getNumber("population_count")$round()
  
  # -- Cobertura del suelo MODIS (500 m) --
  LC <- ee$ImageCollection("MODIS/061/MCD12Q1")$
    filterDate("2014-01-01", after_end)$
    sort("system:index", FALSE)$
    select("LC_Type1")$
    first()$
    clip(aoi)
  MODISprojection <- LC$projection()
  flooded_res     <- flooded$reproject(crs = MODISprojection)
  
  # Cultivos (clases 12 y 14)
  cropmask <- LC$eq(12L)$Or(LC$eq(14L))
  cropland <- LC$updateMask(cropmask)
  cropland_affected <- flooded_res$updateMask(cropland)
  crop_stats <- cropland_affected$multiply(ee$Image$pixelArea())$
    reduceRegion(reducer = ee$Reducer$sum(), geometry = aoi, scale = 500L,
                 maxPixels = 1e13, bestEffort = TRUE)
  crop_area_ha <- crop_stats$getNumber(polarization)$divide(10000)$round()
  
  # Urbano (clase 13)
  urban <- LC$updateMask(LC$eq(13L))
  urban_affected <- urban$updateMask(flooded_res)$updateMask(urban)
  urban_stats <- urban_affected$multiply(ee$Image$pixelArea())$
    reduceRegion(reducer = ee$Reducer$sum(), geometry = aoi, scale = 500L,
                 maxPixels = 1e13, bestEffort = TRUE)
  urban_area_ha <- urban_stats$getNumber("LC_Type1")$divide(10000)$round()
  
  MODIS_date <- ee$String(LC$get("system:index"))$slice(0L, 4L)
  
  # -- Todos los resultados en una sola consulta al servidor --
  res <- ee$Dictionary(list(
    fechas_antes   = dates(before_collection),
    fechas_despues = dates(after_collection),
    inundacion_ha  = flood_area_ha,
    poblacion      = number_pp_exposed,
    cultivos_ha    = crop_area_ha,
    urbano_ha      = urban_area_ha,
    anio_modis     = MODIS_date
  ))$getInfo()
  
  num <- function(x) if (is.null(x)) 0 else as.numeric(x)
  
  fila <- data.frame(
    region              = nombre,
    polarizacion        = polarization,
    orbita              = orbita,
    imagenes_antes      = n_before,
    imagenes_despues    = n_after,
    fechas_antes        = res$fechas_antes,
    fechas_despues      = res$fechas_despues,
    inundacion_ha       = num(res$inundacion_ha),
    poblacion_expuesta  = num(res$poblacion),
    cultivos_afect_ha   = num(res$cultivos_ha),
    urbano_afect_ha     = num(res$urbano_ha),
    anio_modis          = res$anio_modis,
    stringsAsFactors    = FALSE
  )
  print(fila[, c("region", "inundacion_ha", "poblacion_expuesta",
                 "cultivos_afect_ha", "urbano_afect_ha")], row.names = FALSE)
  
  # -- Exportaciones a Google Drive --
  tareas <- list()
  if (exportar) {
    pref <- paste0(slug(nombre), "_", polarization, "_", gsub("-", "", after_start))
    
    tareas$raster <- ee$batch$Export$image$toDrive(
      image = flooded$toByte(), description = paste0("Flood_raster_", pref),
      folder = carpeta_drive, fileNamePrefix = paste0("flooded_", pref),
      region = aoi, scale = 10L, crs = crs_salida, maxPixels = 1e13)
    
    tareas$poblacion <- ee$batch$Export$image$toDrive(
      image = population_exposed, description = paste0("Exposed_pop_", pref),
      folder = carpeta_drive, fileNamePrefix = paste0("population_exposed_", pref),
      region = aoi, scale = 100L, crs = crs_salida, maxPixels = 1e13)
    
    if (exportar_vector) {
      flooded_vec <- flooded$reduceToVectors(
        scale = escala_vector, geometryType = "polygon", geometry = aoi,
        eightConnected = FALSE, bestEffort = TRUE, tileScale = 4L)
      tareas$vector <- ee$batch$Export$table$toDrive(
        collection = flooded_vec, description = paste0("Flood_vector_", pref),
        folder = carpeta_drive, fileFormat = "SHP",
        fileNamePrefix = paste0("flooded_vec_", pref))
    }
    for (t in tareas) t$start()
    message("Exportaciones iniciadas: ", length(tareas), " (carpeta Drive: ", carpeta_drive, ")")
  }
  
  list(
    resultados = fila,
    tareas     = tareas,
    aoi        = aoi,
    capas      = list(
      before = before_filtered, after = after_filtered, difference = difference,
      flooded = flooded, population = population_count,
      population_exposed = population_exposed, LC = LC,
      cropland = cropland, cropland_affected = cropland_affected,
      urban = urban, urban_affected = urban_affected
    )
  )
}


# ---- Mapa interactivo (reemplaza Map.addLayer, panel y leyenda de GEE) -------
ver_mapa <- function(salida) {
  cp <- salida$capas
  b  <- salida$aoi$bounds(1)$coordinates()$getInfo()[[1]]
  lons <- sapply(b, `[[`, 1); lats <- sapply(b, `[[`, 2)
  
  pop_vis <- list(min = 0, max = 50, palette = list("yellow", "orange", "red"))
  lc_vis  <- list(min = 1, max = 17, palette = list(
    "05450a", "086a10", "54a708", "78d203", "009900", "c6b044", "dcd159",
    "dade48", "fbff13", "b6ff05", "27ff87", "c24f44", "a5a5a5", "ff6d4c",
    "69fff8", "f9ffa4", "1c0dff"))
  
  m <- leaflet() %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satélite") %>%
    addProviderTiles(providers$OpenStreetMap,     group = "OSM") %>%
    fitBounds(min(lons), min(lats), max(lons), max(lats))
  
  m <- capa_ee(m, cp$before,     list(min = -25, max = 0), "Antes")
  m <- capa_ee(m, cp$after,      list(min = -25, max = 0), "Después")
  m <- capa_ee(m, cp$difference, list(min = 0, max = 2),   "Diferencia")
  m <- capa_ee(m, cp$population, list(min = 0, max = 50,
                                      palette = list("060606", "337663", "337663", "ffffff")),
               "Densidad de población")
  m <- capa_ee(m, cp$LC,       lc_vis, "Cobertura MODIS")
  m <- capa_ee(m, cp$cropland, list(min = 0, max = 14, palette = list("30b21c")), "Cultivos")
  m <- capa_ee(m, cp$urban,    list(min = 0, max = 13, palette = list("grey")),   "Urbano")
  m <- capa_ee(m, cp$flooded,  list(palette = list("0000FF")), "Inundación")
  m <- capa_ee(m, cp$population_exposed, pop_vis, "Población expuesta")
  m <- capa_ee(m, cp$cropland_affected,
               list(min = 0, max = 14, palette = list("30b21c")), "Cultivos afectados")
  m <- capa_ee(m, cp$urban_affected,
               list(min = 0, max = 13, palette = list("grey")), "Urbano afectado")
  
  grupos <- c("Antes", "Después", "Diferencia", "Densidad de población",
              "Cobertura MODIS", "Cultivos", "Urbano", "Inundación",
              "Población expuesta", "Cultivos afectados", "Urbano afectado")
  
  r <- salida$resultados
  m %>%
    addLayersControl(baseGroups = c("Satélite", "OSM"), overlayGroups = grupos,
                     options = layersControlOptions(collapsed = TRUE)) %>%
    hideGroup(c("Antes", "Diferencia", "Densidad de población",
                "Cobertura MODIS", "Cultivos", "Urbano")) %>%
    addLegend("bottomright", colors = c("#0000FF", "#30b21c", "grey"),
              labels = c("Áreas potencialmente inundadas", "Cultivos afectados",
                         "Urbano afectado"), title = "Leyenda") %>%
    addLegend("bottomright",
              pal = colorNumeric(c("yellow", "orange", "red"), domain = c(0, 50)),
              values = c(0, 50), title = "Población expuesta (hab/píxel)") %>%
    addControl(html = paste0(
      "<b>", r$region, " (", r$polarizacion, ", ", r$orbita, ")</b><br>",
      "Después: ", r$fechas_despues, "<br>",
      "Inundación: <b>", format(r$inundacion_ha, big.mark = ","), " ha</b><br>",
      "Población expuesta: <b>", format(r$poblacion_expuesta, big.mark = ","), "</b><br>",
      "Cultivos afectados: <b>", format(r$cultivos_afect_ha, big.mark = ","), " ha</b><br>",
      "Urbano afectado: <b>", format(r$urbano_afect_ha, big.mark = ","), " ha</b><br>",
      "<small>MODIS ", r$anio_modis, " · GHSL 2020 · sin validación de campo</small>"),
      position = "bottomleft")
}


# =============================================================================
# ---- 4. EJECUCIÓN -----------------------------------------------------------
# =============================================================================

lista_regiones <- construir_regiones()

salidas <- lapply(names(lista_regiones), function(n) {
  tryCatch(mapear_inundacion(n, lista_regiones[[n]]),
           error = function(e) { warning(n, ": ", conditionMessage(e)); NULL })
})
names(salidas) <- names(lista_regiones)
salidas <- Filter(Negate(is.null), salidas)

if (length(salidas) > 0) {
  tabla <- do.call(rbind, lapply(salidas, `[[`, "resultados"))
  rownames(tabla) <- NULL
  print(tabla)
  write.csv(tabla, archivo_csv, row.names = FALSE, fileEncoding = "UTF-8")
  message("Resultados guardados en: ", normalizePath(archivo_csv))
  
  # Mapa de la primera región (para otra: ver_mapa(salidas[["Los Rios"]]))
  print(ver_mapa(salidas[[1]]))
}

# Estado de las exportaciones (vuelve a correr esta línea para actualizar):
# for (s in salidas) for (t in s$tareas) print(t$status()[c("description", "state")])
