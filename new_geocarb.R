library(pacman)

p_load(tidyverse)
p_load(stringr)
p_load(scales)
p_load(xml2)
p_load(reticulate)

use_condaenv("python3")

column_names = c(
  year = "year",
  tco2 = "co2.total",
  pCO2 = "co2.atmos",
  alk = "alkalinity.ocean",
  d13Cocn = "delta.13C.ocean",
  d13Catm = "delta.13C.atmos",
  CO3 = "carbonate.ocean",
  WeatC = "carbonate.weathering",
  WeatS = "silicate.weathering",
  TotW = "total.weathering",
  BurC = "carbon.burial",
  Degas = "degassing.rate",
  Emiss = "emissions",
  Tatm = "temp.atmos",
  Tocn = "temp.ocean"
)

column_descr = c(
  year = "year",
  tco2 = "Total CO2",
  pCO2 = "Atmospheric CO2",
  alk = "Ocean Alkalinity",
  d13Cocn = "Delta 13C (ocean)",
  d13Catm = "Delta 13C (atmosphere)",
  CO3 = "Ocean carbonate concentration",
  WeatC = "Carbonate Weathering Rate",
  WeatS = "Silicate Weathering Rate",
  TotW = "Total Weathering Rate",
  BurC = "Carbon Burial Rate",
  Degas = "Degassing Rate",
  Emiss = "CO2 Emissions",
  Tatm = "Atmospheric Temperature",
  Tocn = "Ocean Temperature"
)

columns = tibble(
  index = names(column_names),
  name = column_names) %>%
  full_join(
    tibble(index = names(column_descr),
           description = column_descr),
    by = "index"
  )


load_geocarb = function(python_script = "geocarb_varco2") {
  path = dirname(python_script)
  module = basename(python_script) %>% str_split("\\.", n = 2) %>%
    simplify() %>% head(1)
  geocarb_module <- import_from_path(module, path)
  assign("geocarb_module", geocarb_module, envir = globalenv())
}

run_geocarb = function(filename,
                       co2_spike = 0,
                       co2_emissions = list(
                         0,
                         seq(1, 100, 1), # Go from 1 to 100 GT in steps of 1.
                         0
                       ),
                       periods = c(5E6, 100, 2E6),
                       time_steps = c(50, 1, 50),
                       degas = 7.5,
                       plants = TRUE,
                       land_area = 1,
                       delta_t2x = 3.0,
                       million_years_ago = 0,
                       mean_latitude_continents = 30,
                       start_recording = -2E6) {
  if (! exists("geocarb_module", envir = globalenv()))
    load_geocarb()

  gc <- geocarb_module$geocarb(co2_spike, degas,
                               time_steps, periods,
                               delta_t2x,
                               million_years_ago, mean_latitude_continents,
                               plants, land_area,
                               (periods[1] + start_recording - time_steps[1]),
                               co2_emissions)

  geocarb_module$save(gc, filename)
}

read_geocarb <- function(filename) {
  f <- file(filename,"r")
  lines <- readLines(f, warn=F)
  close(f)
  lines %>% str_trim() %>% str_replace_all('[ \t]+', ',') %>%
    str_c(collapse = "\n") %>% read_csv() -> df
  names(df) <- column_names[names(df)]
  invisible(df)
}
