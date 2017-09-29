# New Geocarb Model

These files provide a new version of the GEOCARB model that allows you to simulate a CO2 emissions trajectory
(e.g., from burning fossil fuels).

To run this, you need to have python version 3.5 or higher installed, and you need the python `numpy` and `pandas` packages installed.
You also need the R `reticulate` package, but this should install automatically the first time you run the script.

To use this package:

* First execute `source("new_geocarb.R")` to load the script.
* Then execute `load_geocarb()` to load the python GEOCARB library.
* Then execute `run_geocarb` (details below).
    * You can specify two or more periods with different emissions profiles.
    * The first period is always the spinup period. By default it's 5 million years long, and we keep the last 2 million years.
    * The second and subsequent periods are the "simulation" periods, where we change the degassing rate, human CO2 emissions, etc.
* For any parameter that you don't want to change, you can just specify a number.
* For any parameter that you do want to change from one period to another, you pass a vector or list whose length is the number of
  periods.
* For instance, to have three periods, with a degassing rate of 7.5 in the spinup, 10.0 in period 1, and 12.5 in period 2, you would
  specify `degassing = c(7.5, 10.0, 12.5)`.
* You can specify human emissions of CO2 that changes over time with the `co2_emissions` parameter. This should be a `list`
  with one element for each period.
    * For any period with no human emissions, that list element should just be the number 0.
    * For any period with constant emissions, that list element should just be a number representing the emissions, in
      billion tons per year.
    * For any period with changing emissions, that list element should be a vector or list whose length is the number of
      time steps in that period (for a 100 year period with a 5 year time step, this would be length of 20).
    * To specify three periods with no emissions in the spinup, a steady increase in the first simulation period,
      from 1 billion tons per year in the first year, up to 100 billion tons in the 100th year,
      and no emissions in the third period, you'd write,
      `co2_emissions = list(0, seq(1, 100, 1), 0)`
* You specify the periods with two parameters:
    * `periods` is a list or vector of numbers representing the number of years in each period. If you want to use
      three periods, with a 5 million year spin-up, a 100-year first simulation period, and a 2 million year second
      simulation period, you'd say, `periods = c(5E6, 100, 2E6)`.
    * `time_steps` is a list or vector of numbers representing the number of years per time-step for each period.
      If you want to use three periods with a 50-year time step for the spinup,
      a 1 year time step for the first simulation period, and a 50 year period for the second simulation period,
      you'd say, `time_steps = c(50, 1, 50)`. 50 years is a good time step when you're simulating millions of years,
      but if you're simulating a short period with rapid changes, you might want to use a shorter time-step.

```
source("new_geocarb.R")

load_geocarb() # This loads the python geocarb library.

run_geocarb("geocarb_test_run.csv", # file to save the results in.
            co2_spike = 0,    # to put a single spike of 1000 GT CO2
                              # at the beginning of the first simulation
                              # period, you'd say
                              # co2_spike = c(0, 1000, 0)
            co2_emissions = list(
                0,                # no emissions in spinup
                seq(1, 100, 1),   # ramp up emissions for 100 years
                0                 # emissions stop abruptly at
                ),                # beginning of 2nd sim period.
             periods = c(5E6, 100, 2E6),  # 3 periods
             time_steps = c(50, 1, 50),   # time steps for the three periods.
             degas = 7.5,        # constant degassing rate for all three periods.
             plants = TRUE, land_area = 1,
             delta_t_2x = 3.0,   # Climate sensitivity, in degrees C per doubling of CO2.
             million_years_ago = 0,
             mean_latitude_continents = 30,
             start_recording = -2E6)  # Start recording data when
                                      # the spinup gets to -2 million years.

gc = load_geocarb("geocarb_test_run.csv")
```
