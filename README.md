# Noise-induced cycles in the Togashi-Kaneko model with species-dependent degradation

Contains a package and scripts written in Julia used to generate the figures for the article located at:


## Structure

If you want to just see how the figures have been created at a high level, look at the `scripts` directory to see the scripts used to create each figure. This also contains formatting presets and functions for the figure aesthetics.

To see how the Gillespie simulations have been done, look at the package `TK2.jl` in `src`.

The timeseries figure is made using `Catalyst.jl`'s implementation of the Gillespie algorithm.

For longer simulations a specific Gillespie was written for this system in `custom_gillespie.jl` with the aggregated values collected in`outputhist.jl`. This was done to greatly save on memory and time over the `Catalyst` implementation by only saving binned histogram information of interest on the molecular counts.

`parameters` contains the simulation specifications for each figure. The simulations will likely take a long time so `T_end` can be reduced to obtain results more quickly.

![Noise-induced cycles figure](https://raw.githubusercontent.com/jeremyworsfold/Noise-Induced-cycles-TK2/main/figures/02-timeseries-and-cycles.png)