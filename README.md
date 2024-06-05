[![CI](https://github.com/ProjectTorreyPines/ADAS.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/ProjectTorreyPines/ADAS.jl/actions/workflows/CI.yml)

[![runtests_FUSE](https://github.com/ProjectTorreyPines/ADAS.jl/actions/workflows/runtest_FUSE.yml/badge.svg?branch=master)](https://github.com/ProjectTorreyPines/ADAS.jl/actions/workflows/runtest_FUSE.yml)

# ADAS.jl

This package provides ADAS ionization, recombination, radiation, ... rates and some functions to retrieve cooling rates, Zeff and effective charge state.

# Usage

Basic usage:
```julia
using ADAS
# retrieve data for carbon
data = retrieve_ADAS_data("C"; year="latest", type="scd", metastable=false)
# show  available ADAS data
show_ADAS_data()
# show  available ADAS data for C
show_ADAS_data(:C)
# show adf11 format available
show_adf11_types()
```

Additional features
```julia
using ADAS
using Plots
# retrieve effective charge state for Kr
zeff = ADAS.get_Zeff(:Kr)
plot(zeff)
# retrieve cooling rate for W
cr = ADAS.get_cooling_rates(:W)
plot(cr)
# retrieve abundance fraction
af = ADAS.get_abundance_fraction(:Ne)
plot(af, Te=1:1.0:10000.0)
```



