using ADAS, Plots

em_cr = ADAS.get_emission_rate(:cr; Z=0, model="llu");
plot(em_cr; ne=[1e19, 1e20], box = :on, palette=:viridis, legendfontsize=8, legend=:bottomright, legendcolumns=2, size=(800, 800))