using ADAS, Plots

em_cr = ADAS.get_emission_rate(:cr; Z=0, model="llu");
plot(em_cr; ne=[1e19, 1e20], box=:on, palette=:viridis, legendfontsize=8, legend=:bottomright, legendcolumns=2, size=(800, 800))

em_c = ADAS.get_emission_rate(:c; Z=1, model="");
plot!(em_c; iwl = 25, color=:red, xscale=:linear, ylims=[1e-24,1e-18], xlims=[10.0,100.0], xticks=[10.0,20.0,0,40.0,60.0,80.0,100.0], ne=[1e19, 1e20], box=:on, legendfontsize=8, legend=:bottomright, legendcolumns=2, size=(800, 800))

findall(wl -> (wl > 5120 && wl < 5180), em_c.wls)
em_h = ADAS.get_emission_rate(:h; Z=0, model="pju");
plot!(em_h; iwl = 19, color=:blue, xscale=:linear, ylims=[1e-24,1e-18], xlims=[10.0,100.0], xticks=[10.0,20.0,0,40.0,60.0,80.0,100.0], ne=[1e19, 1e20], box=:on, legendfontsize=8, legend=:bottomright, legendcolumns=2, size=(800, 800))

findall(wl -> (wl > 4250 && wl < 4350), em_h.wls)