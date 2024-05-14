#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
Company: General Atomics
ADAS.jl (c) 2024
=#

using RecipesBase
using Colors
@recipe function rplot_dp(af::AbundanceFraction; ne = 1e20, Te = [1:1:1000]..., ylims=[-0.1,1.2])
    colors = Colors.colormap("blues", length(af.Z))
    colors = distinguishable_colors(length(af.Z), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    #colors = map(col -> (red(col), green(col), blue(col)), cols)
    ne_ = ne .+ 0 .* Te
    for (i,Z_) in enumerate(af.Z)
    Z =  float.(Z_).+ 0 .* Te
    fZ = af.fZ.(Z,ne_, Te) 
    title := "Abundance fraction $(af.imp) nₑ = $(ne) m⁻³"
    @series begin 
        seriestype  := :line
        color := colors[i]
        linestyle := :solid
        linewidth := 2.0
        xlabel := "Te [eV]"
        xscale := :log10
        ylabel := "fraction"
        label := missing
        
        Te,fZ
    end
    @series begin 
        seriestype  := :path
        linestyle := :dash
        color := colors[i]
        linewidth := 1.0
        xlabel := "Te [eV]"
        xscale := :log10
        ylabel := "fraction"
        label := missing
        
        Temax = Te[argmax(fZ)]
        [Temax,Temax],[0.0,1.0]
    end
    @series begin 
        seriestype  := :scatter
        if iseven(i)
        series_annotations := Main.Plots.series_annotations(Main.Plots.text("$Z_+",:bottom,10, color=colors[i]))
        ymax=1.0
        else
        series_annotations := Main.Plots.series_annotations(Main.Plots.text("$Z_+",:top,10, color=colors[i]))  
        ymax = 0.0 
        end
        color := colors[i]
        xlabel := "Te [eV]"
        xscale := :log10
        ylabel := "fraction"
        label := missing
        markersize := 0
        Temax = Te[argmax(fZ)]
        [Temax],[ymax]
    end
    
end
end

@recipe function rplot_dp(rr::RadiationRates; ne = 1e20, Zmin = 0, Zmax = 1000,Te = [1:1:1000]... , ylims=[1e-35,1e-30]) 

    colors = distinguishable_colors(length(rr.Z), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    #colors = map(col -> (red(col), green(col), blue(col)), cols)
    ne_ = ne .+ 0 .* Te
    ZZ = filter(s->(s>=Zmin && s <= Zmax),rr.Z) 

    for (i,Z_) in enumerate(ZZ)
    Z =  float.(Z_).+ 0 .* Te
    rates = rr.rates.(Z,ne_, Te) 
    @series begin 
        seriestype  := :line
        color := colors[i]
        linestyle := :solid
        linewidth := 2.0
        xlabel := "Te [eV]"
        label := "$(rr.imp)[Z=$Z_,nₑ = $(ne) m⁻³]"
        xscale := :log10
        yscale := :log10
        ylabel := "radiation rate [W.cm^3]"
        #label := "$(af.imp) Z=$Z_"
        Te,rates
    end

    # @series begin 
    #     seriestype  := :scatter
    #     if iseven(i)
    #     series_annotations := Main.Plots.series_annotations(Main.Plots.text("$Z_+",:bottom,10, color=colors[i]))
    #     ymax=1.0
    #     else
    #     series_annotations := Main.Plots.series_annotations(Main.Plots.text("$Z_+",:top,10, color=colors[i]))  
    #     ymax = 0.0 
    #     end
    #     color := colors[i]
    #     xlabel := "Te [eV]"
    #     xscale := :log10
    #     ylabel := "fraction"
    #     #label := "$(af.imp) Z=$Z_"
    #     markersize := 0
    #     Temax = Te[argmax(fZ)]
    #     [Temax],[ymax]
    # end
    
end
end

@recipe function rplot_dp(Lz_data::Lz_ADAS; ne = 1e20, Te = [1:1:1000]..., ylims=[1e-35,1e-30])
    #colors = map(col -> (red(col), green(col), blue(col)), cols)
    Lz = Lz_data.Lztot(ne,Te)

    @series begin 
        seriestype  := :line
        linestyle := :solid
        linewidth := 2.0
        xlabel := "Te [eV]"
        xscale := :log10
        yscale := :log10
        title := "ADAS cooling rates"
        ylabel := "cooling rates [W.m^3]"
        label := "$(Lz_data.imp)"
        Te,Lz
    end
    
end

    @recipe function rplot_dp(zeff::Zeff; ne = 1e20, fraction = [0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02], Te = collect(1:1:1000), ylims=[1.0,8.00])
        #colors = map(col -> (red(col), green(col), blue(col)), cols)

        for f in fraction
        @series begin 
            seriestype  := :line
            linestyle := :solid
            linewidth := 2.0
            xlabel := "Te [eV]"
            xscale := :log10
            ylim := ylims
            title := "Z_{eff} C-E for $(zeff.af.imp)"
            ylabel := "Z_{eff}"
            label := "xi_{$(zeff.af.imp)} = $f"
            Te,zeff(f,ne,Te)
        end
    end
        
    end
    
