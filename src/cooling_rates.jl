struct Lz_ADAS{U,T,V,F,S}
    ne :: U
    Te :: U
    Lz :: T
    Lbrem :: T
    Lztot :: T
    Lz_grid :: V
    Lbrem_grid :: V
    Lztot_grid :: V
    ff :: F
    imp :: S
end


import Interpolations
function get_cooling_rates(imp::Union{String,Symbol})
    plt = retrieve_ADAS_data(imp;type="plt")
    prb = retrieve_ADAS_data(imp;type="prb")
    scd = retrieve_ADAS_data(imp;type="scd")
    acd = retrieve_ADAS_data(imp;type="acd")

    ndens = length(acd.data.axis.ne)
    ntemp = length(acd.data.axis.Te)
    ntrans = length(plt.data.rates)
    f = zeros(ntrans+1, ndens, ntemp)
    ff = zeros(ntrans+1,ndens, ntemp)


        f[1,:,:] .= 1.   #Initial neutral fraction - everthing scaled to this
        for itrans in 1:ntrans
            f[itrans+1,:,:] .= f[itrans,:,:] .*  scd.data.rates[itrans].values[:,:] ./ acd.data.rates[itrans].values[:,:]
        end

    for i in axes(ff,1)
        ff[i,:,:] .= f[i,:,:]./sum(f,dims=1)[1,:,:]
    end
    
    Lz = zeros(ndens,ntemp)
    Lz_brem = zeros(ndens,ntemp)

        for itrans in 1:ntrans
           #Radiation coefficients are for charge states 0-Zimp-1; bremmstrahlung rates are for charge states 1-Zimp
            Lz[:,:] .= Lz[:,:] .+ ff[itrans,:,:] .* plt.data.rates[itrans].values[:,:]
            Lz_brem[:,:] .= Lz_brem[:,:] .+ ff[itrans+1,:,:] .* prb.data.rates[itrans].values[:,:]
        end

    Lz_tot = Lz .+ Lz_brem
    ne = scd.data.axis.ne
    Te = scd.data.axis.Te
    Lz_tot_ = Interpolations.linear_interpolation((ne,Te),Lz_tot, extrapolation_bc = Interpolations.Flat())
    Lz_brem_ = Interpolations.linear_interpolation((ne,Te),Lz_brem, extrapolation_bc = Interpolations.Flat())
    Lz_ = Interpolations.linear_interpolation((ne,Te),Lz, extrapolation_bc = Interpolations.Flat())
    return  Lz_ADAS(plt.data.axis.ne,plt.data.axis.Te,Lz_,Lz_brem_,Lz_tot_,Lz,Lz_brem,Lz_tot ,ff,imp)
    end

struct AbundanceFraction{U,I,A,S}
    fZ :: U
    fZ_grid :: Array{Float64,3}
    Te :: Vector{Float64}
    ne :: Vector{Float64}
    Z :: Vector{Int64}
    scd :: I
    acd :: A
    imp :: S
end

struct RadiationRates{U,I,S}
    rates :: U
    rates_grid :: Array{Float64,3}
    Te :: Vector{Float64}
    ne :: Vector{Float64}
    Z :: Vector{Int64}
    plt :: I
    imp :: S
end

function get_radiation_rates(imp::Union{String, Symbol}; kw...)
    plt = retrieve_ADAS_data(imp;type="plt", kw...)
    rates = zeros(length(plt.data.rates),length(plt.data.axis.ne),length(plt.data.axis.Te))
    for Z in 1:length(plt.data.rates)
        rates[Z,:,:] =plt.data.rates[Z].values[:,:]
    end
    Te = plt.data.axis.Te
    ne = plt.data.axis.ne
    Z = collect(1:length(plt.data.rates))
    rates_ = Interpolations.linear_interpolation((float.(Z), ne,Te),rates, extrapolation_bc = Interpolations.Flat())
    return RadiationRates(rates_, rates, Te, ne, Z, plt, imp)
end
struct EffectiveCharge{Z,T,AF<:AbundanceFraction}
    t :: Z
    t_grid :: T
    af :: AF
end

function get_effective_charge(imp::Union{String, Symbol}; kw...)
    a = get_abundance_fraction(imp; kw...)
    nZ, nne, nTe = size(a.fZ_grid)
    nZ = nZ-1
    t = zeros(nZ,nne,nTe)
    for Z = 1:nZ
        t[Z,:,:] = a.fZ_grid[Z+1,:,:] .* (Z.^2 - Z)
    end
    t_grid = sum(t;dims=1)[1,:,:]
    Te = a.Te
    ne = a.ne
    t_ = Interpolations.linear_interpolation((ne,Te),t_grid, extrapolation_bc = Interpolations.Flat())
    return EffectiveCharge(t_,t_grid,a)
end

get_effective_charge(imps::Vector{<:Union{String, Symbol}}; kw...) = [get_abundance_fraction(imp; kw...) for imp in imps] 

(ec::EffectiveCharge{Z,T,AF})(fraction,ne,Te) where {Z,T,AF<:AbundanceFraction} = 1.0 .+ fraction  .* ec.t(ne,Te)
function (ecs::Vector{<:EffectiveCharge})(fractions::Vector{Float64},ne,Te)
    @assert length(fractions) == length(ecs) "provide a fraction for each species: $([ec.af.imp in ecs])"
    return 1.0 .+ sum([f .* ec.t(ne,Te) for (f,ec) in zip(fractions,ecs)])
end 
function get_abundance_fraction(imp::Union{String, Symbol}; kw...)
    scd = retrieve_ADAS_data(imp;type="scd", kw...)
    acd = retrieve_ADAS_data(imp;type="acd", kw...)

    ndens = length(acd.data.axis.ne)
    ntemp = length(acd.data.axis.Te)
    nZ = length(acd.data.rates)
    fZ = zeros(nZ+1, ndens, ntemp)
    fZ_ = zeros(nZ+1,ndens, ntemp)


        fZ[1,:,:] .= 1.   #Initial neutral fraction - everthing scaled to this
        for Z in 2:nZ+1
            fZ[Z,:,:] .= fZ[Z-1,:,:] .*  scd.data.rates[Z-1].values[:,:] ./ acd.data.rates[Z-1].values[:,:]
        end

    for i in axes(fZ_,1)
        fZ_[i,:,:] .= fZ[i,:,:]./sum(fZ,dims=1)[1,:,:]
    end
    Te = scd.data.axis.Te
    ne = scd.data.axis.ne
    Z = collect(0:nZ)
    fZ = Interpolations.linear_interpolation((float.(Z), ne,Te),fZ_, extrapolation_bc = Interpolations.Flat())
    return AbundanceFraction(fZ, fZ_, Te, ne, Z, scd, acd, imp)
end



# #= ========================= =#
# #  ActorDivertorHeatFluxTarget #
# #= ========================= =#

# import BoundaryPlasmaModels

# #==========#

# mutable struct ActorDivertorHeatFluxTarget <: PlasmaAbstractActor
#     dd::IMAS.dd
#     par::BoundaryPlasmaModels.DivertorHeatFluxTargetModelParameters
#     model::Union{Nothing,BoundaryPlasmaModels.DivertorHeatFluxModel}
# end

# # """
# # ActorDivertorHeatFluxTarget(dd::IMAS.dd, act::ParametersAllActors; kw...)

# # doc
# # """
# # function ActorDivertorHeatFluxTarget(dd::IMAS.dd, act::ParametersAllActors; kw...)
# #     par = act.ActorDivertorHeatFluxTarget(kw...)
# #     actor = ActorDivertorHeatFluxTarget(dd, par)
# #     step(actor)
# #     finalize(actor)
# #     return actor
# # end

# # function ActorDivertorHeatFluxTarget(dd::IMAS.dd, par::BoundaryPlasmaModels.FUSEparameters__ActorPlasmaFacingSurfaces; kw...)
# #     logging_actor_init(ActorPlasmaFacingSurfaces)
# #     par = par(kw...)
# #     return ActorDivertorHeatFluxTarget(dd, par, nothing)
# # end

# # function _step(actor::ActorDivertorHeatFluxTarget)
# #     actor.model = BoundaryPlasmaModels.DivertorHeatFluxModel(actor.dd,actor.par)
# #     return actor
# # end

# # function _finalize(actor::ActorDivertorHeatFluxTarget; kw...)
# #     BoundaryPlasmaModels.export2dd(actor.dd, actor.heatflux_model; kw...)
# #     actor
# # end
