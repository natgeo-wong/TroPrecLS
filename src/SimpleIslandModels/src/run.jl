function run(
    model :: SimpleIslandModel;
    dt :: Real = 5.,
    nsteps :: Int,
    nstats :: Int,
    Ts0 :: Real,
    Ta0 :: Real,
)

    atm  = model.atm
    sfc  = model.sfc
    cld  = model.cld
    vars = generateVariables(nsteps,nstats,dt)
    initializeVars!(vars,Ts0,Ta0,model,sfc,atm)

    istat = 0
    ifin  = 0
    @showprogress "Running SimpleIslandModels.jl over $nsteps steps, dt = $dt s, model elapsed time = $(nsteps*dt/86400) days ..." for it = 1 : nsteps

        ifin += 1
        n0 = mod(it,nstats)
        stepforward!(vars,model,sfc,atm,cld,n0)
        if iszero(n0)
            istat += 1
            savestats!(vars,nstats,istat)
            ifin = 0
        end

    end

    if !iszero(mod(nsteps,nstats))
        istat += 1
        savestats!(vars,ifin,istat)
    end

    return vars

end

function stepforward!(
    vars :: Variables, m :: SimpleIslandModel,
    sfc :: Surface, atm :: Atmosphere, cld :: CloudAlbedo,
    n0 :: Real
)

    t  = vars.temp[1]
    Ts = vars.temp[3]
    Ta = vars.temp[4]
    αₐ = vars.temp[7]
    δt = vars.δt
    t += δt

    if m.do_cloud
        αₐ = cld.calculateαₐ(αₐ,Ts,Ta,δt,cld.params,atm,sfc,)
    end

    S₀      = calculateS₀(t,m) * (1-αₐ)
    Fₐ, δTa = calculateδTa(S₀,Ts,Ta,δt,m,atm,sfc)
    Fₛ, δTs = calculateδTs(S₀,Ts,Ta,δt,m,atm,sfc)

    vars.temp[1] = t
    vars.temp[2] = S₀
    vars.temp[3] = Ts + δTs
    vars.temp[4] = Ta + δTa
    vars.temp[5] = Fₛ
    vars.temp[6] = Fₐ
    vars.temp[7] = αₐ

    if !iszero(n0)
        @views @. vars.stat += vars.temp
    else
        @views @. vars.stat += vars.temp * 0.5
    end

    return nothing

end

function savestats!(vars::Variables,nstats::Int,istat::Int)

    vars.t[istat]  = vars.stat[1] / (nstats)
    vars.S₀[istat] = vars.stat[2] / (nstats)
    vars.Tₛ[istat] = vars.stat[3] / (nstats)
    vars.Tₐ[istat] = vars.stat[4] / (nstats)
    vars.Fₛ[istat] = vars.stat[5] / (nstats)
    vars.Fₐ[istat] = vars.stat[6] / (nstats)
    vars.αₐ[istat] = vars.stat[7] / (nstats)

    @views @. vars.stat = vars.temp * 0.5

end