using NCDatasets
using Printf
using Statistics

function createsndmean(
    sndname::AbstractString;
    exp::AbstractString, config::AbstractString,
    ndays::Integer=100, psfc::Real, extract=false
)
    
    mkpath(projectdir("exp/snd")); fsnd = projectdir("exp/snd/$(sndname)")
    statds = NCDataset(datadir("$exp/$config/OUT_STAT/RCE_TroPrecLS-$exp.nc"))

    nz = statds.dim["z"]; nt = statds.dim["time"]; z = statds["z"][:]
    dt = (statds["time"][end] - statds["time"][1]) / (nt-1)
    dystep = round(Int,1/dt); beg = ndays*dystep-1

    snddata = zeros(nz,6)
    snddata[:,1] .= z; snddata[:,2] .= -999.0
    snddata[:,3] .= dropdims(mean(statds["THETA"][:,(end-beg):end],dims=2),dims=2)
    snddata[:,4] .= dropdims(mean(statds["QT"][:,(end-beg):end],dims=2),dims=2)
    close(statds)

    printsnd(fsnd,snddata,psfc)

    if extract; return z,snddata end

end

function createsndmean(
    sndname::AbstractString,
    snddata::AbstractArray;
    psfc::Real, extract=false
)

    mkpath(projectdir("exp/snd")); fsnd = projectdir("exp/snd/$(sndname)")
    printsnd(fsnd,snddata,psfc)

    if extract; return z,snddata end

end

function createsnddiurnal(
    sndname::AbstractString;
    exp::AbstractString, config::AbstractString,
    ndays::Integer=100, psfc::Real, extract=false
)
    mkpath(projectdir("exp/snd")); fsnd = projectdir("exp/snd/$(sndname)")
    statds = NCDataset(datadir("$exp/$config/OUT_STAT/RCE_TroPrecLS-$exp.nc"))

    nz = statds.dim["z"]; nt = statds.dim["time"]; p = statds["p"][:]
    dt = (statds["time"][end] - statds["time"][1]) / (nt-1)
    dystep = round(Int,1/dt); beg = ndays*dystep-1
    dbeg = Int(floor(statds["time"][1]))
    dend = Int(floor(statds["time"][end]))

    snddata = zeros(nz,6,dystep);
    snddata[:,1,:] .= statds["z"][:]; snddata[:,2,:] .= -999.0

    tmp = mod.(statds["time"][(end-beg):end],1)
    raw = dropdims(mean(reshape(tmp,dystep,ndays),dims=2),dims=2)
    ind = sortperm(raw); t = raw[ind]

    tmp  = zeros(nz,ndays*dystep)
    tmp .= statds["THETA"][:,(end-beg):end]
    raw  = dropdims(mean(reshape(tmp,nz,dystep,ndays),dims=3),dims=3)
    snddata[:,3,:] .= raw[:,ind]

    tmp .= statds["QT"][:,(end-beg):end];
    raw  = dropdims(mean(reshape(tmp,nz,dystep,ndays),dims=3),dims=3)
    snddata[:,4,:] .= raw[:,ind]

    close(statds)

    printsnd(fsnd,snddata,psfc,t,[dbeg,dend])

    if extract; return t,p,snddata end

end

function printsnd(fsnd::AbstractString, snddata::Array{<:Real,2}, p::Real)

    nz = size(snddata,1)

    open(fsnd,"w") do io
        @printf(io,"      z[m]      p[mb]      tp[K]    q[g/kg]     u[m/s]     v[m/s]\n")
    end

    open(fsnd,"a") do io
        for t in [0.0,10000.0]
            @printf(io,"%10.2f, %10d, %10.2f\n",t,nz,p)
            for iz = 1 : nz
                @printf(
                    io,"%16.8f\t%16.8f\t%16.8f\t%16.8f\t%16.8f\t%16.8f\n",
                    snddata[iz,1],snddata[iz,2],snddata[iz,3],
                    snddata[iz,4],snddata[iz,5],snddata[iz,6]
                )
            end
        end
    end

end

function printsnd(
    fsnd::AbstractString, snddata::Array{<:Real,3},
    p::Real, t::AbstractVector, dvec::AbstractVector
)

    nz = size(snddata,1); nt = length(t)

    open(fsnd,"w") do io
        @printf(io,"      z[m]      p[mb]      tp[K]    q[g/kg]     u[m/s]     v[m/s]\n")
    end

    open(fsnd,"a") do io
        for dy in dvec[1] : dvec[2], it = 1 : nt
            @printf(io,"%10.4f, %10d, %10.2f\n",t[it]+dy,nz,p)
            for iz = 1 : nz
                @printf(
                    io,"%16.8f\t%16.8f\t%16.8f\t%16.8f\t%16.8f\t%16.8f\n",
                    snddata[iz,1,it],snddata[iz,2,it],snddata[iz,3,it],
                    snddata[iz,4,it],snddata[iz,5,it],snddata[iz,6,it]
                )
            end
        end
    end

end
