using Glob
using NCDatasets
using NumericalIntegration
using Printf
using Statistics

function outstatname(
    experiment::AbstractString, config::AbstractString,
    istest::Bool=false,
    isensemble::Bool=false, member::Integer=0
)

    if isensemble
    	  expname = "$(experiment)-member$(@sprintf("%02d",member))"
    else; expname = experiment
    end

    if istest
    	fnc = datadir(joinpath(
    		experiment,config,"OUT_STAT",
    		"RCE_TroPrecLS-$(expname)-test.nc"
    	))
    else
    	fnc = datadir(joinpath(
    		experiment,config,"OUT_STAT",
    		"RCE_TroPrecLS-$(expname).nc"
    	))
    end

    return fnc

end

function out2Dname(
    experiment::AbstractString, config::AbstractString,
    isensemble::Bool=false, member::Integer=0
)

    if isensemble
    	  expname = "$(experiment)-member$(@sprintf("%02d",member))"
    else; expname = experiment
    end

	fnclist = glob(
		"RCE_TroPrecLS-$(expname)*.nc",
		joinpath(datadir(experiment,config,"OUT_2D"))
	); fnc = datadir(fnclist[1])

    return fnc

end

function retrievedims(
    experiment::AbstractString, config::AbstractString;
    istest::Bool=false,
    isensemble::Bool=false, member::Integer=0
)

    rce = NCDataset(outstatname(experiment,config,istest,isensemble,member))
    z = rce["z"][:]
    p = rce["p"][:]
    t = rce["time"][:]
    close(rce)

    return z,p,t

end

function retrievedims(fnc::AbstractString)

    rce = NCDataset(fnc)
    z = rce["z"][:]
    p = rce["p"][:]
    t = rce["time"][:]
    close(rce)

    return z,p,t

end

function retrievedims2D(
    experiment::AbstractString, config::AbstractString;
    isensemble::Bool=false, member::Integer=0
)

    rce = NCDataset(out2Dname(experiment,config,isensemble,member))
	x = rce["x"][:]
    t = rce["time"][:]

	if haskey(rce,"y")
		  isy = true; y = rce["y"][:]
	else; isy = false
	end

    close(rce)

	if isy
    	  return x,y,t
	else; return x,t
	end

end

function retrievedims2D(fnc::AbstractString)

    rce = NCDataset(fnc)
    x = rce["x"][:]
    t = rce["time"][:]

	if haskey(rce,"y")
		  isy = true; y = rce["y"][:]
	else; isy = false
	end

    close(rce)

	if isy
    	  return x,y,t
	else; return x,t
	end

end

function retrievevar(
    variable::AbstractString,
    experiment::AbstractString, config::AbstractString;
    istest::Bool=false,
    isensemble::Bool=false, member::Integer=0
)

    rce = NCDataset(outstatname(experiment,config,istest,isensemble,member))
    var = rce[variable][:]
    close(rce)

    return var

end

function retrievevar(variable::AbstractString, fnc::AbstractString)

    rce = NCDataset(fnc)
    var = rce[variable][:]
    close(rce)

    return var

end

function retrievevar2D(
    variable::AbstractString,
    experiment::AbstractString, config::AbstractString;
    isensemble::Bool=false, member::Integer=0
)

    rce = NCDataset(out2Dname(experiment,config,isensemble,member))
    var = rce[variable][:]
    close(rce)

    return var

end

function retrievevar2D(variable::AbstractString, fnc::AbstractString)

    rce = NCDataset(fnc)
    var = rce[variable][:]
    close(rce)

    return var

end

function calcrh(QV,TAIR,P)

    RH = zeros(size(QV)); np = size(RH,1); nt = size(RH,2)

    for it = 1 : nt, ip = 1 : np
    	RH[ip,it] = QV[ip,it] / tair2qsat(TAIR[ip,it],P[ip]*100)
    end

    return RH

end

function tair2qsat(T,P)

    tb = T - 273.15
    if tb <= 0
    	esat = exp(43.494 - 6545.8/(tb+278)) / (tb+868)^2
    else
    	esat = exp(34.494 - 4924.99/(tb+237.1)) / (tb+105)^1.57
    end


    r = 0.622 * esat / max(esat,P-esat)
    return r / (1+r)

end

function calccsf(RH,P)

    pvec = vcat(0,reverse(P)); nt = size(RH,2)
    pint = integrate(pvec,ones(length(pvec)))
    RHtmp = zeros(length(pvec))
    csf = zeros(nt)

    for it = 1 : nt
    	RHtmp[2:end] .= RH[:,it]
        csf[it] = integrate(pvec,RHtmp) / pint
    end

    return csf

end

function t2d(t::Vector{<:Real}, days::Integer)

    tstep = round(Integer,(length(t)-1)/(t[end]-t[1]))
    t = mod.(t[(end-tstep+1):end],1); tmin = argmin(t)
    tshift = tstep-tmin+1; t = circshift(t,tshift)
    t = vcat(t[end]-1,t,t[1]+1)
    beg = days*tstep - 1

    return t*tstep,tstep,tshift,beg

end

function diurnal2D(data::AbstractVector,tstep::Integer,tshift::Integer)

    data = dropdims(mean(reshape(data,tstep,:),dims=2),dims=2);
    data = circshift(data,tshift)

    return vcat(data[end],data,data[1])

end

function diurnal3D(data::AbstractArray{<:Real,2},tstep::Integer,tshift::Integer)

    nz = size(data,1)
    data = dropdims(mean(reshape(data,nz,tstep,:),dims=3),dims=3);
    data = circshift(data,(0,tshift))

    return cat(dims=2,data[:,end],data,data[:,1])

end
