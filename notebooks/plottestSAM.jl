### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 8102297a-50b7-11eb-0430-f79371a66174
begin
	using DrWatson
	using Pkg
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 9340fa4e-50b4-11eb-253e-ab01deb80456
begin
	@quickactivate "TroPrecLS"
	Pkg.instantiate()
	using NCDatasets
	using NumericalIntegration
	using Printf
	using Statistics
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ addc35d6-50b3-11eb-02dc-452ced2a45ef
md"
# X. Plotting Functions for SAM Output

Here I am creating and testing functions that will be used for the extraction and plotting of SAM statistical output data.
"

# ╔═╡ 19174960-50cf-11eb-12a3-cf977e483262
md"
### 1. Retrieving and Manipulating Statistical Output

We first create the function `retrievedims` to extract the dimensional data: (1) height (in both vertical `z` and pressure `p` coordinates) and (2) time `t`.  Next, we define the function `retrievevar` to extract variables of interest (one variable at a time):
"

# ╔═╡ e4b19fb2-74bf-11eb-162e-2f42a9e40fea
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

# ╔═╡ c1489ae0-5114-11eb-3a56-5b75d263ae63
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

# ╔═╡ 891e5992-50cf-11eb-208e-71e78380caf7
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

# ╔═╡ edc78916-5115-11eb-00a0-498face6f531
md"
Relative humidity in SAM is plotted against the liquid saturation point.  However, we are able to calculate relative humidity using the ice-saturation point for air temperatures below 0ºC, using specific humidity and atmospheric temperature data.
"

# ╔═╡ 6b0b697c-5117-11eb-1c0f-13c2840a65d2
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

# ╔═╡ 48f8e67e-5116-11eb-0037-3b6630655917
function calcrh(QV,TAIR,P)

	RH = zeros(size(QV)); np = size(RH,1); nt = size(RH,2)

	for it = 1 : nt, ip = 1 : np
		RH[ip,it] = QV[ip,it] / tair2qsat(TAIR[ip,it],P[ip]*100)
	end

	return RH

end

# ╔═╡ a8331f82-5117-11eb-158b-b3073e67affd
md"
From the relative humidity that we calculate, we then can calculate column saturation fraction `CSF`, or the overall relative humidity of the column
"

# ╔═╡ dacb2ffe-5117-11eb-33bf-ff68ff8c6504
function calccsf(RH,QV,P)
	
	pvec = vcat(0,reverse(P)) * 100; nt = size(RH,2); np = length(pvec)
	QVtmp = zeros(length(pvec))
	QVsat = zeros(length(pvec))
    csf = zeros(nt)
	swp = zeros(nt)

    for it = 1 : nt
		for ip = 1 : (np-1)
			QVtmp[ip+1] = QV[np-ip,it]
			tmp = QV[np-ip,it] ./ RH[np-ip,it]
			if !isnan(tmp)
				  QVsat[ip+1] = tmp
			else; QVsat[ip+1] = 0
			end
		end
        csf[it] = integrate(pvec,QVtmp) / integrate(pvec,QVsat)
		swp[it] = integrate(pvec,QVsat) / 9.81 / 1000
    end

	return csf,swp
	
end

# ╔═╡ 5dc0fdf6-9326-11eb-1009-db847dddab8e
function calccsf2(RH,P)

    pvec = vcat(0,reverse(P)); nt = size(RH,2)
    pint = integrate(pvec,ones(length(pvec)))
    RHtmp = zeros(length(pvec))
    csf = zeros(nt)

    for it = 1 : nt
    	RHtmp[2:end] .= reverse(RH[:,it])
        csf[it] = integrate(pvec,RHtmp) / pint
    end

    return csf

end

# ╔═╡ f188474a-9319-11eb-1bc9-7754318339b1
function calcpw(QV,P)
	
	pvec = vcat(0,reverse(P)); nt = size(QV,2)
	QVtmp = zeros(length(pvec))
    pw = zeros(nt)
	
    for it = 1 : nt
		QVtmp[2:end] .= QV[:,it]
        pw[it] = integrate(pvec*100,reverse(QVtmp)) / 9.81/1000
    end
	
	return pw
	
end

# ╔═╡ 2fb20134-5119-11eb-1d07-db7f8c326dda
md"
Let us now define functions that will plot the time series 2D and 3D variables on a given axes object `axs`, and axes number `ii`, between days `dbeg` and `dend`.  The function for 3D plot will return information that can be used to plot a colorbar.  The 3D timeseries can also accept optional colormap `cmapname` and levels `lvl` input.
"

# ╔═╡ 4324677a-5119-11eb-24f1-8d80c5a15ceb
function plot2Dtimeseries(axs,ii,t,var;dbeg,dend)
	axs[ii].plot(t,var,lw=1)
	axs[ii].format(xlim=(dbeg,dend))
end

# ╔═╡ e006175a-5119-11eb-02ba-b50a41b540cb
function plot3Dtimeseries(axs,ii,t,p,var;lvl=[],cmapname="Fire",dbeg,dend)
	
	if isempty(lvl)
		  c = axs[ii].contourf(t,p,var,cmap=cmapname,extend="both")
	else; c = axs[ii].contourf(t,p,var,cmap=cmapname,levels=lvl,extend="both")
	end
	
	axs[ii].format(xlim=(dbeg,dend),ylim=(maximum(p),10))
	axs[ii].colorbar(c,loc="r")
	
end

# ╔═╡ f7cf698a-511f-11eb-0276-7d4d7ae912e2
md"
Following this, we define functions that allow us to calculate and plot the diurnal cycle of the variables considered based on the last `ndays` of the simulation.
"

# ╔═╡ 23322f5a-5126-11eb-3d55-f3135e8a8bdc
function t2d(t::Vector{<:Real}, days::Integer)

    tstep = round(Integer,(length(t)-1)/(t[end]-t[1]))
    t = mod.(t[(end-tstep+1):end],1); tmin = argmin(t)
    tshift = tstep-tmin+1; t = circshift(t,tshift)
    t = vcat(t[end]-1,t,t[1]+1)
	beg = days*tstep - 1

    return t*tstep,tstep,tshift,beg

end

# ╔═╡ bea41e8c-5127-11eb-0715-e55c3f2a05ab
function diurnal2D(data::AbstractVector,tstep::Integer,tshift::Integer)

    data = dropdims(mean(reshape(data,tstep,:),dims=2),dims=2);
    data = circshift(data,tshift)

    return vcat(data[end],data,data[1])

end

# ╔═╡ 55946b50-5126-11eb-14e4-570a92dee964
function diurnal3D(data::AbstractArray{<:Real,2},tstep::Integer,tshift::Integer)

    nz = size(data,1)
    data = dropdims(mean(reshape(data,nz,tstep,:),dims=3),dims=3);
    data = circshift(data,(0,tshift))

    return cat(dims=2,data[:,end],data,data[:,1])

end

# ╔═╡ e0d83a56-5127-11eb-26a2-9dc0aa1c9822
function plot2Ddiurnal(axs,ii,t,var;subtractm=true)
	
	mvar = mean(var)
	if subtractm
		  axs[ii].plot(t,var .- mvar,lw=1)
	else; axs[ii].plot(t,var,lw=1)
	end
	
	axs[ii].format(xlim=(0,24))
	if !subtractm; axs[ii].format(ylim=(0,2.5)) end
	
end

# ╔═╡ a8d78482-5126-11eb-2dd1-8751cbd884a7
function plot3Ddiurnal(axs,ii,t,p,var;lvl=[],cmapname="Fire")
	
	mvar = dropdims(mean(var,dims=2),dims=2)
	axs[ii].plot(mvar,p,lw=1)
	
	if isempty(lvl)
		  c = axs[ii+1].contourf(t,p,var,cmap=cmapname,extend="both")
	else; c = axs[ii+1].contourf(t,p,var,cmap=cmapname,levels=lvl,extend="both")
	end
	
	axs[ii+1].format(xlim=(0,24),ylim=(maximum(p),20))
	axs[ii+1].colorbar(c,loc="r")
	
end

# ╔═╡ 6bc27074-5115-11eb-29ed-bf79782322d0
md"
### 2. 2D Statistical Output

The following 2D statistical output are of interest:
* Sea-surface temperature | `SST`
* Surface pressure | `Ps`
* Precipitation | `Prec`
* Precipitable Water | `PW`
* Sensible/Latent Heat Flux | `SHF`/`LHF`
* Surface Net Shortwave/Longwave Flux | `SWNS`/`LWNS`
"

# ╔═╡ 504cbace-5128-11eb-1d27-879bffb48098
md"
### 3. Testing functionality
"

# ╔═╡ 5ffd515e-5128-11eb-3b11-815af069d22f
begin
	# exp = "MakePC"; config = "Forcing0000-Slab10d00"; istst = false;
	expi = "Control"; config = "slab00d05"; istst = false;
	expi = "Slab00d05"; config = "forcingn002"; istst = false;
	# exp = "DiAmp064km"; config = "Slab31d6"; istst = true;
	isen = true; mbr=1
end

# ╔═╡ f440d5ca-5128-11eb-2574-d58f2f7d8fc3
begin
	
	pplt.close(); fts,axsts = pplt.subplots(nrows=3,aspect=3,axwidth=4)
	
	lvls = [
		-100,-70.7,-50,-31.6,-20,-14.1,-10,-7.07,
		-5,5,
		7.07,10,14.1,20,31.6,50,70.7,100
	]/10
	
	z,p,t = retrievedims(expi,config,istest=istst,isensemble=isen,member=mbr)
	nt = length(t)
	# var2D = retrievevar("AREAPREC",exp,config,istest=istst,isensemble=isen,member=mbr)
	var2A = retrievevar("SST",expi,config,istest=istst,isensemble=isen,member=mbr)
	var3T = retrievevar("CLD",expi,config,istest=istst,isensemble=isen,member=mbr)
	var3Q = retrievevar("QBIAS",expi,config,istest=istst,isensemble=isen,member=mbr)
	varob = retrievevar("QVOBS",expi,config,istest=istst,isensemble=isen,member=mbr)
	varTA = retrievevar("TABS",expi,config,istest=istst,isensemble=isen,member=mbr)
	varQV = retrievevar("QV",expi,config,istest=istst,isensemble=isen,member=mbr)/1000
	varPW = retrievevar("PW",expi,config,istest=istst,isensemble=isen,member=mbr)
	varrh = retrievevar("RELH",expi,config,istest=istst,isensemble=isen,member=mbr)
	
	rh  = calcrh(varQV,varTA,p)
	csf,swp = calccsf(rh,varQV,p)
	csf2 = calccsf2(rh,p)
	pw  = calcpw(varQV,p)
	# csfa = calccsf(varrh,p)
	
# 	pbin = 0:10:250
# 	test = fit(Histogram,(mod.(t,1)*24,var2A),(0:24,pbin)).weights
	
# 	nt = length(t); nhr = 24; np = length(pbin) - 1
# 	test = test/nt*np*nhr
# 	test = vcat(test,test)
# 	test[iszero.(test)] .= NaN
	
	plot3Dtimeseries(
		axsts,1,t.-80,p,
		varQV*100,
		dbeg=0,dend=300,lvl=lvls,cmapname="RdBu"
	)
	plot3Dtimeseries(
		axsts,2,t.-80,p,
		var3Q./varob*100,
		dbeg=0,dend=300,lvl=lvls*10,cmapname="drywet"
	)
	# plot3Dtimeseries(
	# 	axsts,1,t.-80,p,varob,
	# 	dbeg=250,dend=300,lvl=0:10:100,cmapname="drywet"
	# )
	# plot2Dtimeseries(
	# 	axsts,1,t.-80,var2A,
	# 	dbeg=0,dend=300
	# )
	# plot2Dtimeseries(
	# 	axsts,3,t.-80,swp*1000,
	# 	dbeg=250,dend=300
	# )
	plot2Dtimeseries(
		axsts,3,t.-80,var2A,
		dbeg=00,dend=50
	)
	# plot2Dtimeseries(
	# 	axsts,3,t.-80,varPW,
	# 	dbeg=250,dend=300
	# )
	# plot2Dtimeseries(
	# 	axsts,3,t.-80,csfa,
	# 	dbeg=250,dend=300
	# )
	# c = axsts[1].pcolormesh(
	# 	-24:24,0:10:250,test',cmap="Blues",cmap_kw=Dict("left"=>0.05),
	# 	levels=10. .^(-1:0.2:1),extend="both"
	# )
	# axsts[1].colorbar(c,loc="r")
	# axsts[1].format(xlim=(-12,12),xlocator=-12:6:12)
	# axsts[1].scatter(mod.(t,1),var2A)
	# plot2Dtimeseries(
	# 	axsts,3,1:50,pr2,
	# 	dbeg=0,dend=200
	# )
	# plot2Dtimeseries(
	# 	axsts,3,t2.-80,var2D2,
	# 	dbeg=0,dend=100
	# )
	# plot2Dtimeseries(
	# 	axsts,3,1:24,dropdims(mean(reshape(var2D2[1:(24*24)],24,:),dims=1),dims=1),
	# 	dbeg=0,dend=100
	# )
	
	axsts[1].format(ylim=(1000,20),yscale="log")
	axsts[2].format(ylim=(1000,20),yscale="log")
	# axsts[3].format(ylim=(0,1))
	# axsts[1].format(xlim=(-12,12),xlocator=-24:3:24)
	
	fts.savefig("plots.png",transparent=false,dpi=200)
	load("plots.png")
	
end

# ╔═╡ eb323e5e-92a2-11eb-3596-4ff64f0cea13
mean(var2A)

# ╔═╡ 8676ab94-81ef-11eb-2f9f-7716d1507bc2
md"
### 4. New things ...
"

# ╔═╡ 3a4f1fe4-92a2-11eb-07c6-fdd645db4f3c
mod.(t,1)*24

# ╔═╡ 3530faa4-931b-11eb-3120-5d25e45d493c
p

# ╔═╡ Cell order:
# ╟─addc35d6-50b3-11eb-02dc-452ced2a45ef
# ╟─8102297a-50b7-11eb-0430-f79371a66174
# ╟─9340fa4e-50b4-11eb-253e-ab01deb80456
# ╟─19174960-50cf-11eb-12a3-cf977e483262
# ╠═e4b19fb2-74bf-11eb-162e-2f42a9e40fea
# ╠═c1489ae0-5114-11eb-3a56-5b75d263ae63
# ╠═891e5992-50cf-11eb-208e-71e78380caf7
# ╟─edc78916-5115-11eb-00a0-498face6f531
# ╠═48f8e67e-5116-11eb-0037-3b6630655917
# ╠═6b0b697c-5117-11eb-1c0f-13c2840a65d2
# ╟─a8331f82-5117-11eb-158b-b3073e67affd
# ╠═dacb2ffe-5117-11eb-33bf-ff68ff8c6504
# ╠═5dc0fdf6-9326-11eb-1009-db847dddab8e
# ╠═f188474a-9319-11eb-1bc9-7754318339b1
# ╟─2fb20134-5119-11eb-1d07-db7f8c326dda
# ╠═4324677a-5119-11eb-24f1-8d80c5a15ceb
# ╠═e006175a-5119-11eb-02ba-b50a41b540cb
# ╟─f7cf698a-511f-11eb-0276-7d4d7ae912e2
# ╠═23322f5a-5126-11eb-3d55-f3135e8a8bdc
# ╠═bea41e8c-5127-11eb-0715-e55c3f2a05ab
# ╠═55946b50-5126-11eb-14e4-570a92dee964
# ╠═e0d83a56-5127-11eb-26a2-9dc0aa1c9822
# ╠═a8d78482-5126-11eb-2dd1-8751cbd884a7
# ╟─6bc27074-5115-11eb-29ed-bf79782322d0
# ╟─504cbace-5128-11eb-1d27-879bffb48098
# ╠═5ffd515e-5128-11eb-3b11-815af069d22f
# ╠═f440d5ca-5128-11eb-2574-d58f2f7d8fc3
# ╠═eb323e5e-92a2-11eb-3596-4ff64f0cea13
# ╟─8676ab94-81ef-11eb-2f9f-7716d1507bc2
# ╠═3a4f1fe4-92a2-11eb-07c6-fdd645db4f3c
# ╠═3530faa4-931b-11eb-3120-5d25e45d493c
