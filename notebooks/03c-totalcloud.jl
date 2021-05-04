### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 6f00b8fc-530c-11eb-2242-99d8544f6e14
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 8f30c56c-530c-11eb-2782-33f3c4ed9e89
begin
	@quickactivate "TroPrecLS"
	using Dates
	using DelimitedFiles
	using GeoRegions
	using Interpolations
	using NCDatasets
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("common.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 90fffbc8-524d-11eb-232a-1bada28d5505
md"
# 3c. Total Cloud Cover Data

The first variable that we want to explore among the ERA5 reanalysis variables is Skin Temperature.  More details about the retrieval of reanalysis data can be found in notebook `03-reanalysis.jl`.
"

# ╔═╡ c1fafcba-530f-11eb-1cc2-67d10a3fa606
md"
### A. Modelling the Diurnal Cycle of Skin Temperature

We wish to find the following characteristics of the diurnal cycle of skin temperature:
* The mean $\mu$ skin temperature
* The amplitude $A$ of the diurnal cycle (max-min)/2 of skin temperature
* The hour $\theta$ at which skin temperature is a maximum
"

# ╔═╡ 3565af3c-5311-11eb-34c4-2d228b05b17c
md"
ERA5 reanalysis is defined according to UTC.  Therefore, for every grid point, we must also correct to the local time before fitting the data to the model.
"

# ╔═╡ fd344560-94d3-11eb-2b79-05c905c5953f
longitude2timeshift(longitude::Real) = longitude / 180 * 12

# ╔═╡ a6a688ca-53ab-11eb-2776-b5380ffb26c1
function eradiurnal2model(data,longitude)
	
	nlon,nlat,nt = size(data)
	θ = zeros(nlon,nlat)
	A = zeros(nlon,nlat)
	μ = zeros(nlon,nlat)
	
	idat = zeros(nt+1)
	it   = 0:0.1:24; it = it[1:(end-1)]; nit = length(it)
	ts   = longitude2timeshift.(longitude)
	var  = zeros(nit)
	
	for ilat = 1 : nlat, ilon = 1 : nlon
		
		tl = (0:24) .+ ts[ilon]
		
        idat[1:24] = data[ilon,ilat,:]
		idat[end]  = data[ilon,ilat,1]
		
		itp = interpolate(idat,BSpline(Cubic(Periodic(OnGrid()))))
        stp = scale(itp,tl)
        etp = extrapolate(stp,Periodic())
		var[:] = etp[it]
		
		μ[ilon,ilat] = mean(@view data[ilon,ilat,:])
		A[ilon,ilat] = (maximum(idat) .- minimum(idat))/2
		θ[ilon,ilat] = argmax(var) / nit * 24

    end
	
	return μ,A,θ
	
end

# ╔═╡ aa05317e-530b-11eb-2ec1-93aff65659dd
md"
### B. Retrieving ERA5 Skin Temperature Data
"

# ╔═╡ 103f85e8-530c-11eb-047d-a537aa60075d
function retrieveera()

    ds  = NCDataset(datadir("reanalysis/hourly/era5-TRPx0.25-tcc-sfc-hour.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	var = ds["tcc"][:] * 1
	close(ds)
	
	μ,A,θ = eradiurnal2model(var,lon)

    return lon,lat,μ,A,θ

end

# ╔═╡ 49d13e5c-53af-11eb-29ca-c994a7acd377
sroot = "/n/kuangdss01/lab/"; regID = "TRP"

# ╔═╡ e8141e20-53af-11eb-1a23-81d34293c5eb
begin
	lon,lat,μ,A,θ = retrieveera()
md"Modelling diurnal cycle of surface temperature"
end

# ╔═╡ d82366b0-53b1-11eb-26c1-ff1bb6ccb027
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ bb59b8d6-53b1-11eb-3631-87ef61219c4c
begin
	pplt.close(); f,axs = pplt.subplots(nrows=3,axwidth=6,aspect=6)
	
	c = axs[1].contourf(lon,lat,μ',levels=0.1:0.1:0.9,cmap="Blues",extend="both")
	axs[1].plot(x,y,c="k",lw=0.5)
	axs[1].format(urtitle=L"$\mu$")
	axs[1].colorbar(c,loc="r")
	
	c = axs[2].contourf(lon,lat,A',levels=10. .^(-2:0.2:0),extend="both")
	axs[2].plot(x,y,c="k",lw=0.5)
	axs[2].format(urtitle="A")
	axs[2].colorbar(c,loc="r",ticks=[0.01,0.1,1])
	
	c = axs[3].pcolormesh(lon,lat,θ',cmap="romaO",levels=0:0.5:24,extend="both")
	axs[3].plot(x,y,c="k",lw=0.5)
	axs[3].format(urtitle=L"$\theta$ / Hour of Day")
	axs[3].colorbar(c,loc="r")
	
	for ax in axs
		ax.format(
			xlim=(0,360),ylim=(-30,30),xlocator=0:60:360,
			xlabel=L"Longitude / $\degree$",
			ylabel=L"Latitude / $\degree$",
		)
	end
	
	f.savefig(plotsdir("tccmodel_TRP.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("tccmodel_TRP.png"))
end

# ╔═╡ 5c0e5bae-554e-11eb-3f83-a364ae0a2485
md"
Test
"

# ╔═╡ 68cfc46c-5755-11eb-1702-373942539652
md"
### C. Regional Analysis

We can get quick snapshots of the results for different GeoRegions specified in this project.
"

# ╔═╡ 52b39ff8-9426-11eb-2a86-43f7da15f62e
begin
	N,S,E,W = [10,-10,115,95]
	md"Defining region coordinates ..."
end

# ╔═╡ ea7f0956-575b-11eb-3e3f-a1ba3e08b771
begin
	rlon,rlat,rinfo = regiongridvec([N,S,E,W],lon,lat)
	if maximum(rlon) > 360; rlon .= rlon .- 360 end
	rμ = regionextractgrid(μ,rinfo)
	rA = regionextractgrid(A,rinfo)
	rθ = regionextractgrid(θ,rinfo)
	md"Extracting information for region ..."
end

# ╔═╡ 5714c13c-575c-11eb-06d4-838b4e8dbcd7
begin
	asp = (maximum(rlon)-minimum(rlon))/(maximum(rlat)-minimum(rlat))
	pplt.close()
	if asp > 1.5
		freg,areg = pplt.subplots(nrows=3,axwidth=3,aspect=asp)
	else
		freg,areg = pplt.subplots(ncols=3,axwidth=2,aspect=asp)
	end
	
	creg = areg[1].contourf(rlon,rlat,rμ',levels=(1:9)/10,cmap="Blues",extend="both")
	areg[1].plot(x,y,c="k",lw=0.5)
	areg[1].format(urtitle=L"$\mu$ / K")
	areg[1].colorbar(creg,loc="r")
	
	creg = areg[2].contourf(rlon,rlat,rA',levels=10. .^(-2:0.2:0),extend="both")
	areg[2].plot(x,y,c="k",lw=0.5)
	areg[2].format(urtitle="A / K")
	areg[2].colorbar(creg,loc="r",ticks=[0.01,0.1,1,10,100])
	
	creg = areg[3].pcolormesh(rlon,rlat,rθ',cmap="romaO",levels=0:0.5:24)
	areg[3].plot(x,y,c="k",lw=0.5)
	areg[3].format(urtitle=L"$\theta$ / Hour of Day")
	areg[3].colorbar(creg,loc="r",ticks=0:3:24)
	
	for ax in areg
		ax.format(
			xlim=(minimum(rlon),maximum(rlon)),
			ylim=(minimum(rlat),maximum(rlat))
		)
	end
	
	freg.savefig(plotsdir("tccmodel_reg.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("tccmodel_reg.png"))
end

# ╔═╡ c4792bf2-5552-11eb-3b52-997f59fd42f3
md"
### D. Binning of Precipitation Data

We now wish to bin the modelled precipitation data in order to determine the relationship between precipitation rate and the diurnal cycle over land and sea.  Is it different?
"

# ╔═╡ 1fadf4ca-5755-11eb-1ece-a99313019785
begin
	lds = NCDataset(datadir("reanalysis/era5-TRPx0.25-lsm-sfc.nc"))
	lsm = lds["lsm"][:]*1
	close(lds)
	
md"Loading Land-Sea Mask for ERA5 data ..."
end

# ╔═╡ f752b054-57c1-11eb-117c-ed52464aa25f
md"
#### i. Mean Precipitation Rate $\mu$ (Hourly)
"

# ╔═╡ 4b289fa8-57b9-11eb-0923-116c3d9444bb
begin
	lbins = collect(0:0.005:1); lpbin = (lbins[2:end].+lbins[1:(end-1)])/2
	lbin_SEA,lavg_SEA = bindatasfclnd([20,-15,165,90],lbins,μ,lon,lat,lsm)
	lbin_TRA,lavg_TRA = bindatasfclnd([10,-10,40,-10],lbins,μ,lon,lat,lsm)
	lbin_CRB,lavg_CRB = bindatasfclnd([25,15,-60,-90],lbins,μ,lon,lat,lsm)
	lbin_AMZ,lavg_AMZ = bindatasfclnd([10,-10,-45,-75],lbins,μ,lon,lat,lsm)
	
	sbins = collect(0:0.005:1); spbin = (sbins[2:end].+sbins[1:(end-1)])/2
	sbin_SEA,savg_SEA = bindatasfcsea([20,-15,165,90],sbins,μ,lon,lat,lsm)
	sbin_TRA,savg_TRA = bindatasfcsea([10,-10,40,-10],sbins,μ,lon,lat,lsm)
	sbin_CRB,savg_CRB = bindatasfcsea([25,15,-60,-90],sbins,μ,lon,lat,lsm)
	sbin_DTP,savg_DTP = bindatasfcsea([10,-10,360,0],sbins,μ,lon,lat,lsm)
	
	md"Binning mean skin temperatures for different tropical regions ..."
end

# ╔═╡ e7ff7ec8-57b9-11eb-0115-abbe4aa9a1a9
begin
	pplt.close(); fbin,abin = pplt.subplots(ncols=2,aspect=2);
	
	abin[1].plot(lpbin,lbin_SEA,c="b",lw=0.5)
	abin[1].plot(lpbin,lbin_TRA,c="r",lw=0.5)
	abin[1].plot(lpbin,lbin_CRB,c="blue3",lw=0.5)
	abin[1].plot(lpbin,lbin_AMZ,c="g",lw=0.5)
	
	abin[1].plot([1,1]*lavg_SEA,[0.1,50],c="b")
	abin[1].plot([1,1]*lavg_TRA,[0.1,50],c="r")
	abin[1].plot([1,1]*lavg_CRB,[0.1,50],c="blue3")
	abin[1].plot([1,1]*lavg_AMZ,[0.1,50],c="g")
	
	abin[2].plot(spbin,sbin_SEA,c="b",lw=0.5);
	abin[2].plot(spbin,sbin_TRA,c="r",lw=0.5);
	abin[2].plot(spbin,sbin_CRB,c="blue3",lw=0.5);
	abin[2].plot(spbin,sbin_DTP,c="k",lw=0.5);
	
	abin[2].plot([1,1]*savg_SEA,[0.1,50],c="b")
	abin[2].plot([1,1]*savg_TRA,[0.1,50],c="r")
	abin[2].plot([1,1]*savg_CRB,[0.1,50],c="blue3")
	abin[2].plot([1,1]*savg_DTP,[0.1,50],c="k")
	
	abin[1].format(
		xlim=(minimum(lbins),maximum(lbins)),
		ylim=(0,25),#yscale="log",
		ylabel="Density",
		ltitle="(a) Land"
	)
	
	abin[2].format(
		xlim=(minimum(sbins),maximum(sbins)),
		# ylim=(0.1,30),yscale="log",
		xlabel="Total Cloud Cover Fraction",
		ltitle="(b) Ocean"
	)
	
	fbin.savefig(plotsdir("tccdiurnalmean.png"),transparent=false,dpi=200)
	load(plotsdir("tccdiurnalmean.png"))
end

# ╔═╡ 0fbb0b46-57c2-11eb-365a-a73a2ebda8e4
md"
#### ii. Diurnal Amplitude $A$ vs $\mu$
"

# ╔═╡ 252508a8-57c2-11eb-08b5-8fa673b1ac8a
begin
	lvec = collect(0:0.001:0.2); lAbin = (lvec[2:end].+lvec[1:(end-1)])/2
	lAbin_SEA,lAavg_SEA = bindatasfclnd([20,-15,165,90],lvec,A,lon,lat,lsm)
	lAbin_TRA,lAavg_TRA = bindatasfclnd([10,-10,40,-10],lvec,A,lon,lat,lsm)
	lAbin_CRB,lAavg_CRB = bindatasfclnd([25,15,-60,-90],lvec,A,lon,lat,lsm)
	lAbin_AMZ,lAavg_AMZ = bindatasfclnd([10,-10,-45,-75],lvec,A,lon,lat,lsm)
	
	svec = collect(0:0.001:0.2); sAbin = (svec[2:end].+svec[1:(end-1)])/2
	sAbin_SEA,sAavg_SEA = bindatasfcsea([20,-15,165,90],svec,A,lon,lat,lsm)
	sAbin_TRA,sAavg_TRA = bindatasfcsea([10,-10,40,-10],svec,A,lon,lat,lsm)
	sAbin_CRB,sAavg_CRB = bindatasfcsea([25,15,-60,-90],svec,A,lon,lat,lsm)
	sAbin_DTP,sAavg_DTP = bindatasfcsea([10,-10,360,0],svec,A,lon,lat,lsm)
	
	md"Binning amplitude of the diurnal cycle for skin temperatures in different tropical regions ..."
end

# ╔═╡ 5f58ae9c-57c2-11eb-1f04-2ddbaf2b4f1b
begin
	pplt.close(); fA,aA = pplt.subplots(ncols=2,aspect=2);
	
	aA[1].plot(lAbin,lAbin_SEA,c="b",lw=0.5)
	aA[1].plot(lAbin,lAbin_TRA,c="r",lw=0.5)
	aA[1].plot(lAbin,lAbin_CRB,c="blue3",lw=0.5)
	aA[1].plot(lAbin,lAbin_AMZ,c="g",lw=0.5)
	
	aA[1].plot([1,1]*lAavg_SEA,[0.1,50],c="b")
	aA[1].plot([1,1]*lAavg_TRA,[0.1,50],c="r")
	aA[1].plot([1,1]*lAavg_CRB,[0.1,50],c="blue3")
	aA[1].plot([1,1]*lAavg_AMZ,[0.1,50],c="g")
	
	aA[2].plot(sAbin,sAbin_SEA,c="b",lw=0.5);
	aA[2].plot(sAbin,sAbin_TRA,c="r",lw=0.5);
	aA[2].plot(sAbin,sAbin_CRB,c="blue3",lw=0.5);
	aA[2].plot(sAbin,sAbin_DTP,c="k",lw=0.5);
	
	aA[2].plot([1,1]*sAavg_SEA,[0.1,50],c="b")
	aA[2].plot([1,1]*sAavg_TRA,[0.1,50],c="r")
	aA[2].plot([1,1]*sAavg_CRB,[0.1,50],c="blue3")
	aA[2].plot([1,1]*sAavg_DTP,[0.1,50],c="k")
	
	aA[1].format(
		xlim=(minimum(lvec),maximum(lvec)),
		ylim=(0,10),#yscale="log",
		ylabel="Density",
		ltitle="(a) Land"
	)
	
	aA[2].format(
		xlim=(minimum(svec),maximum(svec)),
		ylim=(0.1,15),#yscale="log",
		xlabel="A",
		ltitle="(b) Ocean"
	)
	
	fA.savefig(plotsdir("tccdiurnalamplitude.png"),transparent=false,dpi=200)
	load(plotsdir("tccdiurnalamplitude.png"))
end

# ╔═╡ 1432fa12-57c7-11eb-0606-7be0389e8fb3
md"
#### iii. Phase Shift $\theta$
"

# ╔═╡ 76627730-57c7-11eb-2037-3f608e085a04
begin
	θvec = collect(0:0.5:24); pθbin = (θvec[2:end].+θvec[1:(end-1)])/24*pi
	pθbin = vcat(pθbin,pθbin[1]+2*pi)
	
	lθbin_SEA,lθavg_SEA = bindatasfclnd([20,-15,165,90],θvec,θ,lon,lat,lsm)
	lθbin_JAV,lθavg_JAV = bindatasfclnd([-2,-12,120,100],θvec,θ,lon,lat,lsm)
	lθbin_TRA,lθavg_TRA = bindatasfclnd([10,-10,40,-10],θvec,θ,lon,lat,lsm)
	lθbin_CRB,lθavg_CRB = bindatasfclnd([25,15,-60,-90],θvec,θ,lon,lat,lsm)
	lθbin_AMZ,lθavg_AMZ = bindatasfclnd([10,-10,-45,-75],θvec,θ,lon,lat,lsm)
	
	sθbin_SEA,sθavg_SEA = bindatasfcsea([20,-15,165,90],θvec,θ,lon,lat,lsm)
	sθbin_JAV,sθavg_JAV = bindatasfcsea([-2,-12,120,100],θvec,θ,lon,lat,lsm)
	sθbin_TRA,sθavg_TRA = bindatasfcsea([10,-10,40,-10],θvec,θ,lon,lat,lsm)
	sθbin_CRB,sθavg_CRB = bindatasfcsea([25,15,-60,-90],θvec,θ,lon,lat,lsm)
	sθbin_DTP,sθavg_DTP = bindatasfcsea([10,-10,360,0],θvec,θ,lon,lat,lsm)
	
	md"Binning hour of maximum skin temperature for different tropical regions ..."
end

# ╔═╡ 8d739d0a-57c7-11eb-16b6-736f595e329e
begin
	pplt.close(); fθ,aθ = pplt.subplots(ncols=2,aspect=2,proj="polar");
	
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_SEA,lθbin_SEA[1])),c="b")
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_TRA,lθbin_TRA[1])),c="r")
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_CRB,lθbin_CRB[1])),c="blue3")
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_AMZ,lθbin_AMZ[1])),c="g")
	
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_SEA,sθbin_SEA[1])),c="b");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_TRA,sθbin_TRA[1])),c="r");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_CRB,sθbin_CRB[1])),c="blue3");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_DTP,sθbin_DTP[1])),c="k");
	
	aθ[1].format(theta0="N",thetaformatter="tau",ltitle="(a) Land")
	aθ[2].format(theta0="N",thetaformatter="tau",ltitle="(b) Ocean")
	aθ[1].format(suptitle=L"$\theta$ / Fraction of Day")
	
	fθ.savefig(plotsdir("tccdiurnalphase.png"),transparent=false,dpi=200)
	load(plotsdir("tccdiurnalphase.png"))
end

# ╔═╡ Cell order:
# ╟─90fffbc8-524d-11eb-232a-1bada28d5505
# ╟─6f00b8fc-530c-11eb-2242-99d8544f6e14
# ╟─8f30c56c-530c-11eb-2782-33f3c4ed9e89
# ╟─c1fafcba-530f-11eb-1cc2-67d10a3fa606
# ╟─3565af3c-5311-11eb-34c4-2d228b05b17c
# ╠═fd344560-94d3-11eb-2b79-05c905c5953f
# ╠═a6a688ca-53ab-11eb-2776-b5380ffb26c1
# ╟─aa05317e-530b-11eb-2ec1-93aff65659dd
# ╟─103f85e8-530c-11eb-047d-a537aa60075d
# ╟─49d13e5c-53af-11eb-29ca-c994a7acd377
# ╟─e8141e20-53af-11eb-1a23-81d34293c5eb
# ╟─d82366b0-53b1-11eb-26c1-ff1bb6ccb027
# ╟─bb59b8d6-53b1-11eb-3631-87ef61219c4c
# ╟─5c0e5bae-554e-11eb-3f83-a364ae0a2485
# ╟─68cfc46c-5755-11eb-1702-373942539652
# ╠═52b39ff8-9426-11eb-2a86-43f7da15f62e
# ╟─ea7f0956-575b-11eb-3e3f-a1ba3e08b771
# ╟─5714c13c-575c-11eb-06d4-838b4e8dbcd7
# ╟─c4792bf2-5552-11eb-3b52-997f59fd42f3
# ╟─1fadf4ca-5755-11eb-1ece-a99313019785
# ╟─f752b054-57c1-11eb-117c-ed52464aa25f
# ╟─4b289fa8-57b9-11eb-0923-116c3d9444bb
# ╟─e7ff7ec8-57b9-11eb-0115-abbe4aa9a1a9
# ╟─0fbb0b46-57c2-11eb-365a-a73a2ebda8e4
# ╟─252508a8-57c2-11eb-08b5-8fa673b1ac8a
# ╟─5f58ae9c-57c2-11eb-1f04-2ddbaf2b4f1b
# ╟─1432fa12-57c7-11eb-0606-7be0389e8fb3
# ╟─76627730-57c7-11eb-2037-3f608e085a04
# ╟─8d739d0a-57c7-11eb-16b6-736f595e329e
