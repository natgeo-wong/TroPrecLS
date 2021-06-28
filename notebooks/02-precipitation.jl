### A Pluto.jl notebook ###
# v0.14.7

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
	using ClimateSatellite
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
# 2. Investigation of Precipitation Data

We first investigate GPM IMERGv6 precipitation data from 2001-2019 over the domains listed in the notebook `01-domains`.  The data was retrieved from the website of NASA's Precipitation Measurement Mission (PMM), and is found in 0.1º x 0.1º horizontal resolution every half-hour.

The data was downloaded using the `ClimateSatellite.jl` package, which downloads the data in monthly batches in the date range specified.  More information can be found in the package documentation, and instructions on how to use the package to download the data can be found there as well.
"

# ╔═╡ c1fafcba-530f-11eb-1cc2-67d10a3fa606
md"
### A. Characteristics of the Diurnal Cycle of Precipitation

We wish to find the following characteristics of the diurnal cycle of precipitation rate:
* The mean $\mu$ rate of precipitation 
* The amplitude $A$ of the diurnal cycle (max-min)/2 of precipitation
* The hour $\theta$ at which the rate of precipitation is a maximum
"

# ╔═╡ 3565af3c-5311-11eb-34c4-2d228b05b17c
md"
GPS IMERGv6 output time is defined according to UTC.  Therefore, for every grid point, we must also correct to the local time before fitting the data to the model.
"

# ╔═╡ 577a59d8-5311-11eb-2e2e-53bbeff12648
longitude2timeshift(longitude::Real) = longitude / 180 * 12

# ╔═╡ a6a688ca-53ab-11eb-2776-b5380ffb26c1
function gpm2model(data,longitude)
	
	nlon,nlat,nt = size(data)
	θ = zeros(nlon,nlat)
	A = zeros(nlon,nlat)
	μ = zeros(nlon,nlat)
	
	p0 = [0,0.5,0]; t = ((1:(nt+1)).-0.5)/2
	ts = longitude2timeshift.(longitude)
	it = 0:0.1:24; it = it[1:(end-1)]; nit = length(it)
	
	idat = zeros(nt+1); var  = zeros(nit)
	
	for ilat = 1 : nlat, ilon = 1 : nlon

        tl = t .+ ts[ilon]
		
        idat[1:nt] = data[ilon,ilat,:]
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
### B. Retrieving GPM Precipitation Data
"

# ╔═╡ bb90be66-554c-11eb-05de-a5db553ad4b1
md"
We fit the model to GPM IMERGv6 data from 2001 to 2018 for each individual spatial point (i.e. at 0.1º resolution).  Our results are plotted below:
"

# ╔═╡ a116023a-53ad-11eb-25e0-d5dfa8338a1b
function savediurnalphase(regID,rlon,rlat,μ,A,θ)
	
	fnc = datadir("gpmimerg-prcp_rate-$regID.nc")
	if isfile(fnc); rm(fnc) end
	ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));
	
	ds.dim["longitude"] = length(rlon)
    ds.dim["latitude"]  = length(rlat)
	
	nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))
	
	ncμ = defVar(ds,"modelmean",Float32,("longitude","latitude",),attrib = Dict(
        "units" => "kg m**-2 s**-1"
    ))
	
	ncA = defVar(ds,"modelamplitude",Float32,("longitude","latitude",),attrib = Dict(
        "units" => "kg m**-2 s**-1"
    ))
	
	ncθ = defVar(ds,"modelphase",Float32,("longitude","latitude",),attrib = Dict(
        "units" => "hr"
    ))
	
	nclongitude[:] = rlon
	nclatitude[:]  = rlat
	ncμ[:] = μ
	ncA[:] = A
	ncθ[:] = θ
	
	close(ds)
	
end

# ╔═╡ 103f85e8-530c-11eb-047d-a537aa60075d
function retrievegpm(
    sroot::AbstractString;
    regID::AbstractString="GLB",
	timeID::Vector
)

    glon,glat = gpmlonlat()
	rlon,rlat,_ = gregiongridvec(regID,glon,glat)
    nlon = length(rlon)
	nlat = length(rlat)

    gpm  = zeros(nlon,nlat,48)

    nt = timeID[2] - timeID[1] + 1
    for yr = timeID[1] : timeID[2]

        rds,rvar = clisatanaread(
            "gpmimerg","prcp_rate","domain_yearly_mean_hourly",
            Date(yr),regID,path=sroot
        );
        gpm += rvar[:]*3600; close(rds);

    end
	
	gpm .= gpm / nt
	μ,A,θ = gpm2model(gpm,rlon)
	
	savediurnalphase(regID,rlon,rlat,μ,A,θ)

    return rlon,rlat,μ,A,θ

end

# ╔═╡ 49d13e5c-53af-11eb-29ca-c994a7acd377
sroot = "/n/kuangdss01/lab/"; regID = "TRP"

# ╔═╡ e8141e20-53af-11eb-1a23-81d34293c5eb
begin
	if isfile(datadir("gpmimerg-prcp_rate-$regID.nc"))
		ds = NCDataset(datadir("gpmimerg-prcp_rate-$regID.nc"))
		lon = ds["longitude"][:]
		lat = ds["latitude"][:]
		μ  = ds["modelmean"][:]
		A  = ds["modelamplitude"][:]
		θ  = ds["modelphase"][:]
		close(ds)
	else lon,lat,μt,μ,A,θ = retrievegpm(sroot,regID=regID,timeID=[2001,2018])
	end
	md"Extracting diurnal information for GPM IMERGv6 precipitation ..."
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
	
	μlvls = [0.5,0.707,1,1.41,2,3.16,5]/10
	Alvls = [0.1,0.2,0.5,0.9,1.11,2,5,10]
	
	c = axs[1].contourf(lon,lat,μ',cmap="Blues",levels=μlvls,extend="both")
	axs[1].plot(x,y,c="k",lw=0.5)
	axs[1].format(urtitle=L"$\mu$ / mm hr$^{-1}$")
	axs[1].colorbar(c,loc="r",ticks=[0.05,0.1,0.2,0.5])
	
	c = axs[2].contourf(lon,lat,(A./μ)',cmap="broc_r",levels=Alvls,extend="both")
	axs[2].plot(x,y,c="k",lw=0.5)
	axs[2].format(urtitle=L"A/$\mu$")
	axs[2].colorbar(c,loc="r",ticks=[0.1,0.5,1,2,10])
	
	c = axs[3].pcolormesh(lon,lat,θ',cmap="romaO",levels=0:24)
	axs[3].plot(x,y,c="k",lw=0.5)
	axs[3].format(urtitle=L"$\theta$ / Hour of Day")
	axs[3].colorbar(c,loc="r",ticks=6)
	
	for ax in axs
		ax.format(
			xlim=(0,360),ylim=(-30,30),xlocator=0:60:360,
			xlabel=L"Longitude / $\degree$",
			ylabel=L"Latitude / $\degree$",
		)
	end
	
	f.savefig(plotsdir("gpmspatial_TRP.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("gpmspatial_TRP.png"))
end

# ╔═╡ 5c0e5bae-554e-11eb-3f83-a364ae0a2485
md"
The precipitation band is concentrated along a narrow band in the tropical oceans that represents the ITCZ and SPCZ, except around the Indo-Pacific warmpool, where we see that regions of high precipitation extend farther in latitude about the equator.

We see that the presence of land also serves to extend the distribution of precipitation over latitude.  Furthermore, the amplitude $A$ of the diurnal cycle is enhanced over land as compared to over the ocean for similar values of mean precipitation rate.

Lastly, we also see that the peak of precipitation rainfall is drastically different over land compared to the ocean.  Over land, the precipitation rate generally peaks from 1800 hours to midnight, while over the tropical ocean the rainfall peaks in the early morning..
"

# ╔═╡ 68cfc46c-5755-11eb-1702-373942539652
md"
### C. Regional Analysis

We can get quick snapshots of the results for different GeoRegions specified in this project.
"

# ╔═╡ 271c15e6-9426-11eb-3c64-912620b1a8ce
begin
	geo = GeoRegion("BRN")
	md"Defining region coordinates ..."
end

# ╔═╡ ea7f0956-575b-11eb-3e3f-a1ba3e08b771
begin
	N,S,E,W = geo.N,geo.S,geo.E,geo.W
	ggrd = RegionGrid(geo,lon,lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
	rμ = zeros(nlon,nlat)
	rA = zeros(nlon,nlat)
	rθ = zeros(nlon,nlat)
	if typeof(ggrd) <: PolyGrid
		mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end
	for glat in 1 : nlat, glon in 1 : nlon
		rμ[glon,glat] = μ[ilon[glon],ilat[glat]] * mask[glon,glat]
		rA[glon,glat] = A[ilon[glon],ilat[glat]] * mask[glon,glat]
		rθ[glon,glat] = θ[ilon[glon],ilat[glat]] * mask[glon,glat]
	end
	md"Extracting information for region ..."
end

# ╔═╡ 5714c13c-575c-11eb-06d4-838b4e8dbcd7
begin
	asp = (E-W+2)/(N-S+2)
	pplt.close()
	if asp > 1.5
		freg,areg = pplt.subplots(nrows=3,axwidth=asp*1.2,aspect=asp)
	else
		freg,areg = pplt.subplots(ncols=3,axwidth=2,aspect=asp)
	end
	
	creg = areg[1].contourf(
		ggrd.glon,ggrd.glat,rμ',
		cmap="Blues",levels=μlvls,extend="both"
	)
	areg[1].plot(x,y,c="k",lw=0.5)
	areg[1].format(rtitle=L"$\mu$ / mm hr$^{-1}$")
	areg[1].colorbar(creg,loc="r")
	
	creg = areg[2].contourf(
		ggrd.glon,ggrd.glat,(rA./rμ)',
		cmap="broc_r",levels=Alvls,extend="both"
	)
	areg[2].plot(x,y,c="k",lw=0.5)
	areg[2].format(rtitle=L"A / $\mu$")
	areg[2].colorbar(creg,loc="r")
	
	creg = areg[3].pcolormesh(ggrd.glon,ggrd.glat,rθ',cmap="romaO",levels=0:0.5:24)
	areg[3].plot(x,y,c="k",lw=0.5)
	areg[3].format(rtitle=L"$\theta$ / Hour of Day")
	areg[3].colorbar(creg,loc="r",ticks=6)
	
	for ax in areg
		ax.format(
			xlim=(ggrd.glon[1].-1,ggrd.glon[end].+1),
			xlabel=L"Longitude / $\degree$",
			ylim=(S-1,N+1),ylabel=L"Latitude / $\degree$",
			grid=true
		)
	end
	
	freg.savefig(plotsdir("gpmspatial_$(geo.regID).png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("gpmspatial_$(geo.regID).png"))
end

# ╔═╡ c4792bf2-5552-11eb-3b52-997f59fd42f3
md"
### D. Binning of Precipitation Data

We now wish to bin the modelled precipitation data in order to determine the relationship between precipitation rate and the diurnal cycle over land and sea.  Is it different?
"

# ╔═╡ 1fadf4ca-5755-11eb-1ece-a99313019785
begin
	lds = NCDataset(datadir("GPM_IMERG_LandSeaMask-TRP.nc"))
	lsm = lds["landseamask"][:]*1
	close(lds)
	
md"Loading Land-Sea Mask for GPM data ..."
end

# ╔═╡ ffeb8103-85a1-43b6-9cfa-689d157b326a
begin
	lsc = pplt.Colors("Delta_r",15)
	md"Colours for different regions ..."
end

# ╔═╡ f752b054-57c1-11eb-117c-ed52464aa25f
md"
#### i. Mean Precipitation Rate $\mu$ (Hourly)
"

# ╔═╡ 4b289fa8-57b9-11eb-0923-116c3d9444bb
begin
	lvec = collect(-3:0.02:0);
	lbins = 10 .^lvec; lpbin = 10 .^((lvec[2:end].+lvec[1:(end-1)])/2)
	
	lbin_DTP,lavg_DTP = bindatasfclnd(GeoRegion("DTP"),lbins,μ,lon,lat,lsm)
	lbin_SEA,lavg_SEA = bindatasfclnd(GeoRegion("SEA"),lbins,μ,lon,lat,lsm)
	lbin_CRB,lavg_CRB = bindatasfclnd(GeoRegion("CRB"),lbins,μ,lon,lat,lsm)
	lbin_TRA,lavg_TRA = bindatasfclnd(GeoRegion("TRA"),lbins,μ,lon,lat,lsm)
	lbin_AMZ,lavg_AMZ = bindatasfclnd(GeoRegion("AMZ"),lbins,μ,lon,lat,lsm)
	
	svec = collect(-3:0.02:0);
	sbins = 10 .^svec; spbin = 10 .^((svec[2:end].+svec[1:(end-1)])/2)
	
	sbin_DTP,savg_DTP = bindatasfcsea(GeoRegion("DTP"),sbins,μ,lon,lat,lsm)
	sbin_SEA,savg_SEA = bindatasfcsea(GeoRegion("SEA"),sbins,μ,lon,lat,lsm)
	sbin_CRB,savg_CRB = bindatasfcsea(GeoRegion("CRB"),sbins,μ,lon,lat,lsm)
	sbin_EPO,savg_EPO = bindatasfcsea(GeoRegion("AR6_EPO"),sbins,μ,lon,lat,lsm)
	sbin_EIO,savg_EIO = bindatasfcsea(GeoRegion("AR6_EIO"),sbins,μ,lon,lat,lsm)
	sbin_EAO,savg_EAO = bindatasfcsea(GeoRegion("AR6_EAO"),sbins,μ,lon,lat,lsm)
	
	md"Binning mean precipitation rates for different tropical regions ..."
end

# ╔═╡ e7ff7ec8-57b9-11eb-0115-abbe4aa9a1a9
begin
	pplt.close(); fbin,abin = pplt.subplots(ncols=2,aspect=2,axwidth=3);
	
	lgd = Dict("ncol"=>1,"frame"=>false)
	abin[1].plot(lpbin,lbin_DTP,c="k",lw=0.5);
	abin[1].plot(lpbin,lbin_CRB,c=lsc[10],lw=0.5)
	abin[1].plot(lpbin,lbin_SEA,c=lsc[5],lw=0.5)
	abin[1].plot(lpbin,lbin_AMZ,c=lsc[4],lw=0.5)
	abin[1].plot(lpbin,lbin_TRA,c=lsc[3],lw=0.5)
	
	abin[1].plot([1,1]*lavg_DTP,[0.1,50],c="k")
	abin[1].plot([1,1]*lavg_CRB,[0.1,50],c=lsc[10])
	abin[1].plot([1,1]*lavg_SEA,[0.1,50],c=lsc[5])
	abin[1].plot([1,1]*lavg_AMZ,[0.1,50],c=lsc[4])
	abin[1].plot([1,1]*lavg_TRA,[0.1,50],c=lsc[3])
	
	abin[2].plot(spbin,sbin_DTP,c="k",lw=0.5);
	abin[2].plot(spbin,sbin_EPO,c=lsc[13],lw=0.5);
	abin[2].plot(spbin,sbin_EIO,c=lsc[12],lw=0.5);
	abin[2].plot(spbin,sbin_EAO,c=lsc[11],lw=0.5);
	abin[2].plot(spbin,sbin_CRB,c=lsc[10],lw=0.5);
	abin[2].plot(spbin,sbin_SEA,c=lsc[5],lw=0.5);
	
	abin[2].plot([1,1]*savg_DTP,[0.1,50],c="k",label="DTP",legend="r",legend_kw=lgd)
	abin[2].plot([1,1]*savg_EPO,[0.1,50],c=lsc[13],label="AR6_EPO",legend="r")
	abin[2].plot([1,1]*savg_EIO,[0.1,50],c=lsc[12],label="AR6_EIO",legend="r")
	abin[2].plot([1,1]*savg_EAO,[0.1,50],c=lsc[11],label="AR6_EAO",legend="r")
	abin[2].plot([1,1]*savg_CRB,[0.1,50],c=lsc[10],label="CRB",legend="r")
	abin[2].plot([1,1]*savg_SEA,[0.1,50],c=lsc[5],label="SEA",legend="r")
	abin[2].plot([1,1]*NaN,[0.1,50],c=lsc[4],label="AMZ",legend="r")
	abin[2].plot([1,1]*NaN,[0.1,50],c=lsc[3],label="TRA",legend="r")
	
	abin[1].format(
		xlim=(minimum(lbins),maximum(lbins)),xscale="log",
		ylim=(0,20),#yscale="log",
		ylabel="Density",
		ltitle="(a) Land"
	)
	
	abin[2].format(
		xlim=(minimum(sbins),maximum(sbins)),xscale="log",
		ylim=(0,20),#yscale="log",
		xlabel=L"Precipitation Rate / mm hr$^{-1}$",
		ltitle="(b) Ocean"
	)
	
	fbin.savefig(plotsdir("gpmmean.png"),transparent=false,dpi=200)
	load(plotsdir("gpmmean.png"))
end

# ╔═╡ 0fbb0b46-57c2-11eb-365a-a73a2ebda8e4
md"
#### ii. Diurnal Amplitude $A$ vs $\mu$
"

# ╔═╡ 252508a8-57c2-11eb-08b5-8fa673b1ac8a
begin
	Avec = collect(0:0.02:2.5); pAbin = (Avec[2:end].+Avec[1:(end-1)])/2
	
	lAbin_DTP,lAavg_DTP = bindatasfclnd(GeoRegion("DTP"),Avec,A./μ,lon,lat,lsm)
	lAbin_SEA,lAavg_SEA = bindatasfclnd(GeoRegion("SEA"),Avec,A./μ,lon,lat,lsm)
	lAbin_CRB,lAavg_CRB = bindatasfclnd(GeoRegion("CRB"),Avec,A./μ,lon,lat,lsm)
	lAbin_TRA,lAavg_TRA = bindatasfclnd(GeoRegion("TRA"),Avec,A./μ,lon,lat,lsm)
	lAbin_AMZ,lAavg_AMZ = bindatasfclnd(GeoRegion("AMZ"),Avec,A./μ,lon,lat,lsm)
	
	sAbin_DTP,sAavg_DTP = bindatasfcsea(GeoRegion("DTP"),Avec,A./μ,lon,lat,lsm)
	sAbin_SEA,sAavg_SEA = bindatasfcsea(GeoRegion("SEA"),Avec,A./μ,lon,lat,lsm)
	sAbin_CRB,sAavg_CRB = bindatasfcsea(GeoRegion("CRB"),Avec,A./μ,lon,lat,lsm)
	sAbin_EPO,sAavg_EPO = bindatasfcsea(GeoRegion("AR6_EPO"),Avec,A./μ,lon,lat,lsm)
	sAbin_EIO,sAavg_EIO = bindatasfcsea(GeoRegion("AR6_EIO"),Avec,A./μ,lon,lat,lsm)
	sAbin_EAO,sAavg_EAO = bindatasfcsea(GeoRegion("AR6_EAO"),Avec,A./μ,lon,lat,lsm)
	
	md"Binning diurnal amplitude precipitation rates for different tropical regions as a ratio of mean precipitation rate ..."
end

# ╔═╡ 5f58ae9c-57c2-11eb-1f04-2ddbaf2b4f1b
begin
	pplt.close(); fA,aA = pplt.subplots(ncols=2,aspect=2,axwidth=3);
	
	aA[1].plot(pAbin,lAbin_DTP,c="k",lw=0.5)
	aA[1].plot(pAbin,lAbin_SEA,c=lsc[10],lw=0.5)
	aA[1].plot(pAbin,lAbin_CRB,c=lsc[5],lw=0.5)
	aA[1].plot(pAbin,lAbin_AMZ,c=lsc[4],lw=0.5)
	aA[1].plot(pAbin,lAbin_TRA,c=lsc[3],lw=0.5)
	
	aA[1].plot([1,1]*lAavg_DTP,[0.1,50],c="k")
	aA[1].plot([1,1]*lAavg_SEA,[0.1,50],c=lsc[10])
	aA[1].plot([1,1]*lAavg_CRB,[0.1,50],c=lsc[5])
	aA[1].plot([1,1]*lAavg_AMZ,[0.1,50],c=lsc[4])
	aA[1].plot([1,1]*lAavg_TRA,[0.1,50],c=lsc[3])
	
	aA[2].plot(pAbin,sAbin_DTP,c="k",lw=0.5);
	aA[2].plot(pAbin,sAbin_EPO,c=lsc[13],lw=0.5);
	aA[2].plot(pAbin,sAbin_EIO,c=lsc[12],lw=0.5);
	aA[2].plot(pAbin,sAbin_EAO,c=lsc[11],lw=0.5);
	aA[2].plot(pAbin,sAbin_CRB,c=lsc[10],lw=0.5);
	aA[2].plot(pAbin,sAbin_SEA,c=lsc[5],lw=0.5);
	
	aA[2].plot([1,1]*sAavg_DTP,[0.1,50],c="k",label="DTP",legend="r",legend_kw=lgd)
	aA[2].plot([1,1]*sAavg_EPO,[0.1,50],c=lsc[13],label="AR6_EPO",legend="r")
	aA[2].plot([1,1]*sAavg_EIO,[0.1,50],c=lsc[12],label="AR6_EIO",legend="r")
	aA[2].plot([1,1]*sAavg_EAO,[0.1,50],c=lsc[11],label="AR6_EAO",legend="r")
	aA[2].plot([1,1]*sAavg_CRB,[0.1,50],c=lsc[10],label="CRB",legend="r")
	aA[2].plot([1,1]*sAavg_SEA,[0.1,50],c=lsc[5],label="SEA",legend="r")
	aA[2].plot([1,1]*NaN,[0.1,50],c=lsc[4],label="TRA",legend="r")
	aA[2].plot([1,1]*NaN,[0.1,50],c=lsc[3],label="AMZ",legend="r")
	
	aA[1].format(
		xlim=(minimum(Avec),maximum(Avec)),#xscale="log",
		ylim=(0,10),#yscale="log",
		ylabel="Density",
		ltitle="(a) Land"
	)
	
	aA[2].format(
		xlim=(minimum(Avec),maximum(Avec)),
		ylim=(0,15),#yscale="log",
		xlabel=L"A/$\mu$",
		ltitle="(b) Ocean"
	)
	
	fA.savefig(plotsdir("gpmdiurnalamplitude.png"),transparent=false,dpi=200)
	load(plotsdir("gpmdiurnalamplitude.png"))
end

# ╔═╡ 1432fa12-57c7-11eb-0606-7be0389e8fb3
md"
#### iii. Phase Shift $\theta$
"

# ╔═╡ 76627730-57c7-11eb-2037-3f608e085a04
begin
	θvec = collect(0:0.5:24); pθbin = (θvec[2:end].+θvec[1:(end-1)])/24*pi
	pθbin = vcat(pθbin,pθbin[1]+2*pi)
	
	lθbin_DTP,lθavg_DTP = bindatasfclnd(GeoRegion("DTP"),θvec,θ,lon,lat,lsm)
	lθbin_CRB,lθavg_CRB = bindatasfclnd(GeoRegion("CRB"),θvec,θ,lon,lat,lsm)
	lθbin_SEA,lθavg_SEA = bindatasfclnd(GeoRegion("SEA"),θvec,θ,lon,lat,lsm)
	lθbin_TRA,lθavg_TRA = bindatasfclnd(GeoRegion("TRA"),θvec,θ,lon,lat,lsm)
	lθbin_AMZ,lθavg_AMZ = bindatasfclnd(GeoRegion("AMZ"),θvec,θ,lon,lat,lsm)
	
	sθbin_DTP,sθavg_DTP = bindatasfcsea(GeoRegion("DTP"),θvec,θ,lon,lat,lsm)
	sθbin_EPO,sθavg_EPO = bindatasfcsea(GeoRegion("AR6_EPO"),θvec,θ,lon,lat,lsm)
	sθbin_EIO,sθavg_EIO = bindatasfcsea(GeoRegion("AR6_EIO"),θvec,θ,lon,lat,lsm)
	sθbin_EAO,sθavg_EAO = bindatasfcsea(GeoRegion("AR6_EAO"),θvec,θ,lon,lat,lsm)
	sθbin_CRB,sθavg_CRB = bindatasfcsea(GeoRegion("CRB"),θvec,θ,lon,lat,lsm)
	sθbin_SEA,sθavg_SEA = bindatasfcsea(GeoRegion("SEA"),θvec,θ,lon,lat,lsm)
	
	md"Binning hour of maximum skin temperature for different tropical regions ..."
end

# ╔═╡ 8d739d0a-57c7-11eb-16b6-736f595e329e
begin
	pplt.close(); fθ,aθ = pplt.subplots(ncols=2,aspect=2,proj="polar");
	
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_DTP,lθbin_DTP[1])),c="k")
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_CRB,lθbin_CRB[1])),c=lsc[10])
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_SEA,lθbin_SEA[1])),c=lsc[5])
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_TRA,lθbin_TRA[1])),c=lsc[4])
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_AMZ,lθbin_AMZ[1])),c=lsc[3])
	
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_DTP,sθbin_DTP[1])),c="k",label="DTP",legend="r");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_EPO,sθbin_EPO[1])),c=lsc[13],label="AR6_EPO",legend="r");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_EIO,sθbin_EIO[1])),c=lsc[12],label="AR6_EIO",legend="r");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_EAO,sθbin_EAO[1])),c=lsc[11],label="AR6_EAO",legend="r");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_CRB,sθbin_CRB[1])),c=lsc[10],label="CRB",legend="r");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_SEA,sθbin_SEA[1])),c=lsc[5],label="SEA",legend="r");
	aθ[2].plot(pθbin,pθbin*NaN,c=lsc[4],label="TRA",legend="r",legend_kw=lgd)
	aθ[2].plot(pθbin,pθbin*NaN,c=lsc[3],label="AMZ",legend="r")
	
	aθ[1].format(ltitle="(a) Land")
	aθ[2].format(ltitle="(b) Ocean")
	aθ[1].format(suptitle=L"$\theta$ / Fraction of Day")
	
	for ax in aθ
		ax.format(
			theta0="N",thetaformatter="tau",
			rlim=(0,4),rlocator=0:4
		)
	end
	
	fθ.savefig(plotsdir("gpmdiurnalphase.png"),transparent=false,dpi=200)
	load(plotsdir("gpmdiurnalphase.png"))
end

# ╔═╡ Cell order:
# ╟─90fffbc8-524d-11eb-232a-1bada28d5505
# ╟─6f00b8fc-530c-11eb-2242-99d8544f6e14
# ╟─8f30c56c-530c-11eb-2782-33f3c4ed9e89
# ╟─c1fafcba-530f-11eb-1cc2-67d10a3fa606
# ╟─3565af3c-5311-11eb-34c4-2d228b05b17c
# ╠═577a59d8-5311-11eb-2e2e-53bbeff12648
# ╠═a6a688ca-53ab-11eb-2776-b5380ffb26c1
# ╟─aa05317e-530b-11eb-2ec1-93aff65659dd
# ╟─bb90be66-554c-11eb-05de-a5db553ad4b1
# ╠═103f85e8-530c-11eb-047d-a537aa60075d
# ╟─a116023a-53ad-11eb-25e0-d5dfa8338a1b
# ╠═49d13e5c-53af-11eb-29ca-c994a7acd377
# ╟─e8141e20-53af-11eb-1a23-81d34293c5eb
# ╟─d82366b0-53b1-11eb-26c1-ff1bb6ccb027
# ╟─bb59b8d6-53b1-11eb-3631-87ef61219c4c
# ╟─5c0e5bae-554e-11eb-3f83-a364ae0a2485
# ╟─68cfc46c-5755-11eb-1702-373942539652
# ╠═271c15e6-9426-11eb-3c64-912620b1a8ce
# ╟─ea7f0956-575b-11eb-3e3f-a1ba3e08b771
# ╟─5714c13c-575c-11eb-06d4-838b4e8dbcd7
# ╟─c4792bf2-5552-11eb-3b52-997f59fd42f3
# ╟─1fadf4ca-5755-11eb-1ece-a99313019785
# ╟─ffeb8103-85a1-43b6-9cfa-689d157b326a
# ╟─f752b054-57c1-11eb-117c-ed52464aa25f
# ╟─4b289fa8-57b9-11eb-0923-116c3d9444bb
# ╟─e7ff7ec8-57b9-11eb-0115-abbe4aa9a1a9
# ╟─0fbb0b46-57c2-11eb-365a-a73a2ebda8e4
# ╟─252508a8-57c2-11eb-08b5-8fa673b1ac8a
# ╟─5f58ae9c-57c2-11eb-1f04-2ddbaf2b4f1b
# ╟─1432fa12-57c7-11eb-0606-7be0389e8fb3
# ╟─76627730-57c7-11eb-2037-3f608e085a04
# ╟─8d739d0a-57c7-11eb-16b6-736f595e329e
