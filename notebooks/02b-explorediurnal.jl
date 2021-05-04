### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 712c777e-f782-48d6-87be-ed08a74e333a
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ e04d025e-6a3f-4a9f-8fe1-6034f6505070
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

# ╔═╡ 0a0cc842-964c-11eb-2123-67de7995978d
md"
# 2b. Exploratory Analysis of the Diurnal Cycle of Precipitation

Here, we have a look at some samples of the diurnal cycle of precipitation over different regions.
"

# ╔═╡ 99e672ac-2492-4d06-92dd-d407b44b21ff
md"
### A. Standardizing the Diurnal Cycle
"

# ╔═╡ b2a84c27-fa3a-4cd5-8df0-c4f646696fae
longitude2timeshift(longitude::Real) = longitude / 180 * 12

# ╔═╡ 31b2610b-6191-4336-9825-ceeaaaaa8ff0
function gpm2diurnal(data,longitude)
	
	nlon,nlat,nt = size(data)
	drn = zeros(nlon,nlat,nt)
	
	p0 = [0,0.5,0]; t = ((1:(nt+1)).-0.5)/2
	ts = longitude2timeshift.(longitude)
	it = 0:0.5:24; it = it[1:(end-1)]; nit = length(it)
	
	idat = zeros(nt+1); var  = zeros(nit)
	
	for ilat = 1 : nlat, ilon = 1 : nlon

        tl = t .+ ts[ilon]
		
        idat[1:nt] = data[ilon,ilat,:]
		idat[end]  = data[ilon,ilat,1]
		
		itp = interpolate(idat,BSpline(Cubic(Periodic(OnGrid()))))
        stp = scale(itp,tl)
        etp = extrapolate(stp,Periodic())
		drn[ilon,ilat,:] = etp[it]

    end
	
	return drn
	
end

# ╔═╡ 2c2e948c-1a2d-45c4-9ef8-3810c669a265
md"
### B. Retrieving GPM Precipitation Data
"

# ╔═╡ 654c0f36-a149-4a08-be81-55650a2bcb84
function retrievegpmdiurnal(
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
	drn = gpm2diurnal(gpm,rlon)

    return rlon,rlat,drn

end

# ╔═╡ 26f3bf40-60c7-46d5-8e60-1b72cffd5afd
sroot = "/n/kuangdss01/lab/"; regID = "TRP"

# ╔═╡ 1f7c5b9d-f496-4547-9445-d4de202a9a1f
begin
	ds = NCDataset(datadir("gpmimerg-prcp_rate-$regID.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	μ  = ds["modelmean"][:]
	A  = ds["modelamplitude"][:]
	θ  = ds["modelphase"][:]
	close(ds)
	
	_,_,drn = retrievegpmdiurnal(sroot,regID=regID,timeID=[2001,2018])
	md"Extracting diurnal information for GPM IMERGv6 precipitation ..."
end

# ╔═╡ 5c7bde7f-5f9b-4348-a5fa-ae8e189518a7
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ 6f532b3c-4653-4fb6-86c1-f18b240b636b
prect(N::Real,S::Real,W::Real,E::Real) = [W,E,E,W,W],[S,S,N,N,S]

# ╔═╡ e2398f04-828b-4185-b181-d6240f68ebdb
begin
	N,S,E,W = [15,-10,120,95]
	rbox = prect(N,S,W,E)
end

# ╔═╡ da5aa5b4-563d-4487-834f-aefd6289b58c
begin
	pplt.close(); f,axs = pplt.subplots(nrows=3,axwidth=6,aspect=6)
	
	μlvls = [0.5,0.707,1,1.41,2,3.16,5]/10
	Alvls = [0.1,0.2,0.5,0.9,1.11,2,5,10]
	
	c = axs[1].contourf(lon,lat,μ',cmap="Blues",levels=μlvls,extend="both")
	axs[1].plot(x,y,c="k",lw=0.5)
	axs[1].plot(rbox[1],rbox[2],c="k")
	axs[1].format(urtitle=L"$\mu$ / mm hr$^{-1}$")
	axs[1].colorbar(c,loc="r",ticks=[0.05,0.1,0.2,0.5])
	
	c = axs[2].contourf(lon,lat,(A./μ)',cmap="broc_r",levels=Alvls,extend="both")
	axs[2].plot(x,y,c="k",lw=0.5)
	axs[2].plot(rbox[1],rbox[2],c="k")
	axs[2].format(urtitle=L"A/$\mu$")
	axs[2].colorbar(c,loc="r",ticks=[0.1,0.5,1,2,10])
	
	c = axs[3].pcolormesh(lon,lat,θ',cmap="romaO",levels=0:24)
	axs[3].plot(x,y,c="k",lw=0.5)
	axs[3].plot(rbox[1],rbox[2],c="k")
	axs[3].format(urtitle=L"$\theta$ / Hour of Day")
	axs[3].colorbar(c,loc="r",ticks=6)
	
	for ax in axs
		ax.format(xlim=(0,360),ylim=(-30,30),xlocator=0:60:360)
	end
	
	f.savefig("test.png",transparent=false,dpi=200)
	PNGFiles.load("test.png")
end

# ╔═╡ edf19d47-5ee4-4293-9953-6b34c45113b6
begin
	lds = NCDataset(datadir("GPM_IMERG_LandSeaMask-TRP.nc"))
	lsm = lds["landseamask"][:]*1
	close(lds); rm("test.png")
	
md"Loading Land-Sea Mask for GPM data ..."
end

# ╔═╡ c75539a9-bfe0-499f-a5dd-3172da4fea04
function extractdiurnalreg(coords,data,lsm,lon,lat)
	
	rlon,rlat,rinfo = regiongridvec(coords,lon,lat)
	if maximum(rlon) > 360; rlon .= rlon .- 360 end
	rlsm = regionextractgrid(lsm,rinfo)
	rdrn = regionextractgrid(drn,rinfo)
	
	ldrn = zeros(48); sdrn = zeros(48)
	
	for it = 1 : 48
		rdrnit = @view rdrn[:,:,it]
		ldrn[it] = mean(rdrnit[rlsm.>0.5])
		sdrn[it] = mean(rdrnit[rlsm.<0.5])
	end
	
	return ldrn,sdrn
end

# ╔═╡ bfb34bf8-77c3-4506-8d4c-8111a6d2d127
begin
	ldrn_SEA,sdrn_SEA = extractdiurnalreg([20,-15,165,90],drn,lsm,lon,lat)
	ldrn_CRB,sdrn_CRB = extractdiurnalreg([25,15,-60,-90],drn,lsm,lon,lat)
	ldrn_TRA,sdrn_TRA = extractdiurnalreg([10,-10,40,-10],drn,lsm,lon,lat)
	ldrn_AMZ,sdrn_AMZ = extractdiurnalreg([10,-10,-45,-75],drn,lsm,lon,lat)
	ldrn_DTP,sdrn_DTP = extractdiurnalreg([10,-10,360,0],drn,lsm,lon,lat)
	md"Extracting information for various tropical regions ..."
end

# ╔═╡ b1a43930-9ac8-4200-ad47-d8b09fa31d5b
begin
	pplt.close(); fts,ats = pplt.subplots(ncols=2,axwidth=2.5,aspect=2)
	
	ats[1].plot(0:0.5:24,vcat(ldrn_SEA,ldrn_SEA[1]),c="b")
	ats[1].plot(0:0.5:24,vcat(ldrn_CRB,ldrn_CRB[1]),c="Blue3")
	ats[1].plot(0:0.5:24,vcat(ldrn_TRA,ldrn_TRA[1]),c="r")
	ats[1].plot(0:0.5:24,vcat(ldrn_AMZ,ldrn_AMZ[1]),c="g")
	ats[1].plot(0:0.5:24,vcat(ldrn_DTP,ldrn_DTP[1]),c="k")
	ats[1].format(xlim=(0,24),ylim=(0,1),xlocator=0:3:24,ultitle="(a) Land")
	p1 = ats[1].panel(side="r",space=0,width="2em")
	p1.scatter(0.5,mean(ldrn_SEA),c="b",s=5)
	p1.scatter(0.5,mean(ldrn_CRB),c="Blue3",s=5)
	p1.scatter(0.5,mean(ldrn_TRA),c="r",s=5)
	p1.scatter(0.5,mean(ldrn_AMZ),c="g",s=5)
	p1.scatter(0.5,mean(ldrn_DTP),c="k",s=5)
	p1.format(xticks=[])
	
	ats[2].plot(0:0.5:24,vcat(sdrn_SEA,sdrn_SEA[1]),c="b")
	ats[2].plot(0:0.5:24,vcat(sdrn_CRB,sdrn_CRB[1]),c="Blue3")
	ats[2].plot(0:0.5:24,vcat(sdrn_TRA,sdrn_TRA[1]),c="r")
	ats[2].plot(0:0.5:24,vcat(sdrn_AMZ,sdrn_AMZ[1]),c="g")
	ats[2].plot(0:0.5:24,vcat(sdrn_DTP,sdrn_DTP[1]),c="k")
	ats[2].format(xlim=(0,24),ylim=(0,1),xlocator=0:3:24,ultitle="(b) Ocean")
	p2 = ats[2].panel(side="r",space=0,width="2em")
	p2.scatter(0.5,mean(sdrn_SEA),c="b",s=5)
	p2.scatter(0.5,mean(sdrn_CRB),c="Blue3",s=5)
	p2.scatter(0.5,mean(sdrn_TRA),c="r",s=5)
	p2.scatter(0.5,mean(sdrn_AMZ),c="g",s=5)
	p2.scatter(0.5,mean(sdrn_DTP),c="k",s=5)
	p2.format(xticks=[])
	
	for ax in ats
		ax.format(
			xlim=(0,24),ylim=(0,0.75),
			ylabel=L"Rain Rate / mm hr$^{-1}$",
			xlabel="Hour of Day",
			suptitle="Diurnal Cycle of Rainfall"
		)
	end
	
	fts.savefig("testdrn.png",transparent=false,dpi=200)
	PNGFiles.load("testdrn.png")
end

# ╔═╡ Cell order:
# ╟─0a0cc842-964c-11eb-2123-67de7995978d
# ╟─712c777e-f782-48d6-87be-ed08a74e333a
# ╟─e04d025e-6a3f-4a9f-8fe1-6034f6505070
# ╟─99e672ac-2492-4d06-92dd-d407b44b21ff
# ╠═b2a84c27-fa3a-4cd5-8df0-c4f646696fae
# ╠═31b2610b-6191-4336-9825-ceeaaaaa8ff0
# ╟─2c2e948c-1a2d-45c4-9ef8-3810c669a265
# ╠═654c0f36-a149-4a08-be81-55650a2bcb84
# ╠═26f3bf40-60c7-46d5-8e60-1b72cffd5afd
# ╠═1f7c5b9d-f496-4547-9445-d4de202a9a1f
# ╟─5c7bde7f-5f9b-4348-a5fa-ae8e189518a7
# ╠═6f532b3c-4653-4fb6-86c1-f18b240b636b
# ╠═e2398f04-828b-4185-b181-d6240f68ebdb
# ╟─da5aa5b4-563d-4487-834f-aefd6289b58c
# ╟─edf19d47-5ee4-4293-9953-6b34c45113b6
# ╠═c75539a9-bfe0-499f-a5dd-3172da4fea04
# ╟─bfb34bf8-77c3-4506-8d4c-8111a6d2d127
# ╟─b1a43930-9ac8-4200-ad47-d8b09fa31d5b
