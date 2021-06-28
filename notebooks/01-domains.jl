### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ bcfd5bd8-51f5-11eb-1d79-ab69c65febaf
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ c1fdcad0-51f5-11eb-2df9-b1f4c8dca09d
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using GeoRegions
	using NCDatasets
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ f9606c7e-51f4-11eb-2e24-d998e4a91b9a
md"
# 1. Tropical Domains to be Investigated

Here, we plot the different domains to be investigated in our study and exploration of reanalysis data.
"

# ╔═╡ e58cf4ca-51f8-11eb-218d-8d3141dc8289
md"
### A. Regions of Interest

We consider the following domains / regions to be of interest in our study with the Tropics:
* Deep Tropics (DTP)
* Southeast Asia / Maritime Continent (SEA)
* Tropical Rainforest Africa (TRA)
* Amazon Rainforest (AMZ)
* Caribbean Islands (CRB)

These regions were defined as new `GeoRegions` in the files \"src/addgeorect.txt\" and \"src/addgeopoly.txt\".
"

# ╔═╡ 3f3a7695-ce0b-444b-a8b2-29421050d828
begin
	addGeoRegions(srcdir("addgeorect.txt"))
	addGeoRegions(srcdir("addgeopoly.txt"))
	
md"Adding new custom GeoRegions of interest from text files ..."
end

# ╔═╡ a02d93e9-c53e-4a82-bdf8-34d38377ff99
begin
	blnSEA,bltSEA,slnSEA,sltSEA = coordGeoRegion(GeoRegion("SEA"))
	blnTRA,bltTRA,slnTRA,sltTRA = coordGeoRegion(GeoRegion("TRA"))
	blnAMZ,bltAMZ,slnAMZ,sltAMZ = coordGeoRegion(GeoRegion("AMZ"))
	blnCRB,bltCRB = coordGeoRegion(GeoRegion("CRB"))
	blnEPO,bltEPO,slnEPO,sltEPO = coordGeoRegion(GeoRegion("AR6_EPO"))
	blnEIO,bltEIO,slnEIO,sltEIO = coordGeoRegion(GeoRegion("AR6_EIO"))
	blnEAO,bltEAO,slnEAO,sltEAO = coordGeoRegion(GeoRegion("AR6_EAO"))
	md"Loading shapes and bounds of GeoRegions of interest ..."
end

# ╔═╡ f3c8ffe6-51f5-11eb-362d-293e727b14c7
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ 1af88ca2-6a55-481e-a789-705e348d53da
begin
	lsc = pplt.Colors("Delta_r",15)
	md"Colours for different regions ..."
end

# ╔═╡ a4458db2-5248-11eb-3ecf-b9bcdde2ec37
begin
	ds  = NCDataset(srcdir("koppen.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	kctrop = ds["koppenclass_tropics"][:]
	kcarid = ds["koppenclass_arid"][:]
	kctemph = ds["koppenclass_temperatehumid"][:]
	kctemps = ds["koppenclass_temperatedrysummer"][:]
	kctempw = ds["koppenclass_temperatedrywinter"][:]
	close(ds)
	
md"
We also load the Koeppen Climate classification (transformed into gridded data), and plot it along with the domains shown below.
"
end

# ╔═╡ 170ae0f0-51f6-11eb-153c-e59511b6a82d
begin
	pplt.close(); f,axs = pplt.subplots(aspect=6,axwidth=6);
	
	axs[1].pcolormesh(lon,lat,kctrop',levels=0:6,cmap="Reds_r")
	axs[1].pcolormesh(lon,lat,kcarid',levels=5:9,cmap="Yellow3_r")
	axs[1].pcolormesh(lon,lat,kctemph',levels=6:12,cmap="Green2_r")
	axs[1].pcolormesh(lon,lat,kctemps',levels=11:15,cmap="Brown1_r")
	axs[1].pcolormesh(lon,lat,kctempw',levels=12:18,cmap="Blue3_r")
	
	axs[1].plot(x,y,c="k",lw=0.2)
	axs[1].plot([-150,210,210,-150,-150],[-10,-10,10,10,-10],c="k",lw=1,linestyle="--")
	axs[1].plot(slnEPO.-360,sltEPO,c=lsc[13],lw=1,linestyle="--")
	axs[1].plot(slnEPO,sltEPO,c=lsc[13],lw=1,linestyle="--")
	axs[1].plot(slnEIO,sltEIO,c=lsc[12],lw=1,linestyle="--")
	axs[1].plot(slnEAO,sltEAO,c=lsc[11],lw=1,linestyle="--")
	axs[1].plot(blnCRB,bltCRB,c=lsc[10],lw=1,linestyle="--")
	axs[1].plot(slnSEA,sltSEA,c=lsc[5],lw=1,linestyle="--")
	axs[1].plot(slnAMZ,sltAMZ,c=lsc[4],lw=1,linestyle="--")
	axs[1].plot(slnTRA,sltTRA,c=lsc[3],lw=1,linestyle="--")
	
	axs[1].text(-144,10,"DTP",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(-100,10,"CRB",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(-80,-12,"AMZ",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(10,15,"TRA",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(130,-12,"SEA",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(-140,-12,"AR6_EPO",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(170,10,"AR6_EPO",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(60,-12,"AR6_EIO",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(-30,-12,"AR6_EAO",verticalalignment="center",backgroundcolor="gray2")
	
	axs[1].format(
		xlim=(-150,210),xlocator=-180:60:180,xlabel=L"Longitude / $\degree$",
		ylim=(-30,30),ylocator=-30:10:30,ylabel=L"Latitude / $\degree$",grid=true
	)
	
	f.savefig(plotsdir("domain.png"),transparent=false,dpi=200)
	load(plotsdir("domain.png"))
end

# ╔═╡ ce969b7e-5248-11eb-36c2-dd1900221e34
md"
Red shades are tropical regions, yellow colours denote dry/arid regions, green denotes temperature regions, blue denotes temperature regions that are dry in winter.
"

# ╔═╡ f8ee207c-5700-11eb-1a7a-0b7f629d74b9
md"
### B. Define Regions as GeoRegions

These regions are defined in the `gregionsadd.txt` in the `src` directory, and have been added as custom GeoRegions that can be called by the `GeoRegions.jl` package.
"

# ╔═╡ cfae2866-a99d-4da9-823a-e069c0fc6148
md"
### C. Southeast Asian Subregions
"

# ╔═╡ 8cf9af9f-653f-4610-8aaf-df7c598d32ca
begin
	geo = GeoRegion("SEA")
	md"Define Southeast Asia GeoRegion ..."
end

# ╔═╡ d573cfb4-6eb5-48e7-9476-e8b0a8da9910
begin
	addGeoRegions(srcdir("addSEArect.txt"))
	addGeoRegions(srcdir("addSEApoly.txt"))
	md"Adding Southeast Asian subregions ..."
end

# ╔═╡ 02d84e44-993b-4f09-be7a-f5ccb3d06705
begin
	ggrd = RegionGrid(geo,lon,lat)
	N,S,E,W = geo.N,geo.S,geo.E,geo.W
	md"Loading RegionGrid information ..."
end

# ╔═╡ ea1495f2-edcc-4606-a6b3-f739bff104a3
begin
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
	kcrtrop = zeros(nlon,nlat)
	kcrarid = zeros(nlon,nlat)
	kcrtemph = zeros(nlon,nlat)
	kcrtemps = zeros(nlon,nlat)
	kcrtempw = zeros(nlon,nlat)
	if typeof(ggrd) <: PolyGrid
		mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end
	for glat in 1 : nlat, glon in 1 : nlon
		kcrtrop[glon,glat]  = kctrop[ilon[glon],ilat[glat]]  * mask[glon,glat]
		kcrarid[glon,glat]  = kcarid[ilon[glon],ilat[glat]]  * mask[glon,glat]
		kcrtemph[glon,glat] = kctemph[ilon[glon],ilat[glat]] * mask[glon,glat]
		kcrtemps[glon,glat] = kctemps[ilon[glon],ilat[glat]] * mask[glon,glat]
		kcrtempw[glon,glat] = kctempw[ilon[glon],ilat[glat]] * mask[glon,glat]
	end
	md"Using Grid to extract Regional Data ..."
end

# ╔═╡ a5f1afd5-4b14-4f27-b414-191f947f73fc
begin
	blnSMT,bltSMT,slnSMT,sltSMT = coordGeoRegion(GeoRegion("SMT"))
	blnBRN,bltBRN,slnBRN,sltBRN = coordGeoRegion(GeoRegion("BRN"))
	slnJAV,sltJAV = coordGeoRegion(GeoRegion("JAV"))
	blnSLW,bltSLW,slnSLW,sltSLW = coordGeoRegion(GeoRegion("SLW"))
	blnPHL,bltPHL,slnPHL,sltPHL = coordGeoRegion(GeoRegion("PHL"))
	blnPNG,bltPNG,slnPNG,sltPNG = coordGeoRegion(GeoRegion("PNG"))
	md"Loading subregion shapes of interest within Southeast Asia ..."
end

# ╔═╡ fa5dbc76-84e4-4889-9e28-2b1053a3d6c9
begin
	asp = (E-W+2)/(N-S+2)
	pplt.close(); freg,areg = pplt.subplots(aspect=asp,axwidth=asp*1.5);
	
	c = areg[1].pcolormesh(ggrd.glon,ggrd.glat,kcrtrop',levels=0:6,cmap="Reds_r")
	areg[1].pcolormesh(ggrd.glon,ggrd.glat,kcrarid',levels=5:9,cmap="Yellow3_r")
	areg[1].pcolormesh(ggrd.glon,ggrd.glat,kcrtemph',levels=6:12,cmap="Green2_r")
	areg[1].pcolormesh(ggrd.glon,ggrd.glat,kcrtemps',levels=11:15,cmap="Brown1_r")
	areg[1].pcolormesh(ggrd.glon,ggrd.glat,kcrtempw',levels=12:18,cmap="Blue3_r")
	areg[1].plot(slnSMT,sltSMT,lw=1,linestyle="--")
	areg[1].plot(slnBRN,sltBRN,lw=1,linestyle="--")
	areg[1].plot(slnJAV,sltJAV,lw=1,linestyle="--")
	areg[1].plot(slnSLW,sltSLW,lw=1,linestyle="--")
	areg[1].plot(slnPHL,sltPHL,lw=1,linestyle="--")
	areg[1].plot(slnPNG,sltPNG,lw=1,linestyle="--")
	areg[1].plot(slnSEA,sltSEA,c="k",lw=1)
	areg[1].plot(x,y,c="k",lw=0.5)
	
	areg[1].format(
		xlim=(W-1,E+1),xlabel=L"Longitude / $\degree$",xlocator=-180:20:360,
		ylim=(S-1,N+1),ylabel=L"Latitude / $\degree$",grid=true
	)
	
	freg.savefig(plotsdir("domain_SEA.png"),transparent=false,dpi=200)
	load(plotsdir("domain_SEA.png"))
end

# ╔═╡ Cell order:
# ╟─f9606c7e-51f4-11eb-2e24-d998e4a91b9a
# ╟─bcfd5bd8-51f5-11eb-1d79-ab69c65febaf
# ╟─c1fdcad0-51f5-11eb-2df9-b1f4c8dca09d
# ╟─e58cf4ca-51f8-11eb-218d-8d3141dc8289
# ╟─3f3a7695-ce0b-444b-a8b2-29421050d828
# ╟─a02d93e9-c53e-4a82-bdf8-34d38377ff99
# ╟─f3c8ffe6-51f5-11eb-362d-293e727b14c7
# ╟─1af88ca2-6a55-481e-a789-705e348d53da
# ╟─a4458db2-5248-11eb-3ecf-b9bcdde2ec37
# ╟─170ae0f0-51f6-11eb-153c-e59511b6a82d
# ╟─ce969b7e-5248-11eb-36c2-dd1900221e34
# ╟─f8ee207c-5700-11eb-1a7a-0b7f629d74b9
# ╟─cfae2866-a99d-4da9-823a-e069c0fc6148
# ╟─8cf9af9f-653f-4610-8aaf-df7c598d32ca
# ╟─d573cfb4-6eb5-48e7-9476-e8b0a8da9910
# ╟─02d84e44-993b-4f09-be7a-f5ccb3d06705
# ╟─ea1495f2-edcc-4606-a6b3-f739bff104a3
# ╟─a5f1afd5-4b14-4f27-b414-191f947f73fc
# ╟─fa5dbc76-84e4-4889-9e28-2b1053a3d6c9
