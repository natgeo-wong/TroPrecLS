### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ f23d2756-d792-11eb-23dd-ab8bf43c5720
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ d09fb352-5b34-47d8-976c-5660cae9a3be
begin
	@quickactivate "TroPrecLS"
	using GeoRegions
	using NCDatasets
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 54079d86-7cd0-405d-b82d-dcaf6f90274b
md"
# 09c. Skin Temperature vs Column Relative Humidity
"

# ╔═╡ 63a4fc36-b7ca-4880-8bd3-27a09f53a7e2
md"
### A. Loading ERA5 Data
"

# ╔═╡ d3b32933-1ca7-4a34-93d9-94eb65761b8b
begin
	ds  = NCDataset(datadir("reanalysis/era5-TRPx0.25-sstvcsf.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	skt = ds["skt"][:]; pskt = (skt[2:end] + skt[1:(end-1)]) /2
	csf = ds["csf"][:]
	frq = ds["bin_frq"][:]
	close(ds)
	md"Loading column relative humidity and skin temperature bin data ..."
end

# ╔═╡ 017e64a2-3a8c-49f9-9ca6-afca658bcb0b
begin
	dsl = NCDataset(datadir("reanalysis","era5-TRPx0.25-lsm-sfc.nc"))
	lsm = dsl["lsm"][:] * 1
	close(dsl)
	md"Loading land-sea mask data ..."
end

# ╔═╡ 940316b0-ebe1-4bdc-aec0-a1edae7bd71d
md"
### B. Extraction of Regional Data
"

# ╔═╡ 19e09271-e54a-48f2-afac-461cb801cfd1
function extractocnlnd(geo,frq,lon,lat,lsm)
	
	ggrd = RegionGrid(geo,lon,lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
	nskt = size(frq,3)
	ncsf = size(frq,4)
    rfrq = zeros(Int32,nlon,nlat)
    rlsm = zeros(nlon,nlat)
	ofrq = zeros(Int32,nskt,ncsf)
	lfrq = zeros(Int32,nskt,ncsf)

	if typeof(ggrd) <: PolyGrid
		  mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end

	for icsf = 1 : ncsf, iskt = 1 : nskt
		for glat in 1 : nlat, glon in 1 : nlon
			rfrq[glon,glat] = frq[ilon[glon],ilat[glat],iskt,icsf]
			rlsm[glon,glat] = lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
		end
		ofrq[iskt,icsf] = sum(rfrq[rlsm.<0.5])
		lfrq[iskt,icsf] = sum(rfrq[rlsm.>0.5])
	end
	
	return ofrq,lfrq
	
end

# ╔═╡ 6acdcde3-2190-4cec-8a4f-9bf096602b52
begin
	geo = "SEA"
	md"Defining GeoRegion ..."
end

# ╔═╡ d2c09607-2d8c-4744-bc9c-ff5366dd906d
begin
	ofrq,lfrq = extractocnlnd(GeoRegion(geo),frq,lon,lat,lsm)
	tofrq = sum(ofrq,dims=2); nofrq = ofrq ./ tofrq * 100
	tlfrq = sum(lfrq,dims=2); nlfrq = lfrq ./ tlfrq * 100
	md"Extracting and binning regional data ..."
end

# ╔═╡ b2acf38d-01ec-4b9c-a8da-ce28a2c4656b
begin
	pplt.close(); f,a = pplt.subplots(nrows=2,aspect=3,axwidth=4)
	
	lvls = [0.1,0.141,0.2,0.316,0.5,0.707,1,1.41,2,3.16,5,7.07,10]
	
	c = a[1].pcolormesh(csf,skt,nofrq,levels=lvls,extend="both",cmap=geo)
	a[1].format(grid=true,ylocator=285:5:310,xlocator=0:10:100)
	
	p1 = a[1].panel("l",space="1em",width="4em")
	p1.plot(dropdims(tofrq,dims=2)/sum(tofrq)*25,pskt)
	p1.format(xscale="log")
	
	a[2].pcolormesh(csf,skt,nlfrq,levels=lvls,extend="both",cmap=geo)
	a[2].format(grid=true,ylocator=285:5:310,xlocator=0:10:100)
	f.colorbar(c,loc="r",label=L"Conditional Density f$_{r|T_{sk}}$")
	
	p2 = a[2].panel("l",space="1em",width="4em")
	p2.plot(dropdims(tlfrq,dims=2)/sum(tlfrq)*25,pskt)
	p2.format(xscale="log",xlim=(100,0.001),xlocator=[0.01,1,100])
	
	for ax in a
		ax.format(
			ylabel=L"Skin Temperature (T$_{sk}$) / K",
			xlabel="Column Relative Humidity (r) / %",
			suptitle="$(GeoRegion(geo).name)"
		)
	end
	
	fimg = "sktvcsf_$(GeoRegion(geo).regID).png"
	f.savefig(plotsdir(fimg),dpi=200,transparent=false)
	PNGFiles.load(plotsdir(fimg))
end

# ╔═╡ Cell order:
# ╟─54079d86-7cd0-405d-b82d-dcaf6f90274b
# ╟─f23d2756-d792-11eb-23dd-ab8bf43c5720
# ╟─d09fb352-5b34-47d8-976c-5660cae9a3be
# ╟─63a4fc36-b7ca-4880-8bd3-27a09f53a7e2
# ╟─d3b32933-1ca7-4a34-93d9-94eb65761b8b
# ╟─017e64a2-3a8c-49f9-9ca6-afca658bcb0b
# ╟─940316b0-ebe1-4bdc-aec0-a1edae7bd71d
# ╠═19e09271-e54a-48f2-afac-461cb801cfd1
# ╠═6acdcde3-2190-4cec-8a4f-9bf096602b52
# ╟─d2c09607-2d8c-4744-bc9c-ff5366dd906d
# ╟─b2acf38d-01ec-4b9c-a8da-ce28a2c4656b
