### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ abdaeba6-56a5-11eb-222b-a12819fce07c
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ ae7e6c26-56a5-11eb-2444-694b368de12c
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using GeoRegions
	using NCDatasets
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ f838a206-56a1-11eb-35dc-f1489317cd8a
md"
# 2a. GPM land-sea mask data

For GPM IMERGv6 precipitation, NASA has provided a land-sea mask denoting the percentage of the grid-box that is considered to be water.  In this notebook, we do a quick exploration of this dataset and extract the data for certain regions.  We also convert the land-sea mask from percentage of water, into the 0-1 range similar to that found in ECMWF reanalysis land-sea mask datasets.
"

# ╔═╡ a907ca40-56dd-11eb-07e6-dd780174001e
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ b28b2cbe-56a7-11eb-3a3c-09b80ad4e0f2
md"
### A. Loading and presenting GPM IMERG land-sea mask

We start by loading and plotting the GPM IMERG land-sea mask ...
"

# ╔═╡ c43c087a-56a7-11eb-2c12-7db44d0a7eae
begin
	ds  = NCDataset(datadir("GPM_IMERG_LandSeaMask.2.nc4"))
	lon = ds["lon"][:]
	lat = ds["lat"][:]
	lsm = ds["landseamask"][:]'*1; lsm[iszero.(lsm)] .= 0.00001
	close(ds)
end

# ╔═╡ ed63a2e2-56a7-11eb-10f2-2b291133bb4e
begin
	pplt.close(); f,axs = pplt.subplots(aspect=6,axwidth=6)

	c = axs[1].contourf(
		lon,lat,lsm',levels=0:10:100,
		cmap="Delta_r",extend="both"
	)
	axs[1].colorbar(c,loc="r")

	for ax in axs
		ax.format(xlim=(-180,180),ylim=(-30,30),xlocator=-180:60:180)
	end

	f.savefig(plotsdir("GPMlsm-TRP.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("GPMlsm-TRP.png"))
end

# ╔═╡ b5841fc6-56ab-11eb-0c0d-cdf45571f6d7
md"
### B. Extracting and Regridding the ETOPO1 data

Here, we define the functions that are used to extract and regrid the land-sea mask, using the `GeoRegions.jl` package.  The land-sea mask, is recalculated from percentage water, into fraction which is land, similar to ECMWF land-sea mask data.
"

# ╔═╡ b4613c36-56f2-11eb-1865-e33bc1e1b75f
thr = 70

# ╔═╡ e82fb766-56f1-11eb-34aa-f35e027df299
begin
	lsm01 = lsm / 100
	lsm01 = 1 .- lsm01
	
md"Binning land-sea mask into 0s and 1s"
end

# ╔═╡ 3443392c-56ad-11eb-3a44-137601990567
md"Using the functions `regiongridvec` and `regionextractgrid`, we are able to extract the land-sea mask for a given region."

# ╔═╡ 0b25e87e-56f0-11eb-1cb9-4574948fd72a
regID = "SEA"

# ╔═╡ 0b9d45f8-56dd-11eb-3928-2b11e9d7bdf7
begin
	rlon,rlat,rinfo = regiongridvec(gregionbounds(regID),lon,lat)
	# rlon,rlat,rinfo = regiongridvec([6,-6,107,95],lon,lat)
	rlsm = regionextractgrid(lsm01,rinfo)
	md"Extracting Land-Sea Mask for the region ..."
end

# ╔═╡ 4536cd18-56dd-11eb-0a9b-a1530a862491
begin
	asp = (maximum(rlon)-minimum(rlon))/(maximum(rlat)-minimum(rlat))
	pplt.close(); freg,areg = pplt.subplots(aspect=asp,axwidth=asp*2)
	
	creg = areg[1].contourf(
		rlon,rlat,rlsm',levels=(0:10)/10,
		cmap="Delta",extend="both"
	)
	areg[1].plot(x,y,c="k",lw=0.5)
	areg[1].colorbar(creg,loc="r")
	
	for ax in areg
		ax.format(
			xlim=(90,165),xlocator=90:15:165,xminorlocator=90:5:165,
			ylim=(-15,20)
			# xlim=(minimum(rlon),maximum(rlon)),
			# ylim=(minimum(rlat),maximum(rlat))
		)
	end
	
	freg.savefig(plotsdir("GPMlsmreg.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("GPMlsmreg.png"))
end

# ╔═╡ 54903304-57ae-11eb-0654-c7c351bcae9f
function savelsm(regID,lon,lat,lsm)
	
	fnc = datadir("GPM_IMERG_LandSeaMask-$regID.nc")
	if isfile(fnc); rm(fnc) end
	ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));
	
	ds.dim["longitude"] = length(lon)
    ds.dim["latitude"]  = length(lat)
	
	nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))
	
	nclsm = defVar(ds,"landseamask",Float32,("longitude","latitude",),attrib = Dict(
        "units" => "0 - 1"
    ))
	
	nclongitude[:] = lon
	nclatitude[:]  = lat
	nclsm[:] = lsm
	
	close(ds)
	
end

# ╔═╡ c6e405de-57ae-11eb-20b5-dbe72ee915fd
if regID == "TRP"
	savelsm(regID,rlon,rlat,rlsm)
end

# ╔═╡ Cell order:
# ╟─f838a206-56a1-11eb-35dc-f1489317cd8a
# ╟─abdaeba6-56a5-11eb-222b-a12819fce07c
# ╟─ae7e6c26-56a5-11eb-2444-694b368de12c
# ╠═a907ca40-56dd-11eb-07e6-dd780174001e
# ╟─b28b2cbe-56a7-11eb-3a3c-09b80ad4e0f2
# ╠═c43c087a-56a7-11eb-2c12-7db44d0a7eae
# ╟─ed63a2e2-56a7-11eb-10f2-2b291133bb4e
# ╟─b5841fc6-56ab-11eb-0c0d-cdf45571f6d7
# ╠═b4613c36-56f2-11eb-1865-e33bc1e1b75f
# ╠═e82fb766-56f1-11eb-34aa-f35e027df299
# ╟─3443392c-56ad-11eb-3a44-137601990567
# ╠═0b25e87e-56f0-11eb-1cb9-4574948fd72a
# ╟─0b9d45f8-56dd-11eb-3928-2b11e9d7bdf7
# ╠═4536cd18-56dd-11eb-0a9b-a1530a862491
# ╟─54903304-57ae-11eb-0654-c7c351bcae9f
# ╟─c6e405de-57ae-11eb-20b5-dbe72ee915fd
