### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ a675084e-f638-11eb-2930-b70006b9445c
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 4e08778c-cd46-4a79-9b6a-96ca6447e975
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using GeoRegions
	using ImageFiltering
	using NCDatasets

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ faa41f45-5a58-4bcb-a48c-dcd17e8bad3a
md"
# 01c. Land-Sea Mask Filtering
"

# ╔═╡ 75ea588a-64df-4b31-9dae-f99f425ee55a
md"
### A. Loading Land-Sea Mask Data
"

# ╔═╡ c557eeb4-fab4-4751-af4b-bea2115e85a9
begin
	ds  = NCDataset(datadir("reanalysis","era5-TRPx0.25-lsm-sfc.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	lsm = ds["lsm"][:]
	close(ds)
	md"Loading land-sea mask for the Tropical Region ..."
end

# ╔═╡ dc6929b0-c95d-4c71-9a4c-cb9fa4b01e34
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]; y = coast[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ a9bae1a3-c994-4b36-b9f8-740aa14a3929
begin
	pplt.close(); f1,a1 = pplt.subplots(aspect=1,axwidth=3)
	
	c = a1[1].pcolormesh(lon,lat,lsm',levels=(0:10)/10,extend="both",cmap="Delta")
	# a1[1].plot(x,y,c="k",lw=0.5)
	a1[1].plot([95.5,106,105.5,95,95.5],[6.5,-6,-6.5,6,6.5],c="r",lw=1)
	a1[1].plot([95,105.5,105,94.5,95],[6,-6.5,-7,5.5,6],c="r",lw=1)
	a1[1].plot([94.5,105,104.5,94,94.5],[5.5,-7,-7.5,5,5.5],c="r",lw=1)
	a1[1].plot([94,104.5,104,93.5,94],[5,-7.5,-8,4.5,5],c="r",lw=1)
	a1[1].plot([93.5,104,103.5,93,93.5],[4.5,-8,-8.5,4,4.5],c="r",lw=1)
	a1[1].plot([93,103.5,103,92.5,93],[4,-8.5,-9,3.5,4],c="r",lw=1)
	a1[1].format(xlim=(90,110),ylim=(-10,10))
	a1[1].colorbar(c,loc="r")
	
	f1.savefig("testlsm.png",transparent=false,dpi=200)
	PNGFiles.load("testlsm.png")
end

# ╔═╡ 95ffc474-8416-451c-ba9a-4452dc582626
md"
### B. Smoothing the Land-Sea Mask?
"

# ╔═╡ a66e0d88-7f48-4d7d-b716-3a424d93e90b
function filterlsm(olsm,iterations=1)
	
	it = 0
	nlsm = deepcopy(olsm)
	while it < iterations
		nlsm = log10.(imfilter(10. .^nlsm, Kernel.gaussian(5),"circular"));
		nlsm = (nlsm.+olsm)/2
		it += 1
	end
	
	nlsm[nlsm.<0] .= minimum(nlsm[nlsm.>0])
	
	return nlsm
	
end

# ╔═╡ d3864cc7-284b-4278-9499-0b6f4d274608
begin
	nlsm = filterlsm(lsm,3)
	md"Performing gaussian filtering/smoothing on land-sea mask ..."
end

# ╔═╡ 05c13f7e-58f8-4136-830d-fb20a9dc0e6f
begin
	geo = GeoRegion("AMZ")
	md"Defining region coordinates ..."
end

# ╔═╡ 8471606c-fead-450a-80f0-fd4cf308c3e5
begin
	N,S,E,W = geo.N,geo.S,geo.E,geo.W
	ggrd = RegionGrid(geo,lon,lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
	rlsm = zeros(nlon,nlat)
	flsm = zeros(nlon,nlat)
	if typeof(ggrd) <: PolyGrid
		mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end
	for glat in 1 : nlat, glon in 1 : nlon
		rlsm[glon,glat] =  lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
		flsm[glon,glat] = nlsm[ilon[glon],ilat[glat]] * mask[glon,glat]
	end
	md"Extracting information for region ..."
end

# ╔═╡ 51d55158-32df-4135-ab12-186e295ce499
begin
	asp = (E-W+2)/(N-S+2)
	if asp > 1.5
		freg,areg = pplt.subplots(nrows=2,axwidth=asp*1.2,aspect=asp)
	else
		freg,areg = pplt.subplots(ncols=2,axwidth=2.5,aspect=asp)
	end
	
	lvls = vcat(10. .^(-6:-1),0.2,0.5,0.8,0.9)

	creg = areg[1].pcolormesh(
		ggrd.glon,ggrd.glat,rlsm',
		levels=lvls,cmap="delta",extend="both"
	)
	areg[1].plot(x,y,c="k",lw=0.5)
	areg[1].format(urtitle="Raw")
	
	areg[2].pcolormesh(
		ggrd.glon,ggrd.glat,flsm',
		levels=lvls,cmap="delta",extend="both"
	)
	areg[2].plot(x,y,c="k",lw=0.5)
	areg[2].format(urtitle="Filtered")
	
	freg.colorbar(creg,loc="r")

	for ax in areg
		ax.format(
			xlim=(ggrd.glon[1].-1,ggrd.glon[end].+1),
			xlabel=L"Longitude / $\degree$",
			ylim=(S-1,N+1),ylabel=L"Latitude / $\degree$",
			suptitle="Land-Sea Mask",
			grid=true
		)
	end

	freg.savefig(plotsdir("01c-flsm_$(geo.regID).png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("01c-flsm_$(geo.regID).png"))
end

# ╔═╡ Cell order:
# ╟─faa41f45-5a58-4bcb-a48c-dcd17e8bad3a
# ╟─a675084e-f638-11eb-2930-b70006b9445c
# ╟─4e08778c-cd46-4a79-9b6a-96ca6447e975
# ╟─75ea588a-64df-4b31-9dae-f99f425ee55a
# ╟─c557eeb4-fab4-4751-af4b-bea2115e85a9
# ╟─dc6929b0-c95d-4c71-9a4c-cb9fa4b01e34
# ╟─a9bae1a3-c994-4b36-b9f8-740aa14a3929
# ╟─95ffc474-8416-451c-ba9a-4452dc582626
# ╠═a66e0d88-7f48-4d7d-b716-3a424d93e90b
# ╟─d3864cc7-284b-4278-9499-0b6f4d274608
# ╟─05c13f7e-58f8-4136-830d-fb20a9dc0e6f
# ╟─8471606c-fead-450a-80f0-fd4cf308c3e5
# ╟─51d55158-32df-4135-ab12-186e295ce499