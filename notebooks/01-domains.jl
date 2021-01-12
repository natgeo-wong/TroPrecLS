### A Pluto.jl notebook ###
# v0.12.18

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
* Equatorial Africa (EQA)
* Tropical Rainforest Africa (TRA)
* Amazon Rainforest (AMZ)
* Caribbean Islands (CRB)
"

# ╔═╡ 5afacbc6-524c-11eb-08fb-8b6116cc35de
prect(N::Real,S::Real,W::Real,E::Real) = [W,E,E,W,W],[S,S,N,N,S]

# ╔═╡ d9101792-51f6-11eb-055c-051e9c48c5b3
begin
	DTP = prect(15,-15,-150,210)
	SEA = prect(15,-15,90,165)
	EQA = prect(12,-15,-20,40)
	TRA = prect(5,-5,5,40)
	AMZ = prect(10,-10,-75,-45)
	CRB = prect(25,15,-90,-60)
	
md"Defining bounds of regions of interest ..."
end

# ╔═╡ f3c8ffe6-51f5-11eb-362d-293e727b14c7
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ a4458db2-5248-11eb-3ecf-b9bcdde2ec37
begin
	ds  = NCDataset(datadir("koppen.nc"))
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
	axs[1].plot(DTP[1],DTP[2],c="k",lw=1,linestyle="--")
	axs[1].plot(SEA[1],SEA[2],c="b",lw=1,linestyle="--")
	axs[1].plot(TRA[1],TRA[2],c="r",lw=1,linestyle="--")
	axs[1].plot(EQA[1],EQA[2],c="orange",lw=1,linestyle="--")
	axs[1].plot(AMZ[1],AMZ[2],c="g",lw=1,linestyle="--")
	axs[1].plot(CRB[1],CRB[2],c="blue3",lw=1,linestyle="--")
	
	axs[1].text(-144,5,"DTP",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(-100,10,"CRB",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(-80,-15,"AMZ",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(35,15,"EQA",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(-10,-5,"TRA",verticalalignment="center",backgroundcolor="gray2")
	axs[1].text(150,5,"SEA",verticalalignment="center",backgroundcolor="gray2")
	
	axs[1].format(
		xlim=(-150,210),xlocator=-180:60:180,
		ylim=(-30,30),ylocator=-30:10:30,grid=true
	)
	
	f.savefig(plotsdir("domain.png"),transparent=false,dpi=200)
	load(plotsdir("domain.png"))
end

# ╔═╡ ce969b7e-5248-11eb-36c2-dd1900221e34
md"
Red shades are tropical regions, yellow colours denote dry/arid regions, green denotes temperature regions, blue denotes temperature regions that are dry in winter.
"

# ╔═╡ Cell order:
# ╟─f9606c7e-51f4-11eb-2e24-d998e4a91b9a
# ╟─bcfd5bd8-51f5-11eb-1d79-ab69c65febaf
# ╟─c1fdcad0-51f5-11eb-2df9-b1f4c8dca09d
# ╟─e58cf4ca-51f8-11eb-218d-8d3141dc8289
# ╟─d9101792-51f6-11eb-055c-051e9c48c5b3
# ╟─5afacbc6-524c-11eb-08fb-8b6116cc35de
# ╟─f3c8ffe6-51f5-11eb-362d-293e727b14c7
# ╟─a4458db2-5248-11eb-3ecf-b9bcdde2ec37
# ╟─170ae0f0-51f6-11eb-153c-e59511b6a82d
# ╟─ce969b7e-5248-11eb-36c2-dd1900221e34