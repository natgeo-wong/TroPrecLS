### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 36a4e004-5a26-11eb-1ba4-db83d18c7845
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 3a8bbab2-5a26-11eb-20e0-d551c2a68858
begin
	@quickactivate "TroPrecLS"
	using NCDatasets
	using PlutoUI
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("common.jl"))
	include(srcdir("samlsf.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 582bcb26-5a22-11eb-3b52-09e3ec08b3f9
md"
# 5b. Moisture-Flux Profiles

In this notebook, we calculate and find the vertical profile of mean moisture-fluxes into our predefined domains.  These will then be used to force the mean-state of our SAM models into wetter/drier states, such that we can reconstruct the P-r curve in Bretherton et al. [2006].

We use mean values of ERA5 specific humidity, zonal and meridional wind to calculate the moisture flux vertical profile at each level, since this is not given as raw data in the ERA5 CDS data store.  By using only the mean-values calculated, and not the raw hourly output, we assume the following

$$\overline{qv}
= \overline{q}\overline{v} + \overline{q'}\overline{v'}
\approx \overline{q}\overline{v}$$

i.e., that $\overline{q'}\overline{v'} \ll \overline{q}\overline{v}$.
"

# ╔═╡ abde634c-5a83-11eb-00d6-69b233e0c00f
md"
### A. Extracting Moisture Fluxes for Tropical Domains
"

# ╔═╡ e3c519fc-5a35-11eb-0d34-cd4ca8b50bf9
begin
	qfcds = NCDataset(datadir("reanalysis/era5-TRPx0.25-qfc_air.nc"))
	lon = qfcds["longitude"][:]
	lat = qfcds["latitude"][:]
	lvl = qfcds["level"][:]
	qfc = qfcds["qfc_air"][:] * 1
	close(qfcds)
md"Retrieving calculated MFC for the Tropical Region ..."
end

# ╔═╡ e0251f16-5a91-11eb-1bc8-a70691df76e2
begin
	lds = NCDataset(datadir("reanalysis/era5-TRPx0.25-lsm-sfc.nc"))
	lsm = lds["lsm"][:]*1
	close(lds)
	
md"Loading Land-Sea Mask for ERA5 data ..."
end

# ╔═╡ 05f8a0be-5a92-11eb-0187-e19a897d9bf9
begin
	lqfc_SEA,sqfc_SEA = getmean([10,-10,135,90],qfc,lon,lat,length(lvl),lsm)
	lqfc_TRA,sqfc_TRA = getmean([10,-10,40,-10],qfc,lon,lat,length(lvl),lsm)
	lqfc_CRB,sqfc_CRB = getmean([25,15,-60,-90],qfc,lon,lat,length(lvl),lsm)
	lqfc_AMZ,sqfc_AMZ = getmean([10,-10,-45,-75],qfc,lon,lat,length(lvl),lsm)
	lqfc_DTP,sqfc_DTP = getmean([10,-10,360,0],qfc,lon,lat,length(lvl),lsm)
end

# ╔═╡ 46671d58-5a92-11eb-2a04-975238148a1b
begin
	pplt.close(); fqfc,aqfc = pplt.subplots(ncols=2,aspect=0.5,axwidth=1.5);
	
	lgd = Dict("frame"=>false,"ncols"=>1)
	
	aqfc[1].plot(lqfc_SEA*1e8,lvl,c="b")
	aqfc[1].plot(lqfc_TRA*1e8,lvl,c="r")
	aqfc[1].plot(lqfc_CRB*1e8,lvl,c="blue3")
	aqfc[1].plot(lqfc_AMZ*1e8,lvl,c="g")
	
	aqfc[2].plot(sqfc_DTP*1e8,lvl,c="k",label="DTP",legend="r")
	aqfc[2].plot(sqfc_SEA*1e8,lvl,c="b",label="SEA",legend="r",legend_kw=lgd)
	aqfc[2].plot(sqfc_TRA*1e8,lvl,c="r",label="TRA",legend="r")
	aqfc[2].plot(sqfc_CRB*1e8,lvl,c="blue3",label="CRB",legend="r")
	aqfc[2].plot(lqfc_AMZ*1e8*NaN,lvl,c="g",label="AMZ",legend="r")
	
	aqfc[1].format(
		xlim=(-12.5,12.5),ylim=(1000,50),urtitle="Land",
		ylabel="Pressure / hPa",
		xlabel=L"$-\nabla\cdot(q\vec{u})$ / 10$^{-8}$ kg kg$^{-1}$"
	)
	aqfc[2].format(xlim=(-4,4),ylim=(1000,50),urtitle="Ocean")
	
	fqfc.savefig(plotsdir("qfc.png"),transparent=false,dpi=200)
	load(plotsdir("qfc.png"))
end

# ╔═╡ 8774c3bc-5aa0-11eb-2dd9-d1f3daf9d815
md"We then proceed to replot this by height"

# ╔═╡ 06eebd8c-5aa1-11eb-12d9-7314d15a98c9
begin
	zds = NCDataset(datadir("reanalysis/era5-TRPx0.25-z_air.nc"))
	zlon = zds["longitude"][:]
	zlat = zds["latitude"][:]
	zlvl = zds["level"][8:end]
	zair = zds["z_air"][:,:,8:end] /1e3/9.81
	close(zds)
end

# ╔═╡ c4795296-5aa0-11eb-2038-1be17d8ce5cd
begin
	lzair_SEA,szair_SEA = getmean([10,-10,135,90],zair,zlon,zlat,length(zlvl),lsm)
	lzair_TRA,szair_TRA = getmean([10,-10,40,-10],zair,zlon,zlat,length(zlvl),lsm)
	lzair_CRB,szair_CRB = getmean([25,15,-60,-90],zair,zlon,zlat,length(zlvl),lsm)
	lzair_AMZ,szair_AMZ = getmean([10,-10,-45,-75],zair,zlon,zlat,length(zlvl),lsm)
	lzair_DTP,szair_DTP = getmean([10,-10,360,0],zair,zlon,zlat,length(zlvl),lsm)
end

# ╔═╡ f3bd1c4a-5aa0-11eb-049c-8d6783190671
begin
	pplt.close(); fzair,azair = pplt.subplots(ncols=2,aspect=0.5,axwidth=1.5);
	
	azair[1].plot(lqfc_SEA*1e8,lzair_SEA,c="b")
	azair[1].plot(lqfc_TRA*1e8,lzair_TRA,c="r")
	azair[1].plot(lqfc_CRB*1e8,lzair_CRB,c="blue3")
	azair[1].plot(lqfc_AMZ*1e8,lzair_AMZ,c="g")
	
	azair[2].plot(sqfc_DTP*1e8,szair_DTP,c="k",label="DTP",legend="r")
	azair[2].plot(sqfc_SEA*1e8,szair_SEA,c="b",label="SEA",legend="r",legend_kw=lgd)
	azair[2].plot(sqfc_TRA*1e8,szair_TRA,c="r",label="TRA",legend="r")
	azair[2].plot(sqfc_CRB*1e8,szair_CRB,c="blue3",label="CRB",legend="r")
	azair[2].plot(sqfc_AMZ*1e8*NaN,szair_AMZ/1e3,c="g",label="AMZ",legend="r")
	
	azair[1].format(
		xlim=(-12.5,12.5),ylim=(0.1,20),urtitle="Land",
		ylabel="Height / km",
		xlabel=L"$-\nabla\cdot(q\vec{u})$ / 10$^{-8}$ kg kg$^{-1}$"
	)
	azair[2].format(xlim=(-4,4),urtitle="Ocean",)#yscale="log")
	
	fzair.savefig(plotsdir("zair.png"),transparent=false,dpi=200)
	load(plotsdir("zair.png"))
end

# ╔═╡ 13d785a6-5d03-11eb-290b-0dec1275019e
md"Create LSF files? $(@bind dolsf PlutoUI.Slider(0:1))"

# ╔═╡ eae6df2a-5d02-11eb-3d49-1b2b62d806d2
if isone(dolsf)
	
	lsfl = lsfinit(length(zlvl)); lsfl[:,2] .= -999.0; lsfl[:,1] .= reverse(lzair_SEA)
	lsfs = lsfinit(length(zlvl)); lsfs[:,2] .= -999.0; lsfs[:,1] .= reverse(szair_SEA)
	mvec = vcat(-20:-1,1:20)
	
	if !isdir(projectdir("exp/lsf/land/")); mkpath(projectdir("exp/lsf/land/")) end
	if !isdir(projectdir("exp/lsf/ocean/")); mkpath(projectdir("exp/lsf/ocean/")) end
	
	for mul in mvec
		
		lsfl[:,4] .= reverse(lqfc_SEA) * mul
		lsfs[:,4] .= reverse(sqfc_SEA) * mul
		
		if mul > 0
			  mstr = @sprintf("p%02d",abs(mul))
		else; mstr = @sprintf("n%02d",abs(mul))
		end

		lsfprint(projectdir("exp/lsf/land/$(mstr)"),lsfl,1009.32)
		lsfprint(projectdir("exp/lsf/ocean/$(mstr)"),lsfs,1009.32)
		
	end
	
md"Based on these profiles, we create large-scale forcing profiles for moisture flux convergence to be used in our WTG simulations ..."
else
md"We have decided not to override any preexisting large-scale forcing profiles for moisture flux convergence."
end

# ╔═╡ Cell order:
# ╟─582bcb26-5a22-11eb-3b52-09e3ec08b3f9
# ╟─36a4e004-5a26-11eb-1ba4-db83d18c7845
# ╠═3a8bbab2-5a26-11eb-20e0-d551c2a68858
# ╟─abde634c-5a83-11eb-00d6-69b233e0c00f
# ╟─e3c519fc-5a35-11eb-0d34-cd4ca8b50bf9
# ╟─e0251f16-5a91-11eb-1bc8-a70691df76e2
# ╠═05f8a0be-5a92-11eb-0187-e19a897d9bf9
# ╟─46671d58-5a92-11eb-2a04-975238148a1b
# ╟─8774c3bc-5aa0-11eb-2dd9-d1f3daf9d815
# ╠═06eebd8c-5aa1-11eb-12d9-7314d15a98c9
# ╠═c4795296-5aa0-11eb-2038-1be17d8ce5cd
# ╟─f3bd1c4a-5aa0-11eb-049c-8d6783190671
# ╠═13d785a6-5d03-11eb-290b-0dec1275019e
# ╟─eae6df2a-5d02-11eb-3d49-1b2b62d806d2
