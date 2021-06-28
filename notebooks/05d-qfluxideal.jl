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

# ╔═╡ 6140aa6f-73f3-4e5d-baf8-6f6adfeaa162
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 8a7dab37-3bc4-454c-a97b-09d000fd1d46
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using Dierckx
	using GeoRegions
	using NCDatasets
	using PlutoUI
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("samlsf.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 33beb74c-d6f5-11eb-12f8-679ce305002d
md"
# 05d. Idealized QFlux Profile

Here, we create idealized qflux to be input as large-scale forcing in our SAM runs, based on previous profiles of qflux that we have constructed.
"

# ╔═╡ 28e54fe4-5777-439c-8e5a-63598d5039d1
md"
### A. Reanalysis QFlux Profile

We have previously calculated the qflux profiles for different regions in the Deep Tropical Region:
"

# ╔═╡ c00e1e75-376c-4350-8d62-7e146c35783c
PNGFiles.load(plotsdir("qfluxprofile.png"))

# ╔═╡ 387c8f5a-8b87-47f9-b179-d9ee14b0a8ae
md"And we see that over the ocean basins in the Deep Tropical (DTP) region, the qflux seems to roughly follow a tangent profile"

# ╔═╡ 90c7a533-f20c-410c-9ab8-8ccff1938b91
begin
	pre = [
		1,2,3,5,10,20,30,50,70,
		100,125,150,175,200,225,250,300,350,
		400,450,500,550,600,650,700,750,775,
		800,825,850,875,900,925,950,975,1000
	]*1.0; np = length(pre)
	pre2 = 600:10:1000; np2 = length(pre2)
	md"ERA5 reanalysis pressure levels ..."
end

# ╔═╡ c3e65359-14ff-417e-83dc-12c2ac09c1dd
begin
	qflux = zeros(np)
	qflux2 = zeros(np2)
	expnt = 2
	for ip = 1 : np
		if pre[ip] > 800
			if pre[ip] < 900
				  qflux[ip] = ((pre[ip]-800)/100).^expnt * 0.5
			else; qflux[ip] = 1 - ((1000-pre[ip])/100).^expnt * 0.5
			end
		end
	end
	for ip = 1 : np2
		if pre2[ip] > 800
			if pre2[ip] < 900
				  qflux2[ip] = ((pre2[ip]-800)/100).^expnt * 0.5
			else; qflux2[ip] = 1 - ((1000-pre2[ip])/100).^expnt * 0.5
			end
		end
	end
	md"Constructing idealized QFlux profile ..."
end

# ╔═╡ be8839ab-7431-41dd-bdbb-7967ab43bba1
md"
### B. Orographic Height
"

# ╔═╡ 896cf5d2-cdc6-4f71-83ad-b25dd1051483
begin
	dsz  = NCDataset(datadir("reanalysis/era5-TRPx0.25-z_air.nc"))
	lon  = dsz["longitude"][:]
	lat  = dsz["latitude"][:]
	zair = dsz["z_air"][:] / 9.81
	close(dsz)
	md"Loading orographic height data ..."
end

# ╔═╡ 86136b90-eb26-4db7-8e5b-35fab173966d
function extract3d(var,ggrd)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
	nlvl = size(var,3)
	rvar = zeros(nlon,nlat,nlvl)
	if typeof(ggrd) <: PolyGrid
		mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end
	for ilvl = 1 : nlvl, glat in 1 : nlat, glon in 1 : nlon
		rvar[glon,glat,ilvl] = var[ilon[glon],ilat[glat],ilvl] * mask[glon,glat]
	end
	return rvar
end

# ╔═╡ 8d8a0b53-9c07-41bf-a7b0-8bd98f46d664
begin
	ggrd_DTP = RegionGrid(GeoRegion("DTP"),lon,lat);
	zair_DTP = extract3d(zair,ggrd_DTP)
	mza_DTP  = dropdims(mean(zair_DTP,dims=(1,2)),dims=(1,2))
	md"Extracting regional 3D z_air data"
end

# ╔═╡ 482c7636-75f6-4a8d-af94-549f82c67d3f
begin
	iza_DTP = evaluate(Spline1D(pre,mza_DTP),pre2)
	md"Interpolating to finer grids ..."
end

# ╔═╡ 655d61b0-089d-4e33-a87b-69b4698c34e9
md"
### C. Idealized vs Reanalysis
"

# ╔═╡ e28f10d6-31d8-44c0-8ba7-c40a1fddd3ea
begin
	pplt.close(); f,a = pplt.subplots(ncols=2,aspect=0.5,axwidth=1.5,sharex=0)
	
	a[1].plot(qflux,pre)
	a[1].plot(qflux2,pre2)
	a[1].format(
		ylim=(1000,50),xlim=(-1.5,1.5),ylabel="Pressure / hPa",
		xlabel=L"$\nabla\cdot(q\vec{u})$ / 10$^{-8}$ kg kg$^{-1}$"
	)
	
	a[2].plot(mza_DTP/1000,pre)
	a[2].plot(iza_DTP/1000,pre2)
	a[2].format(xlim=(0,20))
	
	f.savefig(plotsdir("qflux-ideal.png"),dpi=150,transparent=false)
	PNGFiles.load(plotsdir("qflux-ideal.png"))
end

# ╔═╡ ab3edf45-01bd-4b38-9e4b-09bc78fe197a
md"Create LSF files? $(@bind dolsf PlutoUI.Slider(0:1))"

# ╔═╡ 1ebb803d-fc84-4ced-bff6-dd34d538ec88
if isone(dolsf)
	
	lsf = lsfinit(np2); lsf[:,2] .= -999.0; lsf[:,1] .= reverse(iza_DTP)
	mvec = -5 : 0.5 : 5
	
	if !isdir(projectdir("exp/lsf/ideal/")); mkpath(projectdir("exp/lsf/ideal/")) end
	
	for mul in mvec
		
		lsf[:,4] .= reverse(qflux2) * mul * 1e-8
		
		if mul > 0
			mstr = @sprintf("p%03.1f",abs(mul))
		elseif mul < 0
			mstr = @sprintf("n%03.1f",abs(mul))
		end
		
		if !iszero(mul)
			mstr = replace(mstr,"."=>"d")
			lsfprint(projectdir("exp/lsf/ideal/$(mstr)"),lsf,1009.32)
		end
		
	end
	
md"Based on these profiles, we create large-scale forcing profiles for moisture flux convergence to be used in our WTG simulations ..."
else
md"We have decided not to override any preexisting large-scale forcing profiles for moisture flux convergence."
end

# ╔═╡ Cell order:
# ╟─33beb74c-d6f5-11eb-12f8-679ce305002d
# ╟─6140aa6f-73f3-4e5d-baf8-6f6adfeaa162
# ╟─8a7dab37-3bc4-454c-a97b-09d000fd1d46
# ╟─28e54fe4-5777-439c-8e5a-63598d5039d1
# ╟─c00e1e75-376c-4350-8d62-7e146c35783c
# ╟─387c8f5a-8b87-47f9-b179-d9ee14b0a8ae
# ╟─90c7a533-f20c-410c-9ab8-8ccff1938b91
# ╟─c3e65359-14ff-417e-83dc-12c2ac09c1dd
# ╟─be8839ab-7431-41dd-bdbb-7967ab43bba1
# ╟─896cf5d2-cdc6-4f71-83ad-b25dd1051483
# ╟─86136b90-eb26-4db7-8e5b-35fab173966d
# ╟─8d8a0b53-9c07-41bf-a7b0-8bd98f46d664
# ╟─482c7636-75f6-4a8d-af94-549f82c67d3f
# ╟─655d61b0-089d-4e33-a87b-69b4698c34e9
# ╠═e28f10d6-31d8-44c0-8ba7-c40a1fddd3ea
# ╟─ab3edf45-01bd-4b38-9e4b-09bc78fe197a
# ╠═1ebb803d-fc84-4ced-bff6-dd34d538ec88
