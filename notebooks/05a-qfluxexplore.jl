### A Pluto.jl notebook ###
# v0.12.20

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
	using GeoRegions
	using DelimitedFiles
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
# 5a. Reanalysis Moisture-Fluxes

In this notebook, we calculate and find the vertical profile of mean moisture-fluxes into our predefined domains.  These will then be used to force the mean-state of our SAM models into wetter/drier states, such that we can reconstruct the P-r curve in Bretherton et al. [2006].

We use mean values of ERA5 specific humidity, zonal and meridional wind to calculate the moisture flux vertical profile at each level, since this is not given as raw data in the ERA5 CDS data store.  By using only the mean-values calculated, and not the raw hourly output, we assume the following

$$\overline{qv}
= \overline{q}\overline{v} + \overline{q'}\overline{v'}
\approx \overline{q}\overline{v}$$

i.e., that $\overline{q'}\overline{v'} \ll \overline{q}\overline{v}$.
"

# ╔═╡ 2b736522-5a33-11eb-0eea-473a5b8566f3
md"
### A. Calculating Moisture Flux Convergence

We calculate the gradients and divergence using finite difference.  Then, we calculate the moisture flux convergence via the following:

$$\text{MFC} = - \vec{u}\cdot\nabla q - q\nabla\cdot\vec{u}$$
"

# ╔═╡ 28ffa404-5a33-11eb-3dc7-7d716d6cf06f
function ddx(
	data::AbstractArray{<:Real,3},
	lon::AbstractVector{<:Real},
	lat::AbstractVector{<:Real}
)
	
	nx = length(lon); ny = length(lat); nlvl = size(data,3)
	rx = 6378e3; dx = rx * (lon[2]-lon[1]) * pi / 180 * cosd.(lat)
	ddatadx = zeros(nx,ny,nlvl)
	
	for ilvl = 1 : nlvl, iy = 1 : ny, ix = 2 : (nx-1)

        ddatadx[ix,iy,ilvl] = (data[ix+1,iy,ilvl] - data[ix-1,iy,ilvl]) / (2 * dx[iy])

    end

    for iy = 1 : ny, ilvl = 1 : nlvl

        ddatadx[1,iy,ilvl] = (data[2,iy,ilvl] - data[1,iy,ilvl]) / dx[iy]
        ddatadx[nx,iy,ilvl] = (data[nx,iy,ilvl] - data[nx-1,iy,ilvl]) / dx[iy]

    end
	
	return ddatadx
	
end

# ╔═╡ 9e01c7ea-5a34-11eb-1c41-8b991a399584
function ddy(
	data::AbstractArray{<:Real,3},
	lat::AbstractVector{<:Real}
)
	
	nx = size(data,1); ny = length(lat); nlvl = size(data,3)
	ry = 6357e3; dy = ry * (lat[2]-lat[1]) * pi / 180
	ddatady = zeros(nx,ny,nlvl)
	
	for ilvl = 1 : nlvl, iy = 2 : (ny-1), ix = 1 : nx

        ddatady[ix,iy,ilvl] = (data[ix,iy+1,ilvl] - data[ix,iy-1,ilvl]) / (2 * dy)

    end

    for ilvl = 1 : nlvl, ix = 1 : nx

        ddatady[ix,1,ilvl] = (data[ix,1,ilvl] - data[ix,2,ilvl]) / dy
        ddatady[ix,ny,ilvl] = (data[ix,ny-1,ilvl] - data[ix,ny,ilvl]) / dy

    end
	
	return ddatady
	
end

# ╔═╡ abde634c-5a83-11eb-00d6-69b233e0c00f
md"
### B. Extracting Moisture Fluxes for Tropical Region
"

# ╔═╡ 4927b4c2-5a8a-11eb-37e4-53efb9445b23
function saveqfc(qfc,rlon,rlat,lvl,regID)
	
	fnc = datadir("reanalysis/era5-$(regID)x0.25-qfc_air.nc")
	if isfile(fnc); rm(fnc) end
	ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));
	
	ds.dim["longitude"] = length(rlon)
    ds.dim["latitude"]  = length(rlat)
    ds.dim["level"]  = length(lvl)
	
	scale,offset = ncoffsetscale(qfc)
	
	nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))
	
	nclevel = defVar(ds,"level",Int32,("level",),attrib = Dict(
        "units"     => "millibars",
        "long_name" => "pressure_level",
    ))
	
	ncqfc = defVar(ds,"qfc_air",Int16,("longitude","latitude","level"),attrib = Dict(
        "scale_factor"  => scale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
        "units" 	    => "kg kg**-1 s**-1",
		"long_name"     => "moisture_flux_convergence",
    ))
	
	nclongitude[:] = rlon
	nclatitude[:]  = rlat
	nclevel[:]     = lvl
	ncqfc[:]       = qfc
	
	close(ds)
	
end

# ╔═╡ 87a146d4-5a29-11eb-1571-b9d9266ee54b
function calcqfc(regID::AbstractString)
	
	ds = NCDataset(datadir("reanalysis/era5-$(regID)x0.25-q_air.nc"))
    lon = ds["longitude"][:]; nlon = length(lon)
    lat = ds["latitude"][:];  nlat = length(lat)
	lvl = ds["level"][:];	  nlvl = length(lvl)
    q_air = ds["q_air"][:]*1
    close(ds)
	
	ds = NCDataset(datadir("reanalysis/era5-$(regID)x0.25-u_air.nc"))
    u_air = ds["u_air"][:]*1
    close(ds)
	
	ds = NCDataset(datadir("reanalysis/era5-$(regID)x0.25-v_air.nc"))
    v_air = ds["v_air"][:]*1
    close(ds)
	
	tmp = zeros(nlon,nlat,nlvl)
	rlon,rlat,rinfo = gregiongridvec(regID,lon,lat)
	
	ru = regionextractgrid!(u_air,rinfo,tmp)
	rv = regionextractgrid!(v_air,rinfo,tmp)
	rq = regionextractgrid!(q_air,rinfo,tmp)
	
	ru_x = ddx(ru,rlon,rlat); rv_y = ddy(rv,rlat)
	rq_x = ddx(rq,rlon,rlat); rq_y = ddy(rq,rlat)
	
	qfc = -((ru_x + rv_y) .* rq .+ ru .* rq_x .+ rv .* rq_y)
	
	saveqfc(qfc,rlon,rlat,lvl,regID)
	
	return rlon,rlat,lvl,qfc
	
end

# ╔═╡ 795a822a-5a29-11eb-37bb-3b6d62cec6be
regID = "TRP"

# ╔═╡ e3c519fc-5a35-11eb-0d34-cd4ca8b50bf9
if isfile(datadir("reanalysis/era5-$(regID)x0.25-qfc_air.nc"))
	qfcds = NCDataset(datadir("reanalysis/era5-$(regID)x0.25-qfc_air.nc"))
	lon = qfcds["longitude"][:]
	lat = qfcds["latitude"][:]
	lvl = qfcds["level"][:]; nlvl = length(lvl)
	qfc = qfcds["qfc_air"][:] * 1
	close(qfcds)
else lon,lat,lvl,qfc = calcqfc(regID)
end

# ╔═╡ e0251f16-5a91-11eb-1bc8-a70691df76e2
begin
	lds = NCDataset(datadir("reanalysis/era5-TRPx0.25-lsm-sfc.nc"))
	lsm = lds["lsm"][:]*1
	close(lds)
	
md"Loading Land-Sea Mask for ERA5 data ..."
end

# ╔═╡ 8fc1c5e2-5d07-11eb-1b2b-9f48919a115f
md"Pressure Level: $(@bind ip PlutoUI.Slider(10:nlvl))"

# ╔═╡ 9d32ee54-5d07-11eb-0561-a16bf1f4bbec
md"
### C. Exploratory analysis for low-level Moisture Flux

We explore low-level moisture flux at the $(lvl[ip]) hPa level.  We also explore the moisture flux over land where the average surface pressure is higher than $(lvl[ip]) hPa.
"

# ╔═╡ 3f21c35c-5d08-11eb-18c9-d90f1cb9ba2d
begin
	spds = NCDataset(datadir("reanalysis/era5-$(regID)x0.25-sp-sfc.nc"))
	p_sfc = spds["sp"][:,:,argmin(abs.(lvl.-lvl[ip]))] / 100
	qfc_p = qfc[:,:,argmin(abs.(lvl.-lvl[ip]))]
	close(spds)
md"Loading averaged surface pressures ..."
end

# ╔═╡ 368e8334-5d0f-11eb-3c66-c1ff6923e90c
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ 57df47c0-5d24-11eb-16aa-834dd327a306
begin
	DTP = prect(10,-10,0,360)
	SEA = prect(20,-15,90,165)
	TRA = prect(10,-10,0,40)
	TR2 = prect(10,-10,345,360)
	AMZ = prect(10,-10,285,315)
	CRB = prect(25,15,270,300)
	
md"Defining bounds of regions of interest ..."
end

# ╔═╡ bac9bf96-5d0d-11eb-21f7-e3a1dda7af36
begin
	qplt = deepcopy(qfc_p)
	qplt[p_sfc.<lvl[ip]] .= NaN
md"If surface pressure is lower than pressure level, set data to NaN for plotting."
end

# ╔═╡ 872b3c7c-5d0e-11eb-18be-0dc521676c0d
begin
	pplt.close(); f,axs = pplt.subplots(ncols=1,aspect=6,axwidth=6);
	
	c = axs[1].contourf(
		lon,lat,qplt',
		levels=vcat(
			-100,-50,-20,-10,-5,-2,-1,-0.5,-0.2,
			-0.1,0,0.1,0.2,0.5,1,2,5,10,20,50,100
			)*1e-8,
		cmap="drywet",cmap_kw=Dict("cut"=>0.1),extend="both"
	)
	axs[1].plot(x,y,c="k",lw=0.5)
	axs[1].plot(DTP[1],DTP[2],c="k",lw=1,linestyle="--")
	axs[1].plot(SEA[1],SEA[2],c="b",lw=1,linestyle="--")
	axs[1].plot(TRA[1],TRA[2],c="r",lw=1,linestyle="--")
	axs[1].plot(TR2[1],TR2[2],c="r",lw=1,linestyle="--")
	axs[1].plot(AMZ[1],AMZ[2],c="g",lw=1,linestyle="--")
	axs[1].plot(CRB[1],CRB[2],c="blue3",lw=1,linestyle="--")
	axs[1].format(xlim=(0,360),ylim=(-30,30))
	
	f.colorbar(c,loc="b",label=L"Moisture Flux Convergence / kg kg$^{-1}$")
	f.savefig(plotsdir("qfc_spatial-$(lvl[ip])hPa.png"),transparent=false,dpi=200)
	load(plotsdir("qfc_spatial-$(lvl[ip])hPa.png"))
end

# ╔═╡ Cell order:
# ╟─582bcb26-5a22-11eb-3b52-09e3ec08b3f9
# ╟─36a4e004-5a26-11eb-1ba4-db83d18c7845
# ╟─3a8bbab2-5a26-11eb-20e0-d551c2a68858
# ╟─2b736522-5a33-11eb-0eea-473a5b8566f3
# ╠═28ffa404-5a33-11eb-3dc7-7d716d6cf06f
# ╠═9e01c7ea-5a34-11eb-1c41-8b991a399584
# ╠═87a146d4-5a29-11eb-1571-b9d9266ee54b
# ╟─abde634c-5a83-11eb-00d6-69b233e0c00f
# ╟─4927b4c2-5a8a-11eb-37e4-53efb9445b23
# ╟─795a822a-5a29-11eb-37bb-3b6d62cec6be
# ╠═e3c519fc-5a35-11eb-0d34-cd4ca8b50bf9
# ╟─e0251f16-5a91-11eb-1bc8-a70691df76e2
# ╟─9d32ee54-5d07-11eb-0561-a16bf1f4bbec
# ╠═8fc1c5e2-5d07-11eb-1b2b-9f48919a115f
# ╟─3f21c35c-5d08-11eb-18c9-d90f1cb9ba2d
# ╟─368e8334-5d0f-11eb-3c66-c1ff6923e90c
# ╟─57df47c0-5d24-11eb-16aa-834dd327a306
# ╟─bac9bf96-5d0d-11eb-21f7-e3a1dda7af36
# ╠═872b3c7c-5d0e-11eb-18be-0dc521676c0d
