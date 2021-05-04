### A Pluto.jl notebook ###
# v0.14.2

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

# ╔═╡ 24abbe3a-5bb7-11eb-160b-1323efad463b
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 24601ef8-5bb7-11eb-1cd8-198dac960d3a
begin
	@quickactivate "TroPrecLS"
	using Statistics
	using StatsBase
	using PlutoUI
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 420e093c-5ba4-11eb-07da-c9a80044c8f1
md"
# 8b. Binning the Diurnal Cycle

In this notebook, we do some basic investigation into the relationship between slab depth and the amplitude of the diurnal cycle when the model is run in RCE state.

Relevant model parameters:
* SST = 301.7 K
* Insolation Peak = 1354.23 W m$^{-2}$
* Momentum Damping $a_m$ = 2
* Momentum Damping $a_m$ exponent = 0
"

# ╔═╡ 2b81f1ca-07eb-494d-b959-2c94632b6d5b
config = "Slab00d7"

# ╔═╡ b41d77a2-420e-4e10-86b7-8869c1b8ba2d
nmem = 5

# ╔═╡ 950f6d71-4f36-4ad2-bbd7-bb2c12ea7ce5
begin
	_,_,t = retrievedims("DiAmp064km",config,isensemble=true,member=1)
	t = t[9601:end]; nt = length(t)
	t_sst = zeros(nt,nmem)
	prcp  = zeros(nt,nmem)
	
	for imem = 1 : nmem
		fnc = outstatname("DiAmp064km",config,false,true,imem)
		t_sst[:,imem] = retrievevar("SST",fnc)[9601:end]
		prcp[:,imem]  = retrievevar("PREC",fnc)[9601:end]/24
	end
	
	t = repeat(t,outer=[nmem])
	t_sst = t_sst[:]
	prcp  = prcp[:]
	
	md"Loading sea surface temperature and precipitation data for $(config) configuration ..."
end

# ╔═╡ 7d6a3ec5-5fef-4cb8-962d-c573aa46bd2a
sstmax = 315

# ╔═╡ 9e5b9009-a18a-4f20-b169-f4abc020fe8f
sstmin = 290

# ╔═╡ 79b1debd-3365-4709-8fbe-4088608c9675
md"The minimum and maximum SSTs are $sstmin K and $sstmax K respectively"

# ╔═╡ 94a9f7c8-a7f1-4fcc-a0ef-7cf6bdc62c86
begin
	time_bins = 0:0.5:24;								nh = length(time_bins) - 1
	prcp_bins = 0:0.25:12.5;							np = length(prcp_bins) - 1
	tsst_bins = range(sstmin,stop=sstmax,length=51);	ns = length(tsst_bins) - 1
	prcp_hist = fit(Histogram,(mod.(t,1)*24,prcp),(time_bins,prcp_bins)).weights
	tsst_hist = fit(Histogram,(mod.(t,1)*24,t_sst),(time_bins,tsst_bins)).weights
	prcp_hist = prcp_hist/nt*np*nh; prcp_hist[iszero.(prcp_hist)] .= NaN
	tsst_hist = tsst_hist/nt*ns*nh; tsst_hist[iszero.(tsst_hist)] .= NaN
	mldepth = replace(config,"Slab"=>"")
	mldepth = replace(mldepth,"d"=>".")
	mldepth = parse(Float64,mldepth)
	md"Binning SST and Precipitation Rate data ..."
end

# ╔═╡ 016955c3-aa3a-42d3-9894-0b2c85221c5e
md"Make Figure? $(@bind mkfig PlutoUI.Slider(0:1))"

# ╔═╡ f48c072b-6a95-4b48-9123-cc722bb3e44a
begin
	if isone(mkfig)
		pplt.close(); f1,a1 = pplt.subplots(ncols=2,aspect=1.5,axwidth=2.5,sharey=0)
		
		lvls = [0.2,0.5,1,2,5,10,20,50]
		# lvls = 2:2:18

		c = a1[1].pcolormesh(
			0:0.5:24,prcp_bins,prcp_hist',
			cmap="Blues",levels=lvls,extend="both"
		)
		a1[1].pcolormesh(
			-24:0.5:0,prcp_bins,prcp_hist',
			cmap="Blues",levels=lvls,extend="both"
		)
		a1[1].format(ltitle=L"(a) Precipitation Rate / mm hr$^{-1}$")
		f1.colorbar(c,loc="r",ticks=[])

		c = a1[2].pcolormesh(
			0:0.5:24,tsst_bins,tsst_hist',
			cmap="Fire",levels=lvls,extend="both"
		)
		a1[2].pcolormesh(
			-24:0.5:0,tsst_bins,tsst_hist',
			cmap="Fire",levels=lvls,extend="both"
		)
		a1[2].format(ltitle="(b) Sea Surface Temperature / K")
		f1.colorbar(c,loc="r",ticks=[0.2,0.5,1,2,5,10,20,50],label="Density")

		for ax in a1
			ax.format(
				xlim=(-12,12),xlocator=-24:6:24,xlabel="Hour of Day",
				suptitle="Mixed Layer Depth: $(mldepth) m"
			)
		end

		f1.savefig(plotsdir("diurnalbin-$(config).png"),transparent=false,dpi=200)
	end
	PNGFiles.load(plotsdir("diurnalbin-$(config).png"))
end

# ╔═╡ Cell order:
# ╟─420e093c-5ba4-11eb-07da-c9a80044c8f1
# ╟─24abbe3a-5bb7-11eb-160b-1323efad463b
# ╟─24601ef8-5bb7-11eb-1cd8-198dac960d3a
# ╠═2b81f1ca-07eb-494d-b959-2c94632b6d5b
# ╠═b41d77a2-420e-4e10-86b7-8869c1b8ba2d
# ╟─950f6d71-4f36-4ad2-bbd7-bb2c12ea7ce5
# ╠═7d6a3ec5-5fef-4cb8-962d-c573aa46bd2a
# ╠═9e5b9009-a18a-4f20-b169-f4abc020fe8f
# ╟─79b1debd-3365-4709-8fbe-4088608c9675
# ╟─94a9f7c8-a7f1-4fcc-a0ef-7cf6bdc62c86
# ╟─016955c3-aa3a-42d3-9894-0b2c85221c5e
# ╠═f48c072b-6a95-4b48-9123-cc722bb3e44a
