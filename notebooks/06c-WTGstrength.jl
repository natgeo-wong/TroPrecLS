### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 681658b0-5914-11eb-0d65-bbace277d145
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 6dce35fc-5914-11eb-0ce2-0d4e164e1898
begin
	@quickactivate "TroPrecLS"
	using SpecialFunctions
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 6c. Momentum Damping Strength

In this notebook, we investigate and develop a way to implement the WTG forcing gradually in the System of Atmospheric Modelling.  The sudden introduction of WTG large-scale forcing often causes a model to enter a \"shocked\" state that unnaturally forces the model into a different state.  Here, we develop a method that gradually increases the strength of the WTG momentum-damping parameter from a near/pseudo-RCE state.
"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
exp = "3P"

# ╔═╡ bdfa9872-5b35-11eb-059e-ad1ac171d295
am_wtg = 0

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	expname = "$(exp)WTGamExp$(am_wtg)"
	configvec = [
		"damping001",
		"damping002",
		"damping004",
		"damping008",
		"damping016",
		"damping032",
		"damping064",
		"damping128",
	]
	ncon = length(configvec)
	blues = pplt.Colors("Blues",(ncon+2))
	reds  = pplt.Colors("Reds",(ncon+2))
	lgd = Dict("frame"=>false,"ncols"=>4)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 55230f4a-7661-11eb-1c37-8b022b95e08e
begin
	pplt.close()
	fts,ats = pplt.subplots(nrows=2,aspect=3,axwidth=4,hspace=0.2,sharey=0)
	
	for ic in 1 : ncon
		config = configvec[ic]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		imem = 0
		
		while imem < 100; imem += 1
			fnc = outstatname(expname,configvec[ic],false,true,imem)
			if isfile(fnc)
				_,_,t = retrievedims(fnc); t = t .- 80
				pr = retrievevar("PREC",fnc)
				sw = retrievevar("SWNS",fnc); lw = retrievevar("LWNS",fnc)
				sh = retrievevar("SHF",fnc);  lh = retrievevar("LHF",fnc)
				seb = sw .- lw .- sh .- lh
				ats[1].plot(t,pr,lw=1,color=blues[ic+1])
				if imem == 1
					ats[2].plot(
						t,seb,lw=1,color=reds[ic+1],
						label=(L"$a_m =$" * " $config"),
						legend="b",legend_kw=lgd
					)
				else
					ats[2].plot(t,seb,lw=1,color=reds[ic+1])
				end
			end
		end
		
	end
	
	for imem = 1 : 5
		fnc = outstatname(expname,"dampingInf",false,true,imem)
		if isfile(fnc)
			_,_,t = retrievedims(fnc); t = t .- 80
			pr = retrievevar("PREC",fnc)
			sw = retrievevar("SWNS",fnc); lw = retrievevar("LWNS",fnc)
			sh = retrievevar("SHF",fnc);  lh = retrievevar("LHF",fnc)
			seb = sw .- lw .- sh .- lh
			ats[1].plot(t,pr,lw=1,color="k")
			if imem == 1
				ats[2].plot(t,seb,lw=1,color="k",label="RCE",legend="b")
			else
				ats[2].plot(t,seb,lw=1,color="k")
			end
		end
	end
	
	ats[1].format(
		xlim=(0,400),xlabel="Time / Days",
		ylim=(0.01,200),ylabel=L"Rainfall Rate / mm day$^{-1}$",yscale="log",
		suptitle=expname
	)
	
	ats[2].format(
		xlim=(0,400),xlabel="Time / Days",
		ylim=(-350,350),ylabel=L"SEB / W m$^{-2}$",ylocator=(-3:3)*100,
		suptitle=expname
	)
	
	fts.savefig(plotsdir(
		"rce2wtg-$(expname).png"),
		transparent=false,dpi=200
	)
	load(plotsdir("rce2wtg-$(expname).png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╠═d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╠═bdfa9872-5b35-11eb-059e-ad1ac171d295
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╠═55230f4a-7661-11eb-1c37-8b022b95e08e
