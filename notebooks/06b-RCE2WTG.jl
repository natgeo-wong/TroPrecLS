### A Pluto.jl notebook ###
# v0.12.21

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
# 6b. Transitioning from RCE to WTG

In this notebook, we investigate and develop a way to implement the WTG forcing gradually in the System of Atmospheric Modelling.  The sudden introduction of WTG large-scale forcing often causes a model to enter a \"shocked\" state that unnaturally forces the model into a different state.  Here, we develop a method that gradually increases the strength of the WTG momentum-damping parameter from a near/pseudo-RCE state.
"

# ╔═╡ c22f840c-5af6-11eb-084d-d30d8891f460
md"
### A. Implementing the WTG Approximation in Models

Blossey et al. [2009] were the first to implement to WTG approximation in SAM.  The strength of the WTG is affected by the momentum-damping coefficient $a_m$ in the below formula:

$$\frac{\partial}{\partial p}
\left( \frac{f^2+a_m^2}{a_m}\frac{\partial\omega'}{\partial p} \right)
\approx \frac{k^2R_d}{p} T_v'$$

In our CRM simulations, we assume that $f$, the coriolis parameter, is 0, within the tropics (and therefore negligible).  Thus, the formula simplifies itself to

$$\frac{\partial}{\partial p}
\left( a_m\frac{\partial\omega'}{\partial p} \right)
\approx \frac{k^2R_d}{p} T_v'$$

If we do a further simplification by assuming that $a_m$ is constant in height (which was noted in Blossey et al. [2009], though they favoured that $a_m \propto p/p_\text{ref}$), this further simplifies the equation down to

$$\frac{\partial^2\omega'}{\partial p^2} \approx \frac{k^2}{a_m} \frac{R_d}{p}T_v'$$

Therefore, we see that for a bigger $a_m$, the induced $\omega'$ (or the WTG-induced vertical adjustment) is smaller.  Thus, bigger $a_m$ would result in less deviation from an RCE solution.  In other words, the magnitude of the WTG forcing is proportional to $a_m^{-1}$.
"

# ╔═╡ 205c5c9a-5914-11eb-09af-4313417b0df2
wtgstrength(am::Real) = 1/am

# ╔═╡ 75ed9ab0-5af7-11eb-1cba-ad51d005bd0d
md"
### B. Implementing a smooth transition

Implementing the Weak-Temperature Gradient approximation suddenly can cause a \"shock\" to the model, which will then transfer into a different regime state.  Therefore, we implement a smooth transition from a pseudo-RCE state ($a_m\gg1$) to a WTG state ($a_m = a_{m.0}$).  This was by varying $a_m$ in the form of an error function.  This will allow the increase in the strength of the WTG approximation to taper off gently, instead of increasing more and more rapidly, as it approaches $a_{m.0}$.

The error function is given by

$$\text{erf}\> t = \frac{2}{\sqrt{t}} \int_0^t e^{-x^2} \>\text{d}{x}$$

We use the error function to vary the momentum damping parameter $a_m$, or more specifically, the parameter $x$ where $x = \log_2a_m$.

The formula for $x$ as a function of time $t$, is given by

$$x(t) = \left( \frac{10-\log_2a_{m.0}}{2} \right)
\text{erf}\> \left( \frac{t_\text{max}/2-t}{t_\text{max}/5} \right)
+ 5 + \frac{1}{2}\log_2a_{m.0}$$
$$a_m(t) = 2^{x(t)}$$

Where $t_\text{max}$ is a scaling parameter, that defines the time when $a_m \approx a_{m.0}$ as a fraction of the total time that the WTG approximation is applied.  So if the WTG approximation is applied over 100 days, then if $t_\text{max} = 0.25$, $a_m \approx a_{m.0}$ occurs at day 25.
"

# ╔═╡ e7bd22f2-5ae9-11eb-1f82-09333889a5b2
twtg_scale = 0.25

# ╔═╡ 310fd792-5914-11eb-0d81-bd1eb22fab05
begin
	pplt.close(); f,axs = pplt.subplots(nrows=3,aspect=3,axwidth=3,sharey=0);
	
	for am in 2 .^(1:9)
		tmax = 1; twtg_max = twtg_scale * tmax
		t = 0:0.001:tmax; 
		test = erf.((twtg_max/2 .-t)/(twtg_max/5)) * (10-log(2,am))/2 .+ 5 .+ log(2,am)/2
		amc = 2 .^(test)
		wtg = wtgstrength.(amc)

		axs[1].plot(t,test,c="k")
		axs[2].plot(t,amc,c="k")
		axs[3].plot(t,wtg,c="k")
	end
	
	axs[1].plot([1,1]*twtg_scale,[0,10],c="r")
	axs[2].plot([1,1]*twtg_scale,[0,1024],c="r")
	axs[3].plot([1,1]*twtg_scale,[0,1],c="r")
	
	axs[1].format(ylim=(0,10),ylabel=L"$\log_2 (a_m)$")
	axs[2].format(ylim=(0,1024),ylabel=L"a_m")
	axs[3].format(ylim=(0,1),ylabel=L"WTG $\propto$ $a_m^{-1}$",#yscale="log",
		xlim=(0,1),xlabel="Nondimensionalized time",
		suptitle="twtg_scale = $(twtg_scale)")
	
	f.savefig(plotsdir("wtgcoeff.png"),transparent=false,dpi=200)
	load(plotsdir("wtgcoeff.png"))
end

# ╔═╡ 2cc46a1a-5b35-11eb-2b84-d536d6f7f7d0
md"
### C. Implementation in SAM
"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
exp = "D"

# ╔═╡ bdfa9872-5b35-11eb-059e-ad1ac171d295
am_wtg = 0

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	expname = "WTGam$(exp)"
	configvec = [
		"damping001",
		# "damping002",
		# "damping004",
		# "damping008",
		# "damping016",
		# "damping032",
		# "damping064",
		# "damping128",
	]
	ncon = length(configvec)
	blues = pplt.Colors("Blues",(ncon+2))
	reds  = pplt.Colors("Reds",(ncon+2))
	teals = pplt.Colors("Teal",(ncon+2))
	lgd = Dict("frame"=>false,"ncols"=>4)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 223b4286-8811-11eb-0e67-4da65e1999a5
function daymean(data)
	
	return dropdims(mean(reshape(data,24,:),dims=1),dims=1)
	
end

# ╔═╡ 55230f4a-7661-11eb-1c37-8b022b95e08e
begin
	pplt.close()
	fts,ats = pplt.subplots(nrows=3,aspect=3,axwidth=3,hspace=0.2,sharey=0)
	
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
				td = daymean(t)
				rn = retrievevar("PREC",fnc)
				pr = retrievevar("AREAPREC",fnc)
				rn = rn ./ pr
				sw = daymean(retrievevar("SWNS",fnc))
				lw = daymean(retrievevar("LWNS",fnc))
				sh = daymean(retrievevar("SHF",fnc))
				lh = daymean(retrievevar("LHF",fnc))
				pw = daymean(retrievevar("PW",fnc))
				seb = sw .- lw .- sh .- lh
				ats[1].plot(t,rn,lw=1,color=blues[ic+1])
				ats[2].plot(td,seb,lw=1,color=reds[ic+1])
				if imem == 1
					ats[3].plot(
						td,pw,lw=1,color=teals[ic+1],
						label=(L"$a_m =$" * " $config"),
						legend="b",legend_kw=lgd
					)
				else
					ats[3].plot(td,pw,lw=1,color=teals[ic+1])
				end
			end
		end
		
	end
	
	for imem = 1 : 5
		fnc = outstatname(expname,"dampingInf",false,true,imem)
		if isfile(fnc)
			_,_,t = retrievedims(fnc); t = t .- 80
				t  = daymean(t)
				pr = daymean(retrievevar("PREC",fnc))
				sw = daymean(retrievevar("SWNS",fnc))
				lw = daymean(retrievevar("LWNS",fnc))
				sh = daymean(retrievevar("SHF",fnc))
				lh = daymean(retrievevar("LHF",fnc))
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
		xlim=(0,250),xlabel="Time / Days",
		# ylocator=10. .^(-3:3),ylabel=L"Rainfall / mm day$^{-1}$",yscale="log",
		# ylim=(5,200),
		suptitle=expname
	)
	
	ats[2].format(
		xlim=(0,250),xlabel="Time / Days",
		ylim=(-100,100),ylabel=L"SEB / W m$^{-2}$",ylocator=(-3:3)*100,
		suptitle=expname
	)
	
	ats[3].format(
		xlim=(225,250),xlabel="Time / Days",
		ylim=(0,80),ylabel="PW / mm",
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
# ╟─c22f840c-5af6-11eb-084d-d30d8891f460
# ╟─205c5c9a-5914-11eb-09af-4313417b0df2
# ╟─75ed9ab0-5af7-11eb-1cba-ad51d005bd0d
# ╠═e7bd22f2-5ae9-11eb-1f82-09333889a5b2
# ╟─310fd792-5914-11eb-0d81-bd1eb22fab05
# ╟─2cc46a1a-5b35-11eb-2b84-d536d6f7f7d0
# ╠═d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╠═bdfa9872-5b35-11eb-059e-ad1ac171d295
# ╠═a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─223b4286-8811-11eb-0e67-4da65e1999a5
# ╠═55230f4a-7661-11eb-1c37-8b022b95e08e
