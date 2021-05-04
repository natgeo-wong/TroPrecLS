### A Pluto.jl notebook ###
# v0.14.0

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

From notebook `06b-RCE2WTG.jl`, we see that adding the WTG approximation into a limited-area domain forces the initial RCE state into one or two regimes: a wet regime and a dry regime.  These regimes are analogues to the wet and dry regimes found in a large-area domain, where self-aggregation of convection naturally occurs in RCE.

In this notebook, we explore the characteristics of these wet and dry regimes under different $a_m$ strength.
"

# ╔═╡ b6892634-9199-11eb-38d5-8dda8da16ed7
md"
### A. The Wet and Dry Regimes in WTG Ensemble Runs

We recall the momentum damping equation (assuming $a_m$ is constant with height):

$$\frac{\partial^2\omega'}{\partial p^2} \approx \frac{k^2}{a_m} \frac{R_d}{p}T_v'$$

Recall that $k$ is the wavenumber.  Therefore, increasing $a_m$ increases the wavenumber required for the WTG response to be the same.  Or in other words:

$$k' = \frac{k}{\sqrt{a_m}}$$

The pseudo-wavenumber $k'$ is smaller than the actual wavenumber $k$, which implies that the wavelength is much, much larger.  This is an analogue to the domain being farther away from the baseline RCE domain, which is taken to be a large-scale domain average.  So as $a_m$ increases, we should see the dry and wet states converge back into the initial RCE state.
"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
expi = "P"

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	expname = "WTGam$(expi)"
	configvec = [
		"damping001",
		"damping002",
		"damping004",
		"damping008",
		"damping016",
		"damping032",
		"damping064",
		"damping128",
		"damping256",
		"damping512",
	]
	ncon = length(configvec)
	blues = pplt.Colors("Blues",(ncon+2))
	grns  = pplt.Colors("Teal",(ncon+2))
	brwns = pplt.Colors("Brown",(ncon+2))
	lgd = Dict("frame"=>false,"ncols"=>4)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 55230f4a-7661-11eb-1c37-8b022b95e08e
begin
	pplt.close()
	fts,ats = pplt.subplots(ncols=4,aspect=0.5,axwidth=1.2,sharex=0)
	
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
				pr = retrievevar("PREC",fnc)./24
				pa = retrievevar("AREAPREC",fnc)
				pw = retrievevar("PW",fnc)
				pra = pr./pa; pra[isnan.(pra)] .= 0; pra[pra.==Inf] .= 0
				ats[1].scatter(mean(pr[(end-99):end]),config,lw=1,color=blues[ic+1])
				ats[2].scatter(mean(pa[(end-99):end]),config,lw=1,color=blues[ic+1])
				ats[3].scatter(mean(pra[(end-99):end]),config,lw=1,color=blues[ic+1])
				ats[4].scatter(mean(pw[(end-99):end]),config,lw=1,color=blues[ic+1])
			end
		end
		
	end
	
	for imem = 1 : 10
		fnc = outstatname("Control","$(expi)INSOL",false,true,imem)
		if isfile(fnc)
			_,_,t = retrievedims(fnc); t = t .- 80
			pr = retrievevar("PREC",fnc)./24
			pa = retrievevar("AREAPREC",fnc)
			pw = retrievevar("PW",fnc)
			pra = pr./pa; pra[isnan.(pra)] .= 0; pra[pra.==Inf] .= 0
			ats[1].plot([1,1]*mean(pr[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[2].plot([1,1]*mean(pa[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[3].plot([1,1]*mean(pra[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[4].plot([1,1]*mean(pw[(end-499):end]),[0,2000],c="grey",lw=1)
		end
	end
	
	ats[1].format(
		ylim=(0.5,2000),ylabel=L"$a_m$ / day$^{-1}$",yscale="log",
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.1),
		xlim=(0,10),xlabel=L"Domain $P$ / mm hr$^{-1}$",
		suptitle=L"Sensitivity to $a_m$ | " * "$(expname)",
		ultitle="(a)"
	)
	
	ats[2].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.1),
		xlim=(0,1),xlabel="Rain Area Fraction",
		ultitle="(b)"
	)
	
	ats[3].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>1),
		xlim=(0,10),xlabel=L"Rain Area $P$ / mm hr$^{-1}$",
		ultitle="(c)"
	)
	
	ats[4].format(
		xlim=(0,75),xlabel="PWV / mm",
		ultitle="(d)"
	)
	
	for ax in ats
		ax.format(lrtitle="Wet",lltitle="Dry")
	end
	
	fts.savefig(plotsdir(
		"wtgstrength-$(expname).png"),
		transparent=false,dpi=200
	)
	load(plotsdir("wtgstrength-$(expname).png"))
end

# ╔═╡ 9cf4fa56-91a8-11eb-2710-955eefd10142
md"
We see the following:
* As $a_m$ increases, both the wet and dry model regimes converge back into the initial RCE state.
* When $a_m$ becomes too small, all model states collapse into a dry regime.
* As $a_m$ decreases, the change in rain-area fraction in the domain, rather than rain-rate in the rain-area, which is responsible for changes in the domain-averaged precipitation.

In conclusion, the overall transition from RCE to WTG forcing (decreasing $a_m$) is as follows:
1. Domain becomes more moist / heavier precipitation / wetter
2. Eventually, a dry regime separates out
3. Wet regime begins to show drier characteristics
4. Eventually, wet regime crosses a threshold and model fully enters into a dry-regime.

We note that the simulations that end up in the wet regime are more numerous than those that end up in the dry regime.  Overall, assuming complete randomness this seems to indicate that the a wet regime is favoured.

I have decided on using precipitable water as the prognostic variable for determining if the model is in a dry or wet regime.
"

# ╔═╡ 364a1ce8-91ba-11eb-29a8-b948110e6125
md"
### B. Exploring some Variables
"

# ╔═╡ 489b5bea-91b4-11eb-358b-3fe61c900511
begin
	pplt.close()
	feb,aeb = pplt.subplots(ncols=5,aspect=0.5,axwidth=1.2,sharex=0); pw = zeros(10)
	
	for imem = 1 : 10
		fnc = outstatname("Control","$(expi)INSOL",false,true,imem)
		if isfile(fnc)
			_,_,t = retrievedims(fnc); t = t .- 80
			sw = retrievevar("SWNS",fnc)
			lw = -retrievevar("LWNS",fnc)
			sh = -retrievevar("SHF",fnc)
			lh = -retrievevar("LHF",fnc)
			sb = sw .+ lw .+ sh .+ lh
			pw[imem] = mean(retrievevar("PW",fnc)[(end-499):end])
			aeb[1].plot([1,1]*mean(sw[(end-499):end]),[0,2000],c="grey",lw=1)
			aeb[2].plot([1,1]*mean(lw[(end-499):end]),[0,2000],c="grey",lw=1)
			aeb[3].plot([1,1]*mean(sh[(end-499):end]),[0,2000],c="grey",lw=1)
			aeb[4].plot([1,1]*mean(lh[(end-499):end]),[0,2000],c="grey",lw=1)
			aeb[5].plot([1,1]*mean(sb[(end-499):end]),[0,2000],c="grey",lw=1)
		end
	end
	
	pw = mean(pw)
	
	for ic in 1 : ncon
		icon = configvec[ic]
		icon = replace(icon,"damping"=>"")
		icon = replace(icon,"d"=>".")
		icon = parse(Float64,icon)
		imem = 0
		
		while imem < 100; imem += 1
			fnc = outstatname(expname,configvec[ic],false,true,imem)
			if isfile(fnc)
				_,_,t = retrievedims(fnc); t = t .- 80
				sw = retrievevar("SWNS",fnc)
				lw = -retrievevar("LWNS",fnc)
				sh = -retrievevar("SHF",fnc)
				lh = -retrievevar("LHF",fnc)
				sb = sw .+ lw .+ sh .+ lh
				pwi = mean(retrievevar("PW",fnc)[(end-99):end])
				if pwi < pw
					aeb[1].scatter(mean(sw[(end-99):end]),icon,lw=1,color=brwns[ic+1])
					aeb[2].scatter(mean(lw[(end-99):end]),icon,lw=1,color=brwns[ic+1])
					aeb[3].scatter(mean(sh[(end-99):end]),icon,lw=1,color=brwns[ic+1])
					aeb[4].scatter(mean(lh[(end-99):end]),icon,lw=1,color=brwns[ic+1])
					aeb[5].scatter(mean(sb[(end-99):end]),icon,lw=1,color=brwns[ic+1])
				else
					aeb[1].scatter(mean(sw[(end-99):end]),icon,lw=1,color=grns[ic+1])
					aeb[2].scatter(mean(lw[(end-99):end]),icon,lw=1,color=grns[ic+1])
					aeb[3].scatter(mean(sh[(end-99):end]),icon,lw=1,color=grns[ic+1])
					aeb[4].scatter(mean(lh[(end-99):end]),icon,lw=1,color=grns[ic+1])
					aeb[5].scatter(mean(sb[(end-99):end]),icon,lw=1,color=grns[ic+1])
				end
			end
		end
		
	end
	
	aeb[1].format(
		ylim=(0.5,2000),ylabel=L"$a_m$ / day$^{-1}$",yscale="log",
		xlim=(0,400),xlabel="Net Shortwave",
		suptitle=L"Energy Balance / W m$^{-2}$ | " * "$(expname)",
		ultitle="(a)"
	)
	
	aeb[2].format(xlim=(-150,0),xlabel="Net Longwave",ultitle="(b)")
	aeb[3].format(xlim=(-75,0),xlabel="Sensible Heat",ultitle="(c)")
	aeb[4].format(xlim=(-300,0),xlabel="Latent Heat",ultitle="(d)")
	aeb[5].format(xlim=(-350,350),xlabel="Surface Balance",ultitle="(e)")
	
	feb.savefig(plotsdir(
		"wtgstrength-$(expname)-seb.png"),
		transparent=false,dpi=200
	)
	load(plotsdir("wtgstrength-$(expname)-seb.png"))
end

# ╔═╡ e967eb5c-91c4-11eb-3066-05ccaa40bd11
begin
	pplt.close()
	f3D,a3D = pplt.subplots(ncols=4,aspect=0.5,axwidth=1.2,sharex=0)
	clc = zeros(64)
	tab = zeros(64)
	
	for ic in 1 : ncon
		icon = configvec[ic]
		icon = replace(icon,"damping"=>"")
		icon = replace(icon,"d"=>".")
		icon = parse(Float64,icon)
		imem = 0
		
		while imem < 100; imem += 1
			fnc = outstatname(expname,configvec[ic],false,true,imem)
			if isfile(fnc)
				_,p,t = retrievedims(fnc); t = t .- 80
				clci = mean(retrievevar("CLD",fnc)[:,(end-99):end],dims=2)*100
				tabi = mean(retrievevar("TABS",fnc)[:,(end-99):end],dims=2)
				qvi  = mean(retrievevar("QV",fnc)[:,(end-99):end],dims=2) / 10
				rhi  = calcrh(qvi,tabi,p)
				wwtg = mean(retrievevar("WWTG",fnc)[:,(end-99):end],dims=2) * 3.6
				pwi  = mean(retrievevar("PW",fnc)[(end-99):end])
				relh = mean(retrievevar("RELH",fnc)[:,(end-99):end],dims=2)
				if pwi < pw
					a3D[1].plot(dropdims(clci,dims=2),p,lw=1,c=brwns[ic+1])
					a3D[2].plot(dropdims(tabi,dims=2),p,lw=1,c=brwns[ic+1])
					a3D[3].plot(dropdims(wwtg,dims=2),p,lw=1,c=brwns[ic+1])
					a3D[4].plot(dropdims(rhi,dims=2),p,lw=1,c=brwns[ic+1])
					# a3D[4].plot(dropdims(relh,dims=2),p,lw=1,c=brwns[ic+1])
				else
					a3D[1].plot(dropdims(clci,dims=2),p,lw=1,c=grns[ic+1])
					a3D[2].plot(dropdims(tabi,dims=2),p,lw=1,c=grns[ic+1])
					a3D[3].plot(dropdims(wwtg,dims=2),p,lw=1,c=grns[ic+1])
					a3D[4].plot(dropdims(rhi,dims=2),p,lw=1,c=grns[ic+1])
					# a3D[4].plot(dropdims(relh,dims=2),p,lw=1,c=grns[ic+1])
				end
			end
		end
		
	end
	
	for imem = 1 : 10
		fnc = outstatname("Control","$(expi)INSOL",false,true,imem)
		if isfile(fnc)
			_,p,_ = retrievedims(fnc);
			clc[:] = mean(retrievevar("CLD",fnc)[:,(end-499):end],dims=2) * 100
			tab[:] = mean(retrievevar("TABS",fnc)[:,(end-499):end],dims=2)
			qv = mean(retrievevar("QV",fnc)[:,(end-99):end],dims=2) / 10
			rh = calcrh(qv,tab,p)
			a3D[1].plot(clc,p,color="k")
			a3D[2].plot(tab,p,color="k")
			a3D[4].plot(dropdims(rh,dims=2),p,color="k")
		end
	end
	
	a3D[1].format(
		ylim=(1000,20),ylabel="Pressure / hPa",yscale="log",
		xlim=(0,100),xlabel="Cloud Fraction / %",
		suptitle="3D Vertical Profiles | $(expname)",ultitle="(a)"
	)
	
	a3D[2].format(xlim=(150,325),xlabel="Temperature / K",ultitle="(b)")
	a3D[3].format(xlim=(-2,2),xlabel=L"$w_{WTG}$ / km hr$^{-1}$",xscale="symlog",
	xscale_kw=Dict("linthresh"=>0.1),ultitle="(c)")
	a3D[4].format(xlim=(0,110),xlabel="Relative Humidity / %",ultitle="(d)")
	# a3D[5].format(xlim=(-350,350),xlabel="Surface Balance",)
	
	f3D.savefig(plotsdir(
		"wtgstrength-$(expname)-3D.png"),
		transparent=false,dpi=200
	)
	load(plotsdir("wtgstrength-$(expname)-3D.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─b6892634-9199-11eb-38d5-8dda8da16ed7
# ╠═d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─55230f4a-7661-11eb-1c37-8b022b95e08e
# ╟─9cf4fa56-91a8-11eb-2710-955eefd10142
# ╟─364a1ce8-91ba-11eb-29a8-b948110e6125
# ╟─489b5bea-91b4-11eb-358b-3fe61c900511
# ╟─e967eb5c-91c4-11eb-3066-05ccaa40bd11
