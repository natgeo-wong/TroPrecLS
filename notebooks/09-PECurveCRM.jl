### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 247edb6d-65a6-4fb7-beab-46604b78cfe9
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 5ebccf8b-f701-436c-b654-9acb6856473d
begin
	@quickactivate "TroPrecLS"
	using Statistics
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 440fd8da-8f4a-11eb-108e-993e7f81183e
md"
# 9. Constructing a PE Curve using SAM

We have explored the impacts of varying the slab depth and the resulting response of the diurnal cycle of temperature and precipitation in the notebook `08b-diurnal.jl`.

Relevant model parameters:
* SST = 301.7 K
* Insolation Peak = 1354.23 W m$^{-2}$
* Momentum Damping $a_m$ = 2
* Momentum Damping $a_m$ exponent = 0
"

# ╔═╡ d8248390-1416-4322-bebe-9b913e4cdc1c
md"
### A. Loading Data for WTG Control Experiments
"

# ╔═╡ 41ed6c57-c7fc-4bbe-bcb7-08414d9f84b3
function retrieve2D(variable,configlist)
	
	nconfig = length(configlist);
	d2D  = zeros(4800,1,nconfig)
	
	for ic = 1 : nconfig
		for ii = 1 : 1
			fnc = outstatname("Control",configlist[ic],false,true,ii)
			d2D[:,ii,ic] .= retrievevar(variable,fnc)[9601:end]
		end
	end
	
	d2D = reshape(d2D,4800*1,nconfig)
	
	return d2D
	
end

# ╔═╡ eed5fd10-e1c8-448e-80c7-44eea38eca01
function retrievecsf(exp,configlist)
	
	nconfig = length(configlist);
	var  = zeros(4800,1,nconfig)
	
	for ic = 1 : nconfig
		for ii = 1
			fnc = outstatname(exp,configlist[ic],false,true,ii)
			tcw = retrievevar("PW",fnc)[9601:end] / 1000
			ta  = retrievevar("TABS",fnc)[:,9601:end]
			qv  = retrievevar("QV",fnc)[:,9601:end] / 1000
			pp  = retrievevar("p",fnc)
			rh  = calcrh(qv,ta,pp)
			_,swp = calccsf(rh,qv,pp)
			var[:,ii,ic] += tcw ./ swp
		end
	end
	
	var = reshape(var,4800*1,nconfig)
	
	return var
	
end

# ╔═╡ 4a1ad4f6-5a83-467f-823e-b3f5d92591b7
begin
	configs = [
		"slab00d05",#"slab00d07",
		#"slab00d10","slab00d14","slab00d20",
		#"slab00d32",
		"slab00d50",#"slab00d71","slab01d00","slab01d41",
		#"slab02d00","slab03d16",
		"slab05d00",#"slab07d07","slab10d00",
		#"slab14d14","slab20d00","slab31d62",
		"slab50d00",
	]
	ncon = length(configs)
	colors = pplt.Colors("blues",ncon+4)
	lndocn = pplt.Colors("Delta_r",ncon+4)
	md"Defining experimental configurations for WTG experiments ..."
end

# ╔═╡ 375f5fe4-cd08-4361-8976-577d17fd0ee3
begin
	prcp = retrieve2D("PREC",configs)
	csf  = retrievecsf("Control",configs) * 100
	ptst = retrievevar("PREC","Slab00d05","forcingn001",isensemble=true,member=1)[9601:end]
	ctst = retrievecsf("Slab00d05",["forcingn001"]) * 100
	ptst = reshape(ptst,4800,1)
	ctst = reshape(ctst,4800,1)
md"Loading results from the WTG Slab-depth experiments ..."
end

# ╔═╡ 73fcc28f-8b94-4d34-ab99-0c12a1b948ac
csf

# ╔═╡ 8d552f5e-42e5-4042-bc3c-a7600a1c8ed6
begin
	mean(ctst),mean(ptst)/24
end

# ╔═╡ 9f580ce2-4dc2-4d05-ad26-0382a5e17b54
md"
### B. Binning Precipitation Rate against Column Saturation Fraction
"

# ╔═╡ 9892a59f-6a61-48a6-aeb2-8aecd66416a4
function csfvprcpbin(prcp,csf,cvec,csep)
	pvec = zeros(length(cvec)); jj = 0;
    pfrq = zeros(Int64,length(cvec))
    for cii in cvec
        pii = @view prcp[ (csf.>(cii-csep)) .& (csf.<=(cii+csep)) ]
        jj = jj + 1; pvec[jj] = mean(pii); pfrq[jj] = length(pii)
    end

    return pvec,pfrq
end

# ╔═╡ 355247fb-e98f-4a20-b7b9-3ce6f227253a
begin
	cbin = collect(0:0.5:100); cstep = (cbin[2]-cbin[1])/2; nbins = length(cbin)
	pvec = zeros(nbins,ncon)
	pfrq = zeros(nbins,ncon)
	
	for icon = 1 : ncon
		pvec[:,icon],pfrq[:,icon] = csfvprcpbin(prcp[:,icon],csf[:,icon],cbin,cstep)
	end
	
	pvtt,pftt = csfvprcpbin(ptst,ctst,cbin,cstep)
end

# ╔═╡ 31c30b01-33e1-42f3-ae30-9eaf4a7fc8c3
begin
	pvec[isnan.(pvec)] .= 0
	pvtt[isnan.(pvtt)] .= 0
end

# ╔═╡ 84637231-78d9-41cc-bb6d-f12d12d42621
pfrq[:,1]

# ╔═╡ 3091c995-be77-4fa5-8517-86e4e4157482
pnew = (pvec[:,1] .+ pvtt) ./ (pfrq[:,1] .+ pftt)

# ╔═╡ 5fb1540c-d4bb-46bd-b3a6-a1253e799a79
begin
	pplt.close(); f,a = pplt.subplots(aspect=2,axwidth=3)
	
	for icon = 1 : 1
		a[1].plot(cbin,pvec[:,icon]./pfrq[:,icon],c=lndocn[icon+2],lw=1)
	end
	a[1].plot(cbin,pvtt./pftt,c="k",lw=1,linestyle=":")
	a[1].plot(cbin,pnew,c="k",lw=1)
	a[1].format(xlim=(00,100),ylim=10. .^(-4,2))
	a[1].format(yscale="log",ylocator=10. .^(-4:2))
	
	f.savefig("test.png",transparent=false,dpi=200)
	PNGFiles.load("test.png")
end

# ╔═╡ 55babf48-524d-42cf-b814-73be9361b727
begin
	pplt.close(); ff,af = pplt.subplots(aspect=2,axwidth=3)
	
	for icon = 1 : 4
		af[1].plot(cbin,pfrq[:,icon]./sum(pfrq[:,icon])*101,c=lndocn[icon+2],lw=1)
	end
	af[1].plot(cbin,pftt./sum(pftt)*101,c="k",lw=1,linestyle=":")
	# a[1].format(xlim=(00,100),ylim=10. .^(-4,2))
	# af[1].format(yscale="log",ylocator=10. .^(-4:2))
	
	ff.savefig("test.png",transparent=false,dpi=200)
	PNGFiles.load("test.png")
end

# ╔═╡ Cell order:
# ╟─440fd8da-8f4a-11eb-108e-993e7f81183e
# ╟─247edb6d-65a6-4fb7-beab-46604b78cfe9
# ╟─5ebccf8b-f701-436c-b654-9acb6856473d
# ╟─d8248390-1416-4322-bebe-9b913e4cdc1c
# ╠═41ed6c57-c7fc-4bbe-bcb7-08414d9f84b3
# ╠═eed5fd10-e1c8-448e-80c7-44eea38eca01
# ╠═4a1ad4f6-5a83-467f-823e-b3f5d92591b7
# ╠═375f5fe4-cd08-4361-8976-577d17fd0ee3
# ╠═73fcc28f-8b94-4d34-ab99-0c12a1b948ac
# ╠═8d552f5e-42e5-4042-bc3c-a7600a1c8ed6
# ╟─9f580ce2-4dc2-4d05-ad26-0382a5e17b54
# ╟─9892a59f-6a61-48a6-aeb2-8aecd66416a4
# ╟─355247fb-e98f-4a20-b7b9-3ce6f227253a
# ╟─31c30b01-33e1-42f3-ae30-9eaf4a7fc8c3
# ╠═84637231-78d9-41cc-bb6d-f12d12d42621
# ╠═3091c995-be77-4fa5-8517-86e4e4157482
# ╠═5fb1540c-d4bb-46bd-b3a6-a1253e799a79
# ╠═55babf48-524d-42cf-b814-73be9361b727
