### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

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
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 420e093c-5ba4-11eb-07da-c9a80044c8f1
md"
# 8b. Amplitude of the Diurnal Cycle

In this notebook, we do some basic investigation into the relationship between slab depth and the amplitude of the diurnal cycle when the model is run in RCE state.

Relevant model parameters:
* SST = 301.7 K
* Insolation Peak = 1354.23 W m$^{-2}$
* Momentum Damping $a_m$ = 2
* Momentum Damping $a_m$ exponent = 0
"

# ╔═╡ 57f52568-5bb9-11eb-1e7f-c34b6efe0bac
function retrieve2D(variable,configlist,beg,tstep,tshift)
	
	nconfig = length(configlist);
	d2D  = zeros(tstep+2,nconfig)
	
	for iconfig = 1 : nconfig
		var = retrievevar2D(variable,"DiAmp064km",configlist[iconfig])
		var = dropdims(mean(var,dims=(1,2)),dims=(1,2))
		d2D[:,iconfig] = diurnal2D(var[(end-beg):end],tstep,tshift);
	end
	
	return d2D
	
end

# ╔═╡ c9859e5e-8feb-11eb-0190-79f138ffc1ab
function retrieve2Dbin(variable,configlist)
	
	nconfig = length(configlist); pbin = 0:5:250
	d2D  = zeros(48,length(pbin)-1,nconfig)
	
	for icon = 1 : nconfig
		x,y,t = retrievedims2D("DiAmp064km",configlist[icon])
		nx = length(x); ny = length(y); nt = length(t)
		var = retrievevar2D(variable,"DiAmp064km",configlist[icon])
		var = dropdims(mean(var,dims=(1,2)),dims=(1,2))
		d2D[:,:,icon] += fit(Histogram,(mod.(t.+1/96,1)*24,var),(0:0.5:24,pbin)).weights
	end
	
	return d2D
	
end

# ╔═╡ 4f40bd7e-5bb9-11eb-34f2-91e1f959c59a
function plot2DWTG(axs,axsii,td,prcp,configlist,colors;islegend::Bool=false)
	
	nconfig = length(configlist)
	
	for iconfig = 1 : nconfig
		
		config = configlist[iconfig]
		config = replace(config,"Slab"=>"")
		config = replace(config,"d"=>".")
		if config != "Inf."
			  config = parse(Float64,config)
		else; config = L"\infty"
		end
		
		if islegend
			axs[axsii].plot(
				td,prcp[:,iconfig],color=colors[iconfig+2],lw=1,
				label="$(config) m",
				legend="b",legend_kw=Dict("frame"=>false,"ncols"=>5)
			)
		else
			axs[axsii].plot(td,prcp[:,iconfig],color=colors[iconfig+2],lw=1)
		end
		
	end
	
end

# ╔═╡ 1ab3ad70-5c2a-11eb-187b-75e524a4581f
begin
	configs = [
		"Slab00d1","Slab00d14","Slab00d2","Slab00d3","Slab00d5",
		"Slab00d7","Slab01d0","Slab01d4","Slab02d0","Slab03d2",
		"Slab05d0","Slab07d1","Slab10d0","Slab14d1","Slab20d0",
		"Slab31d6","Slab50d0","Slab70d7"
	]
	ncon = length(configs)
	colors = pplt.Colors("blues",ncon+4)
end

# ╔═╡ 6e1170f4-5bb9-11eb-0d38-61befdd2ad88
begin
	ndy = 50
md"Loading dimensions from Control Experiment ..."
end

# ╔═╡ d13eeaf8-5e8e-11eb-006e-2d45aa9253fd
begin
	_R,_R,tRCE = retrievedims("Control","DINSOL",isensemble=true,member=1)
	tdRCE,tstepRCE,tshiftRCE,begRCE = t2d(tRCE,ndy); nt = length(tRCE)
	
	prcpRCE = zeros(nt)
	sstRCE  = zeros(nt)
	
	for imem = 1 : 10
		global prcpRCE += retrievevar("PREC","Control","DINSOL",isensemble=true,member=1)
		global sstRCE  += retrievevar("SST","Control","DINSOL",isensemble=true,member=1)
	end
	
	prcpRCE = prcpRCE / 10; prcpRCE = diurnal2D(prcpRCE,tstepRCE,tshiftRCE)
	sstRCE  = sstRCE  / 10; sstRCE  = diurnal2D(sstRCE,tstepRCE,tshiftRCE)
md"Loading results from the RCE Control ..."
end

# ╔═╡ 7b7a5806-8f20-11eb-1e95-3d1834134375
function retrievecsf(configlist,beg,tstep,tshift)
	
	nconfig = length(configlist);
	var  = zeros(nconfig)
	
	for iconfig = 1 : nconfig
		tcw = retrievevar("PW","DiAmp064km",configlist[iconfig]) / 1000
		ta  = retrievevar("TABS","DiAmp064km",configlist[iconfig])
		qv  = retrievevar("QV","DiAmp064km",configlist[iconfig]) / 1000
		pp  = retrievevar("p","DiAmp064km",configlist[iconfig])
		rh  = calcrh(qv,ta,pp)
		_,swp = calccsf(rh,qv,pp)
		var[iconfig] = mean(tcw[(end-beg):end]./swp[(end-beg):end])
	end
	
	return var
	
end

# ╔═╡ 1a31f88e-5c2a-11eb-0c1f-257863a43cf5
begin
	_W,_W,tWTG = retrievedims2D("DiAmp064km","Slab01d0")
	tdWTG,tstepWTG,tshiftWTG,begWTG = t2d(tWTG,ndy);
	prcpWTG = retrieve2D("Prec",configs,begWTG,tstepWTG,tshiftWTG)
	sstWTG  = retrieve2D("SST",configs,begWTG,tstepWTG,tshiftWTG)
	tcwWTG  = retrieve2D("PW",configs,begWTG,tstepWTG,tshiftWTG)
	csfWTG  = retrievecsf(configs,begWTG,tstepWTG,tshiftWTG)
md"Loading results from the WTG Slab-depth experiments ..."
end

# ╔═╡ 7d401cd6-5c2e-11eb-3842-917545e546ef
begin
	
	pplt.close(); f,axs = pplt.subplots(nrows=2,axwidth=4,aspect=2,sharey=0)
	
	plot2DWTG(axs,1,tdWTG/2,sstWTG,configs,colors)
	plot2DWTG(axs,1,tdWTG/2 .-24,sstWTG,configs,colors)
	plot2DWTG(axs,2,tdWTG/2,prcpWTG,configs,colors)
	plot2DWTG(axs,2,tdWTG/2 .-24,prcpWTG,configs,colors,islegend=true)
	
	axs[1].plot(tdRCE,sstRCE,c="k")
	axs[2].plot(tdRCE,prcpRCE,c="k",label="RCE",legend="b")
	
	axs[1].format(ylim=(290,315),urtitle="SST / K")
	axs[2].format(ylim=(0,75),urtitle=L"Rainfall Rate / mm day$^{-1}$")
	
	for ax in axs
		ax.format(xlim=(-12,12),xlabel="Hour of Day",xlocator=-24:3:24)
	end
	
	f.savefig(plotsdir("DiurnalAmp.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("DiurnalAmp.png"))
	
end

# ╔═╡ d7735d46-5e9a-11eb-0fb3-91be12c3508b
begin
	sstWTG_mean  = dropdims(mean(sstWTG[1:(end-1),:],dims=1),dims=1)
	sstWTG_max   = dropdims(maximum(sstWTG[1:(end-1),:],dims=1),dims=1)
	sstWTG_min   = dropdims(minimum(sstWTG[1:(end-1),:],dims=1),dims=1)
	prcpWTG_mean = dropdims(mean(prcpWTG[1:(end-1),:],dims=1),dims=1)
	prcpWTG_max  = dropdims(maximum(prcpWTG[1:(end-1),:],dims=1),dims=1)
	prcpWTG_min  = dropdims(minimum(prcpWTG[1:(end-1),:],dims=1),dims=1)
	tcwWTG_mean  = dropdims(mean(tcwWTG[1:(end-1),:],dims=1),dims=1)
end

# ╔═╡ bc023fc2-5ebf-11eb-39ad-795365f8b162
begin
	
	pplt.close(); f1,a1 = pplt.subplots(ncols=3,nrows=2,sharey=0,axwidth=2)
	
	slabs = [
		0.1,0.141,0.2,0.316,0.5,0.707,
		1,1.41,2,3.16,5,7.07,
		10,14.1,20,31.6,50,70.7
	]
	nslab = length(slabs)
	a1[1].scatter(slabs,sstWTG_mean,c=colors[3:(nslab+2)])
	a1[2].scatter(slabs,prcpWTG_mean,c=colors[3:(nslab+2)])
	a1[3].scatter(slabs,tcwWTG_mean,c=colors[3:(nslab+2)])
	a1[4].scatter(slabs,sstWTG_max.-sstWTG_min,c=colors[3:(nslab+2)])
	a1[5].scatter(slabs,prcpWTG_max.-prcpWTG_min,c=colors[3:(nslab+2)])
	a1[6].scatter(slabs,csfWTG,c=colors[3:(nslab+2)])
	
	a1[1].format(ylim=(295,305))
	a1[2].format(ylim=(0,20))
	a1[3].format(ylim=(0,75))
	a1[4].format(ylim=(0.02,50),yscale="log")
	a1[5].format(ylim=(0,100))
	a1[6].format(ylim=(0,1))
	
	a1[1].format(ultitle="SST / K")
	a1[2].format(ultitle=L"Rain / mm day$^{-1}$")
	a1[3].format(ultitle="Precipitable Water / mm")
	a1[4].format(ultitle="SST (Max - Min) / K")
	a1[5].format(ultitle=L"Rain (Max - Min) / mm day$^{-1}$")
	a1[6].format(ultitle="Column Saturation Fraction")
	
	for ax in a1
		ax.format(xscale="log",xlabel="Slab Depth / m")
	end
	
	f1.savefig(plotsdir("DiurnalAmp_slabmeans.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("DiurnalAmp_slabmeans.png"))
	
end

# ╔═╡ Cell order:
# ╟─420e093c-5ba4-11eb-07da-c9a80044c8f1
# ╟─24abbe3a-5bb7-11eb-160b-1323efad463b
# ╟─24601ef8-5bb7-11eb-1cd8-198dac960d3a
# ╟─57f52568-5bb9-11eb-1e7f-c34b6efe0bac
# ╟─c9859e5e-8feb-11eb-0190-79f138ffc1ab
# ╟─4f40bd7e-5bb9-11eb-34f2-91e1f959c59a
# ╟─1ab3ad70-5c2a-11eb-187b-75e524a4581f
# ╟─6e1170f4-5bb9-11eb-0d38-61befdd2ad88
# ╟─d13eeaf8-5e8e-11eb-006e-2d45aa9253fd
# ╟─1a31f88e-5c2a-11eb-0c1f-257863a43cf5
# ╟─7d401cd6-5c2e-11eb-3842-917545e546ef
# ╟─d7735d46-5e9a-11eb-0fb3-91be12c3508b
# ╟─7b7a5806-8f20-11eb-1e95-3d1834134375
# ╟─bc023fc2-5ebf-11eb-39ad-795365f8b162
