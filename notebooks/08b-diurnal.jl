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
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 420e093c-5ba4-11eb-07da-c9a80044c8f1
md"
# 7b. Amplitude of the Diurnal Cycle

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
		var = retrievevar(variable,"DiurnalAmp",configlist[iconfig])
		d2D[:,iconfig] = diurnal2D(var[(end-beg):end],tstep,tshift);
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
				td,prcp[:,iconfig],color=colors[iconfig+2],label="$(config) m",
				legend="b",legend_kw=Dict("frame"=>false,"ncols"=>5)
			)
		else
			axs[axsii].plot(td,prcp[:,iconfig],color=colors[iconfig+2])
		end
		
	end
	
end

# ╔═╡ 1ab3ad70-5c2a-11eb-187b-75e524a4581f
begin
	configs = [
		"Slab00d1",
		"Slab00d2",
		"Slab00d3","Slab00d5","Slab00d7",
		"Slab01d0","Slab01d4","Slab02d0","Slab03d2","Slab05d0",
		# "Slab07d1",
		"Slab10d0","Slab14d1",
		"Slab20d0",
		"Slab31d6",
		"Slab50d0"
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
	_R,_R,tRCE = retrievedims("Control","3SRCE")
	tdRCE,tstepRCE,tshiftRCE,begRCE = t2d(tRCE,ndy);
	prcpRCE = retrievevar("PREC","Control","3SRCE")/24
	prcpRCE = diurnal2D(prcpRCE,tstepRCE,tshiftRCE)
	sstRCE  = retrievevar("SST","Control","3SRCE")
	sstRCE  = diurnal2D(sstRCE,tstepRCE,tshiftRCE)
md"Loading results from the RCE Control ..."
end

# ╔═╡ 1a31f88e-5c2a-11eb-0c1f-257863a43cf5
begin
	_W,_W,tWTG = retrievedims("DiurnalAmp","Slab05d0")
	tdWTG,tstepWTG,tshiftWTG,begWTG = t2d(tWTG,ndy);
	prcpWTG = retrieve2D("PREC",configs,begWTG,tstepWTG,tshiftWTG)
	sstWTG  = retrieve2D("SST",configs,begWTG,tstepWTG,tshiftWTG)
md"Loading results from the WTG Slab-depth experiments ..."
end

# ╔═╡ 7d401cd6-5c2e-11eb-3842-917545e546ef
begin
	
	pplt.close(); f,axs = pplt.subplots(nrows=2,axwidth=3.5,aspect=2,sharey=0)
	
	plot2DWTG(axs,1,tdWTG,sstWTG,configs,colors)
	plot2DWTG(axs,2,tdWTG,prcpWTG/24,configs,colors,islegend=true)
	
	# axs[1].plot(tdRCE,sstRCE,c="k")
	# axs[2].plot(tdRCE,prcpRCE,c="k",label="RCE",legend="b")
	
	axs[1].format(ylim=(292,312),ylabel=L"SST / K")
	axs[2].format(ylim=(0,2.5),ylabel=L"Rainfall Rate / mm day$^{-1}$")
	
	for ax in axs
		ax.format(xlim=(0,24),xlabel="Hour of Day",xlocator=0:4:24)
	end
	
	f.savefig(plotsdir("DiurnalAmp.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("DiurnalAmp.png"))
	
end

# ╔═╡ d7735d46-5e9a-11eb-0fb3-91be12c3508b
sstWTG_mean = dropdims(mean(sstWTG,dims=1),dims=1)

# ╔═╡ bc023fc2-5ebf-11eb-39ad-795365f8b162
begin
	
	pplt.close(); f1,a1 = pplt.subplots(ncols=2)
	
	slabs = [0.1,0.2,0.316,0.5,0.707,1,1.41,2,3.16,5,7.07,10,14.1,20]
	a1[1].scatter(slabs,sstWTG_mean)
	
	for ax in a1
		ax.format(xscale="log")
	end
	
	f1.savefig(plotsdir("DiurnalAmp_slabmeans.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("DiurnalAmp_slabmeans.png"))
	
end

# ╔═╡ fb578f82-5e9a-11eb-1914-715966e13216
sstWTG_avg = mean(sstWTG_mean[1:(end-1)])

# ╔═╡ de71b57c-5e9b-11eb-374d-b1713311bf4e
sstWTG_std = std(sstWTG_mean[1:(end-1)])

# ╔═╡ Cell order:
# ╟─420e093c-5ba4-11eb-07da-c9a80044c8f1
# ╟─24abbe3a-5bb7-11eb-160b-1323efad463b
# ╟─24601ef8-5bb7-11eb-1cd8-198dac960d3a
# ╠═57f52568-5bb9-11eb-1e7f-c34b6efe0bac
# ╠═4f40bd7e-5bb9-11eb-34f2-91e1f959c59a
# ╠═1ab3ad70-5c2a-11eb-187b-75e524a4581f
# ╠═6e1170f4-5bb9-11eb-0d38-61befdd2ad88
# ╠═d13eeaf8-5e8e-11eb-006e-2d45aa9253fd
# ╠═1a31f88e-5c2a-11eb-0c1f-257863a43cf5
# ╠═7d401cd6-5c2e-11eb-3842-917545e546ef
# ╠═d7735d46-5e9a-11eb-0fb3-91be12c3508b
# ╠═bc023fc2-5ebf-11eb-39ad-795365f8b162
# ╠═fb578f82-5e9a-11eb-1914-715966e13216
# ╠═de71b57c-5e9b-11eb-374d-b1713311bf4e
