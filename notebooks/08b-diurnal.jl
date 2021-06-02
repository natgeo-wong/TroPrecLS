### A Pluto.jl notebook ###
# v0.14.5

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
# 8b. The Diurnal Cycle in SAM

In this notebook, we do some basic investigation into the relationship between slab depth and the amplitude of the diurnal cycle when the model is run in RCE state.

Relevant model parameters:
* SST = 301.7 K
* Insolation Peak = 1354.23 W m$^{-2}$
* Momentum Damping $a_m$ = 2
* Momentum Damping $a_m$ exponent = 0
"

# ╔═╡ 00872bde-94ce-11eb-3ad0-3dafcb957c19
md"
### A. Exploratory Analysis of the Diurnal Cycle

We first load the 2D data of surface temperature and precipitation, bin it according to the diurnal cycle, and then plot the average diurnal cycle of surface temperature and precipitation for some preliminary analysis.
"

# ╔═╡ 57f52568-5bb9-11eb-1e7f-c34b6efe0bac
function retrieve2D(variable,configlist,beg,tstep,tshift)
	
	nconfig = length(configlist);
	d2D  = zeros(48,nconfig)
	
	for ic = 1 : nconfig
		for ii = 1 : 5
			fnc = outstatname("DiAmp064km",configlist[ic],false,true,ii)
			var = retrievevar(variable,fnc)[9601:end]
			d2D[:,ic] += dropdims(mean(reshape(var,48,:),dims=2),dims=2)
		end
	end
	
	return d2D / 5
	
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
				legend="b",legend_kw=Dict("frame"=>false,"ncols"=>6)
			)
		else
			axs[axsii].plot(td,prcp[:,iconfig],color=colors[iconfig+2],lw=1)
		end
		
	end
	
end

# ╔═╡ 1ab3ad70-5c2a-11eb-187b-75e524a4581f
begin
	configs = [
		"Slab00d05","Slab00d07",
		"Slab00d1","Slab00d14","Slab00d2","Slab00d3","Slab00d5",
		"Slab00d7","Slab01d0","Slab01d4","Slab02d0","Slab03d2",
		"Slab05d0","Slab07d1","Slab10d0","Slab14d1","Slab20d0",
		"Slab31d6",#"Slab50d0","Slab70d7"
	]
	ncon = length(configs)
	colors = pplt.Colors("blues",ncon+4)
	lndocn = pplt.Colors("Delta_r",ncon+4)
	md"Defining experimental configurations for WTG experiments ..."
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

# ╔═╡ 1de751c6-94cf-11eb-045c-97560536c7ee
md"
### B. Basic Statistics and comparison to Observations

Now that we've done the exploratory analysis above, we do some basic statistics on the domain-averaged output and compare them to observational satellite and reanalysis data over different tropical regions.
"

# ╔═╡ 7b7a5806-8f20-11eb-1e95-3d1834134375
function retrievecsf(configlist,beg,tstep,tshift)
	
	nconfig = length(configlist);
	var  = zeros(nconfig)
	
	for ic = 1 : nconfig
		for ii = 1 : 5
			fnc = outstatname("DiAmp064km",configlist[ic],false,true,ii)
			tcw = retrievevar("PW",fnc) / 1000
			ta  = retrievevar("TABS",fnc)
			qv  = retrievevar("QV",fnc) / 1000
			pp  = retrievevar("p",fnc)
			rh  = calcrh(qv,ta,pp)
			_,swp = calccsf(rh,qv,pp)
			var[ic] += mean(tcw[(end-beg):end]./swp[(end-beg):end])
		end
	end
	
	var = var / 5
	
	return var
	
end

# ╔═╡ 1a31f88e-5c2a-11eb-0c1f-257863a43cf5
begin
	_W,_W,tWTG = retrievedims("DiAmp064km","Slab01d0")
	tdWTG,tstepWTG,tshiftWTG,begWTG = t2d(tWTG,ndy);
	prcpWTG = retrieve2D("PREC",configs,begWTG,tstepWTG,tshiftWTG)
	sstWTG  = retrieve2D("SST",configs,begWTG,tstepWTG,tshiftWTG)
	tcwWTG  = retrieve2D("PW",configs,begWTG,tstepWTG,tshiftWTG)
	capeWTG = retrieve2D("CAPE",configs,begWTG,tstepWTG,tshiftWTG)
	cinWTG  = retrieve2D("CIN",configs,begWTG,tstepWTG,tshiftWTG)
	sshfWTG = retrieve2D("SHF",configs,begWTG,tstepWTG,tshiftWTG)
	slhfWTG = retrieve2D("LHF",configs,begWTG,tstepWTG,tshiftWTG)
	cldWTG  = retrieve2D("CLDSHD",configs,begWTG,tstepWTG,tshiftWTG)
	csfWTG  = retrievecsf(configs,begWTG,tstepWTG,tshiftWTG)
md"Loading results from the WTG Slab-depth experiments ..."
end

# ╔═╡ 7d401cd6-5c2e-11eb-3842-917545e546ef
begin
	
	pplt.close(); f,axs = pplt.subplots(ncols=3,nrows=3,axwidth=2.5,aspect=1.5,sharey=0)
	
	plot2DWTG(axs,1,-11.75:0.5:12,sstWTG,configs,colors)
	plot2DWTG(axs,2,-11.75:0.5:12,prcpWTG/24,configs,colors)
	plot2DWTG(axs,3,-11.75:0.5:12,cldWTG*100,configs,colors)
	plot2DWTG(axs,4,-11.75:0.5:12,capeWTG,configs,colors)
	plot2DWTG(axs,5,-11.75:0.5:12,cinWTG,configs,colors)
	plot2DWTG(axs,6,-11.75:0.5:12,capeWTG.-cinWTG,configs,colors)
	plot2DWTG(axs,7,-11.75:0.5:12,sshfWTG,configs,colors)
	plot2DWTG(axs,8,-11.75:0.5:12,slhfWTG,configs,colors,islegend=true)
	plot2DWTG(axs,9,-11.75:0.5:12,sshfWTG.+slhfWTG,configs,colors)
	
	axs[1].format(ylim=(290,315),ultitle="(a) SST / K")
	axs[2].format(ylim=(0,3.5),ultitle=L"(b) Rainfall Rate / mm hr$^{-1}$")
	axs[3].format(ylim=(0,100),ultitle="(c) Cloud Cover / %")
	axs[4].format(ylim=(-5,35).*100,ultitle=L"(d) CAPE / J kg$^{-1}$")
	axs[5].format(ylim=(0,200),ultitle=L"(e) CIN / J kg$^{-1}$",yscale="symlog",yscale_kw=Dict("linthresh"=>0.1))
	axs[6].format(ylim=(-5,35).*100,ultitle=L"(f) CAPE - CIN / J kg$^{-1}$")
	axs[7].format(ylim=(0,200),ultitle=L"(g) Sensible Heat Flux / W m$^{-2}$")
	axs[8].format(ylim=(0,800),ultitle=L"(h) Latent Heat Flux / W m$^{-2}$")
	axs[9].format(ylim=(0,800),ultitle=L"(i) Surface Fluxes / W m$^{-2}$")
	
	for ax in axs
		ax.format(xlim=(-12,12),xlabel="Hour of Day",xlocator=-24:3:24)
	end
	
	f.savefig(plotsdir("DiurnalAmp.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("DiurnalAmp.png"))
	
end

# ╔═╡ d7735d46-5e9a-11eb-0fb3-91be12c3508b
begin
	sstWTG_mean  = dropdims(mean(sstWTG,dims=1),dims=1)
	sstWTG_max   = dropdims(maximum(sstWTG,dims=1),dims=1)
	sstWTG_min   = dropdims(minimum(sstWTG,dims=1),dims=1)
	sstWTG_dicy  = (sstWTG_max.-sstWTG_min)./2
	prcpWTG_mean = dropdims(mean(prcpWTG,dims=1),dims=1)
	prcpWTG_max  = dropdims(maximum(prcpWTG,dims=1),dims=1)
	prcpWTG_min  = dropdims(minimum(prcpWTG,dims=1),dims=1)
	prcpWTG_dicy  = (prcpWTG_max.-prcpWTG_min)./2
	tcwWTG_mean  = dropdims(mean(tcwWTG,dims=1),dims=1)
	md"Doing some basic statistics ..."
end

# ╔═╡ bc023fc2-5ebf-11eb-39ad-795365f8b162
begin
	
	pplt.close(); f1,a1 = pplt.subplots(ncols=3,nrows=2,sharey=0,axwidth=2)
	
	slabs = [
		0.05,0.0707,
		0.1,0.141,0.2,0.316,0.5,0.707,
		1,1.41,2,3.16,5,7.07,
		10,14.1,20,31.6
	]
	nslab = length(slabs)
	a1[1].scatter(slabs,sstWTG_mean,s=20)
	a1[1].fill_between([0,100],298,298.5,color="g",alpha=0.3)
	a1[1].fill_between([0,100],300.5,302.5,color="b",alpha=0.3)
	a1[1].plot([0,100],[1,1]*301.7,color="k")
	a1[1].plot([0,100],[1,1]*301.545,color="b",linestyle="--")
	a1[1].plot([0,100],[1,1]*300.683,color="k",linestyle="--")
	a1[1].plot([0,100],[1,1]*300.382,color="blue3",linestyle="--")
	a1[1].plot([0,100],[1,1]*298.169,color="b",linestyle="-.")
	a1[1].plot([0,100],[1,1]*298.363,color="blue3",linestyle="-.")
	
	a1[2].scatter(slabs,prcpWTG_mean/24,s=20)
	a1[2].plot([0,100],[1,1]*0.126,color="k")
	a1[2].plot([0,100],[1,1]*0.284,color="b",linestyle="--")
	a1[2].plot([0,100],[1,1]*0.213,color="k",linestyle="--")
	a1[2].plot([0,100],[1,1]*0.133,color="blue3",linestyle="--")
	a1[2].plot([0,100],[1,1]*0.297,color="b",linestyle="-.")
	a1[2].plot([0,100],[1,1]*0.183,color="blue3",linestyle="-.")
	
	a1[3].scatter(slabs,tcwWTG_mean,s=20)
	
	a1[4].scatter(slabs,sstWTG_dicy,s=20)
	a1[4].fill_between([0,100],4,5,color="g",alpha=0.3)
	a1[4].fill_between([0,100],0.15,0.25,color="b",alpha=0.3)
	a1[4].plot([0,100],[1,1]*0.218,color="b",linestyle="--")
	a1[4].plot([0,100],[1,1]*0.170,color="k",linestyle="--")
	a1[4].plot([0,100],[1,1]*0.160,color="blue3",linestyle="--")
	a1[4].plot([0,100],[1,1]*4.461,color="b",linestyle="-.")
	a1[4].plot([0,100],[1,1]*4.659,color="blue3",linestyle="-.")
	
	a1[5].scatter(slabs,prcpWTG_dicy./prcpWTG_mean,s=20)
	a1[5].fill_between([0,100],0.3,2,color="g",alpha=0.3)
	a1[5].fill_between([0,100],0.1,1,color="b",alpha=0.3)
	a1[5].plot([0,100],[1,1]*0.358,color="b",linestyle="--")
	a1[5].plot([0,100],[1,1]*0.397,color="k",linestyle="--")
	a1[5].plot([0,100],[1,1]*0.350,color="blue3",linestyle="--")
	a1[5].plot([0,100],[1,1]*0.889,color="b",linestyle="-.")
	a1[5].plot([0,100],[1,1]*1.240,color="blue3",linestyle="-.")
	
	a1[6].scatter(slabs,csfWTG,s=20)
	
	a1[1].format(ylim=(295,305))
	a1[2].format(ylim=(0,0.8))
	a1[3].format(ylim=(0,75))
	a1[4].format(ylim=(0.01,50),yscale="log")
	a1[5].format(ylim=(0.01,3))
	a1[6].format(ylim=(0,1))
	
	a1[1].format(ultitle="(a) SST / K")
	a1[2].format(ultitle=L"(b) Rain / mm hr$^{-1}$")
	a1[3].format(ultitle="(c) Precipitable Water / mm")
	a1[4].format(ultitle="(d) SST (Diurnal) / K")
	a1[5].format(ultitle="(e) Rain (Diurnal) / Mean Rain")
	a1[6].format(ultitle="(f) Column Saturation Fraction")
	
	for ax in a1
		ax.format(xscale="log",xlabel="Slab Depth / m",xlim=(0.05,100))
	end
	
	f1.savefig(plotsdir("DiurnalAmp_slabmeans.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("DiurnalAmp_slabmeans.png"))
	
end

# ╔═╡ 5e4ebf22-94e0-11eb-3c58-1ff954f15556
md"We see that the diurnal variability in SST for our WTG experiments roughly follows an inverse relationship with slab-depth.  However, this relationship begins to break down at a slab depth of around 0.2m, where we see that the diurnal variability no longer shows significant increase as slab depth decreases.

It also seems that at <0.2 mixed-layer ocean depth, the model transitions into a regime where the surface temperature is cooler that expected.  In fact, the surface temperature for such regimes coincides very closely with the average temperature observed over land in the tropical regions."

# ╔═╡ 6724ddbc-94da-11eb-2610-39d5c0ca6096
md"
### C. Different Diurnal Regimes

However, in contrast to what the statistical averages of the diurnal cycle may imply, there are several different regimes that the model fluctuates between over the course of the simulation.  This is best seen when we do time-series plots.
"

# ╔═╡ 6ddbb5c8-94f2-11eb-2198-41cc2dc26102
lcon = "Slab10d0"

# ╔═╡ e13e3d69-54bf-48b4-819b-5ab45dcb05bd
begin
	wtglvls = 10. .^(-2:0.5:0)
	wtglvls = vcat(reverse(-wtglvls),wtglvls)
	md"Defining WTG plotting levels"
end

# ╔═╡ cf060c26-94da-11eb-3899-818d2bc6d9d2
begin
	_,p01d0,ts01d0 = retrievedims("DiAmp064km",lcon)
	prcp_ts01d0 = retrievevar("PREC","DiAmp064km",lcon) / 24
	sst_ts01d0  = retrievevar("SST","DiAmp064km",lcon)
	tbo_ts01d0  = retrievevar("TBIAS","DiAmp064km",lcon)
	wtg_ts01d0  = retrievevar("WWTG","DiAmp064km",lcon)
	tab_ts01d0  = retrievevar("TABS","DiAmp064km",lcon)
	qv_ts01d0   = retrievevar("QV","DiAmp064km",lcon) / 1000
	relh_ts01d0 = calcrh(qv_ts01d0,tab_ts01d0,p01d0) * 100
	cld_ts01d0  = retrievevar("CLD","DiAmp064km",lcon) * 100
	ts01d0 = ts01d0 .- 80
	tsend = floor(ts01d0[end])
md"Loading results from the $(lcon) experiment ..."
end

# ╔═╡ fa25ce46-94da-11eb-36b0-09c7ea40fd33
begin
	
	pplt.close()
	fts1,ats1 = pplt.subplots(nrows=3,ncols=2,sharey=0,aspect=2.5,axwidth=3)
	
	ats1[1].plot(ts01d0,prcp_ts01d0)
	ats1[1].format(ylim=(0,12.5),ltitle=L"(a) Rain / mm hr$^{-1}$")
	
	ats1[2].plot(ts01d0,sst_ts01d0)
	ats1[2].format(ltitle="(b) SST / K")
	
	cts1 = ats1[3].contourf(
		ts01d0,p01d0,tbo_ts01d0,
		cmap="RdBu_r",levels=wtglvls*10,extend="both"
	)
	ats1[3].colorbar(cts1,loc="b",ticks=[-1,-0.1,-0.01,0.01,0.1,1]*10)
	ats1[3].format(ltitle="(c) T - TOBS / K")
	
	cts1 = ats1[4].contourf(
		ts01d0,p01d0,wtg_ts01d0,
		cmap="RdBu_r",levels=wtglvls,extend="both"
	)
	ats1[4].colorbar(cts1,loc="b",ticks=[-1,-0.1,-0.01,0.01,0.1,1])
	ats1[4].format(ltitle=L"(d) w$_{WTG}$ / m s$^{-1}$")
	
	cts1 = ats1[5].contourf(
		ts01d0,p01d0,relh_ts01d0,
		cmap="Blues",levels=10:10:100,extend="both"
	)
	ats1[5].colorbar(cts1,loc="b",ticks=10:10:100)
	ats1[5].format(ltitle="(e) Relative Humidity / %")
	
	cts1 = ats1[6].contourf(
		ts01d0,p01d0,cld_ts01d0,
		cmap="Blues",levels=10:10:90,extend="both"
	)
	ats1[6].colorbar(cts1,loc="b",ticks=10:10:90)
	ats1[6].format(ltitle="(f) Cloud Fraction / %")
	
	for ax in ats1
		ax.format(xlim=(tsend-25,tsend))#,xlocator=250:5:300)
	end
	
	for iax in 3:6
		ats1[iax].format(yscale="log")
	end
	
	fts1.savefig(plotsdir("DiurnalAmp_ts$(lcon).png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("DiurnalAmp_ts$(lcon).png"))
	
end

# ╔═╡ 49e3106a-94e0-11eb-3495-7378a30ce4d7
md"We see that for the mixed layer depths that are of 0.2 m or thicker, that the model fluctuates between wet and dry regimes.  There are days when almost no days occur, but there are days when precipitation is extremely heavy.  When daily precipitation is high, then precipitation mostly occurs near midnight.  However, when precipitation is low, then precipitation occurs earlier in the day."

# ╔═╡ a7a9b46c-94e3-11eb-2959-bb8d0aa67644
begin
	_1,p00d1,ts00d1 = retrievedims("DiAmp064km","Slab00d1")
	prcp_ts00d1 = retrievevar("PREC","DiAmp064km","Slab00d1") / 24
	sst_ts00d1  = retrievevar("SST","DiAmp064km","Slab00d1")
	tbo_ts00d1  = retrievevar("TBIAS","DiAmp064km","Slab00d1")
	wtg_ts00d1  = retrievevar("WWTG","DiAmp064km","Slab00d1")
	tab_ts00d1  = retrievevar("TABS","DiAmp064km","Slab00d1")
	qv_ts00d1   = retrievevar("QV","DiAmp064km","Slab00d1") / 1000
	relh_ts00d1 = calcrh(qv_ts00d1,tab_ts00d1,p00d1) * 100
	cld_ts00d1  = retrievevar("CLD","DiAmp064km","Slab00d1") * 100
	ts00d1 = ts00d1 .- 80
md"Loading results from the Slab00d1 experiment ..."
end

# ╔═╡ e07d0e30-94e3-11eb-36b6-69db4d2ba24f
begin
	
	pplt.close()
	fts2,ats2 = pplt.subplots(nrows=3,ncols=2,sharey=0,aspect=2.5,axwidth=3)
	
	ats2[1].plot(ts00d1,prcp_ts00d1)
	ats2[1].format(ylim=(0,12.5),ltitle=L"(a) Rain / mm hr$^{-1}$")
	
	ats2[2].plot(ts00d1,sst_ts00d1)
	ats2[2].format(ltitle="(b) SST / K")
	
	cts2 = ats2[3].contourf(
		ts00d1,p00d1,tbo_ts00d1,
		cmap="RdBu_r",levels=wtglvls*10,extend="both"
	)
	ats2[3].colorbar(cts2,loc="b",ticks=[-1,-0.1,-0.01,0.01,0.1,1]*10)
	ats2[3].format(ltitle="(c) T - TOBS / K")
	
	cts2 = ats2[4].contourf(
		ts00d1,p00d1,wtg_ts00d1,
		cmap="RdBu_r",levels=wtglvls,extend="both"
	)
	ats2[4].colorbar(cts2,loc="b",ticks=[-1,-0.1,-0.01,0.01,0.1,1])
	ats2[4].format(ltitle=L"(d) w$_{WTG}$ / m s$^{-1}$")
	
	cts2 = ats2[5].contourf(
		ts00d1,p00d1,relh_ts00d1,
		cmap="Blues",levels=10:10:100,extend="both"
	)
	ats2[5].colorbar(cts2,loc="b",ticks=10:10:100)
	ats2[5].format(ltitle="(e) Relative Humidity / %")
	
	cts2 = ats2[6].contourf(
		ts00d1,p00d1,cld_ts00d1,
		cmap="Blues",levels=10:10:90,extend="both"
	)
	ats2[6].colorbar(cts2,loc="b",ticks=10:10:90)
	ats2[6].format(ltitle="(f) Cloud Fraction / %")
	
	for ax in ats2
		ax.format(xlim=(275,300),xlocator=250:5:300)
	end
	
	for iax in 3:6
		ats2[iax].format(yscale="log")
	end
	
	fts2.savefig(plotsdir("DiurnalAmp_tsSlab00d1.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("DiurnalAmp_tsSlab00d1.png"))
	
end

# ╔═╡ f4519a32-94fc-11eb-1401-c3f85279a602
md"However, when the mixed layer depth is thinner than 0.2 m, then the model also fluctuates between dry and wet states, but the dry states have significant daily precipitation compared to model runs when the mixed layer depth is thicker than 0.2 m.  It is also of note that wet states occur every 2 days, whereas for thicker ocean mixed-layer depths, it takes some days for dry regimes to build up enough moisture to cause a very wet state."

# ╔═╡ 3e6b4913-6ff7-4298-86e4-dbe49367d584
md"
### D. Hour of Maximum Rainfall
"

# ╔═╡ 9bb67a0d-8d35-449a-b7a6-d4c5b01cac64
function retrievemaxtime(variable,configlist)
	
	nconfig = length(configlist);
	d2D  = zeros(48,nconfig)
	vtmp = zeros(48)
	
	for ic = 1 : nconfig
		for ii = 1 : 5
			fnc = outstatname("DiAmp064km",configlist[ic],false,true,ii)
			var = retrievevar(variable,fnc)[9601:end]
			var = reshape(var,48,:); ndy = size(var,2)
			for idy = 1 : ndy
				vtmp .= var[:,idy]; imax = argmax(vtmp)
				if vtmp[imax] > 0.25*24
					d2D[imax,ic] += 1
				end
			end
		end
		d2D[:,ic] .= d2D[:,ic] ./ sum(d2D[:,ic]) * 48
	end
	
	d2D = 0.5*d2D .+ 0.25 *(circshift(d2D,(1,0)) .+ circshift(d2D,(-1,0)))
	
	return d2D
	
end

# ╔═╡ c9acf590-7c4d-4691-ba10-0ec455acd3e0
begin
	prcphr = retrievemaxtime("PREC",configs)
	ssthr = retrievemaxtime("SST",configs)
	md"Binning the hour of maximum rainfall"
end

# ╔═╡ b99c20a5-5f88-49c5-94aa-e2e1bedebb5c
begin
	pplt.close(); fθ,aθ = pplt.subplots(ncols=2,proj="polar");
	
	θvec = collect(-12:0.5:12)/24*2*pi
	θvec = (θvec[1:(end-1)].+θvec[2:end])/2
	θvec = vcat(θvec,θvec[1]+2*pi)
	
	for islab = 1 : 8
		aθ[1].plot(
			θvec,sqrt.(vcat(prcphr[:,islab],prcphr[1,islab])),
			c=lndocn[islab+2],lw=1,
			label="$(slabs[islab]) m",
			legend="r",legend_kw=Dict("frame"=>false,"ncols"=>1)
		)
	end
	for islab = 11 : 18
		aθ[2].plot(
			θvec,sqrt.(vcat(prcphr[:,islab],prcphr[1,islab])),
			c=lndocn[islab+2],lw=1,
			label="$(slabs[islab]) m",
			legend="r",legend_kw=Dict("frame"=>false,"ncols"=>1)
		)
	end
	
	aθ[1].format(ltitle="(a) Land Analogue")
	aθ[2].format(ltitle="(b) Ocean Analogue")
	
	for ax in aθ
		ax.format(
			theta0="N",thetaformatter="tau",
			#rlim=(0,4),rlocator=0:4,
			suptitle=L"$\theta$ / Fraction of Day"
		)
	end
	
	fθ.savefig(plotsdir("samdiurnalphase.png"),transparent=false,dpi=200)
	load(plotsdir("samdiurnalphase.png"))
end

# ╔═╡ Cell order:
# ╟─420e093c-5ba4-11eb-07da-c9a80044c8f1
# ╟─24abbe3a-5bb7-11eb-160b-1323efad463b
# ╟─24601ef8-5bb7-11eb-1cd8-198dac960d3a
# ╟─00872bde-94ce-11eb-3ad0-3dafcb957c19
# ╟─57f52568-5bb9-11eb-1e7f-c34b6efe0bac
# ╟─c9859e5e-8feb-11eb-0190-79f138ffc1ab
# ╠═4f40bd7e-5bb9-11eb-34f2-91e1f959c59a
# ╟─1ab3ad70-5c2a-11eb-187b-75e524a4581f
# ╟─6e1170f4-5bb9-11eb-0d38-61befdd2ad88
# ╟─d13eeaf8-5e8e-11eb-006e-2d45aa9253fd
# ╠═1a31f88e-5c2a-11eb-0c1f-257863a43cf5
# ╟─7d401cd6-5c2e-11eb-3842-917545e546ef
# ╟─1de751c6-94cf-11eb-045c-97560536c7ee
# ╟─d7735d46-5e9a-11eb-0fb3-91be12c3508b
# ╟─7b7a5806-8f20-11eb-1e95-3d1834134375
# ╟─bc023fc2-5ebf-11eb-39ad-795365f8b162
# ╟─5e4ebf22-94e0-11eb-3c58-1ff954f15556
# ╟─6724ddbc-94da-11eb-2610-39d5c0ca6096
# ╠═6ddbb5c8-94f2-11eb-2198-41cc2dc26102
# ╟─e13e3d69-54bf-48b4-819b-5ab45dcb05bd
# ╟─cf060c26-94da-11eb-3899-818d2bc6d9d2
# ╟─fa25ce46-94da-11eb-36b0-09c7ea40fd33
# ╟─49e3106a-94e0-11eb-3495-7378a30ce4d7
# ╟─a7a9b46c-94e3-11eb-2959-bb8d0aa67644
# ╟─e07d0e30-94e3-11eb-36b6-69db4d2ba24f
# ╟─f4519a32-94fc-11eb-1401-c3f85279a602
# ╟─3e6b4913-6ff7-4298-86e4-dbe49367d584
# ╟─9bb67a0d-8d35-449a-b7a6-d4c5b01cac64
# ╠═c9acf590-7c4d-4691-ba10-0ec455acd3e0
# ╠═b99c20a5-5f88-49c5-94aa-e2e1bedebb5c
