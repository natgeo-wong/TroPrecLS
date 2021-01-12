### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 8102297a-50b7-11eb-0430-f79371a66174
begin
	using DrWatson
	using Pkg
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 9340fa4e-50b4-11eb-253e-ab01deb80456
begin
	@quickactivate "TroPrecLS"
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	include(srcdir("plotSAM.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ addc35d6-50b3-11eb-02dc-452ced2a45ef
md"
# X. Weak-Temperature Gradient Sensitivity Tests in SAM

In this notebook, we do sanity checks on the Weak-Temperature Gradient module in SAM by varying the momentum-damping parameter.  As the momentum-damping rate increases towards infinity, we should be see the model start to look more and more like the RCE state.
"

# ╔═╡ de44e1dc-542d-11eb-33f0-19a085e56922
rceexp = "3P"

# ╔═╡ b9f684ce-5163-11eb-3e19-3313ef47f0f8
begin
	ndy = 50
	z,p,t = retrievedims("Control","$(rceexp)RCE")
	td,tstep,tshift,beg = t2d(t,ndy);
	np = length(p)
end

# ╔═╡ fdfe980a-5163-11eb-3bb1-15a9cfd428e3
md"As an initial setup, we extract the dimensions of the data, and the timeshift that allows us to plot a diurnal cycle based on the last $(ndy) days of statistical output."

# ╔═╡ 19174960-50cf-11eb-12a3-cf977e483262
md"
### 1. Control RCE State

The control RCE state we used for the senstivity test was the `sst301d70` experiment, where the sea-surface temperature was fixed to 301.7ºC, which we have taken to be most representative of the Maritime Continent.

We first begin by plotting the diurnal cycle of precipitation, which is our main variable of interest.
"

# ╔═╡ 235ddba6-5164-11eb-3f62-4b983615f2a4
begin
	prcp_RCE = retrievevar("PREC","Control","$(rceexp)RCE")
	prcp_RCE_d = diurnal2D(prcp_RCE[(end-beg):end],tstep,tshift);
	
	pplt.close(); frce1,arce1 = pplt.subplots(axwidth=3,aspect=2)
	plot2Ddiurnal(arce1,1,td,prcp_RCE_d,subtractm=false)
	arce1[1].format(
		ylim=(0,6),ylabel=L"Rainfall Rate / mm day$^{-1}$",
		xlabel="Hour of Day"
	)
	frce1.savefig(plotsdir("rceprcp_$(rceexp).png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("rceprcp_$(rceexp).png"))
end

# ╔═╡ 891e5992-50cf-11eb-208e-71e78380caf7
begin
	tair_RCE = retrievevar("TABS","Control","$(rceexp)RCE")
	cldf_RCE = retrievevar("CLD","Control","$(rceexp)RCE") * 100
	qvap_RCE = retrievevar("QV","Control","$(rceexp)RCE") / 10
	rhum_RCE = calcrh(qvap_RCE,tair_RCE,p)
	
	tair_RCE_d = diurnal3D(tair_RCE[:,(end-beg):end],tstep,tshift);
	cldf_RCE_d = diurnal3D(cldf_RCE[:,(end-beg):end],tstep,tshift);
	rhum_RCE_d = diurnal3D(rhum_RCE[:,(end-beg):end],tstep,tshift);
	
	arr = [[1,2,2,2,2],[3,4,4,4,4],[5,6,6,6,6]]; pplt.close()
	frce2,arce2 = pplt.subplots(arr,aspect=1/3,axwidth=1/2,wspace=0.15,sharex=1)
	
	plot3Ddiurnal(arce2,1,td,p,tair_RCE_d,lvl=vcat(-5:-1,1:5)/20,cmapname="RdBu_r")
	plot3Ddiurnal(arce2,3,td,p,cldf_RCE_d,lvl=vcat(-5:-1,1:5)*5,cmapname="drywet")
	plot3Ddiurnal(arce2,5,td,p,rhum_RCE_d,lvl=vcat(-5:-1,1:5)/2,cmapname="drywet")
	
	arce2[1].format(xlim=(180,310))
	arce2[3].format(xlim=(0,70))
	arce2[5].format(xlim=(0,120))
	
	frce2.savefig(plotsdir("rce3D_$(rceexp).png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("rce3D_$(rceexp).png"))
end

# ╔═╡ 61ba1e22-5124-11eb-170e-7121c570d476
md"
### List of Configurations
"

# ╔═╡ 785f32a8-5371-11eb-3b1c-cb3d4c9f505c
exper = "$(rceexp)WTGamExp0"

# ╔═╡ f50703e2-5369-11eb-1c4f-a1bd28acce48
configs = [
	"damping02d00","damping08d00",
	"damping16d00",
	"damping32d00","damping64d00",
	"damping128d0","damping256d0",
	"damping512d0","damping2048d",
]

# ╔═╡ 1c9964be-536e-11eb-3da4-77b94fb34bf0
md"### Comparison of WTG Simulations against RCE"

# ╔═╡ fcc8f03a-516c-11eb-21f4-af40a2abee4b
function retrieveprcp(experiment,configlist,beg,tstep,tshift)
	
	nconfig = length(configlist);
	prcp_d  = zeros(tstep+2,nconfig)
	
	for iconfig = 1 : nconfig
		prcp = retrievevar("PREC",experiment,configlist[iconfig])
		prcp_d[:,iconfig] = diurnal2D(prcp[(end-beg):end],tstep,tshift);
	end
	
	return prcp_d
	
end

# ╔═╡ 86f914a0-536c-11eb-0563-b1e67e8c85fc
function plotprcpWTG(axs,axsii,td,prcp,configlist)
	
	nconfig = length(configlist)
	
	for iconfig = 1 : nconfig
		config = configlist[iconfig]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		axs[axsii].plot(td,prcp[:,iconfig],c="blue$iconfig",label=config,legend="r")
	end
	
end

# ╔═╡ 3d052842-5124-11eb-16a4-fd5f78742794
begin
	prcp_WTG = retrieveprcp(exper,configs,beg,tstep,tshift)
	
	lgd = Dict("ncols"=>1,"frame"=>false)
	pplt.close(); fcmp1,acmp1 = pplt.subplots(axwidth=4,aspect=2)
	acmp1[1].plot(td,prcp_RCE_d,c="k",label="RCE",legend="r",legend_kw=lgd)
	plotprcpWTG(acmp1,1,td,prcp_WTG,configs)
	acmp1[1].format(
		ylim=(0,6),ylabel=L"Rainfall Rate / mm day$^{-1}$",
		xlim=(0,24),xlabel="Hour of Day",
		suptitle="WTG - RCE | $exper"
	)
	fcmp1.savefig(plotsdir("prcp_$exper.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("prcp_$exper.png"))
end

# ╔═╡ 2cc1a5d4-5370-11eb-20fb-5d988fc02718
function retrieveWTG3D(variable,experiment,configlist,beg,tstep,tshift,np)
	
	nconfig = length(configlist);
	var_d = zeros(np,nconfig)
	
	for iconfig = 1 : nconfig
		var = retrievevar(variable,experiment,configlist[iconfig])
		var_d[:,iconfig] = dropdims(mean(var[:,(end-beg):end],dims=2),dims=2)
	end
	
	return var_d
	
end

# ╔═╡ 0976bfa4-5372-11eb-367e-a7b8374c3604
function retrieveWTGrhum(experiment,configlist,beg,tstep,tshift,p)
	
	nconfig = length(configlist);
	RH_d = zeros(length(p),nconfig)
	
	for iconfig = 1 : nconfig
		TA = retrievevar("TABS",experiment,configlist[iconfig])
		QV = retrievevar("QV",experiment,configlist[iconfig]) / 10
		RH = calcrh(QV,TA,p)
		RH = diurnal3D(RH[:,(end-beg):end],tstep,tshift);
		RH_d[:,iconfig] = dropdims(mean(RH,dims=2),dims=2)
	end
	
	return RH_d
	
end

# ╔═╡ 64a1f1cc-5373-11eb-27ca-538683b9b9e9
function plotWTG3D(axs,axsii,p,data,configlist;islegend=false)
	
	nconfig = length(configlist)
	lgd = Dict("ncols"=>4,"frame"=>false)
	
	for icon = 1 : nconfig
		config = configlist[icon]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		if islegend
			  axs[axsii].plot(
				data[:,icon],p,c="blue$icon",
				label=config,legend="b",legend_kw=lgd
			  )
		else; axs[axsii].plot(data[:,icon],p,c="blue$icon")
		end
	end
	
end

# ╔═╡ 97e914de-5373-11eb-0495-77717bf39bd1
begin
	tair_CON = dropdims(mean(tair_RCE[:,(end-beg):end],dims=2),dims=2)
	tair_WTG = retrieveWTG3D("TABS",exper,configs,beg,tstep,tshift,np) .- tair_CON
	
	cldf_CON = dropdims(mean(cldf_RCE[:,(end-beg):end],dims=2),dims=2)
	cldf_WTG = retrieveWTG3D("CLD",exper,configs,beg,tstep,tshift,np)*100 .- cldf_CON
	
	rhum_CON = dropdims(mean(rhum_RCE[:,(end-beg):end],dims=2),dims=2)
	rhum_WTG = retrieveWTGrhum(exper,configs,beg,tstep,tshift,p) .- rhum_CON
	
	wtgw_WTG = retrieveWTG3D("WWTG",exper,configs,beg,tstep,tshift,np) * 3600
	dqdt_WTG = retrieveWTG3D("QVTEND",exper,configs,beg,tstep,tshift,np)
	dTdt_WTG = retrieveWTG3D("TTEND",exper,configs,beg,tstep,tshift,np)

	pplt.close()
	fwtg,awtg = pplt.subplots(
		nrows=2,ncols=3,wspace=0.15,
		axwidth=1,aspect=0.5,sharex=0,
	)
	
	plotWTG3D(awtg,1,p,tair_WTG,configs)
	awtg[1].format(
		ylim=(1010,10),xlim=(-15,15),ylabel="Pressure / hPa",
		ultitle="T / K",suptitle="WTG - RCE | $exper"
	)
	
	plotWTG3D(awtg,2,p,cldf_WTG,configs)
	awtg[2].format(xlim=(-30,30),ultitle="CLD / %")
	
	plotWTG3D(awtg,3,p,rhum_WTG,configs)
	awtg[3].format(xlim=(-100,100),ultitle="RH / %")
	
	plotWTG3D(awtg,4,p,wtgw_WTG,configs)
	awtg[4].format(ylim=(1010,10),xlim=(-75,75),ultitle=L"W$_{WTG}$ / m hr$^{-1}$")
	
	plotWTG3D(awtg,5,p,dqdt_WTG,configs,islegend=true)
	awtg[5].format(xlim=(-5,5),ultitle=L"$\frac{dq}{dt}$ / g kg$^{-1}$ s$^{-1}$")
	
	plotWTG3D(awtg,6,p,dTdt_WTG,configs)
	awtg[6].format(xlim=(-5,5),ultitle=L"$\frac{dT}{dt}$ / K s$^{-1}$")
	
	fwtg.savefig(plotsdir("var3D_$exper.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("var3D_$exper.png"))
	
end

# ╔═╡ Cell order:
# ╟─addc35d6-50b3-11eb-02dc-452ced2a45ef
# ╟─8102297a-50b7-11eb-0430-f79371a66174
# ╠═9340fa4e-50b4-11eb-253e-ab01deb80456
# ╟─fdfe980a-5163-11eb-3bb1-15a9cfd428e3
# ╠═b9f684ce-5163-11eb-3e19-3313ef47f0f8
# ╠═de44e1dc-542d-11eb-33f0-19a085e56922
# ╟─19174960-50cf-11eb-12a3-cf977e483262
# ╠═235ddba6-5164-11eb-3f62-4b983615f2a4
# ╟─891e5992-50cf-11eb-208e-71e78380caf7
# ╟─61ba1e22-5124-11eb-170e-7121c570d476
# ╠═785f32a8-5371-11eb-3b1c-cb3d4c9f505c
# ╠═f50703e2-5369-11eb-1c4f-a1bd28acce48
# ╟─1c9964be-536e-11eb-3da4-77b94fb34bf0
# ╟─fcc8f03a-516c-11eb-21f4-af40a2abee4b
# ╟─86f914a0-536c-11eb-0563-b1e67e8c85fc
# ╟─3d052842-5124-11eb-16a4-fd5f78742794
# ╟─2cc1a5d4-5370-11eb-20fb-5d988fc02718
# ╟─0976bfa4-5372-11eb-367e-a7b8374c3604
# ╟─64a1f1cc-5373-11eb-27ca-538683b9b9e9
# ╟─97e914de-5373-11eb-0495-77717bf39bd1
