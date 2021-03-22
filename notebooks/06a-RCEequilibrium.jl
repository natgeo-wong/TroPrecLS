### A Pluto.jl notebook ###
# v0.12.21

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

# ╔═╡ 417ee688-5ade-11eb-2e95-91a301119e88
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 46faa412-5ade-11eb-3c37-23a7e59037a0
begin
	@quickactivate "TroPrecLS"
	using Statistics
	using PlutoUI
	using Printf
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	include(srcdir("samsnd.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 9dd4cd7e-5adb-11eb-2735-a7a4a2bb23b1
md"
# 6a. Finding an Equilibrium RCE State

This notebook will be used to help find an equilibrium RCE state that roughly returns the initial sounding profile used to initialize the model.

SAM reads in \"snd\" (sounding profiles of humidity and potential temperature) files to initialize the model upon spinup.  We run models in RCE mode over 500-2000 day periods and obtain the profile of absolute temperature for the last `ndy` days.  We then compare it to the `TABSOBS` that was fed into the model, and find the temperature difference.  If the difference is too large, then we extract the vertical sounding profile, override the old \"snd\" file and rerun the model.
"

# ╔═╡ b2670c08-81e5-11eb-324e-2b923b289a04
md"
### A. Definining the Model Configuration

There are two broad model configuration categories: (P)erpetual (INSOL)ation, and (D)iurnal (INSOL)ation.
"

# ╔═╡ 810d5a5a-8225-11eb-37b2-6978f37f77c3
md"Is Diurnal? $(@bind isdiurnal PlutoUI.Slider(0:1))"

# ╔═╡ 427511d0-88cb-11eb-2a40-019c91ee1401
md"Toggle Domain Size: $(@bind islarge PlutoUI.Slider(0:1))"

# ╔═╡ ac8b9d4c-5ade-11eb-06f4-33bff063bbde
begin
	if isone(isdiurnal)
		  insol = "D"
	else; insol = "P"
	end
	
	if isone(islarge)
		  domsize = "LARGE"
	else; domsize = "INSOL"
	end
	
	config = "$(insol)$(domsize)"
end

# ╔═╡ 7842e150-822b-11eb-1ded-f35ee4cc6d8c
begin
	arr = [[1,2,2,2,2],[3,4,4,4,4]]
	lvls = vcat(-5,-3.16,-2,-1.41,-1,-0.5,0.5,1,1.41,2,3.16,5)
	md"Defining universal plotting variables"
end

# ╔═╡ 86f0a60a-5ae5-11eb-0b1a-935e8703b842
md"
### B. Initial RCE spinup

Our first step is to do the first spinup for the RCE state from a initial sounding profile.  We run this initial RCE spinup for 500 days, and take the average temperature and humidity profile of the last 100 days to create a new sounding profile.

We extract the temperature and observed temperature (`TABS` and `TABSOBS` respectively) data from the statistical output.  `TABSOBS` is the large-scale temperature profile of the initial sounding profile, which is taken to be constant over the model run.  Then, we plot the difference between the two variables for the last 100 days of the RCE model run.
"

# ╔═╡ ab0909d0-5ade-11eb-2622-0fee6450004b
begin
	z,p,t = retrievedims("Control","$(config)")
	nz = length(z); nt = length(t)
	tem = retrievevar("TABS","Control","$(config)")
	tob = retrievevar("TABSOBS","Control","$(config)")[:,end]
	qvp = retrievevar("QV","Control","$(config)")
	qob = retrievevar("QVOBS","Control","$(config)")[:,end]
md"Loading data from $(config) run ..."
end

# ╔═╡ a198b70a-6d68-11eb-0c84-a1481c131fa3
w = [ 0.0198,0.0659,0.1284,0.1833,0.2052,0.1833,0.1284,0.0659,0.0198 ]

# ╔═╡ f9ed4c2a-6d68-11eb-1b10-ab05cb1adc4b
begin
	temn = zeros(nz,nt-8)
	for iz = 1 : nz, it = 5 : (nt-4)
		temn[iz,it-4] = sum(w .* @view tem[iz,(it-4):(it+4)])
	end
md"Doing some smoothing using Lanzcos filtering for future plotting ..."
end

# ╔═╡ e82cd648-5ade-11eb-2739-ad93c0b6d3af
begin
	tdiff = dropdims(mean(tem[:,(end-100+1):end],dims=2),dims=2).-tob
	tdts  = tem .- tob
	qdiff = dropdims(mean(qvp[:,(end-100+1):end],dims=2),dims=2).-qob
	qdiff = qdiff ./ qob
	qdts  = (qvp .- qob) ./ qob
md"Calculating difference between `TABS` and `TABSOBS` ..."
end

# ╔═╡ 0c8cc7ec-5ae3-11eb-1c52-0fecbb575c43
begin
	trms = sqrt(mean(tdiff.^2)); trms = @sprintf("%.3f",trms)
	qrms = sqrt(mean(qdiff.^2)); qrms = @sprintf("%.3f",qrms)
	
md"Assuming that the tropopause is at 100 hPa, the root-mean-square of the temperature difference between the model temperature `TABS` and the observed temperature `TABSOBS` is $(trms) K.  The profile of the temperature difference is shown below:"
end

# ╔═╡ 978d0442-5b33-11eb-39e5-cdea9c6efa4c
begin
	
	pplt.close(); fts,ats = pplt.subplots(arr,aspect=1/3,axwidth=0.6,sharex=0)
	
	ats[1].plot(tdiff,p)
	ats[1].scatter(tdiff,p,s=7)
	ats[1].format(
		xlim=(-7.5,7.5),xlocator=(-2:2)*5,xlabel=L"T - T$_{OBS}$ / K",
		ylim=(1010,10),yscale="log",ylabel="Pressure / hPa",
		suptitle="RCE Initial Spinup | $config"
	)
	
	ct = ats[2].contourf(
		t.-80,p,tdts,
		# t.-80,p,tem .- mean(tem[:,(end-100+1):end],dims=2),
		# t[5:(end-4)].-80,p,temn .- tob[:,1],
		cmap="RdBu_r",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",
		levels=lvls
	)
	ats[2].format(
		xlim=(00,nt),
		ylim=(1010,20),yscale="log",
		ylabel="Pressure / hPa",
		urtitle=L"T$_{RMS}$" * " = $(trms) K"
	)
	
	ats[3].plot(qdiff,p)
	ats[3].scatter(qdiff,p,s=7)
	ats[3].format(
		xlim=(-0.75,0.75),xlocator=(-2:2)/2,xlabel=L"qr = $\frac{q - q_{OBS}}{q_{OBS}}$",
		ylim=(1010,10),yscale="log",ylabel="Pressure / hPa",
	)
	
	cq = ats[4].contourf(
		t.-80,p,qdts,
		# t.-80,p,qvp .- mean(qvp[:,(end-100+1):end],dims=2),
		cmap="drywet",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",
		levels=lvls/10
	)
	ats[4].format(
		xlim=(00,nt),
		ylim=(1010,20),yscale="log",
		xlabel="time / days",ylabel="Pressure / hPa",
		urtitle=L"$qr_{RMS}$" * " = $(qrms)"# * L" g kg$^{-1}$"
	)
	
	ats[2].colorbar(ct,loc="r",width=0.2)
	ats[4].colorbar(cq,loc="r",width=0.2)
	fts.savefig("rcetdts-$(config).png",transparent=false,dpi=200)
	load("rcetdts-$(config).png")
	
end

# ╔═╡ 39868b46-5ae5-11eb-0bf5-274642855e12
md"
### C. Creating Sounding Profile for initial RCE spinup

From this initial spinup, we create the sounding profile based on the last 100 days of RCE data.  We then have to decide for ourselves whether we want to rerun the RCE model based on this new sounding data, or to move on to using this RCE profile for WTG simulations.
"

# ╔═╡ 4fd6e272-5b32-11eb-2f60-bbd4f8b9fc12
md"Create SND file? $(@bind dosnd PlutoUI.Slider(0:1))"

# ╔═╡ 094999c8-5ae5-11eb-1526-f38c604184cb
if isone(dosnd)
	createsndmean(
		"$(config)",
		exp="Control",config="$(config)",
		psfc=1009.32,ndays=100
	)
md"Creating the SND file $(config) ..."
else
md"We have decided not to create the SND file $(config) yet ..."
end

# ╔═╡ 6fe3de4c-82ec-11eb-023d-7dc6f51fe448
nen = 10

# ╔═╡ 99dca670-81e5-11eb-24b8-cf1b23d8c1f7
md"
### D. Running a $(nen)-member RCE ensemble

From this initial spinup, I ran a $(nen)-member RCE ensemble, with initial conditions perturbed as per SAM's ensemble initialization (which allows for up to 100 members).  We use a model ensemble in order to account for the chaotic nature of convection and how that might perturb the final mean state.  Each member is run for 1000-1500 days, and we take the temperature and humidity profiles of the last 250-500 days of the model run, and create a new sounding profile that overrides the old file.

Should the difference in the final averaged sounding profile be vastly different from the initial sounding profile, the model ensemble is run again, with the final averaged sounding profile as the new initial sounding profile.
"

# ╔═╡ 401a6a3a-8444-11eb-3ee8-594591ed5aa9
nendays = 100

# ╔═╡ 7d905176-81e8-11eb-20d3-b9be287472f5
begin
	z_en,_,t_en = retrievedims("Control","$(config)",isensemble=true,member=1)
	nz_en = length(z_en); nt_en = length(t_en)
	tbi_en = zeros(nz_en,nt_en,nen); tob_en = zeros(nz_en,nen)
	qbi_en = zeros(nz_en,nt_en,nen); qob_en = zeros(nz_en,nen)
	tem_en = zeros(nz_en,nt_en,nen)
	qvp_en = zeros(nz_en,nt_en,nen)
	pre_en = zeros(nz_en,nt_en,nen); plevel = zeros(nz_en,nen)
	for imem = 1 : nen
		tbi_en[:,:,imem] = retrievevar("TBIAS","Control","$(config)",isensemble=true,member=imem)
		qbi_en[:,:,imem] = retrievevar("QBIAS","Control","$(config)",isensemble=true,member=imem)
		tem_en[:,:,imem] = retrievevar("TABS","Control","$(config)",isensemble=true,member=imem)
		qvp_en[:,:,imem] = retrievevar("QV","Control","$(config)",isensemble=true,member=imem)
		pre_en[:,:,imem] = retrievevar("PRES","Control","$(config)",isensemble=true,member=imem)
		tob_en[:,imem] = retrievevar("TABSOBS","Control","$(config)",isensemble=true,member=imem)[:,end]
		qob_en[:,imem] = retrievevar("QVOBS","Control","$(config)",isensemble=true,member=imem)[:,end]
		plevel[:,imem] = retrievevar("p","Control","$(config)",isensemble=true,member=imem)[:,end]
	end
	
	tob_en = dropdims(mean(tob_en,dims=2),dims=2)
	qob_en = dropdims(mean(qob_en,dims=2),dims=2)
	plevel = dropdims(mean(plevel,dims=2),dims=2)
md"Loading data from the $(nen)-member ensemble for the $(config) RCE run ..."
end

# ╔═╡ ad523b4e-81ee-11eb-2d10-a984d5983471
begin
	tdiff_en = dropdims(mean(tbi_en[:,(end-nendays+1):end,:],dims=2),dims=2)
	tdts_en  = dropdims(mean(tbi_en,dims=3),dims=3)
	qdiff_en = dropdims(mean(qbi_en[:,(end-nendays+1):end,:],dims=2),dims=2)
	qdiff_en = qdiff_en ./ qob_en*100
	qdts_en  = dropdims(mean(qbi_en ./ qob_en*100,dims=3),dims=3)
	pts_en   = dropdims(mean(pre_en[:,(end-nendays+1):end,:],dims=2),dims=2)
md"Calculating difference between `TABS` and `TABSOBS` ..."
end

# ╔═╡ c658d650-8614-11eb-252c-c3066fb1d506
ip = plevel .> 70

# ╔═╡ f35ab14c-8226-11eb-29b4-c3b7130f3733
begin
	trms_en = sqrt(mean(tdiff_en[ip,:].^2)); trms_en = @sprintf("%.3f",trms_en)
	qrms_en = sqrt(mean(qdiff_en[ip,:].^2)); qrms_en = @sprintf("%.3f",qrms_en)
	
md"Assuming that the tropopause is at 100 hPa, the root-mean-square of the temperature difference between the model temperature `TABS` and the observed temperature `TABSOBS` is $(trms_en) K.  The profile of the temperature difference is shown below:"
end

# ╔═╡ a3650e8a-822b-11eb-29ca-cd01b90e8099
begin
	
	pplt.close(); fen,aen = pplt.subplots(arr,aspect=1/3,axwidth=0.6,sharex=0)
	
	for im = 1 : nen
		aen[1].scatter(tdiff_en[:,im],pts_en[:,im],s=2,c="gray")
	end
	
	aen[1].plot(
		dropdims(mean(tdiff_en,dims=2),dims=2),
		dropdims(mean(pts_en,dims=2),dims=2),c="k"
	)
	aen[1].format(
		xlim=(-0.15,0.15),xlocator=(-2:2)/10,xlabel=L"T - T$_{OBS}$ / K",
		ylim=(1010,25),yscale="log",ylabel="Pressure / hPa",
		suptitle="Model Ensemble Equilibrium RCE | $config"
	)
	
	cten = aen[2].contourf(
		t_en.-80,plevel,tdts_en,
		# t.-80,p,tem .- mean(tem[:,(end-100+1):end],dims=2),
		# t[5:(end-4)].-80,p,temn .- tob[:,1],
		cmap="RdBu_r",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",
		levels=lvls/10
	)
	aen[2].format(
		xlim=(0,nt_en),
		ylim=(1010,25),yscale="log",
		ylabel="Pressure / hPa",
		urtitle=L"T$_{RMS}$" * " = $(trms_en) K"
	)
	
	for im = 1 : nen
		aen[3].scatter(qdiff_en[:,im],pts_en[:,im],s=2,c="gray")
	end
	
	aen[3].plot(
		dropdims(mean(qdiff_en,dims=2),dims=2),
		dropdims(mean(pts_en,dims=2),dims=2),c="k"
	)
	aen[3].format(
		xlim=(-4,4),xlocator=(-2:2)*2,
		xlabel=L"qr = $\frac{q - q_{OBS}}{q_{OBS}}$ / %",
		ylim=(1010,25),yscale="log",ylabel="Pressure / hPa",
	)
	
	cqen = aen[4].contourf(
		t_en.-80,plevel,qdts_en,
		# t.-80,p,qvp .- mean(qvp[:,(end-100+1):end],dims=2),
		cmap="drywet",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",
		levels=lvls
	)
	aen[4].format(
		xlim=(0,nt_en),
		ylim=(1010,25),yscale="log",
		xlabel="time / days",ylabel="Pressure / hPa",
		urtitle=L"$qr_{RMS}$" * " = $(qrms_en) %"# * L" g kg$^{-1}$"
	)
	
	aen[2].colorbar(cten,loc="r",width=0.2)
	aen[4].colorbar(cqen,loc="r",width=0.2)
	fen.savefig("rcetdts-$(config)-ensemble.png",transparent=false,dpi=200)
	load("rcetdts-$(config)-ensemble.png")
	
end

# ╔═╡ 4db59bf0-82ec-11eb-0374-81a982a74216
md"
### D. Creating Sounding from Ensemble
"

# ╔═╡ 549c7744-82ed-11eb-03e8-7bb72c725856
begin
	qvp_μ = dropdims(mean(qvp_en[:,(end-nendays+1):end,:],dims=(2,3)),dims=(2,3))
	tem_μ = dropdims(mean(tem_en[:,(end-nendays+1):end,:],dims=(2,3)),dims=(2,3))
	pre_μ = dropdims(mean(pre_en[:,(end-nendays+1):end,:],dims=(2,3)),dims=(2,3))
	
	pot_μ = tem_μ .* (1000 ./pre_μ).^(287/1004)
	
	snddata = zeros(nz_en,6)
	snddata[:,1] .= z_en;  snddata[:,2] .= pre_μ
	snddata[:,3] .= pot_μ; snddata[:,4] .= qvp_μ
end

# ╔═╡ eba0fc7a-82eb-11eb-104e-d17de6c6c0de
md"Create SND file from ensemble? $(@bind dosnden PlutoUI.Slider(0:1))"

# ╔═╡ f53da85a-82eb-11eb-0f42-ad74458f849f
if isone(dosnden)
	  createsndmean("$(config)",snddata,psfc=1009.32)
	  md"Creating the SND file $(config) from ensemble simulations ..."
else; md"We have decided not to create the ensemble SND file $(config) yet ..."
end

# ╔═╡ Cell order:
# ╟─9dd4cd7e-5adb-11eb-2735-a7a4a2bb23b1
# ╟─417ee688-5ade-11eb-2e95-91a301119e88
# ╟─46faa412-5ade-11eb-3c37-23a7e59037a0
# ╟─b2670c08-81e5-11eb-324e-2b923b289a04
# ╟─810d5a5a-8225-11eb-37b2-6978f37f77c3
# ╟─427511d0-88cb-11eb-2a40-019c91ee1401
# ╟─ac8b9d4c-5ade-11eb-06f4-33bff063bbde
# ╟─7842e150-822b-11eb-1ded-f35ee4cc6d8c
# ╟─86f0a60a-5ae5-11eb-0b1a-935e8703b842
# ╟─ab0909d0-5ade-11eb-2622-0fee6450004b
# ╟─a198b70a-6d68-11eb-0c84-a1481c131fa3
# ╟─f9ed4c2a-6d68-11eb-1b10-ab05cb1adc4b
# ╟─e82cd648-5ade-11eb-2739-ad93c0b6d3af
# ╟─0c8cc7ec-5ae3-11eb-1c52-0fecbb575c43
# ╟─978d0442-5b33-11eb-39e5-cdea9c6efa4c
# ╟─39868b46-5ae5-11eb-0bf5-274642855e12
# ╟─4fd6e272-5b32-11eb-2f60-bbd4f8b9fc12
# ╟─094999c8-5ae5-11eb-1526-f38c604184cb
# ╟─99dca670-81e5-11eb-24b8-cf1b23d8c1f7
# ╠═6fe3de4c-82ec-11eb-023d-7dc6f51fe448
# ╠═401a6a3a-8444-11eb-3ee8-594591ed5aa9
# ╟─7d905176-81e8-11eb-20d3-b9be287472f5
# ╟─ad523b4e-81ee-11eb-2d10-a984d5983471
# ╟─c658d650-8614-11eb-252c-c3066fb1d506
# ╟─f35ab14c-8226-11eb-29b4-c3b7130f3733
# ╟─a3650e8a-822b-11eb-29ca-cd01b90e8099
# ╟─4db59bf0-82ec-11eb-0374-81a982a74216
# ╟─549c7744-82ed-11eb-03e8-7bb72c725856
# ╟─eba0fc7a-82eb-11eb-104e-d17de6c6c0de
# ╟─f53da85a-82eb-11eb-0f42-ad74458f849f
