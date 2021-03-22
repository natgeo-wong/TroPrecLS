### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ a0c2a92e-666a-11eb-3ad7-fd1b55d90d85
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ a8aacd4e-666a-11eb-2f2d-8b8df6bbf218
begin
	@quickactivate "TroPrecLS"
	using Statistics
	using Printf
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 6debc58c-666a-11eb-08bb-177fb59c634f
md"
# 6d. Plot the Time-Series for RCE state

In this notebook we just plot the timeseries of the default RCE states for runs with and without the diurnal cycle.
"

# ╔═╡ b62de6b8-666a-11eb-13f2-5d0b85c175b8
ndys = 200

# ╔═╡ bbc92150-666a-11eb-0f04-b74946b6cc2d
begin
	z3S,p3S,t3S = retrievedims("Control","3SRCE")
	tem3S = retrievevar("TABS","Control","3SRCE")
	tob3S = retrievevar("TABSOBS","Control","3SRCE")
md"Loading data from 3SRCE run ..."
end

# ╔═╡ d9dda92c-666a-11eb-299d-5558967f904f
begin
	z3P,p3P,t3P = retrievedims("Control","3PRCE")
	tem3P = retrievevar("TABS","Control","3PRCE")
	tob3P = retrievevar("TABSOBS","Control","3PRCE")
md"Loading data from 3PRCE run ..."
end

# ╔═╡ ef617b66-666a-11eb-070a-27e725c6215a
begin
	tdiff3S = dropdims(mean(tem3S[:,(end-ndys+1):end],dims=2),dims=2).-tob3S[:,1]
	tdts3S  = tem3S .- tob3S[:,1]
	trms3S = sqrt(mean(tdiff3S.^2))
	trms3S = @sprintf("%.3f",trms3S)
md"Calculating difference between `TABS` and `TABSOBS` for 3SRCE run ..."
end

# ╔═╡ 1ac8f4be-666b-11eb-1ca4-a7e32cec3195
begin
	tdiff3P = dropdims(mean(tem3P[:,(end-ndys+1):end],dims=2),dims=2).-tob3P[:,1]
	tdts3P  = tem3P .- tob3P[:,1]
	trms3P = sqrt(mean(tdiff3P.^2))
	trms3P = @sprintf("%.3f",trms3P)
md"Calculating difference between `TABS` and `TABSOBS` for 3PRCE run ..."
end

# ╔═╡ 33de6d6c-666b-11eb-0165-cf89e3117f57
begin
	
	pplt.close(); fts,ats = pplt.subplots(nrows=2,aspect=3,axwidth=5,sharex=3)
	
	lvls = [-10,-7.07,-5,-3.16,-2,-1.41,-1,-0.5,0.5,1,1.41,2,3.16,5,7.07,10]/10
	
	c = ats[1].contourf(
		t3P.-80,p3P,tdts3P,
		cmap="RdBu_r",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",levels=lvls,
	)
	ats[1].format(urtitle=L"T$_{RMS}$" * " = $(trms3P) K",ultitle="3PRCE")
	
	ats[2].contourf(
		t3S.-80,p3S,tdts3S,
		cmap="RdBu_r",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",levels=lvls,
	)
	ats[2].format(urtitle=L"T$_{RMS}$" * " = $(trms3S) K",ultitle="3SRCE")
	
	for ax in ats
		ax.format(
			xlim=(600,1000),xlabel="time / days",
			ylim=(1010,20),ylabel="Pressure / hPa",yscale="log",
		)
	end
	
	paxs = ats[1].panel("l",width="4em")
	paxs.plot(tdiff3P,p3P,lw=0.5)
	paxs.scatter(tdiff3P,p3P,s=3)
    paxs.format(xlim=(-0.2,0.2))
	
	paxs = ats[2].panel("l",width="4em")
	paxs.plot(tdiff3S,p3S,lw=0.5)
	paxs.scatter(tdiff3S,p3S,s=3)
    paxs.format(xlim=(-0.2,0.2))

	
	fts.colorbar(c,loc="r",width=0.2,label=L"T - T$_{OBS}$ / K")
	fts.savefig(plotsdir("rcetdts.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("rcetdts.png"))
	
end

# ╔═╡ Cell order:
# ╟─6debc58c-666a-11eb-08bb-177fb59c634f
# ╟─a0c2a92e-666a-11eb-3ad7-fd1b55d90d85
# ╟─a8aacd4e-666a-11eb-2f2d-8b8df6bbf218
# ╠═b62de6b8-666a-11eb-13f2-5d0b85c175b8
# ╟─bbc92150-666a-11eb-0f04-b74946b6cc2d
# ╟─d9dda92c-666a-11eb-299d-5558967f904f
# ╟─ef617b66-666a-11eb-070a-27e725c6215a
# ╟─1ac8f4be-666b-11eb-1ca4-a7e32cec3195
# ╟─33de6d6c-666b-11eb-0165-cf89e3117f57
