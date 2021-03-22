### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ abe749ca-717c-11eb-0f5b-8d4339fd2352
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ abb08b92-717c-11eb-19dd-61af1db1ca9b
begin
	@quickactivate "TroPrecLS"
	using NCDatasets
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 3693f7fe-717c-11eb-1050-1758bf54317f
md"
# 7a. Exploring the RCE 2D Fields

In this notebook, we explore the precipitation fields in the Control RCE experiments.  Key questions include whether there is organized convection / convective aggregation in our simulations, like in many other RCE experiments.
"

# ╔═╡ 8844b87a-7ea8-11eb-3ae8-3fcec98ca351
function retrievedim(
    variable::AbstractString,
	experiment::AbstractString,
	config::AbstractString
)

	rce = NCDataset(datadir(joinpath(
        experiment,config,"OUT_2D",
        "RCE_TroPrecLS-$(experiment)_256.2Dbin_1.nc"
    )))
	
	x   = rce["x"][:]/1000
	y   = rce["y"][:]/1000
	t   = rce["time"][:]
	t   = t .- floor(t[1])

    return x,y,t,rce

end

# ╔═╡ 09d74e22-717d-11eb-2d63-c3eac728c855
begin
	x,y,t,rce  = retrievedim("PW","Control","ODCNV")
	pwv  = rce["PW"][:,:,end-0]
	swp  = rce["SWVP"][:,:,end-0]
	var  = rce["W500"][:,:,end-0]
	csf = pwv./swp*100
	md"Extracting Precipitation Data ..."
end

# ╔═╡ db53919a-7ea7-11eb-005b-b58724933383
begin
	pplt.close(); f2D,a2D = pplt.subplots(axwidth=3)
	hr = @sprintf("%04.1f",mod(t[end],1)*24)
	
	lvl = vcat(-5,-2,-1,-0.5,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.5,1,2,5)
	c = a2D[1].pcolormesh(
		x[1:50],y[1:50],var[1:50,1:50]',
		cmap="RdBu_r",cmap_kw=Dict(),
		extend="both",levels=lvl
	)
	# a2D[1].contour(
	# 	x,y,prcp',lw=1,
	# 	levels=vcat(0,10. .^(-2:2)),norm="segmented",
	# 	cmap="Blues",cmap_kw=Dict("left"=>0.2)
	# )
	# a2D[1].contour(x,y,prcp',lw=1,linestyle="--",color="k",levels=[0])
	a2D[1].format(
		xlabel="X / km",ylabel="Y / km",
		ltitle="time = Day $(floor(t[end])), Hour $(hr)"
	)
	
	f2D.colorbar(c,loc="r",width=0.2)
	f2D.savefig("test2D/end.png",transparent=false,dpi=120)
	PNGFiles.load("test2D/end.png")
end

# ╔═╡ afd094a2-8097-11eb-1311-6331f3ac3a6b
sum(var.<0) / (256*256)

# ╔═╡ 309a78f4-717d-11eb-0c95-8947f0ac1288
# begin
	
# 	for it = 2161 : 2401; hr = @sprintf("%04.1f",mod(t[it],1)*24)
# 		swp  = rce["SWVP"][:,:,it]
# 		tcw  = rce["PW"][:,:,it]
# 		csf  = tcw ./ swp * 100
# 		prcp = rce["Prec"][:,:,it] / 24
		
# 		pplt.close(); f2D,a2D = pplt.subplots(axwidth=3)
# 		c = a2D[1].contourf(x,y,csf',levels=10:5:90,cmap="drywet",extend="both")
# 		a2D[1].contour(
# 			x,y,prcp',lw=1,
# 			levels=vcat(0,10. .^(-2:2)),norm="segmented",
# 			cmap="Blues",cmap_kw=Dict("left"=>0.2)
# 		)
# 		a2D[1].contour(x,y,prcp',lw=1,linestyle="--",color="k",levels=[0])
# 		a2D[1].format(
# 			xlabel="X / km",ylabel="Y / km",
# 			ltitle="time = Day $(floor(t[it])), Hour $(hr)"
# 		)

# 		f2D.colorbar(c,loc="r",width=0.2)
# 		f2D.savefig("test2D/$(it).png",transparent=false,dpi=120)
# 	end
	
# end

# ╔═╡ Cell order:
# ╟─3693f7fe-717c-11eb-1050-1758bf54317f
# ╟─abe749ca-717c-11eb-0f5b-8d4339fd2352
# ╟─abb08b92-717c-11eb-19dd-61af1db1ca9b
# ╠═8844b87a-7ea8-11eb-3ae8-3fcec98ca351
# ╠═09d74e22-717d-11eb-2d63-c3eac728c855
# ╠═db53919a-7ea7-11eb-005b-b58724933383
# ╠═afd094a2-8097-11eb-1311-6331f3ac3a6b
# ╠═309a78f4-717d-11eb-0c95-8947f0ac1288
