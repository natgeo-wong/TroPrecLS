### A Pluto.jl notebook ###
# v0.14.7

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

# ╔═╡ 129f73ae-a99d-4fec-9d0a-2db8926c16f7
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 2e6eeaf3-474f-4883-85f0-558ded3e9fb8
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using PlutoUI
	using Printf
	
	include(srcdir("common.jl"))
	include(srcdir("samlsf.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ ee62bca6-d6dc-11eb-25ef-75831cbc1575
md"
# 05c. QFluxes from an LSF File

We do have some qflux data that we used from lsf files, so in order to not have to redownload data again and again when scratch storage expires, extract data from lsf file.
"

# ╔═╡ 37a33543-16a2-4752-855f-466c855bd8b7
md"
### A. Let's read in the LSF data!
"

# ╔═╡ b4c3910d-8428-4848-bbc9-3ffa52048456
begin
	lsf  = readdlm(projectdir("exp/lsf/ocean/p1d0"),',')
	pre  = lsf[2,3]
	nlvl = lsf[2,2]
	lsf  = Float64.(lsf[3:(2+nlvl),:])
	qflx = lsf[:,4]
	md"Loading large-scale forcing data ..."
end

# ╔═╡ b2b7d523-19e2-4562-951c-b6a17a71ca65
md"
### B. And now we create new forcing files!
"

# ╔═╡ a39e37bb-8fd0-4341-9c5f-0efac04f7758
md"Create LSF files? $(@bind dolsf PlutoUI.Slider(0:1))"

# ╔═╡ d6547c85-0ba7-41cf-980e-1a38a889822d
if isone(dolsf)
	
	mvec = -1.5 : 0.1 : 1.5
	
	for mul in mvec
		
		lsf[:,4] .= qflx * mul
		
		if mul > 0
			  mstr = @sprintf("p%03.1f",abs(mul))
		else; mstr = @sprintf("n%03.1f",abs(mul))
		end
		mstr = replace(mstr,"."=>"d")
		lsfprint(projectdir("exp/lsf/ocean/$(mstr)"),lsf,pre)
		
	end
	
	md"Based on these profiles, we create large-scale forcing profiles for moisture flux convergence to be used in our WTG simulations ..."
else
	md"We have decided not to override any preexisting large-scale forcing profiles for moisture flux convergence."
end

# ╔═╡ Cell order:
# ╟─ee62bca6-d6dc-11eb-25ef-75831cbc1575
# ╟─129f73ae-a99d-4fec-9d0a-2db8926c16f7
# ╟─2e6eeaf3-474f-4883-85f0-558ded3e9fb8
# ╟─37a33543-16a2-4752-855f-466c855bd8b7
# ╟─b4c3910d-8428-4848-bbc9-3ffa52048456
# ╟─b2b7d523-19e2-4562-951c-b6a17a71ca65
# ╟─a39e37bb-8fd0-4341-9c5f-0efac04f7758
# ╠═d6547c85-0ba7-41cf-980e-1a38a889822d
