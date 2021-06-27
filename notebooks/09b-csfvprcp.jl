### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 858e2874-57e8-4e18-b4cc-7f65621ec1f5
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ a39537b9-34d4-463e-904a-4e1b51992e94
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using GeoRegions
	using NCDatasets
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ d9faa384-d455-11eb-32a2-1723c787fbc7
md"
# 09b. Column Saturation Fraction vs Precipitation Rate

Here, we explore the probability distribution of column saturation fraction using ERA5 data.
"

# ╔═╡ 29425799-b14a-4200-829b-619866a2024a
md"
### A. The Basic Breakdown over the Tropics
"

# ╔═╡ 22858627-a7f3-4f40-b343-75fad90c0e6e
begin
	dsg = NCDataset(datadir("reanalysis","gpm.jld2"))
	gfq = dsg["freq"][:] * 1
	gpm = dsg["prcp"][:] * 3600
	csf = dsg["csf"][:]
	close(dsg)
	md"Loading GPM CSF vs Precipitation rate data ..."
end

# ╔═╡ 8414f79b-763c-489c-bfc7-cc04ad5aac41
begin
	dse = NCDataset(datadir("reanalysis","era5.jld2"))
	efq = dse["freq"][:] * 1
	era = dse["prcp"][:] * 1000
	close(dse)
	md"Loading ERA5 CSF vs Precipitation rate data ..."
end

# ╔═╡ bcc90c75-8c85-47da-8b1a-b829bfe2f93b
begin
	dsl = NCDataset(datadir("reanalysis","era5-TRPx0.25-lsm-sfc.nc"))
	lon = dsl["longitude"][:] * 1
	lat = dsl["latitude"][:] * 1
	lsm = dsl["lsm"][:] * 1
	close(dsl)
	md"Loading land-sea mask data ..."
end

# ╔═╡ 73c5847f-9783-49dd-a7cb-736113fe75bf
begin
	lsc = pplt.Colors("Delta_r",15)
	md"Colours for different regions ..."
end

# ╔═╡ 91a53908-8090-4835-ac31-127919a175ac
md"
### B. Binning into different regions

Now, we bin the data into the different sub-regions within the Tropics as was defined in the earlier notebook `01-domain.jl`.
"

# ╔═╡ 38ca34a7-772e-4e8c-b7a9-0ea9cdebae5b
function extractocnlnd(geo,prc,frq,lon,lat,lsm)
	
	ggrd = RegionGrid(geo,lon,lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
    rprc = zeros(nlon,nlat); ncsf = size(frq,3)
    rfrq = zeros(nlon,nlat)
	rlsm = zeros(nlon,nlat)
	ocsf = zeros(ncsf)
	lcsf = zeros(ncsf)

	if typeof(ggrd) <: PolyGrid
		  mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end

	for icsf = 1 : ncsf
		for glat in 1 : nlat, glon in 1 : nlon
			rprc[glon,glat] = prc[ilon[glon],ilat[glat],icsf]
			rfrq[glon,glat] = frq[ilon[glon],ilat[glat],icsf]
			rlsm[glon,glat] = lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
		end
		ocsf[icsf] = sum(rprc[rlsm.<0.5]) ./ sum(rfrq[rlsm.<0.5])
		lcsf[icsf] = sum(rprc[rlsm.>0.5]) ./ sum(rfrq[rlsm.>0.5])
	end
	
	return ocsf,lcsf
	
end

# ╔═╡ 1cd2bee0-0e9a-4c82-afaf-f54dbe4a359a
function normbin(bindata)
	return bindata / sum(bindata) * length(bindata)
end

# ╔═╡ 48e74e01-c1ef-456e-9e6c-5f6da0ce93bf
begin
	ogfrq_TRP,lgfrq_TRP = extractocnlnd(GeoRegion("TRP"),gpm,gfq,lon,lat,lsm)
	ogfrq_DTP,lgfrq_DTP = extractocnlnd(GeoRegion("DTP"),gpm,gfq,lon,lat,lsm)
	ogfrq_SEA,lgfrq_SEA = extractocnlnd(GeoRegion("SEA"),gpm,gfq,lon,lat,lsm)
	ogfrq_CRB,lgfrq_CRB = extractocnlnd(GeoRegion("CRB"),gpm,gfq,lon,lat,lsm)
	ogfrq_TRA,lgfrq_TRA = extractocnlnd(GeoRegion("TRA"),gpm,gfq,lon,lat,lsm)
	ogfrq_AMZ,lgfrq_AMZ = extractocnlnd(GeoRegion("AMZ"),gpm,gfq,lon,lat,lsm)
	ogfrq_EAO,lgfrq_EAO = extractocnlnd(GeoRegion("AR6_EAO"),gpm,gfq,lon,lat,lsm)
	ogfrq_EPO,lgfrq_EPO = extractocnlnd(GeoRegion("AR6_EPO"),gpm,gfq,lon,lat,lsm)
	ogfrq_EIO,lgfrq_EIO = extractocnlnd(GeoRegion("AR6_EIO"),gpm,gfq,lon,lat,lsm)
	ogfrq_SAS,lgfrq_SAS = extractocnlnd(GeoRegion("AR6_SAS"),gpm,gfq,lon,lat,lsm)
	md"Extracting Regional Data for GPM data ..."
end

# ╔═╡ 35c157f3-fdcc-4259-888d-ea6af3aa362a
begin
	oefrq_TRP,lefrq_TRP = extractocnlnd(GeoRegion("TRP"),era,efq,lon,lat,lsm)
	oefrq_DTP,lefrq_DTP = extractocnlnd(GeoRegion("DTP"),era,efq,lon,lat,lsm)
	oefrq_SEA,lefrq_SEA = extractocnlnd(GeoRegion("SEA"),era,efq,lon,lat,lsm)
	oefrq_CRB,lefrq_CRB = extractocnlnd(GeoRegion("CRB"),era,efq,lon,lat,lsm)
	oefrq_TRA,lefrq_TRA = extractocnlnd(GeoRegion("TRA"),era,efq,lon,lat,lsm)
	oefrq_AMZ,lefrq_AMZ = extractocnlnd(GeoRegion("AMZ"),era,efq,lon,lat,lsm)
	oefrq_EAO,lefrq_EAO = extractocnlnd(GeoRegion("AR6_EAO"),era,efq,lon,lat,lsm)
	oefrq_EPO,lefrq_EPO = extractocnlnd(GeoRegion("AR6_EPO"),era,efq,lon,lat,lsm)
	oefrq_EIO,lefrq_EIO = extractocnlnd(GeoRegion("AR6_EIO"),era,efq,lon,lat,lsm)
	md"Extracting Regional Data for ERA5 data ..."
end

# ╔═╡ d0aebfcd-6068-4ff0-9aae-fef936be7a99
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=2,aspect=2,axwidth=3)
	
	lgd = Dict("ncol"=>1,"frame"=>false)
	
	a2[1].plot(csf,ogfrq_TRP,lw=1,c="r")
	a2[1].plot(csf,ogfrq_DTP,lw=1,c="k")
	a2[1].plot(csf,ogfrq_EPO,lw=1,c=lsc[13])
	a2[1].plot(csf,ogfrq_EIO,lw=1,c=lsc[12])
	a2[1].plot(csf,ogfrq_EAO,lw=1,c=lsc[11])
	a2[1].plot(csf,ogfrq_CRB,lw=1,c=lsc[10])
	a2[1].plot(csf,ogfrq_SEA,lw=1,c=lsc[5])
	a2[1].plot(csf,lgfrq_TRP,lw=1,c="r",linestyle="--")
	a2[1].plot(csf,lgfrq_DTP,lw=1,c="k",linestyle="--")
	a2[1].plot(csf,lgfrq_CRB,lw=1,c=lsc[10],linestyle="--")
	a2[1].plot(csf,lgfrq_SEA,lw=1,c=lsc[5],linestyle="--")
	a2[1].plot(csf,lgfrq_AMZ,lw=1,c=lsc[4],linestyle="--")
	a2[1].plot(csf,lgfrq_TRA,lw=1,c=lsc[3],linestyle="--")
	# a2[1].plot(csf,lgfrq_SAS,lw=1,c="r",linestyle="--")
	a2[1].format(ltitle="(a) GPM Precipitation")
	
	a2[2].plot(csf,oefrq_DTP,label="TRP",legend="r",lw=1,c="r")
	a2[2].plot(csf,oefrq_DTP,label="DTP",legend="r",lw=1,c="k")
	a2[2].plot(csf,oefrq_EPO,label="AR6_EPO",legend="r",lw=1,c=lsc[13])
	a2[2].plot(csf,oefrq_EIO,label="AR6_EIO",legend="r",lw=1,c=lsc[12])
	a2[2].plot(csf,oefrq_EAO,label="AR6_EAO",legend="r",lw=1,c=lsc[11])
	a2[2].plot(csf,oefrq_CRB,label="CRB",legend="r",lw=1,c=lsc[10])
	a2[2].plot(csf,oefrq_SEA,label="SEA",legend="r",lw=1,c=lsc[5])
	a2[2].plot(csf,csf*NaN,label="AMZ",legend="r",lw=1,c=lsc[4])
	a2[2].plot(csf,csf*NaN,label="TRA",legend="r",lw=1,c=lsc[3],legend_kw=lgd)
	a2[2].plot(csf,lefrq_TRP,lw=1,c="r",linestyle="--")
	a2[2].plot(csf,lefrq_DTP,lw=1,c="k",linestyle="--")
	a2[2].plot(csf,lefrq_CRB,lw=1,c=lsc[10],linestyle="--")
	a2[2].plot(csf,lefrq_SEA,lw=1,c=lsc[5],linestyle="--")
	a2[2].plot(csf,lefrq_AMZ,lw=1,c=lsc[4],linestyle="--")
	a2[2].plot(csf,lefrq_TRA,lw=1,c=lsc[3],linestyle="--")
	a2[2].format(ltitle="(b) ERA5 Precipitation")
	
	for ax in a2
		ax.format(
			xlim=(0,100),xlabel="Column Relative Humidity",
			ylim=10. .^(-3.5,1.5),ylabel=L"Precipitation Rate / mm hr$^{-1}$",
			yscale="log",ylocator=10. .^(-3:1)
		)
	end
	
	f2.savefig(plotsdir("csffreq.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("csffreq.png"))
end

# ╔═╡ Cell order:
# ╟─d9faa384-d455-11eb-32a2-1723c787fbc7
# ╟─858e2874-57e8-4e18-b4cc-7f65621ec1f5
# ╟─a39537b9-34d4-463e-904a-4e1b51992e94
# ╟─29425799-b14a-4200-829b-619866a2024a
# ╟─22858627-a7f3-4f40-b343-75fad90c0e6e
# ╟─8414f79b-763c-489c-bfc7-cc04ad5aac41
# ╟─bcc90c75-8c85-47da-8b1a-b829bfe2f93b
# ╟─73c5847f-9783-49dd-a7cb-736113fe75bf
# ╟─91a53908-8090-4835-ac31-127919a175ac
# ╠═38ca34a7-772e-4e8c-b7a9-0ea9cdebae5b
# ╠═1cd2bee0-0e9a-4c82-afaf-f54dbe4a359a
# ╟─48e74e01-c1ef-456e-9e6c-5f6da0ce93bf
# ╟─35c157f3-fdcc-4259-888d-ea6af3aa362a
# ╟─d0aebfcd-6068-4ff0-9aae-fef936be7a99
