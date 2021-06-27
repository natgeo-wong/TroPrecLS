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
# 09a. The Basic Distribution of Column Saturation Fraction

Here, we explore the probability distribution of column saturation fraction using ERA5 data.
"

# ╔═╡ 29425799-b14a-4200-829b-619866a2024a
md"
### A. The Basic Breakdown over the Tropics
"

# ╔═╡ 22858627-a7f3-4f40-b343-75fad90c0e6e
begin
	dsc = NCDataset(datadir("reanalysis","era5-TRPx0.25-csf-sfc.nc"))
	lon = dsc["longitude"][:] * 1
	lat = dsc["latitude"][:] * 1
	csf = dsc["average"][:] * 1;
	close(dsc)
	md"Loading compiled CSF data ..."
end

# ╔═╡ 786e66a9-b762-481e-a6ea-c1d4fb3b3e82
begin
	dsf = NCDataset(datadir("reanalysis","era5-TRPx0.25-csffreq-sfc.nc"))
	cbn = dsf["csf"][:] * 1; cbn = cbn[1:(end-1)] .+ (cbn[2]-cbn[1])/2;
	frq = dsf["bin_frq"][:] * 1
	close(dsf)
	md"Loading bin frequency data ..."
end

# ╔═╡ bcc90c75-8c85-47da-8b1a-b829bfe2f93b
begin
	dsl = NCDataset(datadir("reanalysis","era5-TRPx0.25-lsm-sfc.nc"))
	lsm = dsl["lsm"][:] * 1
	close(dsl)
	md"Loading land-sea mask data ..."
end

# ╔═╡ 308c4056-6527-4962-8aba-bd112cae72f5
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ 73c5847f-9783-49dd-a7cb-736113fe75bf
begin
	lsc = pplt.Colors("Delta_r",15)
	md"Colours for different regions ..."
end

# ╔═╡ 025f1fb5-84dc-4c13-a7c8-a513b56e0ace
begin
	pplt.close(); fs,as = pplt.subplots(aspect=6,axwidth=6)
	
	cs = as[1].contourf(lon,lat,csf',cmap="Blues",levels=20:5:80,extend="both")
	as[1].plot(x,y,c="k",lw=0.5)
	as[1].format(
		xlim=(0,360),xlabel=L"Longitude / $\degree$",xlocator=0:60:360,
		ylim=(-30,30),ylabel=L"Latitude / $\degree$",ylocator=-30:10:30,
		suptitle="Column Relative Humidity / %"
	)
	as[1].colorbar(cs,loc="r",locator=20:15:80)
	
	fs.savefig(plotsdir("csf-spatial.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("csf-spatial.png"))
end

# ╔═╡ 91a53908-8090-4835-ac31-127919a175ac
md"
### B. Binning into different regions

Now, we bin the data into the different sub-regions within the Tropics as was defined in the earlier notebook `01-domain.jl`.
"

# ╔═╡ 38ca34a7-772e-4e8c-b7a9-0ea9cdebae5b
function extractocnlnd(geo,frq,lon,lat,lsm)
	
	ggrd = RegionGrid(geo,lon,lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
    rcsf = zeros(nlon,nlat); ncsf = size(frq,3)
    rlsm = zeros(nlon,nlat)
	ocsf = zeros(ncsf)
	lcsf = zeros(ncsf)

	if typeof(ggrd) <: PolyGrid
		  mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end

	for icsf = 1 : ncsf
		for glat in 1 : nlat, glon in 1 : nlon
			rcsf[glon,glat] = frq[ilon[glon],ilat[glat],icsf]
			rlsm[glon,glat] = lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
		end
		ocsf[icsf] = sum(rcsf[rlsm.<0.5])
		lcsf[icsf] = sum(rcsf[rlsm.>0.5])
	end
	
	return ocsf,lcsf
	
end

# ╔═╡ 1cd2bee0-0e9a-4c82-afaf-f54dbe4a359a
function normbin(bindata)
	return bindata / sum(bindata) * length(bindata)
end

# ╔═╡ 48e74e01-c1ef-456e-9e6c-5f6da0ce93bf
begin
	mfrq = dropdims(sum(frq,dims=(1,2)),dims=(1,2))
	ofrq_TRP,lfrq_TRP = extractocnlnd(GeoRegion("TRP"),frq,lon,lat,lsm)
	ofrq_DTP,lfrq_DTP = extractocnlnd(GeoRegion("DTP"),frq,lon,lat,lsm)
	ofrq_SEA,lfrq_SEA = extractocnlnd(GeoRegion("SEA"),frq,lon,lat,lsm)
	ofrq_CRB,lfrq_CRB = extractocnlnd(GeoRegion("CRB"),frq,lon,lat,lsm)
	ofrq_TRA,lfrq_TRA = extractocnlnd(GeoRegion("TRA"),frq,lon,lat,lsm)
	ofrq_AMZ,lfrq_AMZ = extractocnlnd(GeoRegion("AMZ"),frq,lon,lat,lsm)
	ofrq_SAS,lfrq_SAS = extractocnlnd(GeoRegion("AR6_SAS"),frq,lon,lat,lsm)
	ofrq_EAO,lfrq_EAO = extractocnlnd(GeoRegion("AR6_EAO"),frq,lon,lat,lsm)
	ofrq_EPO,lfrq_EPO = extractocnlnd(GeoRegion("AR6_EPO"),frq,lon,lat,lsm)
	ofrq_EIO,lfrq_EIO = extractocnlnd(GeoRegion("AR6_EIO"),frq,lon,lat,lsm)
	md"Extracting Regional Data ..."
end

# ╔═╡ 953d0573-9341-4998-b67a-30939b47dff3
begin
	pplt.close(); f1,a1 = pplt.subplots(ncols=1,aspect=2,axwidth=3)
	
	lgd = Dict("ncol"=>1,"frame"=>false)
	
	a1[1].plot(cbn,normbin(ofrq_TRP),c="b",label="TRP Ocean",legend="r")
	a1[1].plot(cbn,normbin(lfrq_TRP),c="g",label="TRP Land",legend="r",legend_kw=lgd)
	for ax in a1
		ax.format(
			xlim=(0,100),xlabel="Column Relative Humidity",
			ylim=(0,2),ylabel="Normalized Frequency"
		)
	end
	
	f1.savefig(plotsdir("csffreq-TRP.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("csffreq-TRP.png"))
end

# ╔═╡ d0aebfcd-6068-4ff0-9aae-fef936be7a99
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=2,aspect=2,axwidth=3)
	
	a2[1].plot(cbn,normbin(ofrq_DTP),c="k")
	a2[1].plot(cbn,normbin(ofrq_EPO),c=lsc[13])
	a2[1].plot(cbn,normbin(ofrq_EIO),c=lsc[12])
	a2[1].plot(cbn,normbin(ofrq_EAO),c=lsc[11])
	a2[1].plot(cbn,normbin(ofrq_CRB),c=lsc[10])
	a2[1].plot(cbn,normbin(ofrq_SEA),c=lsc[5])
	a2[1].format(ltitle="(a) Ocean")
	
	a2[2].plot(cbn,normbin(lfrq_DTP),label="DTP",legend="r",c="k")
	a2[2].plot(cbn,cbn*NaN,label="AR6_EPO",legend="r",c=lsc[13])
	a2[2].plot(cbn,cbn*NaN,label="AR6_EIO",legend="r",c=lsc[12])
	a2[2].plot(cbn,cbn*NaN,label="AR6_EAO",legend="r",c=lsc[11])
	a2[2].plot(cbn,normbin(lfrq_CRB),label="CRB",legend="r",c=lsc[10])
	a2[2].plot(cbn,normbin(lfrq_SEA),label="SEA",legend="r",c=lsc[5])
	a2[2].plot(cbn,normbin(lfrq_SAS),label="AR6_SAS",legend="r",c=lsc[4])
	a2[2].plot(cbn,normbin(lfrq_AMZ),label="AMZ",legend="r",c=lsc[4])
	a2[2].plot(cbn,normbin(lfrq_TRA),label="TRA",legend="r",c=lsc[3],legend_kw=lgd)
	a2[2].format(ltitle="(b) Land")
	
	for ax in a2
		ax.format(
			xlim=(0,100),xlabel="Column Relative Humidity",
			ylim=(0,5),ylabel="Normalized Frequency"
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
# ╟─786e66a9-b762-481e-a6ea-c1d4fb3b3e82
# ╟─bcc90c75-8c85-47da-8b1a-b829bfe2f93b
# ╟─308c4056-6527-4962-8aba-bd112cae72f5
# ╟─73c5847f-9783-49dd-a7cb-736113fe75bf
# ╟─025f1fb5-84dc-4c13-a7c8-a513b56e0ace
# ╟─953d0573-9341-4998-b67a-30939b47dff3
# ╟─91a53908-8090-4835-ac31-127919a175ac
# ╠═38ca34a7-772e-4e8c-b7a9-0ea9cdebae5b
# ╠═1cd2bee0-0e9a-4c82-afaf-f54dbe4a359a
# ╠═48e74e01-c1ef-456e-9e6c-5f6da0ce93bf
# ╠═d0aebfcd-6068-4ff0-9aae-fef936be7a99
