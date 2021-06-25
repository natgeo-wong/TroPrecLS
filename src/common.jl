using Dates
using Statistics
using GeoRegions

prect(N::Real,S::Real,W::Real,E::Real) = [W,E,E,W,W],[S,S,N,N,S]

yrmo2str(date::TimeType) = Dates.format(date,dateformat"yyyymm")

function ncoffsetscale(data::AbstractArray)

    dmax = maximum(data[.!ismissing.(data)])
    dmin = minimum(data[.!ismissing.(data)])
    scale = (dmax-dmin) / 65533;
    offset = (dmax+dmin-scale) / 2;

    return scale,offset

end

function bindatasfclnd(geo::GeoRegion,bins,var,lon,lat,lsm)

	ggrd = RegionGrid(geo,lon,lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
    rvar = zeros(nlon,nlat)
    rlsm = zeros(nlon,nlat)
    rwgt = ones(nlon,nlat) .* cosd.(reshape(ggrd.glat,1,:))

	if typeof(ggrd) <: PolyGrid
		  mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end

	for glat in 1 : nlat, glon in 1 : nlon
		rvar[glon,glat] = var[ilon[glon],ilat[glat]]
		rlsm[glon,glat] = lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
	end

    lvar = rvar[rlsm.>0.5];
	lvar = lvar[.!ismissing.(lvar)]; lvar = lvar[.!isnan.(lvar)]
    lbin = fit(Histogram,lvar,bins).weights;
	lbin = lbin ./ sum(lbin) * (length(bins) - 1)

    rvar = rvar .* cosd.(reshape(ggrd.glat,1,:))
    lvar = rvar[rlsm.>0.5];
	lvar = lvar[.!ismissing.(lvar)]; lvar = lvar[.!isnan.(lvar)]
    lvar = lvar / mean(rwgt[rlsm.>0.5])

    return lbin,mean(lvar)

end

function bindatasfcsea(geo::GeoRegion,bins,var,lon,lat,lsm)

	ggrd = RegionGrid(geo,lon,lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
    rvar = zeros(nlon,nlat)
    rlsm = zeros(nlon,nlat)
    rwgt = ones(nlon,nlat) .* cosd.(reshape(ggrd.glat,1,:))

	if typeof(ggrd) <: PolyGrid
		  mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end

	for glat in 1 : nlat, glon in 1 : nlon
		rvar[glon,glat] = var[ilon[glon],ilat[glat]]
		rlsm[glon,glat] = lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
	end

    svar = rvar[rlsm.<0.5];
	svar = svar[.!ismissing.(svar)]; svar = svar[.!isnan.(svar)]
    sbin = fit(Histogram,svar,bins).weights;
	sbin = sbin ./ sum(sbin) * (length(bins) - 1)

    rvar = rvar .* cosd.(reshape(ggrd.glat,1,:))
    svar = rvar[rlsm.<0.5];
	svar = svar[.!ismissing.(svar)]; svar = svar[.!isnan.(svar)]
    svar = svar / mean(rwgt[rlsm.<0.5])

    return sbin,mean(svar)

end

function getmean(coords,var,lon,lat,nlvl,lsm)

	ggrd = RegionGrid(geo,lon,lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
    rvar = zeros(nlon,nlat,nlvl)
    rlsm = zeros(nlon,nlat)
    rwgt = ones(nlon,nlat) .* cosd.(reshape(ggrd.glat,1,:))

	if typeof(ggrd) <: PolyGrid
		  mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end

	for ilvl = 1 : nlvl, glat in 1 : nlat, glon in 1 : nlon
		rvar[glon,glat] = var[ilon[glon],ilat[glat]]
		rlsm[glon,glat] = lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
	end

    sprf = zeros(nlvl); lprf = zeros(nlvl)

    for ilvl = 1 : nlvl

        ivar = rvar[:,:,ilvl];  ivar = ivar .* cosd.(reshape(tlat,1,:))
        lavg = ivar[rlsm.>0.5]; lavg = lavg[.!ismissing.(lavg)]; lavg = lavg[.!isnan.(lavg)]
        savg = ivar[rlsm.<0.5]; savg = savg[.!ismissing.(savg)]; savg = savg[.!isnan.(savg)]
        lprf[ilvl] = mean(lavg)
        sprf[ilvl] = mean(savg)

    end

    return lprf,sprf

end
