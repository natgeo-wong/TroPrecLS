using NCDatasets
using Statistics
using StatsBase

include(srcdir("common.jl"))

function compilesavesfc(varname::AbstractString)

    tds = NCDataset(datadir("reanalysis/era5-TRPx0.25-sfc.nc"))
    var = dropdims(mean(tds[varname][:]*1,dims=3),dims=3)
    lon = tds["longitude"][:]*1
    lat = tds["latitude"][:]*1

    fnc = datadir("reanalysis/era5-TRPx0.25-$varname-sfc.nc")
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions"               => "CF-1.6",
        "history"                   => "2020-10-21 23:54:22 GMT by grib_to_netcdf-2.16.0: /opt/ecmwf/eccodes/bin/grib_to_netcdf -S param -o /cache/data4/adaptor.mars.internal-1603323636.0468726-6113-3-fb54412a-f0a5-4783-956d-46233705e403.nc /cache/tmp/fb54412a-f0a5-4783-956d-46233705e403-adaptor.mars.internal-1603323636.0477095-6113-1-tmp.grib",
    ))

    # Dimensions

    ds.dim["longitude"] = length(lon)
    ds.dim["latitude"]  = length(lat)

    # Declare variables

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"                     => "degrees_east",
        "long_name"                 => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"                     => "degrees_north",
        "long_name"                 => "latitude",
    ))

    scale,offset = ncoffsetscale(var)

    varattribs = Dict(
        "scale_factor"  => scale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
        "units"         => tds[varname].attrib["units"],
        "long_name"     => tds[varname].attrib["long_name"],
    )

    if haskey(tds[varname].attrib,"standard_name")
        varattribs["standard_name"] = tds[varname].attrib["standard_name"]
    end

    ncvar = defVar(ds,varname,Int16,("longitude","latitude"),attrib = varattribs)

    nclon[:] = lon
    nclat[:] = lat
    ncvar[:] = var

    close(ds)

    close(tds)

end

function compilesavesfcfeb(varname::AbstractString)

    tds = NCDataset(datadir("reanalysis/era5-TRPx0.25-sfc.nc"))
    var = dropdims(mean(tds[varname][:,:,2:12:end]*1,dims=3),dims=3)
    lon = tds["longitude"][:]*1
    lat = tds["latitude"][:]*1

    fnc = datadir("reanalysis/era5feb-TRPx0.25-$varname-sfc.nc")
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions"               => "CF-1.6",
        "history"                   => "2020-10-21 23:54:22 GMT by grib_to_netcdf-2.16.0: /opt/ecmwf/eccodes/bin/grib_to_netcdf -S param -o /cache/data4/adaptor.mars.internal-1603323636.0468726-6113-3-fb54412a-f0a5-4783-956d-46233705e403.nc /cache/tmp/fb54412a-f0a5-4783-956d-46233705e403-adaptor.mars.internal-1603323636.0477095-6113-1-tmp.grib",
    ))

    # Dimensions

    ds.dim["longitude"] = length(lon)
    ds.dim["latitude"]  = length(lat)

    # Declare variables

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"                     => "degrees_east",
        "long_name"                 => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"                     => "degrees_north",
        "long_name"                 => "latitude",
    ))

    scale,offset = ncoffsetscale(var)

    varattribs = Dict(
        "scale_factor"  => scale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
        "units"         => tds[varname].attrib["units"],
        "long_name"     => tds[varname].attrib["long_name"],
    )

    if haskey(tds[varname].attrib,"standard_name")
        varattribs["standard_name"] = tds[varname].attrib["standard_name"]
    end

    ncvar = defVar(ds,varname,Int16,("longitude","latitude"),attrib = varattribs)

    nclon[:] = lon
    nclat[:] = lat
    ncvar[:] = var

    close(ds)

    close(tds)

end

function compilesavepre(varname::AbstractString)

    tds = NCDataset(datadir("reanalysis/$varname/era5-TRPx0.25-$varname-1979.nc"))
    lon = tds["longitude"][:]*1
    lat = tds["latitude"][:]*1
    lvl = tds["level"][:]*1

    var = zeros(length(lon),length(lat),length(lvl))
    vnc = replace(varname,"_air"=>"")

    for yr = 1979 : 2019
        yds  = NCDataset(datadir("reanalysis/$varname/era5-TRPx0.25-$varname-$(yr).nc"))
        var += dropdims(mean(yds[vnc][:]*1,dims=4),dims=4)
        close(yds)
    end

    var = var / 41

    fnc = datadir("reanalysis/era5-TRPx0.25-$varname.nc")
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions"               => "CF-1.6",
        "history"                   => "2020-10-21 23:54:22 GMT by grib_to_netcdf-2.16.0: /opt/ecmwf/eccodes/bin/grib_to_netcdf -S param -o /cache/data4/adaptor.mars.internal-1603323636.0468726-6113-3-fb54412a-f0a5-4783-956d-46233705e403.nc /cache/tmp/fb54412a-f0a5-4783-956d-46233705e403-adaptor.mars.internal-1603323636.0477095-6113-1-tmp.grib",
    ))

    # Dimensions

    ds.dim["longitude"] = length(lon)
    ds.dim["latitude"]  = length(lat)
    ds.dim["level"]     = length(lvl)

    # Declare variables

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"                     => "degrees_east",
        "long_name"                 => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"                     => "degrees_north",
        "long_name"                 => "latitude",
    ))

    ncpre = defVar(ds,"level",Int32,("level",),attrib = Dict(
        "units"                     => "millibars",
        "long_name"                 => "pressure_level",
    ))

    scale,offset = ncoffsetscale(var)

    varattribs = Dict(
        "scale_factor"  => scale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
        "units"         => tds[vnc].attrib["units"],
        "long_name"     => tds[vnc].attrib["long_name"],
    )

    if haskey(tds[vnc].attrib,"standard_name")
        varattribs["standard_name"] = tds[vnc].attrib["standard_name"]
    end

    ncvar = defVar(ds,varname,Int16,("longitude","latitude","level"),attrib = varattribs)

    nclon[:] = lon
    nclat[:] = lat
    ncpre[:] = lvl
    ncvar[:] = var

    close(ds)

    close(tds)

end

function compilesavesfchour(varname::AbstractString)

    tds = NCDataset(datadir("reanalysis/sfc/era5-TRPx0.25-sfc-1979.nc"))
    lon = tds["longitude"][:]*1; nlon = length(lon)
    lat = tds["latitude"][:]*1;  nlat = length(lat)
    t   = tds["time"][:];        nt   = length(t)

    vart = zeros(nlon,nlat,nt)
    for yr = 1979 : 2019
        fnc   = "era5-TRPx0.25-sfc-$yr.nc"
        yds   = NCDataset(datadir("reanalysis/sfc/$(fnc)"))
        vart += yds[varname][:]*1
        close(yds)
    end

    vart = vart / 41
    vart = dropdims(mean(reshape(vart,nlon,nlat,:,12),dims=4),dims=4)
    vart[ismissing.(vart)] .= 0; vart = vart*1

    nhr = size(vart,3); mhr = nhr - 24; mind = 2*mhr
    var = zeros(nlon,nlat,24)
    var[:,:,(mhr+1):end] .= vart[:,:,(mind+1):end]
    var[:,:,1:mhr] .= dropdims(sum(reshape(vart[:,:,1:mind],nlon,nlat,2,:),dims=3),dims=3)

    fnc = datadir("reanalysis/era5-TRPx0.25-$varname-sfc-hour.nc")
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions"               => "CF-1.6",
        "history"                   => "2020-10-21 23:54:22 GMT by grib_to_netcdf-2.16.0: /opt/ecmwf/eccodes/bin/grib_to_netcdf -S param -o /cache/data4/adaptor.mars.internal-1603323636.0468726-6113-3-fb54412a-f0a5-4783-956d-46233705e403.nc /cache/tmp/fb54412a-f0a5-4783-956d-46233705e403-adaptor.mars.internal-1603323636.0477095-6113-1-tmp.grib",
    ))

    # Dimensions

    ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ds.dim["hour"]      = 24

    # Declare variables

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"                     => "degrees_east",
        "long_name"                 => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"                     => "degrees_north",
        "long_name"                 => "latitude",
    ))

    nchr  = defVar(ds,"hour",Int8,("hour",),attrib = Dict(
        "units"                     => "hours",
        "long_name"                 => "Hours past 00:00 (GMT)",
    ))

    scale,offset = ncoffsetscale(var)

    varattribs = Dict(
        "scale_factor"  => scale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
        "units"         => tds[varname].attrib["units"],
        "long_name"     => tds[varname].attrib["long_name"],
    )

    if haskey(tds[varname].attrib,"standard_name")
        varattribs["standard_name"] = tds[varname].attrib["standard_name"]
    end

    ncvar = defVar(ds,varname,Int16,("longitude","latitude","hour"),attrib = varattribs)

    nclon[:] = lon
    nclat[:] = lat
    nchr[:]  = collect(0:23)
    ncvar[:] = var

    close(ds)

    close(tds)

end

function compilesaveprehour(varname::AbstractString;levels::AbstractVector{<:Real})

    tds = NCDataset(datadir("reanalysis/$varname/era5-TRPx0.25-$varname-197901.nc"))
    lon = tds["longitude"][:]*1; nlon = length(lon)
    lat = tds["latitude"][:]*1;  nlat = length(lat)
                                 nlvl = length(levels)

    var = zeros(nlon,nlat,nlvl,24); ilvl = 0;
    vnc = replace(varname,"_air"=>"")

    for yr = 1979 : 2019, mo = 1 : 12
        fnc  = "era5-TRPx0.25-$varname-$(yrmo2str(Date(yr,mo))).nc"
        yds  = NCDataset(datadir("reanalysis/$(varname)/$(fnc)"))
        var += yds[vnc][:]*1
        close(yds)
    end

    var = var / (41*12)

    fnc = datadir("reanalysis/era5-TRPx0.25-$(varname)-hour.nc")
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions"               => "CF-1.6",
        "history"                   => "2020-10-21 23:54:22 GMT by grib_to_netcdf-2.16.0: /opt/ecmwf/eccodes/bin/grib_to_netcdf -S param -o /cache/data4/adaptor.mars.internal-1603323636.0468726-6113-3-fb54412a-f0a5-4783-956d-46233705e403.nc /cache/tmp/fb54412a-f0a5-4783-956d-46233705e403-adaptor.mars.internal-1603323636.0477095-6113-1-tmp.grib",
    ))

    # Dimensions

    ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ds.dim["level"]     = nlvl
    ds.dim["hour"]      = 24

    # Declare variables

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"                     => "degrees_east",
        "long_name"                 => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"                     => "degrees_north",
        "long_name"                 => "latitude",
    ))

    ncpre = defVar(ds,"level",Int32,("level",),attrib = Dict(
        "units"                     => "millibars",
        "long_name"                 => "pressure_level",
    ))

    nchr  = defVar(ds,"hour",Int8,("hour",),attrib = Dict(
        "units"                     => "hours",
        "long_name"                 => "Hours past 00:00 (GMT)",
    ))

    scale,offset = ncoffsetscale(var)

    varattribs = Dict(
        "scale_factor"  => scale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
        "units"         => tds[vnc].attrib["units"],
        "long_name"     => tds[vnc].attrib["long_name"],
    )

    if haskey(tds[vnc].attrib,"standard_name")
        varattribs["standard_name"] = tds[vnc].attrib["standard_name"]
    end

    ncvar = defVar(
        ds,varname,Int16,
        ("longitude","latitude","level","hour"),
        attrib = varattribs
    )

    nclon[:] = lon
    nclat[:] = lat
    ncpre[:] = levels
    nchr[:]  = collect(0:23)
    ncvar[:] = var

    close(ds)

    close(tds)

end

function compilesfceb()

    tds  = NCDataset(datadir("reanalysis/era5-TRPx0.25-sfc.nc"))
    lon  = tds["longitude"][:]*1
    lat  = tds["latitude"][:]*1
    close(tds)

    tds = NCDataset(datadir("reanalysis/era5-TRPx0.25-ssr-sfc.nc"))
    ssr = tds["ssr"][:]*1
    close(tds)
    tds = NCDataset(datadir("reanalysis/era5-TRPx0.25-str-sfc.nc"))
    str = tds["str"][:]*1
    close(tds)
    tds = NCDataset(datadir("reanalysis/era5-TRPx0.25-sshf-sfc.nc"))
    shf = tds["sshf"][:]*1
    close(tds)
    tds = NCDataset(datadir("reanalysis/era5-TRPx0.25-slhf-sfc.nc"))
    lhf = tds["slhf"][:]*1
    close(tds)

    seb = ssr .+ str .+ shf .+ lhf

    fnc = datadir("reanalysis/era5-TRPx0.25-seb-sfc.nc")
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions"               => "CF-1.6",
        "history"                   => "2020-10-21 23:54:22 GMT by grib_to_netcdf-2.16.0: /opt/ecmwf/eccodes/bin/grib_to_netcdf -S param -o /cache/data4/adaptor.mars.internal-1603323636.0468726-6113-3-fb54412a-f0a5-4783-956d-46233705e403.nc /cache/tmp/fb54412a-f0a5-4783-956d-46233705e403-adaptor.mars.internal-1603323636.0477095-6113-1-tmp.grib",
    ))

    # Dimensions

    ds.dim["longitude"] = length(lon)
    ds.dim["latitude"]  = length(lat)

    # Declare variables

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"                     => "degrees_east",
        "long_name"                 => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"                     => "degrees_north",
        "long_name"                 => "latitude",
    ))

    scale,offset = ncoffsetscale(seb)

    varattribs = Dict(
        "scale_factor"  => scale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
        "units"         => "J m**-2",
        "long_name"     => "Surface net energy balance",
        "standard_name" => "surface_net_downward_energy_flux",
    )

    ncvar = defVar(ds,"seb",Int16,("longitude","latitude"),attrib = varattribs)

    nclon[:] = lon
    nclat[:] = lat
    ncvar[:] = seb

    close(ds)

end

function compilesfcebfeb()

    tds  = NCDataset(datadir("reanalysis/era5-TRPx0.25-sfc.nc"))
    lon  = tds["longitude"][:]*1
    lat  = tds["latitude"][:]*1
    close(tds)

    tds = NCDataset(datadir("reanalysis/era5feb-TRPx0.25-ssr-sfc.nc"))
    ssr = tds["ssr"][:]*1
    close(tds)
    tds = NCDataset(datadir("reanalysis/era5feb-TRPx0.25-str-sfc.nc"))
    str = tds["str"][:]*1
    close(tds)
    tds = NCDataset(datadir("reanalysis/era5feb-TRPx0.25-sshf-sfc.nc"))
    shf = tds["sshf"][:]*1
    close(tds)
    tds = NCDataset(datadir("reanalysis/era5feb-TRPx0.25-slhf-sfc.nc"))
    lhf = tds["slhf"][:]*1
    close(tds)

    seb = ssr .+ str .+ shf .+ lhf

    fnc = datadir("reanalysis/era5feb-TRPx0.25-seb-sfc.nc")
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions"               => "CF-1.6",
        "history"                   => "2020-10-21 23:54:22 GMT by grib_to_netcdf-2.16.0: /opt/ecmwf/eccodes/bin/grib_to_netcdf -S param -o /cache/data4/adaptor.mars.internal-1603323636.0468726-6113-3-fb54412a-f0a5-4783-956d-46233705e403.nc /cache/tmp/fb54412a-f0a5-4783-956d-46233705e403-adaptor.mars.internal-1603323636.0477095-6113-1-tmp.grib",
    ))

    # Dimensions

    ds.dim["longitude"] = length(lon)
    ds.dim["latitude"]  = length(lat)

    # Declare variables

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"                     => "degrees_east",
        "long_name"                 => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"                     => "degrees_north",
        "long_name"                 => "latitude",
    ))

    scale,offset = ncoffsetscale(seb)

    varattribs = Dict(
        "scale_factor"  => scale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
        "units"         => "J m**-2",
        "long_name"     => "Surface net energy balance",
        "standard_name" => "surface_net_downward_energy_flux",
    )

    ncvar = defVar(ds,"seb",Int16,("longitude","latitude"),attrib = varattribs)

    nclon[:] = lon
    nclat[:] = lat
    ncvar[:] = seb

    close(ds)

end

function prereanalysis(variable::AbstractString,coords::Vector{<:Real})

    ds = NCDataset(datadir("reanalysis/era5-TRPx0.25-$(variable).nc"))

    lon = ds["longitude"][:]; nlon = length(lon)
    lat = ds["latitude"][:];  nlat = length(lat)
    lvl = ds["level"][:];     nlvl = length(lvl)
    var = ds[variable][:]*1

    long = ds[variable].attrib["long_name"]
    unit = ds[variable].attrib["units"]

    close(ds)

    ds = NCDataset(datadir("reanalysis/era5-TRPx0.25-lsm-sfc.nc"))
    lsm = ds["lsm"][:]*1
    close(ds)

    lprf,sprf = getmean(coords,var,lon,lat,nlvl,lsm)

    return lvl,lprf,sprf

end
