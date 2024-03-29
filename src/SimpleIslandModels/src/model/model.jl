struct SimpleIslandModel{FT<:Real}
     S0 :: FT
      τ :: FT
     Fo :: FT
    sfc :: Surface
    atm :: Atmosphere
    cld :: CloudAlbedo
    do_wtg     :: Bool
    do_diurnal :: Bool
    do_ocnflux :: Bool
    do_cloud   :: Bool
    cloudscheme :: Int
end

function CreateModel(
    sfc :: Surface,
    atm :: Atmosphere,
    cld :: CloudAlbedo = CloudAlbedo(dummycalculateαₐ,DummyCloud()),
    FT = Float64;
    S0  :: Real = 1361.,
    τ   :: Real = 0.,
    Fo  :: Real = 0.,
    do_wtg     :: Bool = false,
    do_diurnal :: Bool = true,
    do_ocnflux :: Bool = false,
    do_cloud   :: Bool = false,
    cloudscheme :: Int = 1
)

    if do_wtg && iszero(τ)
        error("$(modulelog()) - Weak Temperature Gradient parameterization called, but relaxation timescale τ not specified")
    end

    if !do_diurnal
        S0 = S0/π
    end

    return SimpleIslandModel{FT}(
        S0,τ,Fo,
        sfc,atm,cld,
        do_wtg,do_diurnal,do_ocnflux,do_cloud,cloudscheme
    )

end

function show(io::IO, model::SimpleIslandModel)
    print(
		io,
		"SimpleIslandModel:\n",
		" ├─── Surface Model                (sfc) : ", typeof(model.sfc), '\n',
		" ├─── Atmophere Model              (atm) : ", typeof(model.atm), '\n',
		" ├─── Cloud Model                  (cld) : ", typeof(model.cld), '\n',
		" ├─── Insolation Constant           (S0) : ", model.S0,          '\n',
		" ├─── Diurnal Insolation    (do_diurnal) : ", model.do_diurnal,  '\n',
		" ├─── Weak Temperature Gradient (do_wtg) : ", model.do_wtg,      '\n',
        " ├─── Relaxation Timescale / sec     (τ) : ", model.τ,           '\n',
		" ├─── Do Ocean Transport    (do_ocnflux) : ", model.do_ocnflux,  '\n',
        " ├─── Ocean Flux / W m**-2 s**-1    (Fo) : ", model.Fo,          '\n',
		" └─── Do Cloud                (do_cloud) : ", model.do_cloud,    '\n',
	)
end