using DrWatson
@quickactivate "TroPrecLS"
using Printf

exp = "DiAmp064km"
tprm = projectdir("exp","tmp.prm")

slabvec = [ 0.1, 0.1*sqrt(2), 0.2, 0.2*sqrt(2.5), 0.5, 0.5*sqrt(2) ]
slabvec = vcat(slabvec,slabvec*10,slabvec*100)


for slabii in slabvec
    slabstr = "slab$(@sprintf("%04.1f",slabii))"
    slabstr = replace(slabstr,"."=>"d")
    slabprnt = @sprintf("%04.3f",slabii)
    for imember = 1 : 5
        prm = projectdir("exp","prm",exp,slabstr,"member$(imember).prm")
        open(tprm,"w") do fprm
            open(prm,"r") do oprm
                s = read(oprm,String)
                s = replace(s,"member1"=>"member$(imember)")
                s = replace(s,"nensemble     = 1,"=>"nensemble     = $(imember),")
                write(fprm,s)
            end
        end
        mv(tprm,prm,force=true)
    end
end
