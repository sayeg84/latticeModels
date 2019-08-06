include("structs.jl")
module test1
    import ..SpinLattice
    import ..RandomPosition
    function test(x::SpinLattice)
        pos = RandomPosition(x)
        return pos
    end
end
include("statEnsemble.jl")
StatEnsemble.NormalEnergy