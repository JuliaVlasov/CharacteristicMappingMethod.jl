using CharacteristicMappingMethod
using Test

@testset "Set the parameters" begin

    include(joinpath("..","params","PARAMS_landau_damping.jl"))

    @test true

end
