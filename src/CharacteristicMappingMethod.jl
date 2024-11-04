module CharacteristicMappingMethod

export vlasov_main

include(joinpath("physics", "inicond.jl"))
include(joinpath("lib_cmm", "genMeshPoints.jl"))
include(joinpath("lib_cmm", "genWaveNumbers.jl"))
include(joinpath("lib_cmm", "velocity_periodicfication.jl"))
include(joinpath("lib_cmm", "set_grids.jl"))
include(joinpath("lib_cmm", "vlasov_CMM.jl"))

end # module CharacteristicMappingMethod
