#########################################################################
# main script for CMM
#########################################################################

using CharacteristicMappingMethod
################################################
## select case:
#PARAMS_two_stream;
include(joinpath("params", "PARAMS_landau_damping.jl"));
#PARAMS_non_lin_landau_damping;
################################################

## simulate
@time params1,f1,XB1,YB1 = vlasov_CMM(params);


