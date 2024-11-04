using Documenter
using CharacteristicMappingMethod

makedocs(
    sitename = "CharacteristicMappingMethod",
    format = Documenter.HTML(),
    modules = [CharacteristicMappingMethod],
)

deploydocs(
    repo = "github.com/juliavlasov/CharacteristicMappingMethod.jl",
    devbranch = "develop",
)
