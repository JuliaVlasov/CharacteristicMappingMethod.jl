using Parameters

@with_kw struct Grid

    name::Symbol
    size::Vector{Int}
    L::Vector{Float64}
    dx::Float64
    dv::Float64
    x::Vector{Float64}
    v::Vector{Float64}
    kx::Vector{ComplexF64}
    kv::Vector{ComplexF64}
    v_periodic::Vector{Float64}
    sigma::Vector{Float64}

end

function set_grids(params, grid_size_list, grid_names)

    grid_list = Dict{Symbol,Grid}()
    for ig in eachindex(grid_size_list)

        name = Symbol(grid_names[ig])
        grid_size = grid_size_list[ig]
        x, y = genMeshPoints(params.dom, grid_size)
        kx, kv = genWaveNumbers(params.dom, grid_size)

        grid = Grid(
            name = name,
            size = grid_size,
            L = params.L,
            dx = params.L[1] / grid_size[1],
            dv = params.L[2] / grid_size[2],
            x = collect(x),
            v = collect(y) .- params.Lv,
            kx = kx, 
            kv = kv,
            v_periodic = similar(y),
            sigma = similar(y),
        )

        velocity_periodicfication!(grid; bump_transition_width = params.bump_transition_width)

        grid_list[name] = grid
    end

    grid_list

end
