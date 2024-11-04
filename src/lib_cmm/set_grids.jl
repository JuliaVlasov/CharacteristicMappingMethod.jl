struct Grid

    name::String
    size::Int
    L::Float64
    dx::Float64
    dv::Float64
    x::Vector{Float64}
    v::Vector{Float64}
    kx::Vector{Float64}
    kv::Vector{Float64}
    v_periodic::Vector{Float64}
    sigma::Vector{Float64}

end

function set_grids(params, grid_size_list, grid_names)

    grid_list = Dict{Symbol,Grid}()
    for ig in eachindex(grid_size_list)

        x, y = genMeshPoints(params.dom, grid_size)
        kx, kv = genWaveNumbers(params.dom, grid_size)
        v_periodic, sigma = velocity_periodicfication(grid, params.bump_transition_width)

        grid = Grid(
            name = grid_names[ig],
            size = grid_size_list[ig],
            L = params.L,
            dx = params.L[1] / grid_size[1],
            dv = params.L[2] / grid_size[2],
            x = x,
            v = v .- params.Lv,
            v_periodic = v_periodic,
            sigma = sigma,
        )
        grid_list[name] = grid
    end

    grid_list

end
