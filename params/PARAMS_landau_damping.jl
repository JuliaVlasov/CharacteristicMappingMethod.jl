using Parameters

@with_kw struct Params

    n :: Int = 2^5
    Nmap :: Int = n 
    Nfine :: Int = 2^9 
    Nsampling :: Int = 2^8 
    Nplotting :: Int = 1024 
    nv :: Int = n
    Lv :: Float64 = 4 * pi
    detTol :: Float64 = 1e-2 
    l :: Float64 = 1
    k :: Float64 = 0.5
    eps :: Float64 = 5e-2
    v0 :: Float64 = 3
    Lx :: Float64 = 2pi / k
    L :: Vector{Float64} = [Lx, Lv * 2]   
    dom :: Vector{Float64} = [0, 0, Lx, 2 * Lv]    
    CFL :: Float64 = 0.8
    T_end :: Float64 = 80
    dt :: Float64 = 0.001 / 10
    iplot :: Float64 = 100 
    ihist :: Float64 = 100 
    ilog :: Float64 = 10 
    dt_hist :: Float64 = 1
    dt_log :: Float64 = 0.1
    bump_transition_width :: Float64 = 0.1 * Lv
    case :: String = "landau_damping"
    filter :: String = "gauss"

end

params = Params()
