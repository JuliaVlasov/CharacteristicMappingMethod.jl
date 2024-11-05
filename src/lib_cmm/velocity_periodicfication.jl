using Statistics

function velocity_periodicfication!(grid; bump_transition_width = nothing, type = "exp")

    b = isnothing(bump_transition_width) ? 0.2 * grid.L[2] : bump_transition_width

    ###########################################################################
    ### first calculate "vorticity" such that global vorticity omega = 0 
    # (u_1,u_2) = (g(v), d/dx \phi)
    # omega = d/dx u_2 - d/dv u_1 = d^2/(dx)2 \phi(x) - d/dv g(v) == 0 
    # assume phi(x)=0 
    # if g(v) = v 
    # omega = d/dv g(v) = 1 != 0 ):
    # so we do the following:
    # look for a function h(v) = d/dv g(v) that keeps global omega=0, 
    # but is 1 close to the
    # origin and negative on the boundaries (= heaviside function)
    ###########################################################################

    if type == "tanh"
        a = 0.4
        @. grid.sigma = 0.5 - 0.5 * tanh(2 * pi * (abs(grid.v) - a * grid.L(2)) / (b))
    elseif type == "exp"
        nu = x -> exp(-1 / (1 - abs(x)^2)) * (abs(x) < 1)
        Lv = grid.L[2] / 2
        b = 0.2 * Lv
        @. grid.sigma = 1 - nu((abs(grid.v) - Lv) / b)
    else
        @error "type: $type not known."
    end

    h = copy(grid.sigma) #tanh((-abs(grid.v)+a*grid.L(2)/2)/(b))+0.5;
    h .= h .- mean(h)
    h .= h ./ maximum(h)

    kv2 = grid.kv .^ 2
    kv2[begin] = 1
    intu_hat = (-fft(h) ./ kv2)
    intu = ifft(intu_hat)
    grid.v_periodic .= real(ifft(intu_hat .* grid.kv))

end
