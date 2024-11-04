using Statistics

function velocity_periodicfication(params; bump_transition_width = nothing, type = "exp")

    b = isnothing(bump_transition_width) ? 0.2 * params.L[2] : bump_transition_width

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
        sigma = 0.5 - 0.5 * tanh(2 * pi * (abs(params.v) - a * params.L(2)) / (b))
    elseif type == "exp"
        nu = x -> exp(-1 / (1 - abs(x)^2)) * (abs(x) < 1)
        Lv = params.L[2] / 2
        b = 0.2 * Lv
        @. sigma = 1 - nu((abs(params.v) - Lv) / b)
    else
        @error "type: $type not known."
    end

    h = sigma#tanh((-abs(params.v)+a*params.L(2)/2)/(b))+0.5;
    h = h - mean(h)
    h = h ./ maximum(h)

    kv = params.kv
    kv2 = kv .^ 2
    kv2[1] = 1
    intu_hat = (-fft(h) ./ kv2)
    intu = ifft(intu_hat)
    v_periodic = real(ifft(intu_hat .* kv))

    v_periodic, sigma

end
