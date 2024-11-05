function inicond(params)

    k = params.k
    Lv = params.L[2]
    eps = params.eps

    if params.case == "landau_damping"
        f0 = (x, v) -> (1 + eps * cos(x * k)) / sqrt(2 * pi) * exp(-(v - Lv / 2)^2 / 2)
    elseif params.case == "two_stream"
        v0 = params.v0
        f0 = (x, v) -> (1 + eps * cos(k * x)) / (2*sqrt(2*pi)) * (exp.(-(v-Lv/2-v0)^2/2)+exp(-(v-Lv/2+v0)^2/2));
    else
        @error "Case: $(params.case) does not exist!"
    end

    f0

end
