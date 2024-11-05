using FFTW

function genWaveNumbers(dom, gsz)

    MX = gsz[2]
    MY = gsz[1]

    kxn = fftfreq(MX, MX)
    kyn = fftfreq(MY, MY)

    rx = 2pi * 1im ./ dom[3]
    ry = 2pi * 1im ./ dom[4]

    kx = kxn .* rx
    ky = kyn .* ry

    kn = sqrt.(abs.(kx) .^ 2 .+ abs.(ky) .^ 2)

    kx, ky, kn, rx, ry, kxn, kyn

end
