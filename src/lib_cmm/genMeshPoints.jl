function genMeshPoints(dom, gsz)

    x = LinRange(dom[1], dom[3], gsz[2] + 1)[1:end-1]
    y = LinRange(dom[2], dom[4], gsz[1] + 1)[1:end-1]

    x, y

end
