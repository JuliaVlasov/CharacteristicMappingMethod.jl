function uk = distrib2u(distribution)
    global params
    

    stream = poisson(distribution);
    uk(:,:,1) = -1i*params.Ky.*stream;
    uk(:,:,2) = +1i*params.Kx.*stream;
end