

function [k, EK] = spectrum(image)
    u = image;
    % Calculate the 2D Fourier Transform of the image
    spectrum = fftshift(fft2(image));

    % Get the dimensions of the image
    [M, N] = size(image);

    % Create a grid of frequency coordinates
    [kx, ky] = meshgrid(-floor(N/2):floor((N-1)/2), -floor(M/2):floor((M-1)/2));

    % Calculate the radial distance from the center
    r = sqrt(kx.^2 + ky.^2);

    % Calculate the radial power spectrum
    EK = accumarray(round(r(:))+1, abs(spectrum(:)).^2, [], @sum);

    % Plot the radial power spectrum
    k = 0:numel(EK)-1;
end
