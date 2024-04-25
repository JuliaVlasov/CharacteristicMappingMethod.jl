function  upsampledImage = upsample(image,upsamplingFactor)
% Get the size of the image
[rows, cols] = size(image);

% Calculate the new dimensions
new_rows = upsamplingFactor * rows;
new_cols = upsamplingFactor * cols;

% Perform Fourier transform
F = fftshift(fft2(image));

% Upsample the frequency domain
F_upsampled = padarray(F, [(new_rows - rows)/2, (new_cols - cols)/2], 0, 'both');

% Perform inverse Fourier transform
upsampledImage = real(ifft2(ifftshift(F_upsampled)))*upsamplingFactor^2;


