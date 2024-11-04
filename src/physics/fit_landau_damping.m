    function [damping_rate, frequency, amplitude, peaks, locations, fit_line_fun, fit_damping_fun] = fit_landau_damping(t,Epot,min_peak_threshold)

% find peaks to fit
if nargin<3 || isempty(min_peak_threshold)
    [peaks, locations] = findpeaks(Epot);
else
    [peaks, locations] = findpeaks(Epot,"Threshold",min_peak_threshold);
end

% Define the fitting model function
log_E_fun = @(A, t) A(1) * t + A(2);

% Initial guess for parameters
initialGuess = [-1, 1]; % Damping Rate, Amplitude

% Fit the model to the data
fittedParams = lsqcurvefit(log_E_fun, initialGuess, t(locations), log(peaks));

% Extract the fitted parameters
amplitude = fittedParams(2);
damping_rate = fittedParams(1)/2; % 1/2 because of the square when interating over electrical field
frequency = 2*pi/(2*mean(diff(t(locations)))); %2pi/periodendauer
shift = + frequency*t(locations(1)); % define shift such that first maxima intersects with the peak: cos(freq *t_peak + shift) = 1!

% Generate the fitted curve using the fitted parameters
fit_line_fun = @(time) exp(log_E_fun(fittedParams, time));
fit_damping_fun = @(time) exp(log_E_fun(fittedParams, time)).*cos(frequency*time-shift).^2;

end