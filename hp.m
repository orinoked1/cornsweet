function filter_matrix = hp(order, cutoff)
% hp --> Creates a high pass filter
%	 The filter has size (order x order)
%	 and cutoff specified between 0 and 1

%Create desired frequency response
[f1,f2] = freqspace(order,'meshgrid');
d = find(f1.^2+f2.^2 < cutoff^2);
Hd = zeros(order);
Hd(d) = ones(size(d));
%This inverts the frequency response making it a HPF
Hd = 1 - Hd;

% Design the filter's impulse response
filter_matrix = fsamp2(Hd);