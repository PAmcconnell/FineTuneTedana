function [a,f] = calc_ps(data,TR)
% cjl2007@med.cornell.edu; 
% modelled after MELODIC 

a = abs(fft(data)).^2; 
a(1)=[]; 
b = round(length(a)/2);
a = a(1:b);

% sweep data
for i = 1:length(a)
    f(i)= 1 / ((TR*length(data)) / i);
end

end