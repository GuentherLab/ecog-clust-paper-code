function y = erc_test_signal
rng(0);

shift = 6;
fs = 1000;
T = 4;
y = zeros(3, T*fs); 

r1 = 86:96;
r2 = 105:114;

t = linspace(0, 1, fs);

s1 = zeros(1, fs);
for i = 86:96
    s1 = s1 + sin(i*2*pi*t);
end
s1 = s1 ./ max(s1);

s2 = zeros(1, fs);
for i = 105:114
    s2 = s2 + sin(i*2*pi*t);
end
s2 = s2 ./ max(s2);

function [mu, sigma] = noise_params
    mu = 0.01 * randn;
    sigma = 0.5 + 0.1*randn;
end


[mu, sigma] = noise_params;
y(1, 1:fs) = s1 + (mu + sigma*randn(size(s1)));
y(2, 1:fs) = s1 + (mu + sigma*randn(size(s1)));
y(3, 1:fs) = s1 + (mu + sigma*randn(size(s1)));

y(1, fs+1:2*fs) = y(1, 1:fs);
[mu, sigma] = noise_params;
y(2, fs+1:2*fs) = (mu + sigma*randn(size(s1))) + circshift(s1, shift);
[mu, sigma] = noise_params;
y(3, fs+1:2*fs) = (mu + sigma*randn(size(s1))) + circshift(s1, 2*shift);

[mu, sigma] = noise_params;
y(1, 2*fs+1:3*fs) = s2 + (mu + sigma*randn(size(s2)));
y(2, 2*fs+1:3*fs) = s2 + (mu + sigma*randn(size(s2)));
y(3, 2*fs+1:3*fs) = s2 + (mu + sigma*randn(size(s2)));

[mu, sigma] = noise_params;
y(1, 3*fs+1:4*fs) = (mu + sigma*randn(size(s2))) + circshift(s2, 2*shift);
[mu, sigma] = noise_params;
y(2, 3*fs+1:4*fs) = (mu + sigma*randn(size(s2))) + circshift(s2, shift);
y(3, 3*fs+1:4*fs) = y(3, 2*fs+1:3*fs);

end