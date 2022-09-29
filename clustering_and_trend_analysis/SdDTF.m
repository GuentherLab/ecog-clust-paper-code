function [Z, CHI] = SdDTF(S, H)

% partial coherence
[n1, n2, f] = size(S);
CHI = zeros(n1, n2, f);
for k = 1:f
    C = inv(S(:, :, k));
    d = diag(C);
    CHI(:, :, k) = C ./ sqrt(d * d.'); 
end

h = abs(H);
chi = abs(CHI);

Z = (h .* chi) ./ sqrt(sum(h.^2 .* chi.^2, 'all'));

end