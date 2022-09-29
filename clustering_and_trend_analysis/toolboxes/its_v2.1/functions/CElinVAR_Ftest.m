%% Statistics of difference in conditional entropies estimated through linear regression
% modified 22/8/2018 to get as input the variance of the residuals

% Vu: variance of residuals of unrestricted regression
% Vr: variance of residuals of restricted regression
% Nu: number of coefficients for unrestricted regression
% Nr: number of coefficients for restricted regression
% N: data length

function [p_value] = CElinVAR_Ftest(Vu,Vr,Nu,Nr,N)

    % F-statistic for significance
    RSSu=Vu*(N-Nu+1); %inferred residual sum of squares (unrestricted)
    RSSr=Vr*(N-Nr+1); %inferred residual sum of squares (restricted)
    %RSSu=sum(Upu.^2);
    %RSSr=sum(Upr.^2);
    n1=Nu-Nr; %number of restrictions
    % n2=(length(Upu)-size(Vu,1)); % n. of observations - n. tot of coeffs
    n2=N-Nu; % n. of observations - n. tot of coeffs
    %n2=length(Upr)-Nu; % n. of observations - n. tot of coeffs
    f_value=((RSSr-RSSu)/n1)/(RSSu/n2);
    p_value=1-cdff_fun(f_value,n1,n2); % function of Seth Toolbox

end


%% function of Seth Toolbox: CDF computes the cumulative 'F' distribution
function p = cdff_fun(x,v1,v2)
    p = 0;
    t = (v1 <= 0 | v2 <= 0 | isnan(x) | isnan(v1) | isnan(v2));
    p(t) = NaN;
    s = (x==Inf) & ~t;
    if any(s)
       p(s) = 1;
       t = t | s;
    end

    % Compute P when X > 0.
    k = find(x > 0 & ~t & isfinite(v1) & isfinite(v2));
    if any(k), 
        xx = x(k)./(x(k) + v2(k)./v1(k));
        p(k) = betainc(xx, v1(k)/2, v2(k)/2);
    end

    if any(~isfinite(v1(:)) | ~isfinite(v2(:)))
       k = find(x > 0 & ~t & isfinite(v1) & ~isfinite(v2) & v2>0);
       if any(k)
          p(k) = chi2cdf(v1(k).*x(k),v1(k));
       end
       k = find(x > 0 & ~t & ~isfinite(v1) & v1>0 & isfinite(v2));
       if any(k)
          p(k) = 1 - chi2cdf(v2(k)./x(k),v2(k));
       end
       k = find(x > 0 & ~t & ~isfinite(v1) & v1>0 & ~isfinite(v2) & v2>0);
       if any(k)
          p(k) = (x(k)>=1);
       end
    end
end