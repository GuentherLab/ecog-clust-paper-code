%% Statistics of difference in conditional entropies estimated through linear regression
%%% modification including the variances of the two models instead of the residuals (prediction errors)

% Var_u: variance of residuals of unrestricted regression
% Var_r: variance of residuals of restricted regression
% Nu: number of coefficients for unrestricted regression
% Nr: number of coefficients for restricted regression
% N: numer of data points from which parameters of both regressions are estimated

function [p_value] = its_LinReg_Ftest_var(Var_u,Var_r,Nu,Nr,N)
    % F-statistic for significance
    n1=Nu-Nr; %number of restrictions
    % n2=(length(Upu)-size(Vu,1)); % n. of observations - n. tot of coeffs
    n2=N-Nu; % n. of observations - n. tot of coeffs
    f_value=((Var_r-Var_u)/n1)/(Var_u/n2);
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