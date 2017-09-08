function p = GamInvPdf(x,alpha,beta)

p = gamma(alpha)*(beta^alpha)*x.^(-alpha-1).*exp(-beta./x);