function [beta,sigma] = BVAR(y,lags,b,PSI,lambda,draw);

%data matrix manipulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions
[TT,n]=size(y);
k=n*lags;         % # coefficients for each equation

% constructing the matrix of regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=zeros(TT,k);
for i=1:lags
    x(:,(i-1)*n+1:i*n)=lag(y,i);
end

x=x(lags+1:end,:);
y=y(lags+1:end,:);
y0=mean(y(1:lags,:),1);
[T,n]=size(y);


% starting values for the minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lambda=.2;     % std of MN prior
alpha=2;       % lag-decaying parameter of the MN prior
d=n+2;         % df for the covariance


% priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega=zeros(k,1);
for i=1:lags
    omega((i-1)*n+1:i*n)=(d-n-1)*(lambda^2)*(1/(i^alpha))./PSI;
end

%SOC prior
% miu = 3;       % std of SOC
% ydnoc=(1/miu)*diag(y0)*0;
% xdnoc=[(1/miu)*repmat(diag(y0),1,lags)]; 
% y=[y;ydnoc];
% x=[x;xdnoc];
% T=T+n;
    
    
% output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% posterior mode of the VAR coefficients
betahat=(x'*x+diag(1./omega))\(x'*y+diag(1./omega)*b);

% VAR residuals
epshat=y-x*betahat;

sigmahat=(epshat'*epshat + diag(PSI) + (betahat-b)'*diag(1./omega)*(betahat-b))/(T+d+n+1);

stationary = 0;
if draw==1
    while stationary==0
        [V E]=eig(sigmahat*(T+d+n+1));
        Sinv=V*diag(1./abs(diag(E)))*V';
        eta=mvnrnd(zeros(1,n),Sinv,T+d);
        sigma=(eta'*eta)\eye(n);
        %[cholSIGMA,junk]=chol(drawSIGMA);
        %betadraw=betahat+mvnrnd(zeros(k,1),(x'*x+diag(1./omega))\eye(k),n)'*cholSIGMA;
        cholSIGMA=cholred((sigma+sigma')/2);
        cholZZinv = cholred((x'*x+diag(1./omega))\eye(k));
        beta=betahat + cholZZinv'*randn(size(betahat))*cholSIGMA;
        AA = zeros(n*lags);
        AA(1:n,1:n*lags) = beta';
        AA(n+1:end,1:n*(lags-1))=eye(n*(lags-1));
        stationary = (sum(abs(eig(AA))>1)==0);
    end;
    
else
    beta = betahat;
    sigma = sigmahat;
end;

function z = lag(x,n,v)
% PURPOSE: creates a matrix or vector of lagged values
% -------------------------------------------------------
% USAGE: z = lag(x,n,v)
% where: x = input matrix or vector, (nobs x k)
%        n = order of lag
%        v = (optional) initial values (default=0)
% e.g.
%     z = lag(x) creates a matrix (or vector) of x, lagged 1 observations
%     z = lag(x,n) creates a matrix (or vector) of x, lagged n observations
%     z = lag(x,n,v) creates a matrix (or vector) of x, lagged n observations,
%         with initial values taking a value v.
% ------------------------------------------------------
% RETURNS: z = matrix (or vector) of lags (nobs x k)
% ------------------------------------------------------
% NOTES: if n <= 0, z = [] is returned. While you may find this
%        preverse, it is sometimes useful.
%-------------------------------------------------------
% SEE ALSO: mlag() 
%-------------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

switch(nargin)

case 1
   n = 1; v = 0;
   zt = ones(n,cols(x))*v;
   z = [ zt; trimr(x,0,n)];

case 2
   v = 0;
   if n < 1
   z = [];
   return;
   end;
   zt = ones(n,cols(x))*v;
   z = [ zt; trimr(x,0,n)];

case 3
   if n < 1
   z = [];
   return;
   end;
   zt = ones(n,cols(x))*v;
   z = [ zt; trimr(x,0,n)];

otherwise
error('lag: wrong # of input arguments');
end;

function c = cols(x)
% PURPOSE: return columns in a matrix x
% -----------------------------------------
% USAGE: c = cols(x)
% where: x = input matrix
% -----------------------------------------
% RETURNS: c = # of columns in x
% -----------------------------------------

  [r,c] = size(x);

  function z = trimr(x,n1,n2)
% PURPOSE: return a matrix (or vector) x stripped of the specified rows.
% -----------------------------------------------------
% USAGE: z = trimr(x,n1,n2)
% where: x = input matrix (or vector) (n x k)
%       n1 = first n1 rows to strip
%       n2 = last  n2 rows to strip
% NOTE: modeled after Gauss trimr function
% -----------------------------------------------------
% RETURNS: z = x(n1+1:n-n2,:)
% -----------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com

  [n junk] = size(x);
  if (n1+n2) >= n; 
     error('Attempting to trim too much in trimr');
  end;
  h1 = n1+1;   
  h2 = n-n2;
  z = x(h1:h2,:);
  
  
function C = cholred(S)
[v,d] = eig(S);
dd = diag(d);
dd(dd<1e-6)=0;
C = diag(sqrt(dd))*v';
