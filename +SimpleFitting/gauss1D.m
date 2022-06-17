function [FitPar,Fit,resNorm]=gauss1D(A,domain,guess)

switch nargin
    case 2 
        sigGuess = abs((abs(domain(2))-abs(domain(1))))*3;

        [val,idx] = max(A);
        muGuess = domain(idx);
    case 3
        sigGuess = guess.sig;
        muGuess  = guess.mu;
        minMaxDom = guess.minMaxDomain;
    otherwise
        error('wrong number of arguments');
end

A=A(:);
domain =domain(:);
[val,ind] = max(A);

%                   Sigma                       mu             A        y0           
lb        = [1    min(minMaxDom)     0        0];
ub        = [50     max(minMaxDom)     3*val      val];
initguess = [      sigGuess                  muGuess                          val-min(A) min(A)];
opts = optimset('Display','off');
[FitPar,resNorm] = lsqcurvefit(@SimpleFitting.gaussian,initguess,domain,A,lb,ub,opts);

Fit=FitPar(3)*exp(-((domain-FitPar(2))./(sqrt(2).*FitPar(1))).^2)+FitPar(4);
