function [KL] = KLdiv(q,p, basis)
% KULLBACK Kullback-Leibler Divergence between discrete pdf's
% 
% [KL] = KULLBACK(Q, P, <BASIS>)
% 
% Inputs :
% Q, P : Target and Model distributions
% BASIS : Basis for log(.), <Default : e>
%
% Outputs :
% KL(Q||P)
% 
% Usage Example : KL = kullback([0.1 0.9],[0.5 0.5]);

if nargin<3, 
  C = 1; 
else
  C = 1/log(basis);
end;

q = q(:);
p = p(:);

zq = find(q>0);
zp = find(p>0);

if ~isempty(setdiff(zq,zp)),
  KL = Inf;
else
  KL = q(zq)'*log(q(zq)) - q(zp)'*log(p(zp));
end;

KL = C*KL;

end