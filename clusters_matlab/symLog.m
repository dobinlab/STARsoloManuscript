function s=symLog(x, rev)

if nargin==2 && rev
    s=sign(x).* (2 .^(abs(x))-1); % reverse
else
    s=sign(x).*log2(1+abs(x));
end