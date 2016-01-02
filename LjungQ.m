function [ q ] = LjungQ( x, maxlag )
%When fed a 1xn vector and max lag,
%Qljung will produce Q-stats for each lag

%   Pj's must be in lag order 1 to n

pjq=autocorrel(x, maxlag);
n=length(x);
summer=0;
for j=1:maxlag
for k=1:j
    summer=summer+pjq(k)^2/(n-k);
end;
    q(j)=n*(n+2)*summer;
    summer=0;
end;
end

