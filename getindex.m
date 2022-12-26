function [index] = getindex(l)
%Given a monomial x^\alpha, compute the index of the corresponding moment
%Input: 
%   l: the exponent alpha in the monomial x^\alpha
%Ouput:
%   index: the index (position) of the corresponding moment in the moment
%   sequence

s=size(l);
if s(1)>s(2)
    l=l';
    n=s(1);
else
    n=s(2);
end

if l==zeros(n)
    index=1;
    return;
end


if min(size(l))~=1
    disp('invalid input');
end

dd=sum(l);
s1=nchoosek(n+dd-1,dd-1);
s2=nchoosek(n+dd,dd);

index2=0;

for i=n:-1:2
    di=dd-sum(l(i+1:n));
    index2=index2+nchoosek(i-1+di,i-1)-nchoosek(i-1+di-l(i),i-1);
end


index=s2-index2;


end

