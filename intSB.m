function [ints,intb] = intSB(s)
%Given a exponent list s of a monomial, compute its integral on the unit
%sphere and the unit ball. 
%Output: 
%   ints: integral on the unit sphere
%   intb: integral on the unit ball

n=size(s);
n=max(n(1),n(2));

for i=1:n
    if rem(s(i),2)==1
        ints=0;
        intb=0;
        return;
    end
end

ints=2;

for i=1:n
    ints=ints*gamma((s(i)+1)/2);
end

ints=ints/(gamma((sum(s)+n)/2));

R=1;


intb=ints/(sum(s)+n);

intb=intb*R^(sum(s)+n);

end

