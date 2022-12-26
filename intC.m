function [intc] = intC(s)
%Given a exponent list s of a monomial, compute its integral intc on the
%hypercube [-1,1]^n


n=size(s);
n=max(n(1),n(2));

intc=1;

for i=1:n
    intc=intc*((1^(s(i)+1)-(-1)^(s(i)+1))/(s(i)+1));
end

end

