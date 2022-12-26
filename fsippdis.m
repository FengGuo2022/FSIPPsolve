function [fval] = fsippdis(r,N)
%solve Problem (23) in Example 5.2 by discretization method (10) with the
%regular grid (12)
%Input:
%   r: the random vector a in Problem (23)
%   N: the number N in the regular grid (12)
%Output
%   fval: the optimal value of Problem (23) in Example 5.2 by discretization method (10) with the regular grid (12)
%Please refer the paper:
%Feng Guo and Meijun Zhang, AN SDP METHOD FOR FRACTIONAL SEMI-INFINITE
%PROGRAMMING PROBLEMS WITH SOS-CONVEX POLYNOMIALS, arXiv:2110.04737
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=max(size(r));

s='[X1';
for i=2:n
    s=strcat(s,strcat(', X',num2str(i)));
    s=strcat(s,{32});
end
s=strcat(s,']');
cmd=[strcat(s,'=ndgrid(-1:2/N:1);')];
eval(char(cmd));
X=[];
for i=1:n
    cmd=[strcat('X=[X, vec(', strcat('X',num2str(i)), ')];')];
    eval(cmd);
end

n1=size(X,1);
for i=1:n 
    X(:,i)=ones(n1,1)-(X(:,i)-r(i)*ones(n1,1)).^2/4;
end

options = optimoptions('fmincon','MaxFunctionEvaluations',5000);

[x,fval,exitflag]=fmincon(@(x)sum((x-1).^4)/(sum(x)+1), zeros(n,1), [], [], [], [], [], [], @(x)cons(x,X),options); 

disp(['The optimal value by discretization method with N=', num2str(N), ' is ']);
disp([num2str(fval)]);

end

function [c, ceq] = cons(x,X)

c=[-ones(size(X,1),1)+X*(x.^2);-sum(x)];
ceq=[];
end

