function [rho]=fsippsolve_primal(prob,para,order)
%Compute the primal SDP relaxation (P_k) implemented by Yalmip
%Input:
%   prob: the problem structure
%   para: the parameter structure
%   order: the relaxation order
%Output:
%   optval: the optimal value of the order-th SDP relaxation (P_k)
%Please refer the paper:
%Feng Guo and Meijun Zhang, AN SDP METHOD FOR FRACTIONAL SEMI-INFINITE
%PROGRAMMING PROBLEMS WITH SOS-CONVEX POLYNOMIALS, arXiv:2110.04737
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[coef_X,term_Y]=coefficients(prob.p,prob.Y);
s=length(term_Y);


term_Y_expo=[];
for i=1:s
    term_Y_expo=[term_Y_expo; degree(term_Y(i),prob.Y)];
end

[degs, start_index] = deglist(prob.Ynum, 0, order);
r=size(degs,1);

%compute Lp which represents H(p(x,y)) in (P_k)
P=sdpvar(r,r);
Lp=0;
for i=1:r
    for j=1:i
        temp=0;
        for k=1:s
            switch prob.indset
                case 'hypercube'
                    temp=temp+coef_X(k)*intC(term_Y_expo(k,:)+degs(i,:)+degs(j,:));
                case 'sphere'
                    temp2=intSB(term_Y_expo(k,:)+degs(i,:)+degs(j,:));
                    temp=temp+coef_X(k)*temp2;
                case 'ball'
                    [~,temp2]=intSB(term_Y_expo(k,:)+degs(i,:)+degs(j,:));
                    temp=temp+coef_X(k)*temp2;
                case 'polytope'
                    temp=temp+coef_X(k)*int_mono(prob.A, prob.b, term_Y_expo(k,:)+degs(i,:)+degs(j,:),prob.latte_dir);
            end
        end
        Lp=Lp+P(i,j)*temp;
        if i~=j
            Lp=Lp+P(j,i)*temp;
        end
    end
end

d=max([degree(prob.f), degree(prob.g), degree(replace(prob.p,prob.Y,ones(size(prob.Y)))), max(degree(prob.phi))]);
d=ceil(d/2);

[s1,c1] = polynomial(prob.X,2*d-2);

if ~isempty(prob.phi)
    eta=sdpvar(max(size(prob.phi)),1);
end

sdpvar rho 

[s1,c1] = polynomial(prob.X,2*d-2);

if strcmp(class(prob.g),'double')
    if ~isempty(prob.phi)
        F=[sos(prob.f/prob.g-rho+Lp+prob.phi*eta-s1*((prob.R)^2-sum((prob.X).^2))), sos(s1), eta>=0, P>=0];
    else
        F=[sos(prob.f/prob.g-rho+Lp-s1*((prob.R)^2-sum((prob.X).^2))), sos(s1), P>=0];
    end
    
    ops = sdpsettings('solver',para.solver);
    ops = sdpsettings(ops,'verbose',0);
    solvesos(F,-rho,ops,[c1]);
    
    disp(['The optimal value r^primal_k of the ', num2str(order), '-th primal SDP relaxation (P_k) is ']);
    disp([num2str(value(rho))]);

else 
    [s2,c2] = polynomial(prob.X,2*d-degree(prob.g));
    if ~isempty(prob.phi)
        F=[sos(prob.f-rho*prob.g+Lp+prob.phi*eta-s1*((prob.R)^2-sum((prob.X).^2))...
            -s2*(prob.g-prob.gstar)), sos(s1),sos(s2), eta>=0, P>=0];
    else
        F=[sos(prob.f-rho*prob.g+Lp-s1*((prob.R)^2-sum((prob.X).^2))...
            -s2*(prob.g-prob.gstar)), sos(s1),sos(s2), P>=0];
    end   
    
    ops = sdpsettings('solver',para.solver);
    ops = sdpsettings(ops,'verbose',0);
    solvesos(F,-rho,ops,[c1;c2]);
    
    disp(['The optimal value r^primal_k of the ', num2str(order), '-th primal SDP relaxation (P_k) is ']);
    disp([num2str(value(rho))]);
end


end









