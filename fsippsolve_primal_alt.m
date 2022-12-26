function [optval]=fsippsolve_primal_alt(prob,para,order)
%Compute the optimal value of the relaxation using SOS/SDSOS/DSOS by the
%software spotless_isos
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


x = msspoly('x',prob.Xnum);
prog = spotsosprog;
prog = prog.withIndeterminate(x);

coeffs=[];
for i=1:s
    coeffs=[coeffs; sdpvar2msspoly(coef_X(i),prob.X,x)];
end


term_Y_expo=[];
for i=1:s
    term_Y_expo=[term_Y_expo; degree(term_Y(i),prob.Y)];
end


[degs, start_index] = deglist(prob.Ynum, 0, order);
r=size(degs,1);


[prog,rho] = prog.newFree(1);
[prog,Q] = prog.newPSD(r);


Lp=0;
for i=1:r
    for j=1:i
        temp=0;
        for k=1:s
            switch prob.indset
                case 'hypercube'
                    temp=temp+coeffs(k)*intC(term_Y_expo(k,:)+degs(i,:)+degs(j,:));
                case 'sphere'
                    temp2=intSB(term_Y_expo(k,:)+degs(i,:)+degs(j,:));
                    temp=temp+coeffs(k)*temp2;
                case 'ball'
                    [~,temp2]=intSB(term_Y_expo(k,:)+degs(i,:)+degs(j,:));
                    temp=temp+coeffs(k)*temp2;
                case 'polytope'
                    temp=temp+coeffs(k)*int_mono(prob.A, prob.b, term_Y_expo(k,:)+degs(i,:)+degs(j,:),prob.latte_dir);
            end
        end
        Lp=Lp+Q(i,j)*temp;
        if i~=j
            Lp=Lp+Q(j,i)*temp;
        end
    end
end


f = sdpvar2msspoly(prob.f,prob.X,x);
g = sdpvar2msspoly(prob.g,prob.X,x);
phi=[];

if ~isempty(prob.phi)
    for i=1:max(size(prob.phi))
        phi=[phi, sdpvar2msspoly(prob.phi(i),prob.X,x)];
    end
    s=max(size(phi));
   [prog,eta] = prog.newFree(s);
    prog=prog.withPos(eta);
end


d=max([degree(prob.f), degree(prob.g), degree(replace(prob.p,prob.Y,ones(size(prob.Y)))), max(degree(prob.phi))]);
d=ceil(d/2);


vx = monomials(x,0:2*(d-1));
[prog,coeff] = prog.newFree(length(vx));
P=coeff'*vx;

switch para.relaxtype
    case 'sos'
        prog=prog.withSOS(P);
    case 'sdsos'
        prog=prog.withSDSOS(P);
    case 'dsos'
        prog=prog.withDSOS(P);
end
    


if strcmp(class(prob.g),'double')
    if ~isempty(prob.phi)
        switch para.relaxtype
            case 'sos'
                prog=prog.withSOS((1/prob.g)*f-rho+Lp+phi*eta-P*(prob.R^2-sum(x.^2)));
            case 'sdsos'
                prog=prog.withSDSOS((1/prob.g)*f-rho+Lp+phi*eta-P*(prob.R^2-sum(x.^2)));
            case 'dsos'
                prog=prog.withDSOS((1/prob.g)*f-rho+Lp+phi*eta-P*(prob.R^2-sum(x.^2)));
        end
    else
        switch para.relaxtype
            case 'sos'
                prog=prog.withSOS((1/prob.g)*f-rho+Lp-P*(prob.R^2-sum(x.^2)));
            case 'sdsos'
                prog=prog.withSDSOS((1/prob.g)*f-rho+Lp-P*(prob.R^2-sum(x.^2)));
            case 'dsos'
                prog=prog.withDSOS((1/prob.g)*f-rho+Lp-P*(prob.R^2-sum(x.^2)));
        end
    end
else
    vx2 = monomials(x,0:2*d-deg(g));
    [prog,coeff] = prog.newFree(length(vx2));
    P2=coeff'*vx2;
    switch para.relaxtype
        case 'sos'
            prog=prog.withSOS(P2);
        case 'sdsos'
            prog=prog.withSDSOS(P2);
        case 'dsos'
            prog=prog.withDSOS(P2);
    end
    if ~isempty(prob.phi)
        switch para.relaxtype
            case 'sos'
                prog=prog.withSOS(f-rho*g+Lp+phi*eta-P*(prob.R^2-sum(x.^2))-P2*(g-prob.gstar));
            case 'sdsos'
                prog=prog.withSDSOS(f-rho*g+Lp+phi*eta-P*(prob.R^2-sum(x.^2))-P2*(g-prob.gstar));
            case 'dsos'
                prog=prog.withDSOS(f-rho*g+Lp+phi*eta-P*(prob.R^2-sum(x.^2))-P2*(g-prob.gstar));
        end
    else
        switch para.relaxtype
            case 'sos' 
                prog=prog.withSOS(f-rho*g+Lp-P*(prob.R^2-sum(x.^2))-P2*(g-prob.gstar));
            case 'sdsos'
                prog=prog.withSDSOS(f-rho*g+Lp-P*(prob.R^2-sum(x.^2))-P2*(g-prob.gstar));
            case 'dsos'
                prog=prog.withDSOS(f-rho*g+Lp-P*(prob.R^2-sum(x.^2))-P2*(g-prob.gstar));
        end
    end
end

options = spot_sdp_default_options();
% sol = prog.minimize(-rho, @spot_sedumi, options);

%options = spot_sdp_default_options();
options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

sol = prog.minimize(-rho, @spot_mosek, options);

optval = double(sol.eval(rho));
switch para.relaxtype
    case 'sos' 
        disp(['The optimal value  of the ', num2str(order), '-th primal SDP relaxation (P_k) using SOS is ']);
        disp([num2str(value(optval))]);
    case 'sdsos'
        disp(['The optimal value  of the ', num2str(order), '-th primal SDP relaxation (P_k) using SDSOS is ']);
        disp([num2str(value(optval))]);
    case 'dsos'
        disp(['The optimal value  of the ', num2str(order), '-th primal SDP relaxation (P_k) using DSOS is ']);
        disp([num2str(value(optval))]);
end

end









