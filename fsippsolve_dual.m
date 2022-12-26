function [rho]=fsippsolve_dual(prob,para,order)
%Compute the dual SDP relaxation (D_k) implemented by Yalmip
%Input:
%   prob: the problem structure
%   para: the parameter structure
%   order: the relaxation order
%Output:
%   optval: the optimal value of the order-th SDP relaxation (D_k)
%Please refer the paper:
%Feng Guo and Meijun Zhang, AN SDP METHOD FOR FRACTIONAL SEMI-INFINITE
%PROGRAMMING PROBLEMS WITH SOS-CONVEX POLYNOMIALS, arXiv:2110.04737
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=sdpvar(prob.Xnum,1);

d=max([degree(prob.f), degree(prob.g), degree(replace(prob.p,prob.Y,ones(size(prob.Y)))), max(degree(prob.phi))]);
d=ceil(d/2);

[degs1, start_index] = deglist(prob.Xnum, 0, d);
sdegs1=size(degs1,1);

sdegs2=nchoosek(prob.Xnum+2*d,prob.Xnum);
mm=sdpvar(sdegs2); %the moment sequence


%compute the moment matrix
MM=sdpvar(sdegs1);

for i=1:sdegs1
    for j=1:i 
        index=getindex(degs1(i,:)+degs1(j,:));
        MM(i,j)=mm(index);
        MM(j,i)=MM(i,j);
    end
end


[degs0, start_index] = deglist(prob.Xnum, 0, d-1);
sdegs0=size(degs0,1);

%compute the localizing moment matrix at R^2-||x||^2
MM0=sdpvar(sdegs0);
E=2*eye(prob.Xnum);
for i=1:sdegs0
    for j=1:i
        temp=0;
        index=getindex(degs0(i,:)+degs0(j,:));
        temp=temp+prob.R^2*mm(index);
        for k=1:prob.Xnum
            index=getindex(degs0(i,:)+degs0(j,:)+E(k,:));
            temp=temp-mm(index);
        end
        MM0(i,j)=temp;
        MM0(j,i)=temp;
    end
end

sMM=size(MM);
sMM=sMM(1);
MM=[MM zeros(sMM,sdegs0);zeros(sdegs0,sMM) MM0];


%compute the matrix Msig of moments such that Msig>=0 represents the
%constraint -L(p(x,y))\in P^k(Y)
[coef_X,term_Y]=coefficients(prob.p,prob.Y);
s=length(term_Y);


coef_mm=[];
for i=1:s
    [c_X,t_X]=coefficients(coef_X(i),prob.X);
    temp=0;
    for j=1:length(c_X)        
        index=getindex(degree(t_X(j),prob.X));
        temp=temp+c_X(j)*mm(index);
    end
    coef_mm=[coef_mm; temp];
end


term_Y_expo=[];
for i=1:s
    term_Y_expo=[term_Y_expo; degree(term_Y(i),prob.Y)];
end

[degs, start_index] = deglist(prob.Ynum, 0, order);
r=size(degs,1);

Msig=sdpvar(r,r);
for i=1:r
    for j=1:i
        temp=0;
        for k=1:s
            switch prob.indset
                case 'hypercube'
                    temp=temp+coef_mm(k)*intC(term_Y_expo(k,:)+degs(i,:)+degs(j,:));
                case 'sphere'
                    temp2=intSB(term_Y_expo(k,:)+degs(i,:)+degs(j,:));
                    temp=temp+coef_mm(k)*temp2;
                case 'ball'
                    [~,temp2]=intSB(term_Y_expo(k,:)+degs(i,:)+degs(j,:));
                    temp=temp+coef_mm(k)*temp2;
                case 'polytope'
                    temp=temp+coef_mm(k)*int_mono(prob.A, prob.b, term_Y_expo(k,:)+degs(i,:)+degs(j,:),prob.latte_dir);
            end
        end
        Msig(i,j)=-temp;
        Msig(j,i)=-temp;
    end
end


sMM=size(MM,1);
MM=[MM zeros(sMM,r);zeros(r,sMM) Msig];


%compute the localizing moment matrix at g
if strcmp(class(prob.g),'double')
    sMM=size(MM,1);
    MM=[MM zeros(sMM,1);zeros(1,sMM) mm(1,1)-1];
    sMM=size(MM);
    sMM=sMM(1);
    MM=[MM zeros(sMM,1);zeros(1,sMM) mm(1,1)-1];
else
    [cg,t_g]=coefficients(prob.g,prob.X);
    scg=length(t_g);

    vg=[];
    for i=1:scg
        vg=[vg; degree(t_g(i),prob.X)];
    end

    [degsgg, start_index] = deglist(prob.Xnum, 0, floor(d-degree(prob.g)/2));
    sdegsgg=size(degsgg,1);

    MMgg=sdpvar(sdegsgg);
    for i=1:sdegsgg
        for j=1:i
            temp=0;
            index=getindex(degsgg(i,:)+degsgg(j,:));
            temp=temp-prob.gstar*mm(index);
            for k=1:scg
                index=getindex(degsgg(i,:)+degsgg(j,:)+vg(k,:));
                temp=temp+cg(k)*mm(index);
            end
            MMgg(i,j)=temp;
            MMgg(j,i)=temp;
        end
    end

    sMM=size(MM);
    sMM=sMM(1);
    MM=[MM zeros(sMM,sdegsgg);zeros(sdegsgg,sMM) MMgg];

    g=0;
    scg=size(cg,1);
    for i=1:scg
        index=getindex(vg(i,:));
        g=g+cg(i)*mm(index);
    end

    sMM=size(MM);
    sMM=sMM(1);
    MM=[MM zeros(sMM,1);zeros(1,sMM) g-1];
    sMM=size(MM);
    sMM=sMM(1);
    MM=[MM zeros(sMM,1);zeros(1,sMM) 1-g];
end


for i=1:length(prob.phi)
    [c_phi,t_phi]=coefficients(prob.phi(i),prob.X);
    sphi=length(t_phi);

    temp=0;
    for j=1:sphi
        temp=temp+c_phi(j)*mm(getindex(degree(t_phi(j),prob.X)));
    end
    sMM=size(MM);
    sMM=sMM(1);
    MM=[MM zeros(sMM,1);zeros(1,sMM) -temp];
end


[c_f,t_f]=coefficients(prob.f,prob.X);
s_f=length(t_f);

f=0;
for i=1:s_f
    f=f+c_f(i)*mm(getindex(degree(t_f(i),prob.X)));
end

ops = sdpsettings('solver',para.solver);
ops = sdpsettings(ops,'verbose',0);
diagnostic=optimize([MM>=0],f,ops)

%%%%%%%%%%%%%%%%%%%%%%%
%the order of the variables for the minimizer is x(n),...,x(1), so reverse
%%%%%%%%%%%%%%%%%%%%%%
xx=[];, 
for i=prob.Xnum:-1:1
    xx=[xx value(mm(i+1,1))/value(mm(1,1))];
end
rho=value(f);

disp(['The optimal value r^dual_k of the ', num2str(order), '-th dual SDP relaxation (D_k) is ']);
disp([num2str(value(f))]);

disp(['The approximate minimizer computed by the ', num2str(order), '-th dual SDP relaxation (D_k) is ']);
disp([ '[',num2str(xx(1)),', ', num2str(xx(2)),']']);

  
end


