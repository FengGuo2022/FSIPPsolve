function [prob,para] = P_5_1_III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set the FSIPP problem data in the following form
%           min_{x\in R^m}  f(x)/g(x)
%                      s.t. p(x,y)<=0, \forall y\in Y\subset R^n
%                           phi_1(x)<=0, \ldots, phi_s(x)<=0
%
%Note: Y\subset [-1, 1]^n by scaling if necessary
%Y could be (1) hypercube [-1, 1]^n
%           (2) unit sphere
%           (3) unit ball
%           (4) polytope define by Ay<=b
%
%Please refer the paper:
%Feng Guo and Meijun Zhang, AN SDP METHOD FOR FRACTIONAL SEMI-INFINITE
%PROGRAMMING PROBLEMS WITH SOS-CONVEX POLYNOMIALS, arXiv:2110.04737
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
%Number of the variables x
prob.Xnum=2; 
%Number of the parameters y 
prob.Ynum=2; 
%==========================================================================
%Please use x to denote the variables and y to denote the parameters
%DO NOT change
x=sdpvar(prob.Xnum,1);
y=sdpvar(prob.Ynum,1);
prob.X=x;
prob.Y=y;
%==========================================================================
%Numerator f(x) of the objective
prob.f=(x(1)-1)^2+(x(2)-1)^2;
%==========================================================================
%Denominator of the objective
prob.g=1; 
%==========================================================================
%Semi-infinite constraint
prob.p=(y(1)*x(1)-y(2)*x(2))^2/4+(y(2)*x(1)+y(1)*x(2))^2-1;
%==========================================================================
%Finitely many constraints
prob.phi=[]; 
%==========================================================================
%Type of the index set Y. 
%Options:   'hypercube'; 
%           'sphere'; 
%           'ball'; 
%           'polytope'
prob.indset='sphere'; 
%==========================================================================
%If prob.indset='polytope' and latte is installed, 
%please set the data A and b (should be rational) such that Y is defined by Ay<=b,
% prob.A=[-1 0; 0 1; 1 -1];
% prob.b=[1; 1; 0]; 
%please set the directory where the latte command integrate is located
% prob.latte_dir='~/Downloads/latte-integrale-1.7.5/dest/bin/';
%==========================================================================
%constants R and gstar provided by user such that there exists a minimizer
%u^* satisfying ||u^*||<=R and g(u^*)>=gstar 
prob.R=1;
prob.gstar=1/2;
%==========================================================================
%Choose the primal/dual relaxation: 
%Options:   'primal'; (-> solving (P_k))
%           'dual': (-> solving (D_k))
para.relax='dual';
%==========================================================================
%If para.relax='primal', please choose the relaxation type: 
%Options:   'sos'; (using M_d(Q) in (P_k))
%           'sdsos'; (using M_d^sdsos(Q) in (P_k))
%           'dsos'; (using M_d^dsos(Q) in (P_k))
%If para.relax='dual', para.relaxtype='sos' is the only option
para.relaxtype='sos';
%==========================================================================
%If para.relax='primal' and para.relaxtype='sos', please choose the 
%flatform to impletement the relaxation (P_k): 
%Options:   'yalmip'; 
%           'yalmip+spotless_isos'
%If para.relax='dual', para.relaxtype='sos' is the only option
para.platform='yalmip';
%para.platform='yalmip+spotless_isos';
%==========================================================================
%Choose the SDP solver to solve (P_k) or (D_k)
%Options: 'sedumi' (usually slow but accurate)
%         'mosek' (fast)
para.solver='mosek';
%para.solver='sedumi';





