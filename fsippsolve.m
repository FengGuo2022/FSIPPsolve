function [optval]=fsippsolve(prob,para,order)
%Compute the optimal value of the SDP relaxation of the FSIPP problem
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
%Input:
%   prob: the problem structure
%   para: the parameter structure
%   order: the relaxation order
%Output:
%   optval: the optimal value of the order-th SDP relaxation


switch para.relax
    case 'primal'
        switch para.relaxtype
            case 'sos'
                switch para.platform
                    case 'yalmip'
                        optval=fsippsolve_primal(prob,para,order);
                    case 'yalmip+spotless_isos'
                        optval=fsippsolve_primal_alt(prob,para,order);
                end
            otherwise
                optval=fsippsolve_primal_alt(prob,para,order);
        end
    case 'dual'
        optval=fsippsolve_dual(prob,para,order);
end

end
    


