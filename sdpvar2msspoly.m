function [fun] = sdpvar2msspoly(expr,sdp_vars,mss_vars)
% To convert a sdpvar object used by Yalmip to a msspoly object used by
% spotless_isos
%expr: the object
%sdp_vars: the sdpvar variables
%mss_vars: the msspoly variables

[C,T]=coefficients(expr,sdp_vars);

fun=0;
for i=1:max(size(C))
    fun=fun+C(i)*prod((mss_vars').^(degree(T(i),sdp_vars)));
end


end

