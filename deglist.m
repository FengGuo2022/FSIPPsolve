function [degs, start_index] = deglist(num_var, min_tdeg, max_tdeg)
%Given the number num_var of variables, compute the list of exponets of monomials
%of degree from min_tdeg to max_tdeg


if num_var < 1
    error('Number of variables < 1.');
end

if min_tdeg < 0
    min_tdeg = 0;
end

if max_tdeg < min_tdeg
    error(sprintf('max tdeg %.0f < min tdeg %.0f',max_tdeg,min_tdeg));
end

if min_tdeg < 1
    start_index = 0;
else
    start_index = nchoosek(min_tdeg-1+num_var, num_var);
end

if num_var == 1
    degs = (min_tdeg:max_tdeg)';
    return;
end

degs = zeros(nchoosek(max_tdeg+num_var,num_var)-start_index, num_var);
degs(end-max_tdeg:end,1:2) = [(0:max_tdeg);(max_tdeg:-1:0)]';

s = max_tdeg + 1;
for i=2:num_var
    t = s; % (max_tdeg+i-1)!/(max_tdeg!(i-1)!)
    for j=max_tdeg-1:-1:((i==num_var)*min_tdeg)
        t = t*(j+1)/(i+j); % (i-1+j)!/((i-1)!j!)
        degs(end-s-t+1:end-s,1:i-1) = degs(end-s+1:end-s+t,1:i-1);
        degs(end-s-t+1:end-s,i) = degs(end-s+1:end-s+t,i) - 1;
        if i < num_var
            degs(end-s-t+1:end-s,i+1) = max_tdeg - j;
        end
        s = s + t;
    end
end
