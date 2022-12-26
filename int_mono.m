function [val] = int_monos(A, b, expo, latte_dir)
%If the software latte is installed, given a monomial and a polyopte, 
%compute the integral of the monomial over the polytope
%Inputs: 
%   A, b: represent the polytope Ay<=b
%   expo: the exponent list beta of a monomial y^\beta
%   latte_dir: the directory where the latte command integrate is located
%Output: 
%   val: the integral of y^\beta over the polytope defined by Ay<=b


[m,d]=size(A);

fileID = fopen('A_b_file','w');

fprintf(fileID,'%d %d', m, d+1);
fprintf(fileID,'\n');

for i=1:m
    fprintf(fileID,'%d ', b(i));
    for j=1:d
        fprintf(fileID,'%d ', -A(i,j));
    end
    fprintf(fileID,'\n');
end

fclose(fileID);



fileID = fopen('intdata','w');
fprintf(fileID, '[');
fprintf(fileID,'[%d/%d,[',1,1);
for j=1:length(expo)
    if j<length(expo)
        fprintf(fileID,'%d,',expo(j));
    else
        fprintf(fileID,'%d',expo(j));
    end
end
fprintf(fileID,']]]');
fclose(fileID);

dir=pwd;
dir=strrep(dir,' ','\ ');
A_b_dir=strcat(dir,'/A_b_file');
mono_dir=strcat(dir,'/intdata');
latte_cmd=strcat(latte_dir, 'integrate --valuation=integrate --triangulate --monomials=');

command=strcat(latte_cmd,mono_dir,32,A_b_dir);
[~,cmdout] = system(command);

i1=strfind(cmdout,'Answer'); 
i2=i1+8; 
while(cmdout(1,i2)~=' ')  
    i2=i2+1; 
end
val=str2num(cmdout(1,i1+8:i2-1));

delete('A_b_file');
delete('intdata');
fclose all;

end

