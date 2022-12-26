function [val] = int_polys(A, b, fun, vars)
%

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

[C,T]=coefficients(fun,vars);
s=length(C);

T_expo=[];
for i=1:s
    T_expo=[T_expo; degree(T(i),vars)];
end

intvals=[];
for i=1:s
    fileID = fopen('intdata','w');
    fprintf(fileID, '[');
    fprintf(fileID,'[%d/%d,[',1,1);
    for j=1:length(vars)
        if j<length(vars)
            fprintf(fileID,'%d,',T_expo(i,j));
        else
            fprintf(fileID,'%d',T_expo(i,j));
        end
    end
    fprintf(fileID,']]]');
    fclose(fileID);
    
    dir=pwd;
    dir=strrep(dir,' ','\ ');
    A_b_dir=strcat(dir,'/A_b_file');
    mono_dir=strcat(dir,'/intdata');
    
    command=strcat('~/Downloads/latte-integrale-1.7.5/dest/bin/integrate --valuation=integrate --triangulate --monomials=',mono_dir,32,A_b_dir);
    [~,cmdout] = system(command);

    i1=strfind(cmdout,'Answer'); 
    i2=i1+8; 
    while(cmdout(1,i2)~=' ')  
        i2=i2+1; 
    end
    intvals=[intvals,str2num(cmdout(1,i1+8:i2-1))];
end

val=intvals*C;

delete('A_b_file');
delete('intdata');
fclose all;

end

