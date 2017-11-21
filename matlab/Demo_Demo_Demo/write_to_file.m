function write_to_file(X,filename)
[r,c] = size(X);
fid = fopen(filename,'w');
for i = 1:r
    for j = 1:c
        fprintf(fid,'%5d, ',X(i,j));
        if j == c
           fprintf(fid,'\n');
        end
    end
end
fclose(fid);