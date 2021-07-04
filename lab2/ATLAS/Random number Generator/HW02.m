format long;
filename = 'myfilem3.txt';
fname = 'matlabm3.txt';
fid = fopen(filename,'w');
fi2d = fopen(fname,'w');
%fprintf(fid,'%s\n','Problem 1 (a)');

%k1=[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 ...
%    250 300 350 400 450 500];
k1=[5 6];
ik=1;
k=k1(ik);
while k~=k1(length(k1)),
    k=k1(ik);
    A=zeros(k,k);
    B=zeros(k,k);
    C=zeros(k,k);
    A=rand(k,k);
    B=rand(k,k);
    C=rand(k,k);
    
   % y1=y+A*x;
   %C1 = C+A*B;
   k
        
    fprintf(fid,'%i\n',k);
    %fprintf(fi2d,'%i\n',k);
    
    %for i=1:k,
    %    for j=1:k,          
    %    fprintf(fi2d,'%5.20f ',C1(i,j));
    %    end,
    %    fprintf(fi2d,'\n');
    %end,
    %fprintf(fi2d,'\n');
    
    
    for i=1:k,
        for j=1:k,          
        fprintf(fid,'%5.20f ',A(i,j));
        end,
        fprintf(fid,'\n');
    end,
    fprintf(fid,'\n');
    for i=1:k,
        for j=1:k,          
        fprintf(fid,'%5.20f ',B(i,j));
        end,
        fprintf(fid,'\n');
    end,
    fprintf(fid,'\n');
    for i=1:k,
        for j=1:k,          
        fprintf(fid,'%5.20f ',C(i,j));
        end,
        fprintf(fid,'\n');
    end,
    fprintf(fid,'\n');
    ik=ik+1;
end,
