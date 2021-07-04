A1 =importdata('residualijk.txt');

n1=A1(:,1);
r1=A1(:,2);

figure(1);hold on;
figure(1);plot(n1,r1,'.-');
xlabel('Vector Dimension');
ylabel('Residual');
%%%%%
A2 =importdata('residualikj.txt');

n2=A2(:,1);
r2=A2(:,2);

figure(1);plot(n2,r2,'.-');
%%%%%%
A3 =importdata('residualjik.txt');

n3=A3(:,1);
r3=A3(:,2);

figure(1);plot(n2,r2,'.-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B =importdata('result2.txt');

n=B(:,1);
mflops=B(:,2)*10^-6;

figure(1);subplot(1,2,2);plot(n,mflops,'.-');
xlabel('Vector Dimension');
ylabel('MFLOPS/s');