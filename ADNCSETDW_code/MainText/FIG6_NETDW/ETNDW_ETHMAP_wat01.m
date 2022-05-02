clear all; clc;
C=[1,0,0,0;0,1,0,0];
B=[0.0001;0.0002;0.0100;0.0300];
Data = importdata('TVKF_del00001_ETNDW_nor1.txt');
len=length(Data.data(:,1));
pos=Data.data(:,2);
ang=Data.data(:,3);
gamma=Data.data(:,4);
Im=Data.data(:,5);
res1=Data.data(:,6);
res2=Data.data(:,7);
Psi11=Data.data(:,8);
Psi12=Data.data(:,9);
Psi21=Data.data(:,10);
Psi22=Data.data(:,11);
wat1=Data.data(:,12);
wat2=Data.data(:,13);

Data_a = importdata('TVKF_wat01_del00001_ETNDW_ETHMAP.txt');
len_a=length(Data_a.data(:,1));
pos_a=Data_a.data(:,2);
ang_a=Data_a.data(:,3);
gamma_a=Data_a.data(:,4);
Im_a=Data_a.data(:,5);
res1_a=Data_a.data(:,6);
res2_a=Data_a.data(:,7);
Psi11_a=Data_a.data(:,8);
Psi12_a=Data_a.data(:,9);
Psi21_a=Data_a.data(:,10);
Psi22_a=Data_a.data(:,11);
wat1_a=Data_a.data(:,12);
wat2_a=Data_a.data(:,13);

C=[1,0,0,0;0,1,0,0];
Aa=0.1*eye(4,4);
xa=zeros(4,len_a);
xa(:,400)=[0.00001;0.00001;0.00001;0.00001];
for i=400:(len_a-1)
   xa(:,i+1)=Aa*xa(:,i); 
end
a=zeros(2,len_a);
for i=400:len_a
   a(:,i)=gamma_a(i)*(-1*[pos_a(i);ang_a(i)]+C*xa(:,i)); 
end
Sig_aTa=0;
AP=zeros(1,len_a);
for i=1:len_a
    Sig_aTa=Sig_aTa+(a(:,i))'*(a(:,i));
    AP(i)=Sig_aTa/i;
end

trPsi_a=zeros(1,len_a);
for i=1:len_a
   trPsi_a(i)=trace([Psi11(i),Psi12(i);Psi21(i),Psi22(i)]); 
end

count=0;
for i=1:len
    if(gamma(i)==1)
        count=count+1;
    end
end
tri_rate=count/len*100;
bR=zeros(2,2);
bD=zeros(2,1);
ind1=zeros(1,len);
ind2=zeros(1,len);

count_a=0;
for i=1:len_a
    if(gamma_a(i)==1)
        count_a=count_a+1;
    end
end
tri_rate_a=count_a/len_a*100;
bR_a=zeros(2,2);
bD_a=zeros(2,1);
ind1_a=zeros(1,len_a);
ind2_a=zeros(1,len_a);

for i=2:len
    bR=bR+([res1(i)*res1(i),res1(i)*res2(i);res2(i)*res1(i),res2(i)*res2(i)]-[Psi11(i),Psi12(i);Psi21(i),Psi22(i)]);
    bD=bD+([res1(i);res2(i)]*[wat1(i),wat2(i)]);
    ind1(i)=norm(bR/(i-1));
    ind2(i)=norm(bD/(i-1));
end

for i=2:len_a
    bR_a=bR_a+([res1_a(i)*res1_a(i),res1_a(i)*res2_a(i);res2_a(i)*res1_a(i),res2_a(i)*res2_a(i)]-[Psi11_a(i),Psi12_a(i);Psi21_a(i),Psi22_a(i)]);
    bD_a=bD_a+([res1_a(i);res2_a(i)]*[wat1_a(i),wat2_a(i)]);
    ind1_a(i)=norm(bR_a/(i-1));
    ind2_a(i)=norm(bD_a/(i-1));
end

%%ºÏ≤‚„–÷µ
thd1=zeros(1,len);
thd2=zeros(1,len);
kap1=1e-6;
rho1=1;
kap2=2.3e-6;
rho2=1;
for k=2:len
   thd1(k)=sqrt((1+rho1)*kap1*log((k-1))/(k-1))+1e-3;
   thd2(k)=sqrt((1+rho2)*kap2*log((k-1))/(k-1));
end

Fonts=18;
k=1:len_a;
fig1=figure(1);
set(fig1, 'Position', [0 0 1000 350]); 
set(fig1, 'PaperSize', [29.7000 21.0000]); 
set(fig1,'PaperPosition',[1 1 21 7]);
plot(k,pos(k),'b','LineWidth',2);hold on;
plot(k,pos_a(k),'-.r','LineWidth',2);hold off;
le11=xlabel('$k$');
le12=ylabel('$\alpha$ (m)');
set(gca,'FontSize',Fonts);
set(le11,'Interpreter','latex');
set(le12,'Interpreter','latex');

fig2=figure(2);
set(fig2, 'Position', [0 0 1000 350]); 
set(fig2, 'PaperSize', [29.7000 21.0000]); 
set(fig2,'PaperPosition',[1 1 21 7]);
plot(k,ang(k),'b','LineWidth',2);hold on;
plot(k,ang_a(k),'-.r','LineWidth',2);hold on;
le11=xlabel('$k$');
le12=ylabel('$\theta$ (rad)');
set(gca,'FontSize',Fonts);
set(le11,'Interpreter','latex');
set(le12,'Interpreter','latex');
t21=text(len_a,ang_a(len_a),'OFF');
set(t21,'Interpreter','latex','FontSize',Fonts);hold off;

fig3=figure(3);
set(fig3, 'Position', [0 0 1000 350]); 
set(fig3, 'PaperSize', [29.7000 21.0000]); 
set(fig3,'PaperPosition',[1 1 21 7]);
stairs(k,gamma(k),'b','LineWidth',2);hold on;
stairs(k,gamma_a(k),'-.r','LineWidth',2);hold off;
le11=xlabel('$k$');
le12=ylabel('$\gamma(k)$');
ylim([-1,2]);
set(gca,'FontSize',Fonts);
set(le11,'Interpreter','latex');
set(le12,'Interpreter','latex');

fig4=figure(4);
set(fig4, 'Position', [0 0 1000 350]); 
set(fig4, 'PaperSize', [29.7000 21.0000]); 
set(fig4,'PaperPosition',[1 1 21 7]);
semilogy(k,ind1(k),'b','LineWidth',2);hold on;
semilogy(k,ind1_a(k),'-.r','LineWidth',2);hold on;
semilogy(k,thd1(k),'--k','LineWidth',2);hold off;
le11=xlabel('$k$');
le12=ylabel('$\frac{1}{i}{\tilde \mathcal{R}_{n,i}}$');
set(gca,'FontSize',Fonts);
set(le11,'Interpreter','latex');
set(le12,'Interpreter','latex');

fig5=figure(5);
set(fig5, 'Position', [0 0 1000 350]); 
set(fig5, 'PaperSize', [29.7000 21.0000]); 
set(fig5,'PaperPosition',[1 1 21 7]);
semilogy(k,ind2(k),'b','LineWidth',2);hold on;
semilogy(k,ind2_a(k),'-.r','LineWidth',2);hold on;
semilogy(k,thd2(k),'--k','LineWidth',2);hold off;
le11=xlabel('$k$');
le12=ylabel('$\frac{1}{i}{\mathcal{D}_{n,i}}$');
set(gca,'FontSize',Fonts);
set(le11,'Interpreter','latex');
set(le12,'Interpreter','latex');

fig6=figure(6);
set(fig6, 'Position', [0 0 1000 350]); 
set(fig6, 'PaperSize', [29.7000 21.0000]); 
set(fig6,'PaperPosition',[1 1 21 7]);
plot(k,trPsi_a(k),'b','LineWidth',2);hold on;
plot(k,AP(k),'-.r','LineWidth',2);hold off;
le11=xlabel('$k$');
le12=ylabel('$\mathcal{A}(k)$');
ylim([-0.005,0.015]);
set(gca,'FontSize',Fonts);
set(le11,'Interpreter','latex');
set(le12,'Interpreter','latex');


saveas(fig1,'ETNDW01_ETHMAP_x1','pdf');
saveas(fig2,'ETNDW01_ETHMAP_x2','pdf');
saveas(fig3,'ETNDW01_ETHMAP_ga','pdf');
saveas(fig4,'ETNDW01_ETHMAP_rr','pdf');
saveas(fig5,'ETNDW01_ETHMAP_rd','pdf');
saveas(fig6,'ETNDW01_ETHMAP_ap','pdf');