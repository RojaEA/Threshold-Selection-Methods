clc
clear
DataTrain=xlsread('SSMs-all');
format short g;
IntervalNumber=1000;%IntervalNumber
figure(1)
set(gcf,'units','normalized','position',[0.03,0.03,0.35,0.6]);%gcf position[x,y,length,width]
% %%%%%%%%
% % (Intersection point method)
% %%%%%%%%
fprintf('Intersection point method:')
%%%%%
Crash=DataTrain(find(DataTrain(:,3)==1),:);
NonCrash=DataTrain(find(DataTrain(:,3)==0),:);
X0=1:1:IntervalNumber;
X0=X0./IntervalNumber;
[CrashFrequency,CrashMidpoint]=hist(Crash(:,4),X0);
[NonCrashFrequency,NonCrashMidpoint]=hist(NonCrash(:,4),X0);
temp=0;
temp1=0;
for i=1:1:IntervalNumber
    temp=temp+CrashFrequency(i);
    CrashCumulativeRate(i)=1-temp./(sum(CrashFrequency)+1);
    temp1=temp1+NonCrashFrequency(i);
    NonCrashCumulativeRate(i)=temp1./(sum(NonCrashFrequency)+1);
    Minus(i)=abs(CrashCumulativeRate(i)-NonCrashCumulativeRate(i));
end
A=[X0' Minus' CrashCumulativeRate' NonCrashCumulativeRate'];
Location=find(A(:,2)==min(A(:,2)))
Location0=Location(length(Location),1)
PThreshold0=A(Location0,1)
PCrash=A(Location0,3)
PNonCrash=A(Location0,3)
% figure(1)
subplot(2,2,1);
plot(CrashMidpoint,CrashCumulativeRate,'k-','LineWidth',2)
hold on
plot(NonCrashMidpoint,NonCrashCumulativeRate,'k-.','LineWidth',2)
hold on
plot([PThreshold0 PThreshold0],[0 1],'k:')
hold on
plot(PThreshold0,PCrash,'rp','MarkerFaceColor','r')
text(PThreshold0-0.23,PCrash-0.11,'(0.373, 93.84%)','FontName','Times New Roman','color','r');
%legend('Crash','Non-crash','Threshold',0)
title('a) Intersection point method','FontName','Times New Roman','FontSize',10)
xlabel('Conflict risk','FontName','Times New Roman','FontSize',10)
ylabel('Cumulative proportion','FontName','Times New Roman','FontSize',10,'Rotation',90)
axis([0 1 0 1]);
set(gca,'YTick',0:0.2:1)
set(gca,'XTick',0:0.2:1)
set(gca,'YTickLabel',{'0','20%','40%','60%','80%','100%'}) 
set(gca,'FontName','Times New Roman','FontSize',10)
axis square
legend({'Conflict Rate','Non-Conflict Rate'},'Location','southwest','Orientation','vertical')

%%%%%%%%%
%£¨P-tile method£© %
%%%%%%%%%
fprintf('P-tile method:')
%%%%%
CrashNumber=sum(DataTrain(:,3)==1)
NonCrashNumber=sum(DataTrain(:,3)==0)
PNonCrash=NonCrashNumber/(CrashNumber+NonCrashNumber)%PCrash
[Frequency,Midpoint]=hist(DataTrain(:,4),IntervalNumber);%hist£¬Frequency£¬Midpoint
Rate=Frequency./sum(Frequency);%Rate
temp=0;
for i=1:1:IntervalNumber
    temp=temp+Frequency(i);
    CumulativeRate(i)=temp./(sum(Frequency)+1);
    Th(i)=abs(CumulativeRate(i)-PNonCrash);
end
A=[Midpoint' Rate' CumulativeRate' Th'];
Location=find(A(:,4)==min(A(:,4)));
Location0=Location(1,1)
PThreshold0=A(Location0,1)
% figure(2)
subplot(2,2,2);
plot(Midpoint,CumulativeRate,'k-','LineWidth',2)
hold on
plot([PThreshold0 PThreshold0],[0 1],'k:')
hold on
plot(PThreshold0,PNonCrash,'rp','MarkerFaceColor','r')
text(PThreshold0-0.02,PNonCrash-0.05,'(0.52342, 70.00%)','FontName','Times New Roman','color','r');
%legend('CPC','Threshold',0)
title('b) P-tile method','FontName','Times New Roman','FontSize',10)
xlabel('Conflict risk','FontName','Times New Roman','FontSize',10)
ylabel('Cumulative proportion','FontName','Times New Roman','FontSize',10,'Rotation',90)
axis([0 1 0.5 1]);
set(gca,'YTick',0.5:0.1:1)
set(gca,'XTick',0:0.2:1)
set(gca,'YTickLabel',{'50%','60%','70%','80%','90%','100%'}) 
set(gca,'FontName','Times New Roman','FontSize',10)
axis square
legend({'CPC','Threshold'},'Location','northwest','Orientation','vertical')

% %%%%%%%%%
% %£¨Maximum between-class variance method)%
% %%%%%%%%%
fprintf('Maximum between-class variance method:')
%%%%%
CrashNumber=sum(DataTrain(:,3)==1);
NonCrashNumber=sum(DataTrain(:,3)==0);
PCrash=CrashNumber/(CrashNumber+NonCrashNumber);%PCrash
Tabulate=tabulate(DataTrain(:,4));%Tabulate(:,1)£¬Tabulate(:,2)£¬Tabulate(:,3)£¨%£©
Tabulate1=[Tabulate Tabulate(:,1).*Tabulate(:,3)./100];%Tabulate(:,1).*Tabulate(:,2).*Tabulate(:,3)./100
for i=1:1:IntervalNumber
    Th(i)=i/IntervalNumber;
    CumulativeRate(i)=sum(Tabulate(find(Tabulate(:,1)<=Th(i)),3))./100;
    CumulativeExpectation0(i)=sum(Tabulate1(find(Tabulate1(:,1)<=Th(i)),4))/CumulativeRate(i);
    CumulativeExpectation1(i)=sum(Tabulate1(find(Tabulate1(:,1)>Th(i)),4))/(1-CumulativeRate(i));
end
G=CumulativeRate.*(1-CumulativeRate).*((CumulativeExpectation0-CumulativeExpectation1).^2);
A=[Th' G'];
PThreshold0=A(find(A(:,2)==max(A(:,2))),1)%PCrash
MaxVariance=max(A(:,2))
% figure(3)
subplot(2,2,3);
plot(Th,G,'k-','LineWidth',2)
hold on
plot([0.451 0.451],[0 1],'k:')%
hold on
plot(0.451,MaxVariance,'rp','MarkerFaceColor','r')%
text(0.451+0.008, MaxVariance+0.008, '(0.451, 0.15316)','FontName','Times New Roman','color','r');
%legend('\sigma_B^2','Threshold',0)
title('c) Maximum between-class variance method','FontName','Times New Roman','FontSize',10)
xlabel('Conflict risk','FontName','Times New Roman','FontSize',10)
ylabel('Between-class variance','FontName','Times New Roman','FontSize',10,'Rotation',90)
axis([0 1 0 0.18]);
set(gca,'YTick',0:0.04:0.18)% Y
set(gca,'XTick',0:0.2:1)% X
% set(gca,'YTickLabel',{'0','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'}) %Y
set(gca,'FontName','Times New Roman','FontSize',10)%
axis square
legend({'\sigma_b^2','Threshold'},'Location','southwest','Orientation','vertical')

% % %%%%%%%%
% % 
% %%%%%%%%
%%%fprintf('Maximum entropy method:')
%%%%%
%%%CrashNumber=sum(DataTrain(:,3)==1);%
%%%NonCrashNumber=sum(DataTrain(:,3)==0);%
%%%Tabulate=tabulate(DataTrain(:,4));%£¬Tabulate(:,1)£¬Tabulate(:,2)£¬Tabulate(:,3)£¨%£©
%%%Tabulate1=[Tabulate Try];
%%%H=sum(Try);%
%%%for i=1:1:IntervalNumber
   %%% Th(i)=i/IntervalNumber;
   %%% CumulativeRate(i)=sum(Tabulate(find(Tabulate(:,1)<=Th(i)),3))./100;
   %%% Ht(i)=sum(Tabulate1(find(Tabulate1(:,1)<=Th(i)),4));
   %%% HO(i)=log(CumulativeRate(i))+Ht(i)/CumulativeRate(i);
   %%% if abs(CumulativeRate(i)-1)<1e-6
     %%%   HB(i)=HB(i-1);
   %%% else
    %%%    HB(i)=log(1-CumulativeRate(i))+(H-Ht(i))/(1-CumulativeRate(i));
   %%% end
%%%end
%%%Hx=HO+HB;
%%%A=[Th' Hx'];
%%%PThreshold0=A(find(A(:,2)==max(A(:,2))),1)%PCrash
%%%MaxEntroy=max(A(:,2))
% figure(4)
%%%subplot(3,2,5);
%%%plot(Th,Hx,'k-','LineWidth',2)
%%%hold on
%%%plot([PThreshold0 PThreshold0],[0 18],'k:')%
%%%hold on
%%%plot(PThreshold0,MaxEntroy,'rp','MarkerFaceColor','r')%
%%%text(PThreshold0+0.03,MaxEntroy+0.2,'(0.008,12.348)','FontName','Times New Roman','color','r');
%legend('Entroy','Threshold',0)
%%%title('(e)Maximum entropy method','FontName','Times New Roman','FontSize',10)
%%%xlabel('Crash risk','FontName','Times New Roman','FontSize',10)
%%%ylabel('Entropy','FontName','Times New Roman','FontSize',10,'Rotation',90)
%%%axis([0 1.2 8 14]);
%%%set(gca,'YTick',8:2:14)% Y
%%%set(gca,'XTick',0:0.2:1.2)% X
% set(gca,'YTickLabel',{'0','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'}) %Y
%%%set(gca,'FontName','Times New Roman','FontSize',10)%
%%%axis square


% % %%%%%%%%
% % Minimum cross entropy method
% %%%%%%%%
fprintf('Minimum cross entropy method:')
%%%%%
Tabulate=tabulate(DataTrain(:,4));%£¬Tabulate(:,1)£¬Tabulate(:,2)£¬Tabulate(:,3)£¨%£©
Tabulate1=[Tabulate Tabulate(:,1).*Tabulate(:,3)./100 Tabulate(:,1).*Tabulate(:,3)./100.*log(Tabulate(:,1))];%Tabulate(:,1).*Tabulate(:,3)./100
H=sum(Tabulate1(:,5));%
for i=1:1:IntervalNumber
    Th(i)=i/IntervalNumber;
    CumulativeRate(i)=sum(Tabulate(find(Tabulate(:,1)<=Th(i)),3))./100;
    if CumulativeRate(i)==0
        CumulativeExpectation0(i)=1e-3;
    else
        CumulativeExpectation0(i)=sum(Tabulate1(find(Tabulate1(:,1)<=Th(i)),4))/CumulativeRate(i);
    end
    if abs(CumulativeRate(i)-1)<1e-6
        CumulativeExpectation1(i)=CumulativeExpectation1(i-1);
    else
        CumulativeExpectation1(i)=sum(Tabulate1(find(Tabulate1(:,1)>Th(i)),4))/(1-CumulativeRate(i));
    end
    D0(i)=log(CumulativeExpectation0(i))*sum(Tabulate1(find(Tabulate1(:,1)<=Th(i)),4));
    D1(i)=log(CumulativeExpectation1(i))*sum(Tabulate1(find(Tabulate1(:,1)>Th(i)),4));
end
D0=D0';
D1=D1';
D=H-D0-D1;
A=[Th' D];
location=find(A(:,2)==min(A(:,2)));
location0=location(1,1);
PThreshold0=A(location0,1)%PCrash
MminCrossEntroy=A(location0,2)
% figure(4)
subplot(2,2,4);
plot(Th,D,'k-','LineWidth',2)
hold on
plot([PThreshold0 PThreshold0],[0 0.12],'k:')%
hold on
plot(PThreshold0,MminCrossEntroy,'rp','MarkerFaceColor','r')%
text(PThreshold0+0.06,MminCrossEntroy-0.004,'(0.193, 0.03170)','FontName','Times New Roman','color','r');
%legend('Entropy','Threshold',0)
title('d) Minimum cross-entropy method','FontName','Times New Roman','FontSize',10)
xlabel('Conflict risk','FontName','Times New Roman','FontSize',10)
ylabel('Entropy','FontName','Times New Roman','FontSize',10,'Rotation',90)
axis([0 1 0.02 0.12]);
set(gca,'YTick',0.02:0.02:0.12)% Y
set(gca,'XTick',0:0.2:1)% X
% set(gca,'YTickLabel',{'0','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'}) %Y
set(gca,'FontName','Times New Roman','FontSize',10)%
axis square
legend({'Entropy','Threshold'},'Location','northwest','Orientation','vertical')
