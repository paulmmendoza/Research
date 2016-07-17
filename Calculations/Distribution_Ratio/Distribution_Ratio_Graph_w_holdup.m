%This program makes a graph showing how volume affects the mass in each 
% portion of solutionm, and single Step decontamination factors.

N=30; % Make N = 30 for prettier graphs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Variables/Input/Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Volume of TBP added to stock solution
vtl=logspace(-3,0.176091,N);  %Log Space
vt=linspace(0,1.5*2,N);         %Linear Space

%Volume of stock solution (held constant)
vh=0.5;

%Uranium,Pu, Ru, Zr, Nb, Rare earch distribution coefficents
%estimated from reactor handbook for our experiment.
DR=[20,10,0.1,0.01,0.001,0.0001];
%Mass Percent in TBP
pt=zeros(length(DR),length(vt));
ptl=zeros(length(DR),length(vtl));
%Mass Percent in TBP second conatct?
pt2=zeros(length(DR),length(vt));
pt2l=zeros(length(DR),length(vtl));
%Mass Percent in TBP third contact
pt3=zeros(length(DR),length(vt));
pt3l=zeros(length(DR),length(vtl));
%Mass Percent in TBP fourth contact
pt4=zeros(length(DR),length(vt));
pt4l=zeros(length(DR),length(vtl));


%Sum of mass percent after four contacts
spt4=zeros(length(DR),length(vt));
spt4l=zeros(length(DR),length(vt));


%Decontamination factor with respect to Pu
DF=zeros(length(DR)-2,length(vt)); %No U and Pu
DFl=zeros(length(DR)-2,length(vt)); %No U and Pu
%Decontamination factor with respect to Pu second contact?
DF2=zeros(length(DR)-2,length(vt));
DF2l=zeros(length(DR)-2,length(vt));

DF2C=zeros(length(DR)-2,length(vt));
DF2Cl=zeros(length(DR)-2,length(vt));

DF3Cl=zeros(length(DR)-2,length(vt));

%Decontamination factor with respect to Pu combined fourth contact?
DF4C=zeros(length(DR)-2,length(vt));
DF4Cl=zeros(length(DR)-2,length(vt));


%DFOAl=zeros(1,length(vt)); %A single decontamination factor for the listed elements
ptls=zeros(1,length(vt)); %Sum up all the fractions
ptls2=zeros(1,length(vt)); %Sum up all the fractions for the second contact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Calculate the percent in TBP for each Distribution ratio %%%%%%%%%%

%Multiply the separation by 5/7 because we remove only 0.5 ml not 0.7
Hold_up=1;

for i=1:length(DR)
   
   %First contact (remember l stands for logrithmic scale)
   %pt(i,:)=1-((((DR(i).*vt)./vh)+1).^-1);
   pt(i,:)=1./(1+(1/DR(i)).*(vh./vt)).*Hold_up; %Same calculation just another form
   ptl(i,:)=(1-((((DR(i).*vtl)./vh)+1).^-1)).*Hold_up;
   
   
   %Second Contact (if we keep the volume ratio the same between the first
   %and second contacts)
   pt2(i,:)=(1-pt(i,:)).*pt(i,:);
   pt2l(i,:)=(1-ptl(i,:)).*ptl(i,:);
   
   %Third (I think this is right- I didn't think about this)
   pt3(i,:)=(1-pt2(i,:)-pt(i,:)).*pt(i,:);
   pt3l(i,:)=(1-pt2l(i,:)-ptl(i,:)).*ptl(i,:);
   
   %Fourth (I think this is right - I didn't think about this)
   pt4(i,:)=(1-pt3(i,:)-pt2(i,:)-pt(i,:)).*pt(i,:);
   pt4l(i,:)=(1-pt3l(i,:)-pt2l(i,:)-ptl(i,:)).*ptl(i,:);
   
   %Let us add up all the previous
   spt4=pt(i,:)+pt2(i,:)+pt3(i,:)+pt4(i,:);
   spt4l=ptl(i,:)+pt2l(i,:)+pt3l(i,:)+pt4l(i,:);
   
   
end





%%%%%%%%%%%%%%%% Calculate the DF with respect to Pu %%%%%%%%%%%%%%%%%%%%%%

for i=3:length(DR) %i starts at 3 because FP start at 3
    %Decontamination factors for the first contact
    DF(i-2,:)=pt(2,:).*(pt(i,:).^-1);
    DFl(i-2,:)=ptl(2,:).*(ptl(i,:).^-1);
    
    %Decontamination factors for the second contact
    DF2(i-2,:)=pt2(2,:).*(pt2(i,:).^-1);
    DF2l(i-2,:)=pt2l(2,:).*(pt(i,:).^-1);
    
    %Combined decontamination factors (after second contact)
    DF2C(i-2,:)=(pt(2,:)+pt2(2,:)).*((pt(i,:)+pt2(i,:)).^-1);
    DF2Cl(i-2,:)=(ptl(2,:)+pt2l(2,:)).*((ptl(i,:)+pt2l(i,:)).^-1);
    
    %Third and fourth contacts can be carried out in a similiar way with a
    %more clever code
    
    DF3Cl(i-2,:)=(ptl(2,:)+pt2l(2,:)+pt3l(2,:)).*((ptl(i,:)+pt2l(i,:)+pt3l(i,:)).^-1);
    
    %Combined Decontamination factors (after fourth contact)
    DF4C(i-2,:)=(pt(2,:)+pt2(2,:)+pt3(2,:)+pt4(2,:)).*((pt(i,:)+pt2(i,:)+pt3(i,:)+pt4(i,:)).^-1);
    DF4Cl(i-2,:)=(ptl(2,:)+pt2l(2,:)+pt3l(2,:)+pt4l(2,:)).*((ptl(i,:)+pt2l(i,:)+pt3l(i,:)+pt4l(i,:)).^-1);
    
end

%For the overall decontamination factor lets sum all the mass fractions
%We are assuming that all the species have the same mass in the solution,
%which isn't true, but I do not think it will matter. 

for j=1:length(pt(1,:))
    ptls(1,j)=sum(ptl(3:length(ptl(:,1)),j));
    ptls2(1,j)=sum(ptl(3:length(ptl(:,1)),j)+pt2l(3:length(ptl(:,1)),j));
end
DFOAl=ptl(2,:).*(ptls(1,:).^-1)*4;
DFOA2=(ptl(2,:)+pt2l(2,:)).*(ptls2(1,:).^-1)*4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First Contact %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%% Plot the mass Percent in TBP %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Log Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogx(vtl./vh,100*ptl(1,:),'k--','LineWidth',2); %Uranium
hold on;
semilogx(vtl./vh,100*ptl(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Plutonium
semilogx(vtl./vh,100*ptl(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Ru
semilogx(vtl./vh,100*ptl(4,:),'b--','LineWidth',1.5,'MarkerSize',10); %Zr
semilogx(vtl./vh,100*ptl(5,:),'g.-','LineWidth',1.5,'MarkerSize',20); %Nb
semilogx(vtl./vh,100*ptl(6,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Mass Percent of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('U  DR = 20','Pu DR = 10','Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;

%%%%%%%%%%%%%%%%%%%%%%% Plot the mass Percent in TBP %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Log log %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loglog(vtl./vh,100*ptl(1,:),'k--','LineWidth',2); %Uranium
hold on;
loglog(vtl./vh,100*ptl(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Plutonium
loglog(vtl./vh,100*ptl(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Ru
loglog(vtl./vh,100*ptl(4,:),'b--','LineWidth',1.5,'MarkerSize',10); %Zr
loglog(vtl./vh,100*ptl(5,:),'g.-','LineWidth',1.5,'MarkerSize',20); %Nb
loglog(vtl./vh,100*ptl(6,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Mass Percent of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('U  DR = 20','Pu DR = 10','Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;

%%%%%%%%%%%%%%%%%%%%%%% Plot the mass Percent in TBP %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Linear Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(vt./vh,100*pt(1,:),'k--','LineWidth',2); %Uranium
hold on;
plot(vt./vh,100*pt(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Plutonium
plot(vt./vh,100*pt(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Ru
plot(vt./vh,100*pt(4,:),'b--','LineWidth',1.5,'MarkerSize',10); %Zr
plot(vt./vh,100*pt(5,:),'g.-','LineWidth',1.5,'MarkerSize',20); %Nb
plot(vt./vh,100*pt(6,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Mass Percent of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('U  DR = 20','Pu DR = 10','Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);hold off;

%%%%%%%%%%%%%%%%%%%%%%% Plot the mass Percent in TBP %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Semi log y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogy(vt./vh,100*pt(1,:),'k--','LineWidth',2); %Uranium
hold on;
semilogy(vt./vh,100*pt(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Plutonium
semilogy(vt./vh,100*pt(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Ru
semilogy(vt./vh,100*pt(4,:),'b--','LineWidth',1.5,'MarkerSize',10); %Zr
semilogy(vt./vh,100*pt(5,:),'g.-','LineWidth',1.5,'MarkerSize',20); %Nb
semilogy(vt./vh,100*pt(6,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Mass Percent of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('U  DR = 20','Pu DR = 10','Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;



%%%%%%%%%%%%%%%% Plot the Decontamination Factor in TBP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Log Scale X %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogx(vtl./vh,100*DFl(1,:),'k--','LineWidth',2); %Ru
hold on;
semilogx(vtl./vh,100*DFl(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Zr
semilogx(vtl./vh,100*DFl(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Nb
semilogx(vtl./vh,100*DFl(4,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE


xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Single Pass Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;


%%%%%%%%%%%%%%%% Plot the Decontamination Factor in TBP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Log Log Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loglog(vtl./vh,100*DFl(1,:),'k--','LineWidth',2); %Ru
hold on;
loglog(vtl./vh,100*DFl(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Zr
loglog(vtl./vh,100*DFl(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Nb
loglog(vtl./vh,100*DFl(4,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Single Pass Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;

%%%%%%%%%%%%%%%% Plot the Decontamination Factor in TBP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Lin Scale  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(vtl./vh,100*DF(1,:),'k--','LineWidth',2); %Ru
hold on;
plot(vtl./vh,100*DF(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Zr
plot(vtl./vh,100*DF(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Nb
plot(vtl./vh,100*DF(4,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE


xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Single Pass Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;


%%%%%%%%%%%%%%%% Plot the Decontamination Factor in TBP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  semi log Scale y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogy(vtl./vh,100*DF(1,:),'k--','LineWidth',2); %Ru
hold on;
semilogy(vtl./vh,100*DF(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Zr
semilogy(vtl./vh,100*DF(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Nb
semilogy(vtl./vh,100*DF(4,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE


xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Single Pass Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second Contact %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% Plot the mass Percent in TBP %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Linear Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Just the second Contact %%%%%%%%%%%%%%%%%%%%%%

plot(vt./vh,100*pt2(1,:),'k--','LineWidth',2); %Uranium
hold on;
plot(vt./vh,100*pt2(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Plutonium
plot(vt./vh,100*pt2(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Ru
plot(vt./vh,100*pt2(4,:),'b--','LineWidth',1.5,'MarkerSize',10); %Zr
plot(vt./vh,100*pt2(5,:),'g.-','LineWidth',1.5,'MarkerSize',20); %Nb
plot(vt./vh,100*pt2(6,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Mass Percent of solute in TBP Second Contact','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('U  DR = 20','Pu DR = 10','Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
%ylim([0,100]);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%  Just the second Contact %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Plot the mass Percent in TBP %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Log Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogx(vtl./vh,100*pt2l(1,:),'k--','LineWidth',2); %Uranium
hold on;
semilogx(vtl./vh,100*pt2l(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Plutonium
semilogx(vtl./vh,100*pt2l(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Ru
semilogx(vtl./vh,100*pt2l(4,:),'b--','LineWidth',1.5,'MarkerSize',10); %Zr
semilogx(vtl./vh,100*pt2l(5,:),'g.-','LineWidth',1.5,'MarkerSize',20); %Nb
semilogx(vtl./vh,100*pt2l(6,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Mass Percent of solute in Second Contact','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('U  DR = 20','Pu DR = 10','Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;









%%%%%%%%%%%%%%%%%%%%%%% Plot the mass Percent in TBP %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Linear Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  First + Second Contact %%%%%%%%%%%%%%%%%%%%%%%

plot(vt./vh,100*(pt2(1,:)+pt(1,:)),'k--','LineWidth',2); %Uranium
hold on;
plot(vt./vh,100*(pt2(2,:)+pt(2,:)),'r.-','LineWidth',1.5,'MarkerSize',20); %Plutonium
plot(vt./vh,100*(pt2(3,:)+pt(3,:)),'bo-','LineWidth',1.5,'MarkerSize',8);  %Ru
plot(vt./vh,100*(pt2(4,:)+pt(4,:)),'b--','LineWidth',1.5,'MarkerSize',10); %Zr
plot(vt./vh,100*(pt2(5,:)+pt(5,:)),'g.-','LineWidth',1.5,'MarkerSize',20); %Nb
plot(vt./vh,100*(pt2(6,:)+pt(6,:)),'g*-','LineWidth',1.5,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Mass Percent of solute in TBP Second Contact','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('U  DR = 20','Pu DR = 10','Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Just the second contact %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Plot the Decontamination Factor in TBP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Log Log Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loglog(vtl./vh,DF2l(1,:),'k--','LineWidth',2); %Ru
hold on;
loglog(vtl./vh,DF2l(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Zr
loglog(vtl./vh,DF2l(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Nb
loglog(vtl./vh,DF2l(4,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Second Pass Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%% First and second contact %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Plot the Decontamination Factor in TBP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Log Log Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loglog(vtl./vh,DF2Cl(1,:),'k--','LineWidth',4); %Ru
hold on;
loglog(vtl./vh,DF2Cl(2,:),'r.-','LineWidth',2,'MarkerSize',20); %Zr
loglog(vtl./vh,DF2Cl(3,:),'bo-','LineWidth',2,'MarkerSize',8);  %Nb
loglog(vtl./vh,DF2Cl(4,:),'g*-','LineWidth',2,'MarkerSize',10); %RE
loglog(vtl./vh,DFOAl(1,:),'kx-','LineWidth',2,'MarkerSize',12); %Sum

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Combined Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('Ru DR = 0.2','Zr DR = 0.09','Nb DR = 0.025','Rare Earths DR = 0.01','Previous Added');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%% Plot the first pass again to compare
loglog(vtl./vh,DFl(1,:),'k--','LineWidth',2); %Ru
hold on;
loglog(vtl./vh,DFl(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Zr
loglog(vtl./vh,DFl(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Nb
loglog(vtl./vh,DFl(4,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE
loglog(vtl./vh,DFOA2(1,:),'kx-','LineWidth',2,'MarkerSize',12); %Sum

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%% First and second contact (no add) %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Plot the Decontamination Factor in TBP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Log Log Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loglog(vtl./vh,DF2Cl(1,:),'k--','LineWidth',4); %Ru
hold on;
loglog(vtl./vh,DF2Cl(2,:),'r.-','LineWidth',2,'MarkerSize',20); %Zr
loglog(vtl./vh,DF2Cl(3,:),'bo-','LineWidth',2,'MarkerSize',8);  %Nb
loglog(vtl./vh,DF2Cl(4,:),'g*-','LineWidth',2,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Combined Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('DC = 0.2','DC = 0.09','DC = 0.025','DC = 0.01');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%% Plot the first pass again to compare
loglog(vtl./vh,DFl(1,:),'k--','LineWidth',2); %Ru
hold on;
loglog(vtl./vh,DFl(2,:),'r.-','LineWidth',1.5,'MarkerSize',20); %Zr
loglog(vtl./vh,DFl(3,:),'bo-','LineWidth',1.5,'MarkerSize',8);  %Nb
loglog(vtl./vh,DFl(4,:),'g*-','LineWidth',1.5,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
hold off;



%%%%%%%%%%%%%%%%%%%%%%%%%%% First and fourth contact (no add) %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Plot the Decontamination Factor in TBP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Log Log Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loglog(vtl./vh,DF4Cl(1,:),'k^-','LineWidth',2,'MarkerSize',10); %Ru
hold on;
loglog(vtl./vh,DF4Cl(2,:),'r.-','LineWidth',2,'MarkerSize',20); %Zr
loglog(vtl./vh,DF4Cl(3,:),'bo-','LineWidth',2,'MarkerSize',8);  %Nb
loglog(vtl./vh,DF4Cl(4,:),'g*-','LineWidth',2,'MarkerSize',10); %RE

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Combined Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('DC = 0.1','DC = 0.01','DC = 0.001','DC = 0.0001');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%% Plot the first pass again to compare
loglog(vtl./vh,DFl(1,:),'k^-','LineWidth',1.5,'MarkerSize',8); 
hold on;
loglog(vtl./vh,DFl(2,:),'r.-','LineWidth',1.5,'MarkerSize',15); 
loglog(vtl./vh,DFl(3,:),'bo-','LineWidth',1.5,'MarkerSize',6);  
loglog(vtl./vh,DFl(4,:),'g*-','LineWidth',1.5,'MarkerSize',8); 

xlabel ('V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
ylabel ('Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
hold off;

%%%% Add an extra vert line
hold on;
loglog(ones(1,N)*1.4,logspace(0,5,N),'m','LineWidth',2.5)
xlim([0.0015,4])


%vtl(27)/vh;   27 on the log scale is close to what we want with a volume
%ratio of 5/7

%Actual DF values at the 5/7 ratio
disp('Four contact')
disp(DF4Cl(:,27))
disp('First contact')
disp(DFl(:,27))

%The Ratio between 4th contact and the first
disp('Ratio between 4th and first')
disp(DF4Cl(:,27)./DFl(:,27))




hold off;
p(1)=loglog(vtl./vh,DFl(1,:),'bo-','LineWidth',1.5,'MarkerSize',6); 
hold on;
p(2)=loglog(vtl./vh,DF2Cl(1,:),'g*-','LineWidth',1.5,'MarkerSize',8); 
p(3)=loglog(vtl./vh,DF3Cl(1,:),'r.-','LineWidth',1.5,'MarkerSize',15); 
p(4)=loglog(vtl./vh,DF4Cl(1,:),'k^-','LineWidth',2,'MarkerSize',10);
xlabel ('V_{R} = V_{TBP}/V_{HNO_{3}}','FontSize',12,'Fontname','Garamond');
%grid on;
hline=loglog(0.03,1,'w');


ylabel ('Decontamination Factor of solute in TBP','FontSize',13,'Fontname','Garamond');
set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
leg1=legend(p,'One Contact','Two Contacts','Three Contacts','Four Contacts');

%loglog(ones(1,N)*1.4,logspace(0,2,N),'m','LineWidth',2)
xlim([0.0015,4])

