%This program makes a graph showing how volume affects the mass in each 
% portion of solutionm, and single Step decontamination factors.

N=30; % Make N = 30 for prettier graphs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Variables/Input/Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Volume of TBP added to stock solution (extraction)
vt=0.5;         %Linear Space

%Volume of stock solution 
vh=0.5;

%Distribution coefficents
DR=[logspace(-4,1,N)];

%Mass Percent in TBP first contact
pt=zeros(length(DR),1);
%Mass Percent in TBP second conatct
pt2=zeros(length(DR),1);
%Mass Percent in TBP third contact
pt3=zeros(length(DR),1);
%Mass Percent in TBP fourth contact
pt4=zeros(length(DR),1);

%Mass Percent in product 1 extraction 1 back-extraction
ph1e1b=zeros(length(DR),1);
%Mass Percent in product 4 extraction 1 back-extraction
ph4e1b=zeros(length(DR),1);
%Mass Percent in product 4 extraction 2 back-extraction
ph4e2b=zeros(length(DR),1);
%Mass Percent in product 4 extraction 3 back-extraction
ph4e3b=zeros(length(DR),1);

%Sum of mass percent after four contacts
spt2=zeros(length(DR),1);
spt3=zeros(length(DR),1);
spt4=zeros(length(DR),1);

%Sum of mass percent after four extractions and three back-extractions
sph4e3b=zeros(length(DR),1);

%Decontamination factor with respect to Pu 1 extraction 1 back-extraction
DF=zeros(length(DR)-1,1);
%Decontamination factor with respect to Pu 4 extraction 3 back-extraction
DF4=zeros(length(DR)-1,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calculate the percent in TBP for each Distribution ratio %%%%%%%%%%



%Multiply the separation by 5/7 because we remove only 0.5 ml not 0.7
Hold_up=5/7;
%Hold_up=1;

for i=1:length(DR)
   
   %First contact
   pt(i,1)=1./(1+(1/DR(i)).*(vh./vt)).*Hold_up;
   
   %Second Contact (if we keep the volume ratio the same between the first
   %and second contacts)
   pt2(i,1)=(1-pt(i,1)).*pt(i,1);
   
   %Third (I think this is right)
   pt3(i,1)=(1-pt2(i,1)-pt(i,1)).*pt(i,1);
   
   %Fourth (I think this is right)
   pt4(i,1)=(1-pt3(i,1)-pt2(i,1)-pt(i,1)).*pt(i,1);
   
   %Let us add up all the previous
   spt2(i,1)=pt(i,1)+pt2(i,1);
   spt3(i,1)=pt(i,1)+pt2(i,1)+pt3(i,1);
   spt4(i,1)=pt(i,1)+pt2(i,1)+pt3(i,1)+pt4(i,1); 
end

semilogx(DR,pt,'k^-','LineWidth',2,'MarkerSize',10); %One Contact
hold on;
semilogx(DR,spt2,'r.-','LineWidth',2,'MarkerSize',20); %Two Contacts
semilogx(DR,spt3,'bo-','LineWidth',2,'MarkerSize',8);  %Three Contacts
semilogx(DR,spt4,'g*-','LineWidth',2,'MarkerSize',10); %Four Contacts


xlabel ('Distribution Coefficient','FontSize',12,'Fontname','Garamond');
ylabel ('Combined percentage of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('One Contact','Two Contacts','Three Contacts','Four Contacts');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Back Extraction %%%%%%%%%%%%%%%%%%%%%%

%Let us assume Pu back extracts in Pu(III)

%Mass Percent in product 1 extraction 1 back-extraction
%Assume DR=0.0223
hold_up=0.5/0.6;
hold_up=1;
ph1e1b(1,1)=pt(1,1)*(1-1/(1+(1/0.0223)*(0.6/0.5)))*(hold_up);

DR_Pu_be=0.6;
vh_vt_be=2.1/2.0;
hold_up=2.0/2.1;

%Percent in TBP
TBP=1/(1+(1/DR_Pu_be)*vh_vt_be);

%Mass Percent in product 4 extraction 1 back-extraction
ph4e1b(1,1)=spt4(1,1)*(1-TBP)*hold_up;
%Mass Percent in product 4 extraction 2 back-extraction
ph4e2b(1,1)=(spt4(1,1)-ph4e1b(1,1))*(1-TBP)*hold_up;
%Mass Percent in product 4 extraction 3 back-extraction
ph4e3b(1,1)=(spt4(1,1)-ph4e1b(1,1)-ph4e2b(1,1))*(1-TBP)*hold_up;

%We added them all up
sph4e3b(1,1)=ph4e1b(1,1)+ph4e2b(1,1)+ph4e3b(1,1);

for i=2:length(DR)
    %DR here should get lower because we are lowering HNO3 concentration
    ph1e1b(i,1)=pt(i,1)*(1-1/(1+(1/(DR(i)))*(0.6/0.5)))*(0.5/0.6); 
    
    %Percent in TBP
    TBP=1/(1+(1/((DR(i)*1))*vh_vt_be));

    %Mass Percent in product 4 extraction 1 back-extraction
    ph4e1b(i,1)=spt4(i,1)*(1-TBP)*hold_up;
    %Mass Percent in product 4 extraction 2 back-extraction
    ph4e2b(i,1)=(spt4(i,1)-ph4e1b(i,1))*(1-TBP)*hold_up;
    %Mass Percent in product 4 extraction 3 back-extraction
    ph4e3b(i,1)=(spt4(i,1)-ph4e1b(i,1)-ph4e2b(i,1))*(1-TBP)*hold_up;

    %We added them all up
    sph4e3b(i,1)=ph4e1b(i,1)+ph4e2b(i,1)+ph4e3b(i,1);
    
end





%%%%%%%%%%%%%%%% Calculate the DF with respect to Pu %%%%%%%%%%%%%%%%%%%%%%

for i=2:length(DR) %i starts at 3 because FP start at 3
    %Decontamination factors for 1 extraction 1 back extraction
    DF(i-1,:)=ph1e1b(1,:).*(ph1e1b(i,:).^-1);
    
    %Decontamination factors for 4 extractions 3 back-extractions
    DF4(i-1,:)=sph4e3b(1,:).*(sph4e3b(i,:).^-1);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Plot Ratios
plot(DR(1,2:N+1),DF4./DF,'k^-','LineWidth',2,'MarkerSize',10) %Ratio
xlabel ('Distribution Coefficient','FontSize',12,'Fontname','Garamond');
ylabel ('Exp4_{DF}/Exp1_{DF}','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');

hold off;



%%%%%%%%%%%%%%%%%%%%%%%%%%% First and fourth contact (no add)%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Plot the Decontamination Factor in TBP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Log Log Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loglog(DR(1,2:N+1),DF,'k^-','LineWidth',2,'MarkerSize',10); %First Exp
hold on;
loglog(DR(1,2:N+1),DF4,'r.-','LineWidth',2,'MarkerSize',20); %Second Exp

xlabel ('Distribution Coefficient','FontSize',12,'Fontname','Garamond');
ylabel ('Decontamination Factor','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend('First Experiment','Second Experiment');
set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);
hold off;



disp(mean(DF4./DF))
