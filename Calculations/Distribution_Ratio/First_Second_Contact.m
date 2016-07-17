% This will plot the differences between first and second contacts of TBP

N=30; % Make N = 30 for prettier graphs
M=5;  % We have to look at different volume fractions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Variables/Input/Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%It will vary with different Distribution Ratios
DR=logspace(-2,1.30103,N);  %Log Space

%Uranium,Pu, Ru, Zr, Nb, Rare earch distribution coefficents
%estimated from reactor handbook for our experiment.
h_t=linspace(1/3,3,M);
h_tl=logspace(-0.47712125471967,3,M); %3 Vt/VH to 0.3 Vt/VH
%h_tl=logspace(-1,3,M); %10 Vt/VH to 0.3 Vt/VH

%Linear Scale
pt2=zeros(M,N);
pt=zeros(M,N);

%Log Scale
pt2l=zeros(M,N);
ptl=zeros(M,N);

for i=1:length(h_t)
    f=(1./DR).*(h_t(i)); %Factor to make rest of calcs easier
    fl=(1./DR).*(h_tl(i));
    %Percent in TBP after one contact (same total volume)
    pt(i,:)=1./(1+f.*0.5);
    ptl(i,:)=1./(1+fl.*0.5);
    %Percent in TBP after two separate contacts
    pt2(i,:)=1./(2+f.*(1+1./f.^2));  %This is the percent in second contact
    pt2l(i,:)=1./(2+fl.*(1+1./fl.^2));
    pt2(i,:)=pt2(i,:)+1./(1+f);      %This is adding the first contact
    pt2l(i,:)=pt2l(i,:)+1./(1+fl);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% Plot the mass Percent in TBP %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Linear Scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



leg=cell(1,M*2);
str1='V_{TBP}/V_{HNO_{3}}';
str2='V_{TBP}/V_{HNO_{3}}';
linewidth=linspace(1.5,3,M);
h=cell(1,M*2);
markersize=linspace(15,23,M);

for i=1:length(h_t)
    if i~=2
        h{i*2-1}=semilogx(DR,100*ptl(i,:),'k.-','LineWidth',linewidth(i),'MarkerSize',markersize(i)); %One Contact
        hold on;
        h{i*2}=semilogx(DR,100*pt2l(i,:),'r--','LineWidth',linewidth(i)); %Two Contacts
    end
    if i==2
        h{i*2-1}=semilogx(DR,100*ptl(i,:),'g.-','LineWidth',linewidth(i),'MarkerSize',markersize(i)); %One Contact
        hold on;
        h{i*2}=semilogx(DR,100*pt2l(i,:),'g--','LineWidth',linewidth(i)); %Two Contacts  
    end
    if i==1
        leg{i*2-1}=['1C ',str1,' = ',num2str(1./h_tl(i),'%.0f')];
        leg{i*2}=['2C ',str2,' = ',num2str(1./h_tl(i),'%.0f')];
    end
    if i~=1
        leg{i*2-1}=['1C ',str1,' = ',num2str(1./h_tl(i),'%.0e')];
        leg{i*2}=['2C ',str2,' = ',num2str(1./h_tl(i),'%.0e')]; 
    end
end


xlabel ('Distribution Ratio','FontSize',12,'Fontname','Garamond');
ylabel ('Mass Percent of solute in TBP','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');
h_legend=legend([h{1},h{2},h{length(h_t)*2-1},h{length(h_t)*2}],leg{1},leg{2},leg{length(h_t)*2-1},leg{length(h_t)*2});
set(gca,'fontsize',12,'Fontname','Garamond');
set(h_legend,'string');
set(h_legend,'FontSize',10);hold off;
xlim([10^-2,20])
hold off;

%semilogx(DR,100*pt2(1,:),'r--','LineWidth',2); %Uranium

