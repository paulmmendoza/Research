N=32; %Number of sample numbers
X=linspace(0,N,N+1)';
%x refers to sample number
x =[21,31];

%y refers to instrument response
y =[0.0000127004782781985,0.0000115714228650488];

%linear response function
[m,fit] = polyfit(x,y,1);
b_l=m(2);
m_l=m(1);
Y_l=polyval(m,X);%The fit

%Logrithmic response Function
[m,fit] = polyfit(x,log(y),1);
b_L=exp(m(2));
m_L=-m(1);
Y_L=b_L.*exp(-m_L.*X);

%Base b Function
base=3;
[m,fit] = polyfit(x,logb(base,y),1);
b_t=base^(m(2));
m_t=-m(1);
Y_t=b_L.*base.^(-m_L.*X);

%Ploting
plot(x,y,'r.','MarkerSize',20);hold on;           %Plot Points
plot(X,Y_l,'g--','LineWidth',1.5,'MarkerSize',20) %Plot Linear Fit
plot(X,Y_L,'bo-','LineWidth',1.5,'MarkerSize',8)  %Plot Logrithmic Fit
plot(X,Y_t,'k*-','LineWidth',1.5,'MarkerSize',8)  %Plot Base 10 Function

%Plotting Parameters
xlabel ('Sample Number','FontSize',12,'Fontname','Garamond');
ylabel ('Instrument Response (ppb/cps)','FontSize',13,'Fontname','Garamond');
grid on;set(gca,'FontSize',12,'Fontname','Garamond');
set(findall(gcf,'type','text'),'FontSize',12,'Fontname','Garamond');

h_legend=legend('Points','Linear Fit','Logrithmic Fit','Base Function');

set(gca,'fontsize',12,'Fontname','Garamond')
set(h_legend,'FontSize',10);



