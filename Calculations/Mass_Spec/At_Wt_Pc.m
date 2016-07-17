function [NP] = At_Wt_Pc(x,n)
%This Function will calculate weight percent from atom percent or the
%reverse. 

%Copy and paste from the IsotopicMass_NaturalAbundance.pdf file in your
%drive to x.

OP=x(:,2)/100;
W=x(:,1);
Columns=size(OP,1); %Used for Looping
NP=zeros(Columns,1);

if(n==1) %Atom To Weight
    sums=sum(OP.*W);
    for i=1:Columns
       NP(i,1)=(OP(i)*W(i))/sums;
    end
end

if(n==2) %Weight to Atom
    sums=sum(OP./W);
    for i=1:Columns
       NP(1,i)=(OP(i)/W(i))/sums;
    end
end

