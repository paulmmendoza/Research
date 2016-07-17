function [values]=Error_Calculations(x)
tol=0.001;

    %CHANGE THESE FOR NORMAL CALCULATIONS
time=200;              %count time
tc=1;                  %tc is a correction factor for time error
Uranium_Plutonium=0;   %Change to zero for normal, 1 for U&Pu

rows=size(x,1);values=zeros(rows,4);

for i=1:rows
    %Define Variables
a1=x(i,1);            %atomic abundance for element
a2=x(i,2);            %atomic abundance for interference
m1=x(i,3);            %Mass abundance for element
m2=x(i,4);            %Mass abundance for interference
p1=x(i,7);            %ppb for element
p2=x(i,8);            %ppb for interference
c(1)=x(i,5);          %CPS count one
c(2)=x(i,6);          %CPS count two
pc(1)=x(i,9);         %ppb/cps interference count 1
pc(2)=x(i,10);        %ppb/cps interference count 2
    %Start Calculations
opc=pc.*0.01;         %sigma ppb/cps interference
cp=1./pc;             %cps/ppb interference  
ocp=(opc./pc);        %sigma cps/ppb interference
    %Assumed Percent Error in Values
oa1=a1*0.005;         %assume 0.5% error in abundance (atomic)
om1=m1*0.005;         %assume 0.5% error in abundance (mass)
oa2=a2*0.005;         %assume 0.5% error in abundance (atomic)
om2=m2*0.005;         %assume 0.5% error in abundacne (mass)
op1=0.05;             %assume 0.05 ppb error
op2=0.05;             %assume 0.05 ppb error
oc=c.^0.5;            %Assume guassian distribution for cps

if(Uranium_Plutonium==1)
    oa1=a1*0.0005;
    op1=p1*0.005;     %assume a 0.5% error in ppb
    tc=(1/(time^0.5))*(0.447214); %time correction for lesser counts
end


pc1=zeros(1,2);opc1=pc1;pc2=pc1;opc2=pc1;
opc1m=pc1;pc1m=pc1;pc2m=pc1;opc2m=pc1;
for j=1:2

    %No interference Atomic Calculations
    pc1(j)=(p1*a1)/(c(j));  %ppb/cps
    term1=((a1/c(j))^2)*(op1^2);
    term2=((p1/c(j))^2)*((oa1^2));
    term3=(((a1*p1)/(c(j)^2))^2)*((oc(j)*tc)^2);
    opc1(j)=(term1+term2+term3)^0.5; %Error ppb/cps

    %Interference Atomic Calculations
    pc2(j)=(p1*a1)/(c(j)-p2*a2*cp(j)); %ppb/cps
    term1=((a1/(c(j)-(a2*cp(j)*p2)))^2)*(op1^2);
    term2=((p1/(c(j)-a2*cp(j)*p2))^2)*(oa1^2);
    term3=(((a1*p1)/((c(j)-a2*cp(j)*p2)^2))^2)*((oc(j)*tc)^2);
    term4=(((a1*p1*a2*cp(j))/((p2*a2*cp(j)-c(j))^2))^2)*(op2^2);
    term5=(((a1*p1*cp(j)*p2)/((a2*cp(j)*p2-c(j))^2))^2)*(oa2^2);
    term6=(((a1*p1*a2*p2)/((cp(j)*a2*p2-c(j))^2))^2)*(ocp(j)^2);
    opc2(j)=(term1+term2+term3+term4+term5+term6)^0.5; %Error ppb/cps

    %If values do not add up when they should I want to know
    if ( and( abs(pc1(j)-pc2(j))>tol , or(a2==0,m2==0)|p2==0) )
      msg= 'Error ppb/cps not the same when a2==0';
      error(msg);
    end
    if ( and( abs(opc1(j)-opc2(j))>tol , or(a2==0,m2==0)|p2==0))
        msg= 'Error ppb/cps Error not the same when a2==0';
        error(msg);
    end

    %No interference Mass Calculations
    pc1m(j)=(p1*m1)/(c(j));  %ppb/cps
    term1=((m1/c(j))^2)*(op1^2);
    term2=((p1/c(j))^2)*((om1^2));
    term3=(((m1*p1)/(c(j)^2))^2)*((oc(j)*tc)^2);
    opc1m(j)=(term1+term2+term3)^0.5; %Error ppb/cps

    %Interference Mass Calculations
    pc2m(j)=(p1*m1)/(c(j)-p2*m2*cp(j)); %ppb/cps
    term1=((m1/(c(j)-(m2*cp(j)*p2)))^2)*(op1^2);
    term2=((p1/(c(j)-m2*cp(j)*p2))^2)*(om1^2);
    term3=(((m1*p1)/((c(j)-m2*cp(j)*p2)^2))^2)*((oc(j)*tc)^2);
    term4=(((m1*p1*m2*cp(j))/((p2*m2*cp(j)-c(j))^2))^2)*(op2^2);
    term5=(((m1*p1*cp(j)*p2)/((m2*cp(j)*p2-c(j))^2))^2)*(om2^2);
    term6=(((m1*p1*m2*p2)/((cp(j)*m2*p2-c(j))^2))^2)*(ocp(j)^2);
    opc2m(j)=(term1+term2+term3+term4+term5+term6)^0.5; %Error ppb/cps
   
    
    %If values do not add up when they should I want to know
    if ( and( abs(pc1m(j)-pc2m(j))>tol , or(a2==0,m2==0)|p2==0) )
        msg= 'Error ppb/cps not the same when a2==0';
        error(msg);
    end
    if (and(abs(opc1m(j)-opc2m(j))>tol,or(a2==0,m2==0)|p2==0))
        msg= 'Error ppb/cps not the same when a2==0';
        error(msg);
    end

end
 
    %If we have interference we should print the right numbers
if (and(not(a2==0),not(m2==0))&&not(p2==0))
    values(i,1)=pc2(1);
    values(i,2)=opc2(1);
    values(i,3)=pc2(2);
    values(i,4)=opc2(2);
    values(i,5)=pc2m(1);
    values(i,6)=opc2m(1);
    values(i,7)=pc2m(2);
    values(i,8)=opc2m(2);
else
    values(i,1)=pc1(1);
    values(i,2)=opc1(1);
    values(i,3)=pc1(2);
    values(i,4)=opc1(2);
    values(i,5)=pc1m(1);
    values(i,6)=opc1m(1);
    values(i,7)=pc1m(2);
    values(i,8)=opc1m(2); 
end

end


%disp('    ppb/cps1       Sig1          ppb/cps2        sig2');

%spc1=sprintf('%0.7e',pc1);
%sopc1=sprintf('%0.7e',opc1);
%spc2=sprintf('%0.7e',pc2);
%sopc2=sprintf('%0.7e',opc2);

%output=[spc1,'  ',sopc1,'  ',spc2,'  ',sopc2];
%disp(output);

    %Display for Uranium and Plutonium
if(Uranium_Plutonium==1)
    UorPu=(((1/c)^2)*(op1^2)+((p1/(c^2))^2)*((oc*tc)^2))^0.5;

    spc1=sprintf('%0.7e',pc1/c);
    sUorPu=sprintf('%0.7e',UorPu);

    disp('    ppb/cpsU       Sig1');
    disp([spc1,'  ',sUorPu]);
end
format shortE