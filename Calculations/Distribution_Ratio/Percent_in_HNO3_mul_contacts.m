%This program will walk through some calculations:
%This will change (distribution ratio)
DR=[0.051329
];
DRPu=16.1777;

vt=0.7;  %Volume of TBP
vh=0.5;  %Volume of HNO3
vtr=0.5; %Volume of TBP removed
contacts=4; %Number of TBP contacts
for j=1:length(DR)
f=1/(1+DR(j)*(vt/vh)); %A factor I think I will use alot 
fpu=1/(1+DRPu*(vt/vh));
m=1;        %Mass Percent in Stock solution (Mass percent remaining)
mpu=1;
%This will loop over as many cycles as you like
for i=1:contacts
    p=m*f;%After contact i, Percent of mass in HNO3
    ppu=mpu*fpu;
    %Total mass percent (TBP residual+HNO3) left in original solution
    m=m-(vtr/vt)*m*(1-f); 
    mpu=mpu-(vtr/vt)*mpu*(1-fpu);
    %disp('    Contact %In HNO3  %Removed  %Remains')
    %disp([i,p,1-m,m])
end
disp('    For the Contaminate');
disp('    Contact %In HNO3  %Removed  %Remains');
disp([contacts,p,1-m,m]);
disp('    For the Pu');
disp('    Contact %In HNO3  %Removed  %Remains');
disp([contacts,ppu,1-mpu,mpu]);
disp('    Decontamination Factor Estimate');
disp((1-mpu)/(1-m));
end