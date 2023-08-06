
format short e
%% QUESTION 1: 
disp("QUESTION 1:")
EL=131*10^(3); %GPa
ET=11*10^(3); %GPa
GLT=7.2*10^(3);%GPa
VLT=0.32;




SI5=[1/EL -VLT/EL 0;
    -VLT/EL 1/ET 0;
   0 0 1/GLT ];

Qi =inv(SI5);

fprintf("Reduced compliance matrice (all elements have unit [1/GPa]):\n")
disp(SI5)
fprintf("Reduced stiffness matrice (all elements have unit [GPa]):\n")
disp(Qi)
%%
%QUESTION 2:
disp("QUESTION 2:")
disp(" ")
h=0.256;

c0=cosd(0)
c1=cosd(45)
c2=cosd(-45)
c3=cosd(90)

s0=sind(0)
s1=sind(45)
s2=sind(-45)
s3=sind(90)

Ts0 = zeros (3,3);
Ts0(1,1) = c0.^2;
Ts0(1,2) = s0.^2;
Ts0(1,3) = -2*s0*c0;
Ts0(2,1) = s0.^2;
Ts0(2,2) = c0.^2;
Ts0(2,3) = 2*s0*c0;
Ts0(3,1) = s0*c0;
Ts0(3,2) = -s0*c0;
Ts0(3,3) = (c0.^2)-(s0.^2);

Ts45 = zeros (3,3);
Ts45(1,1) = c1.^2;
Ts45(1,2) = s1.^2;
Ts45(1,3) = -2*s1*c1;
Ts45(2,1) = s1.^2;
Ts45(2,2) = c1.^2;
Ts45(2,3) = 2*s1*c1;
Ts45(3,1) = s1*c1;
Ts45(3,2) = -s1*c1;
Ts45(3,3) = (c1.^2)-(s1.^2);

Tsm45 = zeros (3,3);
Tsm45(1,1) = c2.^2;
Tsm45(1,2) = s2.^2;
Tsm45(1,3) = -2*s2*c2;
Tsm45(2,1) = s2.^2;
Tsm45(2,2) = c2.^2;
Tsm45(2,3) = 2*s2*c2;
Tsm45(3,1) = s2*c2;
Tsm45(3,2) = -s2*c2;
Tsm45(3,3) = (c2.^2)-(s2.^2);

Ts90 = zeros (3,3);
Ts90(1,1) = c3.^2;
Ts90(1,2) = s3.^2;
Ts90(1,3) = -2*s3*c3;
Ts90(2,1) = s3.^2;
Ts90(2,2) = c3.^2;
Ts90(2,3) = 2*s3*c3;
Ts90(3,1) = s3*c3;
Ts90(3,2) = -s3*c3;
Ts90(3,3) = (c3.^2)-(s3.^2);

Teinv0 = zeros (3,3);
Teinv0(1,1) = c0.^2;
Teinv0(1,2) = s0.^2;
Teinv0(1,3) = s0*c0;
Teinv0(2,1) = s0.^2;
Teinv0(2,2) = c0.^2;
Teinv0(2,3) = -s0*c0;
Teinv0(3,1) = -2*s0*c0;
Teinv0(3,2) = 2*s0*c0;
Teinv0(3,3) = (c0.^2)-(s0.^2);

Teinv45 = zeros (3,3);
Teinv45(1,1) = c1.^2;
Teinv45(1,2) = s1.^2;
Teinv45(1,3) = s1*c1;
Teinv45(2,1) = s1.^2;
Teinv45(2,2) = c1.^2;
Teinv45(2,3) = -s1*c1;
Teinv45(3,1) = -2*s1*c1;
Teinv45(3,2) = 2*s1*c1;
Teinv45(3,3) = (c1.^2)-(s1.^2);

Teinvm45 = zeros (3,3);
Teinvm45(1,1) = c2.^2;
Teinvm45(1,2) = s2.^2;
Teinvm45(1,3) = s2*c2;
Teinvm45(2,1) = s2.^2;
Teinvm45(2,2) = c2.^2;
Teinvm45(2,3) = -s2*c2;
Teinvm45(3,1) = -2*s2*c2;
Teinvm45(3,2) = 2*s2*c2;
Teinvm45(3,3) = (c2.^2)-(s2.^2);

Teinv90 = zeros (3,3);
Teinv90(1,1) = c3.^2;
Teinv90(1,2) = s3.^2;
Teinv90(1,3) = s3*c3;
Teinv90(2,1) = s3.^2;
Teinv90(2,2) = c3.^2;
Teinv90(2,3) = -s3*c3;
Teinv90(3,1) = -2*s3*c3;
Teinv90(3,2) = 2*s3*c3;
Teinv90(3,3) = (c3.^2)-(s3.^2);

Q0 = Ts0*Qi*Teinv0
Q45 = Ts45*Qi*Teinv45;
Qm45 = Tsm45*Qi*Teinvm45;
Q90 = Ts90*Qi*Teinv90;

A = 2*h*Q0+2*h*Q45+2*h*Qm45+2*h*Q90
B = zeros(3,3) %Symmetrical laminate
D = (1/3)*(((4*h).^3-(3*h).^3)*Q0 +((3*h).^3-(2*h).^3)*Q45+((2*h).^3-(h).^3)*Qm45+((h).^3-(0).^3)*Q90+((0).^3-(-h).^3)*Q90+((-h).^3-(-2*h).^3)*Qm45+((-2*h).^3-(-3*h).^3)*Q45+((-3*h).^3-(-4*h).^3)*Q0)

%%
%QUESTION 3:

disp("QUESTION 3:")
N = zeros (3,1);
N(1,1) = 100*1000; 

straing= inv(A)*N %strain in global 
strain0=zeros (3,1);
strain0=straing;
strain0(1,1)=straing(1,1);
strain0(2,1)=straing(2,1);
strain0(3,1)=straing(3,1)./2;
disp('strain0 : ')
disp(strain0)

straing45 = Teinv45*strain0; %strain in local 45
strain45 = straing45;
strain45(3,1)=straing45(3,1)./2;
disp('strain45 : ')
disp(strain45)

straingm45 = Teinvm45*strain0; %strain in local min45
strainm45 = straingm45;
strainm45(3,1)=straingm45(3,1)./2;
disp('strain -45 : ')
disp(strainm45)

straing90= Teinv90*strain0; %strain in local 90
strain90 = straing90;
strain90(3,1)=straing90(3,1)./2;
disp('strain90 : ')
disp(strain90)

stress0 = Qi*strain0./100
stress45 = Qi*strain45./100;
stress45(3,1)=stress45(3,1)*2;
disp('stress 45 :')
disp(stress45)

stressm45 = Qi*strainm45./100;
stressm45(3,1)=stressm45(3,1)*2;
disp('stress -45 :')
disp(stressm45)

stress90 = Qi*strain90./100;
stress90(3,1)=stress90(3,1)*2;
disp('stress90 :')
disp(stress90)

%%
%QUESTION 3 suite :
xt=1990 %MPA
xc=-1500 %MPA
yt=38 %MPA
sc=70 %MPA
%%disp("QUESTION 3 suite:")
R0 = zeros (3,1);
R0(1,1) =xt./stress0(1,1);
R0(2,1) =yt./stress0(2,1);
R0(3,1) =0;
M0=zeros(3,1);
M0(1,1)=1./R0(1,1);
M0(2,1)=1./R0(2,1);
M0(3,1)=0;

R45 = zeros (3,1);
R45(1,1) =xt./stress45(1,1);
R45(2,1) =yt./stress45(2,1);
R45(3,1) =sc./stress45(3,1);
M45=zeros(3,1);
M45(1,1)=1./R45(1,1);
M45(2,1)=1./R45(2,1);
M45(3,1)=abs(1./R45(3,1));

Rm45 = zeros (3,1);
Rm45(1,1) =xt./stressm45(1,1);
Rm45(2,1) =yt./stressm45(2,1);
Rm45(3,1) =sc./stressm45(3,1);
Mm45=zeros(3,1);
Mm45(1,1)=1./Rm45(1,1);
Mm45(2,1)=1./Rm45(2,1);
Mm45(3,1)=1./Rm45(3,1);

R90 = zeros (3,1);
R90(1,1) =xc./stress90(1,1);
R90(2,1) =yt./stress90(2,1);
R90(3,1) =sc./stress90(3,1);
M90=zeros(3,1);
M90(1,1)=1./R90(1,1);
M90(2,1)=1./R90(2,1);
M90(3,1)=1./R90(3,1);

disp('R0 : ')
disp(R0)
disp('M0 : ')
disp(M0)
disp('R45 : ')
disp(R45)
disp('M45 : ')
disp(M45)

disp('R -45 : ')
disp(Rm45)
disp('M -45 : ')
disp(Mm45)

disp('R90 : ')
disp(R90)
disp('M90 : ')
disp(M90)
%% the maximum criterion failure is the inverse of the reserve factor 


%%
%QUESTION 4:




%%
%QUESTION 5:

EL=131*10^(3); %GPa
ET5=0.0001; %GPa
GLT5=0.0001;%GPa
VLT5=0.0001;

Sbr=[1/EL -VLT5/EL 0;
    -VLT5/EL 1/ET5 0;
   0 0 1/GLT5 ]

Qbr =inv(Sbr)
%%
%QUESTION 6:
Qbr90 = Ts90*Qbr*Teinv90;

A_5 = 2*h*Q0+2*h*Q45+2*h*Qm45+2*h*Qbr90
B_5 = zeros(3,3) %Symmetrical laminate
D_5 = (1/3)*(((4*h).^3-(3*h).^3)*Q0 +((3*h).^3-(2*h).^3)*Q45+((2*h).^3-(h).^3)*Qm45+((h).^3-(0).^3)*Qbr90+((0).^3-(-h).^3)*Qbr90+((-h).^3-(-2*h).^3)*Qm45+((-2*h).^3-(-3*h).^3)*Q45+((-3*h).^3-(-4*h).^3)*Q0)

%%
%QUESTION 7:
disp("QUESTION 7:")
N = zeros (3,1);
N(1,1) = 100*1000; 

straing= inv(A_5)*N %strain in global 
strain0_5=zeros (3,1);
strain0_5=straing;
strain0_5(1,1)=straing(1,1);
strain0_5(2,1)=straing(2,1);
strain0_5(3,1)=straing(3,1)./2;
disp('strain0 : ')
disp(strain0_5)

straing45_5 = Teinv45*strain0_5; %strain in local 45
strain45_5 = straing45_5;
strain45_5(3,1)=straing45_5(3,1)./2;
disp('strain45 : ')
disp(strain45_5)

straingm45_5 = Teinvm45*strain0_5; %strain in local min45
strainm45_5 = straingm45_5;
strainm45_5(3,1)=straingm45_5(3,1)./2;
disp('strain -45 : ')
disp(strainm45_5)

straing90_5= Teinv90*strain0_5; %strain in local 90
strain90_5 = straing90_5;
strain90_5(3,1)=straing90_5(3,1)./2;
disp('strain90 : ')
disp(strain90_5)

stress0_5 = Qi*strain0_5./100
stress45_5 = Qi*strain45_5./100;
stress45_5(3,1)=stress45_5(3,1)*2;
disp('stress 45 :')
disp(stress45_5)

stressm45_5 = Qi*strainm45_5./100;
stressm45_5(3,1)=stressm45_5(3,1)*2;
disp('stress -45 :')
disp(stressm45_5)

stress90_5 = Qbr*strain90_5./100;
stress90_5(3,1)=stress90_5(3,1)*2;
disp('stress90 :')
disp(stress90_5)
%%
%QUESTION 7 suite: determine the maximum failure criterion 
R0_5 = zeros (3,1);
R0_5(1,1) =xt./stress0_5(1,1);
R0_5(2,1) =yt./stress0_5(2,1);
R0_5(3,1) =0;
M0_5=zeros(3,1);
M0_5(1,1)=1./R0_5(1,1);
M0_5(2,1)=1./R0_5(2,1);
M0_5(3,1)=0;

R45_5 = zeros (3,1);
R45_5(1,1) =xt./stress45_5(1,1);
R45_5(2,1) =yt./stress45_5(2,1);
R45_5(3,1) =sc./stress45_5(3,1);
M45_5=zeros(3,1);
M45_5(1,1)=1./R45_5(1,1);
M45_5(2,1)=1./R45_5(2,1);
M45_5(3,1)=abs(1./R45_5(3,1));

Rm45_5 = zeros (3,1);
Rm45_5(1,1) =xt./stressm45_5(1,1);
Rm45_5(2,1) =yt./stressm45_5(2,1);
Rm45_5(3,1) =sc./stressm45_5(3,1);
Mm45_5=zeros(3,1);
Mm45_5(1,1)=1./Rm45_5(1,1);
Mm45_5(2,1)=1./Rm45_5(2,1);
Mm45_5(3,1)=1./Rm45_5(3,1);

R90_5 = zeros (3,1);
R90_5(1,1) =xc./stress90_5(1,1);
R90_5(2,1) =yt./stress90_5(2,1);
R90_5(3,1) =sc./stress90_5(3,1);
M90_5=zeros(3,1);
M90_5(1,1)=1./R90_5(1,1);
M90_5(2,1)=1./R90_5(2,1);
M90_5(3,1)=1./R90_5(3,1);

disp('R0 : ')
disp(R0_5)
disp('M0 : ')
disp(M0_5)

disp('R45 : ')
disp(R45_5)
disp('M45 : ')
disp(M45_5)

disp('R -45 : ')
disp(Rm45_5)
disp('M -45 : ')
disp(Mm45_5)

disp('R90 : ')
disp(R90_5)
disp('M90 : ')
disp(M90_5)

%%
%QUESTION 8: que du texte et de la logique 


%%
%QUESTION 9:We now degrade the elastic properties of the ±45°-plies broken in matrix mode. 

Qbr45 = Ts45*Qbr*Teinv45;
Qbrm45 = Tsm45*Qbr*Teinvm45;
Qbr90 = Ts90*Qbr*Teinv90;

A_9 = 2*h*Q0+2*h*Qbr45+2*h*Qbrm45+2*h*Qbr90
B_9 = zeros(3,3) %Symmetrical laminate
D_9 = (1/3)*(((4*h).^3-(3*h).^3)*Q0 +((3*h).^3-(2*h).^3)*Qbr45+((2*h).^3-(h).^3)*Qbrm45+((h).^3-(0).^3)*Qbr90+((0).^3-(-h).^3)*Qbr90+((-h).^3-(-2*h).^3)*Qbrm45+((-2*h).^3-(-3*h).^3)*Qbr45+((-3*h).^3-(-4*h).^3)*Q0)

%Reduced stiffness matrix of the laminate [unit: MPa] = A/h 
Sred=A_9./(8*h)

%%
%QUESTION 10: 
disp("QUESTION 10:")
N = zeros (3,1);
N(1,1) = 100*1000; 
xc=-1500
yc=-150
straing_10= inv(A_9)*N %strain in global 
strain0_10=zeros (3,1);
strain0_10=straing;
strain0_10(1,1)=straing_10(1,1);
strain0_10(2,1)=straing_10(2,1);
strain0_10(3,1)=straing_10(3,1)./2;
disp('strain0 : ')
disp(strain0_10)

straing45_10 = Teinv45*strain0_10; %strain in local 45
strain45_10 = straing45_10;
strain45_10(3,1)=straing45_10(3,1)./2;
disp('strain45 : ')
disp(strain45_10)

straingm45_10 = Teinvm45*strain0_10; %strain in local min45
strainm45_10 = straingm45_10;
strainm45_10(3,1)=straingm45_10(3,1)./2;
disp('strain -45 : ')
disp(strainm45_10)

straing90_10= Teinv90*strain0_10; %strain in local 90
strain90_10 = straing90_10;
strain90_10(3,1)=straing90_10(3,1)./2;
disp('strain90 : ')
disp(strain90_10)

stress0_10 = Qi*strain0_10./100
stress45_10 = Qbr*strain45_10./100;
stress45_10(3,1)=stress45_10(3,1)*2;
disp('stress 45 :')
disp(stress45_10)

stressm45_10 = Qbr*strainm45_10./100;
stressm45_10(3,1)=stressm45_10(3,1)*2;
disp('stress -45 :')
disp(stressm45_10)

stress90_10 = Qbr*strain90_10./100;
stress90_10(3,1)=stress90_10(3,1)*2;
disp('stress90 :')
disp(stress90_10)


R0_10 = zeros (3,1);
R0_10(1,1) =xt./stress0_10(1,1);
R0_10(2,1) =yc./stress0_10(2,1);
R0_10(3,1) =0;
M0_10=zeros(3,1);
M0_10(1,1)=1./R0_10(1,1);
M0_10(2,1)=1./R0_10(2,1);
M0_10(3,1)=0;

R45_10 = zeros (3,1);
R45_10(1,1) =xt./stress45_10(1,1);
R45_10(2,1) =yt./stress45_10(2,1);
R45_10(3,1) =sc./stress45_10(3,1);
M45_10=zeros(3,1);
M45_10(1,1)=1./R45_10(1,1);
M45_10(2,1)=1./R45_10(2,1);
M45_10(3,1)=abs(1./R45_10(3,1));

Rm45_10 = zeros (3,1);
Rm45_10(1,1) =xt./stressm45_10(1,1);
Rm45_10(2,1) =yt./stressm45_10(2,1);
Rm45_10(3,1) =sc./stressm45_10(3,1);
Mm45_10=zeros(3,1);
Mm45_10(1,1)=1./Rm45_10(1,1);
Mm45_10(2,1)=1./Rm45_10(2,1);
Mm45_10(3,1)=1./Rm45_10(3,1);

R90_10 = zeros (3,1);
R90_10(1,1) =xc./stress90_10(1,1);
R90_10(2,1) =yt./stress90_10(2,1);
R90_10(3,1) =sc./stress90_10(3,1);
M90_10=zeros(3,1);
M90_10(1,1)=1./R90_10(1,1);
M90_10(2,1)=1./R90_10(2,1);
M90_10(3,1)=1./R90_10(3,1);

disp('R0 : ')
disp(R0_10)
disp('M0 : ')
disp(M0_10)

disp('R45 : ')
disp(R45_10)
disp('M45 : ')
disp(M45_10)

disp('R -45 : ')
disp(Rm45_10)
disp('M -45 : ')
disp(Mm45_10)

disp('R90 : ')
disp(R90_10)
disp('M90 : ')
disp(M90_10)
