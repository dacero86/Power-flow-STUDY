clc
close all
clear

%% parcial analisis

V1=1.2*exp(0*i);
V2=1.0*exp(-1.5*i*pi/180);
V3=0.9*exp(-2*i*pi/180);
V4=1.1*exp(-0.5*i*pi/180);

Dp2=-0.048100711598073409973065182903698;
Dp3=-0.13179141972938305701460694649649;
Dp4=0.10778110329790750504515221403167;

Pg1=(abs(V1)*abs(V2)*sin(angle(V1)-angle(V2)))/(0.2)
Pg4=((abs(V4)*abs(V3)*sin(angle(V4)-angle(V3)))/(0.2))
Pd2=((abs(V2)*abs(V1)*sin(angle(V2)-angle(V1)))/(0.2))+((abs(V2)*abs(V3)*sin(angle(V2)-angle(V3)))/(0.2))+Dp2
Pd3=((abs(V3)*abs(V2)*sin(angle(V3)-angle(V2)))/(0.2))+((abs(V3)*abs(V4)*sin(angle(V3)-angle(V4)))/(0.2))+Dp3
%+Dp4
Pl=0.045331413411659;
Pd=Pg1+Pg4-Pl;
E=(Pd)/(Pg1+Pg4);

%% Punto 2
P12=(abs(V1)*abs(V2)*sin(angle(V1)-angle(V2)))/(0.2);
P23=(abs(V2)*abs(V3)*sin(angle(V2)-angle(V3)))/(0.2);
P34=(abs(V4)*abs(V3)*sin(angle(V4)-angle(V3)))/(0.2);

R12=((abs(V1)-abs(V2))^2)/(P12)
R23=((abs(V2)-abs(V3))^2)/(P23)
R34=((abs(V4)-abs(V3))^2)/(Pd3)




