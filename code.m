clc;
clear;

N=5;            %%Total no of buses
Ng=1;           %%No of PV buses
N_PQ=3;         %%No of PQ buses
z=(2*N-Ng-2);        %%dimensions of Jackobian
mism=[1];        %%initial mismatches


%%1.Bus no.   2.P(G)    3.Q(G)   4.P(L)    5.Q(L)    6.V      7.angle     8.bus type
bustable=[1,    0,       0,      0,      0,      1.00,    0,     1;
    2,     0,       0,      8,    2.8,      1.00,     0,    2;
    3,   5.2,       0,    0.8,    0.4,      1.05,     0,    3;
    4,     0,       0,      0,      0,      1.00,      0,   2;
    5,     0,       0,      0,      0,      1.00,      0,   2];
%%BUS ADMITTANCE MATRIX
Ybus=[3.73-j*49.72,     0,               0,           0,         -3.73+j*49.72;
    0,             2.68-j*28.46,        0,      -0.89+j*9.92,   -1.79+j*19.84;
    0,                 0,        7.46-j*99.44,  -7.46+j*99.44,              0;
    0,         -0.89+j*9.92,   -7.46+j*99.44, 11.92-j*147.96,  -3.57+j*39.68;
   -3.73+j*49.72, -1.79+j*19.84           0,      -3.57+j*39.68,  9.09-j*108.58];

%%Data from bus input data table
busnumber=bustable(:,1);
P_G=bustable(:,2);
Q_G=bustable(:,3);
P_L=bustable(:,4);
Q_L=bustable(:,5);
P_scheduled= (P_G-P_L);
Q_scheduled= (Q_G-Q_L);
V=bustable(:,6);
ANGLE=bustable(:,7);
Angle=ANGLE;
v=V;


output=zeros(N,3);
G=real(Ybus);
B=imag(Ybus);
Ybus_Angle=angle(Ybus);

while sum(abs(mism))>0.001
    %% delta P & delta Q
    for i = 2:N
        P_sum=0;
        Q_sum=0;
        for j= 1:N
            P_sum=P_sum+abs(abs(Ybus(i,j))*V(j)*V(i))*cos((Ybus_Angle(i,j))+bustable(j,7)-bustable(i,7));
            Q_sum=Q_sum+(abs((abs(Ybus(i,j)))*V(j)*V(i))*sin((Ybus_Angle(i,j))+bustable(j,7)-bustable(i,7)));
        end
        p_calculated(i)=P_sum;
        delta_pp(i)=P_scheduled(i)-p_calculated(i);
        Q_calculated(i)=-Q_sum;
        delta_QQ(i)=Q_scheduled(i)-Q_calculated(i);
    end
    delta_QQ(3)=[];
    delta_p=delta_pp';
    delta_Q=delta_QQ';
    %%putting delta P and delta Q into the mismatches matrix
    for w=2:5
        mism(w-1,1)=delta_p(w);
    end
    for u=2:5
        for i=2:4
            if bustable(u,8)== 2
                mism(N-2+i,1)=delta_Q(i);
            else
                i=i-1;
            end
        end
    end
    %%
    [J_11 J_12 J_21 J_22 Jackobian]= jackobian ( Ybus,V, Ybus_Angle,bustable, G, B, p_calculated,N,N_PQ, Ng, Q_calculated);
    %%
    format shortG
    correction=inv(Jackobian)*mism;
    %%ANGLE
    ANGLE(2:5) = ANGLE(2:5) + correction(1:4);
    bustable(:,7)=ANGLE;
    %%Voltage
    V(2) = V(2) + correction(5);
    V(4:5) = V(4:5) + correction(6:7);
    bustable(:,6)=V;
    v=V;
    Angle=ANGLE;
end
%%Voltage and angle at each bus
output(1:N,1)=bustable(:,1);
output(1:N,2)=abs(bustable(:,6));
output(1:N,3)=radtodeg(bustable(:,7));
disp('        Bus No.          V           Angle');
format shortG
disp(output);
%%PV bus power
p_calculated = p_calculated*100;
P_pv=p_calculated(3);
Q_calculated = Q_calculated*100;
Q_pv=Q_calculated(3);

%%current at each branch
for s = 1:5
    VV(s) = V(s)*(cos(ANGLE(s))+1i*sin(ANGLE(s)));
end
V13=15;
V245=345;
I_base_13 = 100/(sqrt(3)*V13);
I_base_234 = 100/(sqrt(3)*V245);
for b1 = 1:5
    for b2 = 1:5
        I(b1,b2) = (VV(b1)-VV(b2))*Ybus(b1,b2);
    end
end
I(:,1) = abs(I(:,1))*I_base_13/I_base_234;
I = abs(I) * I_base_234*1000;  %%A

%%slack bus power
for i = 1
    P_sum=0;
    Q_sum=0;
    for j= 1:5
        P_sum=P_sum+abs(abs(Ybus(i,j))*V(j)*V(i))*cos((Ybus_Angle(i,j))+bustable(j,7)-bustable(i,7));
        Q_sum=Q_sum+(abs((abs(Ybus(i,j)))*V(j)*V(i))*sin((Ybus_Angle(i,j))+bustable(j,7)-bustable(i,7)));
    end
    p_slack=P_sum*100;
    Q_slack=-Q_sum*100;
end
%%system losses
P_losses = sum(p_calculated)+p_slack;
actualvolt=[0 0 0 0 0];
%%actual voltage at buses
actualvolt(1) = V(1)*15;
actualvolt(2) = V(2)*345;
actualvolt(3) = V(3)*15;
actualvolt(4:5) = V(4:5)*345;
%%  Fault calculations
clear j
 Xd_G1 = 0.2*100/400;
 Xd_G2 = 0.2*100/800;
 YdG1 = 1/Xd_G1;
 YdG2 = 1/Xd_G2;
 Ybus(1,1) = Ybus(1,1) + Xd_G1;
 Ybus(3,3) = Ybus(3,3) + Xd_G2;
 Zbus = Ybus^-1;
%%at bus 3
%%buses voltages
V_3_1=V(3)*(1-(Zbus(1,3)/Zbus(3,3)));
V_3_2=V(3)*(1-(Zbus(2,3)/Zbus(3,3)));
V_3_4=V(3)*(1-(Zbus(4,3)/Zbus(3,3)));
V_3_5=V(3)*(1-(Zbus(5,3)/Zbus(3,3)));

%%buses current
I_3_bus1=V_3_1*Ybus(1,1)+V_3_5*Ybus(1,5);
I_3_bus2=V_3_2*Ybus(2,2)+V_3_4*Ybus(2,4)+V_3_5*Ybus(2,5);
I_3_bus3=V(3)*Ybus(3,3)+V_3_4*Ybus(3,4);
I_3_bus4=V_3_4*Ybus(4,4)+V_3_5*Ybus(4,5)+V_3_2*Ybus(4,2)+V(3)*Ybus(4,3);
I_3_bus5=V_3_5*Ybus(5,5)+V_3_1*Ybus(5,1)+V_3_4*Ybus(5,4)+V_3_2*Ybus(5,2);

%% at Bus 2
%%buses voltages
V_2_1=V(2)*(1-(Zbus(1,2)/Zbus(2,2)));
V_2_4=V(2)*(1-(Zbus(4,2)/Zbus(2,2)));
V_2_5=V(2)*(1-(Zbus(5,2)/Zbus(2,2)));
V_2_3=V(2)*(1-(Zbus(3,2)/Zbus(2,2)));

%%buses current
I_2_bus1=V_2_1*Ybus(1,1)+V_2_5*Ybus(1,5);
I_2_bus2=V(2)*Ybus(2,2)+V_2_4*Ybus(2,4)+V_2_5*Ybus(2,5);
I_2_bus3=V_2_3*Ybus(3,3)+V_2_4*Ybus(3,4);
I_2_bus4=V_2_4*Ybus(4,4)+V_2_5*Ybus(4,5)+V(2)*Ybus(4,2);
I_2_bus5=V_2_5*Ybus(5,5)+V_2_1*Ybus(5,1)+V_2_4*Ybus(5,4)+V(2)*Ybus(5,2);