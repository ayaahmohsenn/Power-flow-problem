clc;
clear;

N=5;            %%Total no of buses
Ng=1;           %%No of PV buses
N_PQ=3;         %%No of PQ buses
z=(2*N-Ng-2);    %%dimensions of Jackobian
mism=[1];        %%initial mismatches


%%1.Bus no.   2.P(G)    3.Q(G)   4.P(L)    5.Q(L)    6.V      7.angle     8.bus type
bustable=[1,    0,       0,      0,      0,      1.00,    0,     1;
    2,     0,       0,      8,    2.8,      1.00,     0,    2;
    3,   5.2,       0,    0.8,    0.4,      1.05,     0,    3;
    4,     0,       0,      0,      0,      1.00,      0,   2;
    5,     0,       0,      0,      0,      1.00,      0,   2];
%%BUS ADMITTANCE MATRIX
Ybus1=[3.73-j*49.72,     0,               0,           0,         -3.73+j*49.72;
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


%%OLTC at t=0.99 and t=1.01 at 1, 5 buses
OLTC=zeros(2,2);
Ybus=zeros(5,5);
%t=1.01;
t=0.99;
OLTC=[-(t)^2*Ybus1(1,5)   t*Ybus1(1,5);
       t*Ybus1(1,5)         -Ybus1(1,5)];
Ybus1(1,1)=OLTC(1,1);
Ybus1(1,5) = OLTC(1,2);
Ybus1(5,1)= OLTC(2,1);
Ybus1(5,5)=9.0856-j*108.58;

Ybus=Ybus1;
output=zeros(N,3);
G=real(Ybus);
B=imag(Ybus);
Ybus_value=abs(Ybus);
Ybus_Angle=angle(Ybus);

while sum(abs(mism))>0.1
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
Q_calculated = Q_calculated*100;
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
%%PV bus power
P_pv=p_calculated(3);
Q_pv=Q_calculated(3);
for i = 1:5
    VV(i) = V(i)*(cos(ANGLE(i))+1i*sin(ANGLE(i)));
end
%%current at each branch
I_base_1 = 100/(sqrt(3)*15);
I_base_2 = 100/(sqrt(3)*345);
for b1 = 1:5
    for b2 = 1:5
        I(b1,b2) = (VV(b1)-VV(b2))*Ybus(b1,b2);
    end
end
I(:,1) = abs(I(:,1))*I_base_1/I_base_2;
I = abs(I) * I_base_2*1000;
%%system losses
P_losses = sum(p_calculated)+p_slack;