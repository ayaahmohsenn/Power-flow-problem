function [J_11 J_12 J_21 J_22 Jackobian]= jackobian ( Ybus,V, Ybus_Angle,bustable, G, B, p_calculated, N,N_PQ,Ng, Q_calculated)
z=(2*N-Ng-2);        %%dimensions of Jackobian
Jackobian= zeros(z,z);       %%initial jackobian
J_11=zeros(N-1,N-1);
J_12=zeros(N-1,N_PQ);
J_21=zeros(N_PQ,N-1);
J_22=zeros(N_PQ,N_PQ);
%%  parial_p/partial_delta(J_11)
    for i=2:N
        for j=2:N
            if i==j
                for n=1:N
                    if n~=i
                        J_11(i-1,i-1)=-Q_calculated(i)-((abs(V(i)))^2)*B(i,i);
                    end
                end
            else
                J_11(i-1,j-1)= -(abs(V(i)*V(j)*abs(Ybus(i,j))))*sin((Ybus_Angle(i,j))+bustable(j,7)-bustable(i,7));
            end
        end
    end
    for i=1:N-1
        for j=1:N-1
            Jackobian(i,j)=J_11(i,j);
        end
    end
    %%  partial_p/partial_V(J_12)
    for i=2:N
        for j=[2 4 5]
            if i==j
                J_12(i-1,i-1)=p_calculated(i) +(abs((V(i))))^2*G(i,i);
            else
                J_12(i-1,j-1)= abs(V(i)*V(j)*abs(Ybus(i,j)))*cos((Ybus_Angle(i,j))+bustable(j,7)-bustable(i,7));
            end
        end
    end
    J_12(:,2)=[];
    for i=1:N_PQ
        for j=1:N-1
            Jackobian(j,i+N-1)=J_12(j,i);
        end
    end
    %%    partial_Q/partial_delta(J_21)
    for i=[2 4 5]
        for j=2:N
            if i==j
                for n=1:N
                    if n~=i
                       J_21(i-1,i-1)=p_calculated(i)-(abs(V(i)))^2*G(i,i);
                    end
                end
            else
             J_21(i-1,j-1)= -abs(V(i)*V(j)*abs(Ybus(i,j)))*cos((Ybus_Angle(i,j))+bustable(j,7)-bustable(i,7));  
            end
        end
    end
J_21(2,:)=[];
for i=1:N_PQ
    for j=1:N-1
        Jackobian(i+N-1,j)=J_21(i,j);
    end
end
    %%   Partial_Q/partial_V(J22)
    for i=[2 4 5]
        for j=[2 4 5]
                if i==j
                    J_22(i-1,i-1)=Q_calculated(i)-((abs(V(i)))^2)*B(i,i);
                else
                    J_22(i-1,j-1)= -(abs(V(i)*V(j)*abs(Ybus(i,j)))*sin(angle(Ybus_Angle(i,j))+bustable(j,7)-bustable(i,7)));

                end
        end
    end
    J_22(2,:)=[];
    J_22(:,2)=[];
    for i=1:N_PQ
        for j=1:N_PQ
            Jackobian(i+N-1,j+N-1)=J_22(j,i);
        end
    end