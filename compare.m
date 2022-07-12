%% Function for comparison total energy H

function [A,S,E_old,E_old_delta,E_old_pair,Bal_tri]=compare(A,S,E_old,n,Num_Tri,Num_Pairs,S_old,A_old,alpha,x,y,E_old_delta,E_old_pair,Bal_tri_old)
Edelta_sum=0;
Bal_tri=0;
for i=1:(n-2)
    for j=i+1:(n-1)
        for l=j+1:n
            Etriad=-A(i,j)*A(i,l)*A(j,l);
            if Etriad==-1
                Bal_tri=Bal_tri+1;
            end
            Edelta_sum=Etriad+Edelta_sum;
        end
    end
end
delta_Edelta=Edelta_sum/Num_Tri-E_old_delta;
if A(x,y)~=A_old(x,y) || S_old(x)~=S(x) || S_old(y)~=S(y)
    Exy_old=0;
    Exy_new=0;
    for i=[x,y]
        Ep_sum_old=0;
        Ep_sum_new=0;
        for j=1:n
            if ((i==x && j~=x) || (i==y && j~=y)) && ~(i==x && j==y)
                if mod((S_old(i)+S_old(j)),2)
                    Ep_old=A_old(i,j)*(1-S_old(i)-S_old(j))/2;
                else
                    Ep_old=A_old(i,j)*((S_old(i)-S_old(j))^2)/4;
                end
                Ep_sum_old=Ep_old+Ep_sum_old;
                if mod((S(i)+S(j)),2)
                    Ep_new=A(i,j)*(1-S(i)-S(j))/2;
                else
                    Ep_new=A(i,j)*((S(i)-S(j))^2)/4;
                end
                Ep_sum_new=Ep_new+Ep_sum_new;
            end
        end
        Exy_old=Exy_old+Ep_sum_old;
        Exy_new=Exy_new+Ep_sum_new;
    end
    delta_Ep=(Exy_new-Exy_old)/Num_Pairs;
else
    delta_Ep=0;
end
delta_E=alpha*delta_Edelta+(1-alpha)*delta_Ep; %Total energy
if delta_E>=0
    if delta_E==0
        if rand<0.5
            S=S_old;
            A=A_old;
            Bal_tri=Bal_tri_old;
            %fprintf('old config1\n')
        else
            E_old_delta=E_old_delta+delta_Edelta;
            E_old_pair=E_old_pair+delta_Ep;
            % fprintf('new config2\n')
        end
    else
        S=S_old;
        A=A_old;
        Bal_tri=Bal_tri_old;
        % fprintf('old config\n')
    end
else
    E_old=E_old+delta_E;
    E_old_delta=E_old_delta+delta_Edelta;
    E_old_pair=E_old_pair+delta_Ep;
    %fprintf('NEW\n')
end
end