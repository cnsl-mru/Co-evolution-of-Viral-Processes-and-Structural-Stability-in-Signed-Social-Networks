%% Figure 9

clear all;
clc,%close all;
max_nodes=750;
N=100;  %number of simulations
actual_time=400; %time for simulation (number of steps in sim = actual time*delta_t)
points=20; %how many samples to record
p0=0.15; %initial density of infected nodes
r0=0.25; %ratio of friendly/total edges
b0=8; %infection+alerting rate b=b0-k
k_s=linspace(0,b0,5); %alerting rate
a0=0.3; %coef. for decreased infection rate ba=a0*b
d0=9; %recovery rate
alpha=0.5; %coef. for energy function
delta_t=0.01;
dataname=["soc-sign-Slashdot081106.txt","soc-sign-bitcoinotc.txt","senat.txt"];

p_av_final=zeros(3,5);
a_av_final=zeros(3,5);
s_av_final=zeros(3,5);
energy_av_final=zeros(3,5);
energy_pair_av_final=zeros(3,5);
energy_delta_av_final=zeros(3,5);
sum_friendly_final=zeros(3,5);
Bal_tri_sum_final=zeros(3,5);

for dataset=1:3

    M = dlmread(dataname(dataset));%loading the graph parameters


[nRow,nCol] = size(M);
if nRow == nCol %to convert different types of datasets
    G=graph(M);
else
    G = graph(M(:,1)+1,M(:,2)+1,M(:,3));
end

G2=simplify(G,'min');  %  returns a graph without multiple edges or self-loops

if length(M)>max_nodes
    nodes=1:max_nodes;
else
    nodes=1:length(M);
end
    
%%
figure
H = subgraph(G2,nodes);
H=rmnode(H, find(degree(H)==0));
Neg_Edges=sum(H.Edges.Weight==-1);
H.Edges.Weight(:)=-1;
r0_vector=sort(randperm(numedges(H),floor(r0*numedges(H))));
H.Edges.Weight(r0_vector)=1;
Neg_Edges=sum(H.Edges.Weight==-1);
hin=plot(H,'Layout','force');
highlight(hin,'Edges',find(H.Edges.Weight==-1),'EdgeColor', 'r')
B=full(adjacency(H,'weighted'));%edge matrix
n=length(B);
%%


Num_Tri=0; %Number of triangles in the graph
Num_Pairs=numedges(H); %Number of pairs in the graph
myName= mfilename
time=actual_time/delta_t; %time for evolution
down=time/points; %downsampling
cut=time/down; % from 1 until down*cut (1000*30=3*10^4) (to see whole line cut=time/down)

%data collection
S_save=ones(n,length(k_s));
A_save=ones(n,n,length(k_s));
pdif=ones(N,time/down,length(k_s));
adif=ones(N,time/down,length(k_s));
energy_func=ones(N,time/down,length(k_s));
energy_delta=ones(N,time/down,length(k_s));
energy_pair=ones(N,time/down,length(k_s));
sum_friendly=zeros(length(k_s),time/down);
Bal_tri_sum=zeros(length(k_s),time/down);
%% Setting

for p=1:length(k_s)
    k=k_s(p)*delta_t; %alerting rate
    b=(b0-k_s(p))*delta_t; %infection rate
    ba=b*a0;
    d=d0*delta_t;
    for o=1:N
        fprintf('complete: %.1f \n', (o+(p-1)*(N))/(length(k_s)*N)*100)%Number of simulation
        A=B;
        S = ones(n,1); % Creation of the Identity matrix for agent’s state,
        p_old=0;
        while p_old<p0 %Initial infected nodes
            i=randi(n);
            S(i)=-1;
            p_old=sum(S(:)==-1)/n;
        end
        
        Bal_tri=0;
        Edelta_sum=0;
        Num_Tri=0;
        for i=1:(n-2)
            for j=i+1:(n-1)
                for l=j+1:n
                    Etriad=-A(i,j)*A(i,l)*A(j,l);
                    if Etriad~=0
                        Num_Tri=Num_Tri+1; %Number of triangles in the graph
                    end
                    if Etriad==1
                        Bal_tri=Bal_tri+1;
                    end
                    Edelta_sum=Etriad+Edelta_sum;
                end
            end
        end
        
        Ep_sum=0;
        for i=1:(n-1)
            for j=i+1:n
                if mod((S(i)+S(j)),2)
                    Ep=A(i,j)*(1-S(i)-S(j))/2;
                else
                    Ep=A(i,j)*((S(i)-S(j))^2)/4;
                end
                if Ep==0
                    Ep=0;
                end
                Ep_sum=Ep+Ep_sum;
            end
        end
        E_old_delta=Edelta_sum/Num_Tri;
        E_old_pair=Ep_sum/Num_Pairs;
        E_old=alpha*E_old_delta+(1-alpha)*E_old_pair; %Total energy
        %% Simulation
        t=0; %timestep
        Energy=zeros(1,time);% Energy for 1 iteration
        p_new=zeros(1,time);% infection density for 1 iteration
        a_dif=zeros(1,time);% alerting density for 1 iteration
        Energy_delta=zeros(1,time);% Energy for 1 iteration
        Energy_pair=zeros(1,time);% Energy for 1 iteration
        friendly_edges=zeros(1,time);% alerting density for 1 iteration
        Bal_tri_t=zeros(1,time);
        while t<time % Simulate until time limit
            
            t=t+1; %time step
            S_old=S; %Saving old configuration
            A_old=A; %Saving old configuration
            x=randi(n); %Random selection of node
            y=randi(n);
            while A(x,y)==0 %Preventing selection of the same node
                x=randi(n);
                y=randi(n);
            end
            
            r=rand;
            if sign(S(x)+0.1)==sign(S(y)+0.1) %for S[+-]S A[+-]A I[+-]I S[+-]A
                if S(x)==-1  %% I+I || I-I
                    if  A(x,y)==1 %% I+I
                        if r<d*(1-d)   % S+I || S-I
                            S(x)=1;
                        elseif r<2*d*(1-d) % S+I || I-S
                            S(y)=1;
                        else %I-I || I+I
                            A(x,y)=-A(x,y);
                            A(y,x)=-A(y,x);
                        end
                    else %% I-I
                        if r<d   % S-I
                            S(x)=1;
                        elseif r<2*d % I-S
                            S(y)=1;
                        else % I+I
                            A(x,y)=-A(x,y);
                            A(y,x)=-A(y,x);
                        end
                    end
                    
                else % change relation for S[+-]S A[+-]A S[+-]A
                    if ((S(x)==1  && S(y)==0) || (S(x)==0  && S(y)==1)) && A(x,y)==1  % S+A //(S(x)+S(y))*A(x,y)==1
                        if r<(k*1)*(1-ba*1)   % A+A
                            S(x)=0;
                            S(y)=0;
                        else %S-A
                            A(x,y)=-A(x,y);
                            A(y,x)=-A(y,x);
                        end
                    else %  S[+-]S A[+-]A S-A
                        A(x,y)=-A(x,y);
                        A(y,x)=-A(y,x);
                    end
                end
            else %% S[+-]I A[+-]I
                if A(x,y)==1 %% S+I A+I
                    if S(x)==1 || S(y)==1 % S+I
                        if r<(k*1)*(1-d) % A+I
                            if S(x)==1
                                S(x)=0;
                            else
                                S(y)=0;
                            end
                            %fprintf("here1 \n")
                        elseif r<(b*1)*(1-d)+(k*(1-d))*1 % I+I
                            S(x)=-1;
                            S(y)=-1;
                            %fprintf("here2 \n")
                        elseif r<(d*(1-(b*1)-(k*1)))+(b*1)*(1-d)+(k*1)*(1-d) % S+S
                            S(x)=1;
                            S(y)=1;
                            %fprintf("here3 \n")
                        else % S-I
                            A(x,y)=-A(x,y);
                            A(y,x)=-A(y,x);
                            %fprintf("here4 \n")
                        end
                    else % A+I
                        if r<(ba*1)*(1-d) % I+I
                            S(x)=-1;
                            S(y)=-1;
                        elseif r<d*(1-(ba*1))+(ba*1)*(1-d) % A+S
                            if S(x)==0
                                S(x)=0;
                                S(y)=1;
                            else
                                S(x)=1;
                                S(y)=0;
                            end
                        else %A-I
                            A(x,y)=-A(x,y);
                            A(y,x)=-A(y,x);
                        end
                    end
                else %% S-I || A-I
                    if r<d
                        if S(y)==-1  % S-I || A-I
                            S(y)=1; % S-S || A-S
                        else        % I-S || I-A
                            S(x)=1; % S-S || S-A
                        end
                    else
                        A(x,y)=-A(x,y);
                        A(y,x)=-A(y,x); % S+I || A+I
                    end
                end
            end
            [A,S,E_old,E_old_delta,E_old_pair,Bal_tri]=compare(A,S,E_old,n,Num_Tri,Num_Pairs,S_old,A_old,alpha,x,y,E_old_delta,E_old_pair,Bal_tri); %calculating Energy for new configuration and receive old/new configuration
            p_new(t)=length(find(S==-1))/n; %calculating  infection density
            Energy(t)=E_old; %collecting energy
            Energy_delta(t)=E_old_delta; %collecting energy
            Energy_pair(t)=E_old_pair; %collecting energy
            a_dif(t)=length(find(S==0))/n; %calculating  alerting density
            friendly_edges(t)=length(find(A==1))/2;
            Bal_tri_t(t)=Bal_tri; %balanced triads
        end
        p_new=downsample(p_new,down); %downsampling array of infection desity
        a_dif=downsample(a_dif,down); %downsampling array of alerting desity
        Energy=downsample(Energy,down); %downsampling array of energy
        Energy_delta=downsample(Energy_delta,down); %downsampling array of energy
        Energy_pair=downsample(Energy_pair,down); %downsampling array of energy
        adif(o,:,p)=a_dif; %store array (alert. dens.) for corresponding attempt and k_s
        pdif(o,:,p)=p_new; %store array (infect. dens.) for corresponding attempt and k_s
        energy_func(o,:,p)=Energy; %store array (energy) for corresponding attempt and k_s
        energy_delta(o,:,p)=Energy_delta; %store array (energy) for corresponding attempt and k_s
        energy_pair(o,:,p)=Energy_pair; %store array (energy) for corresponding attempt and k_s
        
        %fprintf('Number of friednly edges: %f \n', friendly_edges(t))
        friendly_edges=downsample(friendly_edges,down);
        Bal_tri_t=downsample(Bal_tri_t,down);
        Bal_tri_sum(p,:)=Bal_tri_sum(p,:)+Bal_tri_t;
        sum_friendly(p,:)=sum_friendly(p,:)+friendly_edges;
    end
    A_save(:,:,p)=A;
    S_save(:,p)=S;
    Bal_tri_sum(p,:)=Bal_tri_sum(p,:)/N;
    sum_friendly(p,:)=sum_friendly(p,:)/N;
end
fprintf('dataset completed: %.0f \n', dataset)%Number of simulation
p_av=sum(pdif)/N;
a_av=sum(adif)/N;
s_av=1-p_av-a_av;
energy_av=sum(energy_func)/N;
energy_pair_av=sum(energy_pair)/N;
energy_delta_av=sum(energy_delta)/N;



p_av_final(dataset,:)=p_av(:,end,:);
a_av_final(dataset,:)=a_av(:,end,:);
s_av_final(dataset,:)=1-p_av_final(dataset,:)-a_av_final(dataset,:);
energy_av_final(dataset,:)=energy_av(:,end,:);
energy_pair_av_final(dataset,:)=energy_pair_av(:,end,:);
energy_delta_av_final(dataset,:)=energy_delta_av(:,end,:);
sum_friendly_final(dataset,:)=sum_friendly(:,end)/Num_Pairs;
Bal_tri_sum_final(dataset,:)=Bal_tri_sum(:,end)/Num_Tri;
end

%%  Plotting
%%r_inf
figure
X = categorical({'(8,0)','(6,2)','(4,4)','(2,6)','(0,8)'});
X = reordercats(X,{'(8,0)','(6,2)','(4,4)','(2,6)','(0,8)'});
bar(X,sum_friendly_final)
xlabel('(\beta, \kappa)');
ylabel('r_\infty');

legend({'SL','BC','CS'}) 
% %%energy_delta
% figure
% bar(X,energy_delta_av_final)
% xlabel('Infection and alerting rates');
% ylabel("Triad's energy");
% legend({'Slashdot','Bitcoin','Senat'}) 
% %%energy_pair
% figure
% bar(X,energy_pair_av_final)
% xlabel('Infection and alerting rates');
% ylabel('Pairwise energy');
% legend({'Slashdot','Bitcoin','Senat'}) 
%%energy_av
figure
bar(X,energy_av_final)
xlabel('(\beta, \kappa)');
ylabel('Network Energy');
legend({'SL','BC','CS'}) 
% %% Susceptible
% figure
% bar(X,s_av_final)
% xlabel('Infection and alerting rates');
% ylabel('Susceptible density');
% legend({'Slashdot','Bitcoin','Senat'}) 
%% Alerting
figure
bar(X,a_av_final)
xlabel('(\beta, \kappa)');
ylabel('Alerting density');
legend({'SL','BC','CS'}) 
%% Infection
figure
bar(X,p_av_final)
xlabel('(\beta, \kappa)');
ylabel('Infected Density (\rho_{\infty})');
legend({'SL','BC','CS'})
%% Balanced Triads
figure
bar(X,Bal_tri_sum_final)
xlabel('(\beta, \kappa)');
ylabel('# of Balanced triads');
legend({'SL','BC','CS'}) 
