%% Figure 5

clear all;
clc,%close all;

n=180; %number of nodes
r0=0.25; %ratio of friendly/total edges
N=100;  %number of simulations
actual_time=400; %time for simulation (number of steps in sim = actual time*delta_t)
points=20; %how many samples to record
p0=0.15; %initial density of infected nodes
b0=8; %infection+alerting rate b=b0-k
k_s=[0,2,4,6,8];  %alerting rate
a0=0.3; %coef. for decreased infection rate ba=a0*b
d0=9; %recovery rate
alpha=0.5; %coef. for energy function
delta_t=0.01;


Num_Tri=nchoosek(n,3); %Number of triangles in the graph
Num_Pairs=nchoosek(n,2); %Number of pairs in the graph
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
%creating Signed network
A = zeros(n,n);
r0_vector=sort(randperm(n*(n-1)/2,floor(r0*n*(n-1)/2))); %positive edges
counter=0;
r0_item=1;
for i=1:n %Initial type of edges
    for j=i:n
        if i~=j
            counter=counter+1;
            if ~isempty(r0_vector)
                if counter==r0_vector(r0_item)
                    A(i,j)=1; %
                    A(j,i)=A(i,j);
                    if r0_item==length(r0_vector)
                    else
                        r0_item=r0_item+1;
                    end
                else
                    A(i,j)=-1;
                    A(j,i)=A(i,j);
                end
            else
                A(i,j)=-1;
                A(j,i)=A(i,j);
            end
        end
    end
end
A_init=A; %save this configuration

for p=1:length(k_s)
    k=k_s(p)*delta_t; %alerting rate
    b=(b0-k_s(p))*delta_t; %infection rate
    ba=b*a0;
    d=d0*delta_t;
    for o=1:N
        fprintf('complete: %.1f \n', (o+(p-1)*(N))/(length(k_s)*N)*100)
        A=A_init; %Number of simulation
        S = ones(n,1); % Creation of the Identity matrix for agent’s state,
        p_old=0;
        while p_old<p0 %Initial infected nodes
            i=randi(n);
            S(i)=-1;
            p_old=sum(S(:)==-1)/n;
        end
        
        Bal_tri=0;
        Edelta_sum=0;
        for i=1:(n-2)
            for j=i+1:(n-1)
                for l=j+1:n
                    Etriad=-A(i,j)*A(i,l)*A(j,l);
                    Edelta_sum=Etriad+Edelta_sum;
                    if Etriad==1
                        Bal_tri=Bal_tri+1;
                    end
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
        Bal_tri_t=zeros(1,time);%#Balanced triads
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
                        if r<k*(1-ba)   % A+A
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
                        if r<k*(1-d) % A+I
                            if S(x)==1
                                S(x)=0;
                            else
                                S(y)=0;
                            end
                            %fprintf("here1\n")
                        elseif r<b*(1-d)+k*(1-d) % I+I
                            S(x)=-1;
                            S(y)=-1;
                            %fprintf("here2 \n")
                        elseif r<(d*(1-(b)-(k)))+(b)*(1-d)+(k)*(1-d) % S+S
                            S(x)=1;
                            S(y)=1;
                            %fprintf("here3 \n")
                        else % S-I
                            A(x,y)=-A(x,y);
                            A(y,x)=-A(y,x);
                            %fprintf("here4 \n")
                        end
                    else % A+I
                        if r<(ba)*(1-d) % I+I
                            S(x)=-1;
                            S(y)=-1;
                        elseif r<d*(1-(ba))+(ba)*(1-d) % A+S
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
            Bal_tri_t(t)=Bal_tri;
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
        Bal_tri_t=downsample(Bal_tri_t,down);
        Bal_tri_sum(p,:)=Bal_tri_sum(p,:)+Bal_tri_t;
        %fprintf('Number of friednly edges: %f \n', friendly_edges(t))
        friendly_edges=downsample(friendly_edges,down);

        sum_friendly(p,:)=sum_friendly(p,:)+friendly_edges;
    end
    A_save(:,:,p)=A;
    S_save(:,p)=S;
    sum_friendly(p,:)=sum_friendly(p,:)/N;
    Bal_tri_sum(p,:)=Bal_tri_sum(p,:)/N;
end
%% Plotting the results
%close all;
figure
sub=1;
%ti=time/down;
for p=1:length(k_s)
    B=A_save(:,:,p); %Matrix for highlighting edges
    for i= 1:n
        for j=1:n
            if B(i,j)==-1
                B(i,j)=0;
            end
        end
    end
    Inodes=find(S_save(:,p)==-1);
    Anodes=find(S_save(:,p)==0);
    Snodes=find(S_save(:,p)==1);
    subplot(2,length(k_s),sub);
    sub=sub+1;
    G=graph(A); %Original Graph
    T=graph(B); %Friendly connections
    h=plot(T,'Layout','circle');
    title(['k=',num2str(k_s(p)),', b=',num2str(b0-k_s(p)),', % of Inf=',num2str(length(find(S_save(:,p)==-1))/n*100),'%']);
    highlight(h,T,'LineStyle','-')
    highlight(h,Inodes,'NodeColor','r')
    highlight(h,Anodes,'NodeColor','g')
    highlight(h,Snodes,'NodeColor','b')
    subplot(2,length(k_s),sub+length(k_s)-1);
    % sub=sub+1;
    h1=plot(T,'Layout','force');
    title(['k=',num2str(k_s(p)),', b=',num2str(b0-k_s(p)),', % of Inf=',num2str(length(find(S_save(:,p)==-1))/n*100),'%']);
    highlight(h1,T,'LineStyle','-')
    highlight(h1,Inodes,'NodeColor','r')
    highlight(h1,Anodes,'NodeColor','g')
    highlight(h1,Snodes,'NodeColor','b')
end
sgtitle(['time=',num2str(actual_time),', d=',num2str(d0),', alpha=',num2str(alpha),', r0=',num2str(r0),', p0=',num2str(p0)]);

%%
figure
t = tiledlayout(2,2);
t.Padding = 'none';
t.TileSpacing = 'none';
%% Infection density
nexttile
% p_av=ones(1,time,length(b));
p_av=sum(pdif)/N;
signedge='random';
t=linspace(0,actual_time,length(p_new));
hold on
%set(gca, 'XScale', 'log')
marker=["--*","-o",":d","-.x","--^","-.",":+","-.o","--d","-x",":^","-.."];
for i=1:length(k_s)
    plot(t(1:cut),(p_av(:,1:cut,i)),marker(i),'DisplayName',['b=',num2str(b0-k_s(i)),' k=',num2str(k_s(i))])
end
%title({['Infection density'],['n = ',num2str(n),' alpha = ',num2str(alpha),'; p0 = ' num2str(p0),'; b0 = ',num2str(b0),'; ba =b*k ','; d = ',num2str(d),'.'],[' Initial edges have a ', signedge,  ' sign'],['Number of simulations: ',num2str(N)]})
xlabel('time')
ylabel('Infection density')
legend
box on
hold off


%% Alerted density
nexttile
% a_av=ones(1,time,length(b));
a_av=sum(adif)/N;
hold on
%set(gca, 'XScale', 'log')
for i=1:length(k_s)
    plot(t(1:cut),(a_av(:,1:cut,i)),marker(i),'DisplayName',['b=',num2str(b0-k_s(i)),' k=',num2str(k_s(i))])
end
%title({['Density of alerted nodes'],['n = ',num2str(n),' alpha = ',num2str(alpha),'; p0 = ' num2str(p0),'; b0 = ',num2str(b0),'; ba =b*k ','; d = ',num2str(d),'.'],[' Initial edges have a ', signedge,  ' sign'],['Number of simulations: ',num2str(N)]})
xlabel('time')
ylabel('Density of alerted nodes')
legend
box on
hold off
%% Susceptible density
nexttile
s_av=1-p_av-a_av;
hold on
%set(gca, 'XScale', 'log')
for i=1:length(k_s)
    plot(t(1:cut),(s_av(:,1:cut,i)),marker(i),'DisplayName',['b=',num2str(b0-k_s(i)),' k=',num2str(k_s(i))])
end
%title({['Density of Susceptible nodes'],['n = ',num2str(n),' alpha = ',num2str(alpha),'; p0 = ' num2str(p0),'; b0 = ',num2str(b0),'; ba =b*k ','; d = ',num2str(d),'.'],[' Initial edges have a ', signedge,  ' sign'],['Number of simulations: ',num2str(N)]})
xlabel('time')
ylabel('Susceptible density')
legend
box on
hold off

%% Friendly edges
nexttile
hold on
for i=1:length(k_s)
    plot(t(1:cut),(sum_friendly(i,1:cut)),marker(i),'DisplayName',['b=',num2str(b0-k_s(i)),' k=',num2str(k_s(i))])
end
%title({['Infection density'],['n = ',num2str(n),' alpha = ',num2str(alpha),'; p0 = ' num2str(p0),'; b0 = ',num2str(b0),'; ba =b*k ','; d = ',num2str(d),'.'],[' Initial edges have a ', signedge,  ' sign'],['Number of simulations: ',num2str(N)]})
xlabel('time')
ylabel('Friendly edges')
legend
box on
hold off


%% Energy
figure
t = tiledlayout(1,3);
t.Padding = 'none';
t.TileSpacing = 'none';
t=linspace(0,actual_time,length(p_new));
nexttile
energy_av=sum(energy_func)/N;
hold on
for i=1:length(k_s)
    plot(t(1:cut),(energy_av(:,1:cut,i)),marker(i),'DisplayName',['b=',num2str(b0-k_s(i)),' k=',num2str(k_s(i))])
end
%title({['Infection density'],['n = ',num2str(n),' alpha = ',num2str(alpha),'; p0 = ' num2str(p0),'; b0 = ',num2str(b0),'; ba =b*k ','; d = ',num2str(d),'.'],[' Initial edges have a ', signedge,  ' sign'],['Number of simulations: ',num2str(N)]})
xlabel('time')
ylabel('Energy')
legend
box on
hold off

%% Energy pair
nexttile
energy_av=sum(energy_pair)/N;
hold on
for i=1:length(k_s)
    plot(t(1:cut),(energy_av(:,1:cut,i)),marker(i),'DisplayName',['b=',num2str(b0-k_s(i)),' k=',num2str(k_s(i))])
end
%title({['Infection density'],['n = ',num2str(n),' alpha = ',num2str(alpha),'; p0 = ' num2str(p0),'; b0 = ',num2str(b0),'; ba =b*k ','; d = ',num2str(d),'.'],[' Initial edges have a ', signedge,  ' sign'],['Number of simulations: ',num2str(N)]})
xlabel('time')
ylabel('Energy pair')
legend
box on
hold off


%% Energy delta
nexttile
energy_av=sum(energy_delta)/N;
hold on
for i=1:length(k_s)
    plot(t(1:cut),(energy_av(:,1:cut,i)),marker(i),'DisplayName',['b=',num2str(b0-k_s(i)),' k=',num2str(k_s(i))])
end
%title({['Infection density'],['n = ',num2str(n),' alpha = ',num2str(alpha),'; p0 = ' num2str(p0),'; b0 = ',num2str(b0),'; ba =b*k ','; d = ',num2str(d),'.'],[' Initial edges have a ', signedge,  ' sign'],['Number of simulations: ',num2str(N)]})
xlabel('time')
ylabel('Energy delta')
legend
box on
hold off
