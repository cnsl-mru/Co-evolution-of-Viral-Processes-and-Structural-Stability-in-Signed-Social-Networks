%%  Code for plotting and saving figures (loading saved workspace)
nameload="n750N10t400p075"
load(nameload+".mat")
X = categorical({'(8,0)','(6,2)','(4,4)','(2,6)','(0,8)'});
X = reordercats(X,{'(8,0)','(6,2)','(4,4)','(2,6)','(0,8)'});
%% Infection
figure
bar(X,p_av_final, 'BarWidth', 1)
xlabel('(\beta, \kappa)');
ylabel('Infected Density (\rho_{\infty})');
legend({'SL','BC','CS'})
axis tight
set(gca,'fontname','Century')  % Set it to times
set(gca,'FontSize',14)
ylim([0 .7])

stringsave="Infected Density"+".fig"
saveas(gcf,stringsave)
%% Alerting
figure
bar(X,a_av_final, 'BarWidth', 1)
xlabel('(\beta, \kappa)');
ylabel('Alerted Density (\alpha_{\infty})');
legend({'SL','BC','CS'}) 
axis tight
set(gca,'fontname','Century')  % Set it to times
set(gca,'FontSize',14)
ylim([0 1])
legend('Location','northwest')

stringsave="Alerted Density"+".fig"
saveas(gcf,stringsave)
%% r_inf
figure
bar(X,sum_friendly_final, 'BarWidth', 1)
xlabel('(\beta, \kappa)');
ylabel('Friendly Links (r_\infty)');
axis tight
legend({'SL','BC','CS'}) 
set(gca,'fontname','Century')  % Set it to times
set(gca,'FontSize',14)
ylim([0 1])

stringsave="Friendly Links"+".fig"
saveas(gcf,stringsave)
%% Balanced Triads
figure
bar(X,Bal_tri_sum_final, 'BarWidth', 1)
xlabel('(\beta, \kappa)');
ylabel('Balanced Triads (|\Delta_b|)');
legend({'SL','BC','CS'}) 
axis tight
set(gca,'fontname','Century')  % Set it to times
set(gca,'FontSize',14)
ylim([0 1])

stringsave="Balanced Triads"+".fig"
saveas(gcf,stringsave)
%% Network energy
figure
bar(X,energy_av_final, 'BarWidth', 1)
xlabel('(\beta, \kappa)');
ylabel('Network Energy (E(G))');
legend({'SL','BC','CS'}) 
axis tight
set(gca,'fontname','Century')  % Set it to times
set(gca,'FontSize',14)
ylim([-1 0])

stringsave="Network Energy"+".fig"
saveas(gcf,stringsave)
