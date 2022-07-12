%%  Code for plotting and saving figures (loading saved workspace)
nameload="n750N10t400p075"
load(nameload+".mat")
X = categorical({'(8,0)','(6,2)','(4,4)','(2,6)','(0,8)'});
X = reordercats(X,{'(8,0)','(6,2)','(4,4)','(2,6)','(0,8)'});
%% Infection
figure
str="\rho_0="+string(p0);
fig=gcf;
fig.Position(3:4)=[4000,600];
tiledlayout(1,5,'TileSpacing','compact','Padding','tight')
nexttile
bar(X,p_av_final, 'BarWidth', 1)
xlabel('(\beta, \kappa)');
ylabel('Infected Density (\rho_{\infty})');
legend({'SL','BC','CS'})
axis tight
set(gca,'fontname','Century')  % Set it to times
set(gca,'FontSize',14)
ylim([0 .7])
annotation('textbox', [0.1, 0.7, 0.2, 0.2], 'String', str,'FitBoxToText','on',"FontSize",14,"BackgroundColor","w")
%% Alerting
nexttile
bar(X,a_av_final, 'BarWidth', 1)
xlabel('(\beta, \kappa)');
ylabel('Alerted Density (\alpha_{\infty})');
legend({'SL','BC','CS'}) 
axis tight
set(gca,'fontname','Century')  % Set it to times
set(gca,'FontSize',14)
ylim([0 1])
legend('Location','northwest')
annotation('textbox', [0.3, 0.7, 0.2, 0.2], 'String', str,'FitBoxToText','on',"FontSize",14,"BackgroundColor","w")
%% r_inf
nexttile
bar(X,sum_friendly_final, 'BarWidth', 1)
xlabel('(\beta, \kappa)');
ylabel('Friendly Links (r_\infty)');
axis tight
legend({'SL','BC','CS'}) 
set(gca,'fontname','Century')  % Set it to times
set(gca,'FontSize',14)
ylim([0 1])
annotation('textbox', [0.5, 0.7, 0.2, 0.2], 'String', str,'FitBoxToText','on',"FontSize",14,"BackgroundColor","w")
%% Balanced Triads
nexttile
bar(X,Bal_tri_sum_final, 'BarWidth', 1)
xlabel('(\beta, \kappa)');
ylabel('Balanced Triads (|\Delta_b|)');
legend({'SL','BC','CS'}) 
axis tight
set(gca,'fontname','Century')  % Set it to times
set(gca,'FontSize',14)
ylim([0 1])
annotation('textbox', [0.7, 0.7, 0.2, 0.2], 'String', str,'FitBoxToText','on',"FontSize",14,"BackgroundColor","w")
%% Network energy
nexttile
bar(X,energy_av_final, 'BarWidth', 1)
xlabel('(\beta, \kappa)');
ylabel('Network Energy (E(G))');
legend({'SL','BC','CS'}) 
axis tight
set(gca,'fontname','Century')  % Set it to times
set(gca,'FontSize',14)
ylim([-1 0])
annotation('textbox', [0.9, 0.7, 0.2, 0.2], 'String', str,'FitBoxToText','on',"FontSize",14,"BackgroundColor","w")
stringsave=nameload+".fig"
saveas(gcf,stringsave)
stringsave=nameload+".png"
saveas(gcf,stringsave)