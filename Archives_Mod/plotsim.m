% 
% 
% % fpath='C:\Users\sebas\Google Drive\Var_RTBC\Final\RBC\Output\'; % path for saving pdf
% 
% ttt=size(oo_.exo_simul(:,1));   %periods
% t=ttt(1);
% tt=linspace(0,t,t);
% 
% % Exogenous variables
% tau=oo_.exo_simul(:,1);
% 
% % Endogenous variables
% v1=oo_.endo_simul(1,:);
% v2=oo_.endo_simul(2,:);
% v3=oo_.endo_simul(3,:);
% v4=oo_.endo_simul(4,:);
% v5=oo_.endo_simul(5,:);
% v6=oo_.endo_simul(6,:);
% v7=oo_.endo_simul(7,:);
% v8=oo_.endo_simul(8,:);
% v9=oo_.endo_simul(9,:);
% v10=oo_.endo_simul(10,:);
% v11=oo_.endo_simul(11,:);
% v12=oo_.endo_simul(12,:);
% name1=M_.endo_names_long{1};
% name2=M_.endo_names_long{2};
% name3=M_.endo_names_long{3};
% name4=M_.endo_names_long{4};
% name5=M_.endo_names_long{5};
% name6=M_.endo_names_long{6};
% name7=M_.endo_names_long{7};
% name8=M_.endo_names_long{8};
% name9=M_.endo_names_long{9};
% name10=M_.endo_names_long{10};
% name11=M_.endo_names_long{11};
% name12=M_.endo_names_long{12};
% 
% 
% hFig = figure(1);
% % set(hFig, 'Position', [0 20 750 750])
% subplot(3,3,1)
% plot(tt,v1,'LineWidth',2,'Color',[0 0 1])
% axis([0 t -inf inf])
% % hline(0,'-k')
% title(name1)
% box on
% subplot(3,3,2)
% plot(tt,v2,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name2)
% axis([0 t -inf inf])
% box on
% subplot(3,3,3)
% plot(tt,v3,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name3)
% axis([0 t -inf inf])
% box on
% subplot(3,3,4)
% plot(tt,v4,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name4)
% axis([0 t -inf inf])
% box on
% subplot(3,3,5)
% plot(tt,v5,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name5)
% axis([0 t -inf inf])
% box on
% subplot(3,3,6)
% plot(tt,v6,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name6)
% axis([0 t -inf inf])
% box on
% subplot(3,3,7)
% plot(tt,v7,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name7)
% axis([0 t -inf inf])
% box on
% subplot(3,3,8)
% plot(tt,v8,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name8)
% axis([0 t -inf inf])
% box on
% subplot(3,3,9)
% plot(tt,v9,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name9)
% axis([0 t -inf inf])
% box on
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure;
subplot(2,2,1)
plot(20:5:100,[oo_.steady_state(148:164,1)],'LineWidth',2,'Color',[0 0 1])
title('Epargne')
axis([20 100 -inf inf])
subplot(2,2,2)
plot(20:5:100,[oo_.steady_state(165:181,1)],'LineWidth',2,'Color',[0 0 1])
title('Consommation')
axis([20 100 -inf inf])
subplot(2,2,3)
plot(20:5:60,[oo_.steady_state(182:190,1)],'LineWidth',2,'Color',[0 0 1])
title('Offre de travail')
axis([20 60 -inf inf])
subplot(2,2,4)
plot(20:5:60,[oo_.steady_state(191:199,1)],'LineWidth',2,'Color',[0 0 1])
title('Salaires')
axis([20 60 -inf inf])

% 
% v1d=100*(oo_.endo_simul(1,end)-oo_.endo_simul(1,1))/abs(oo_.endo_simul(1,1));
% v2d=100*(oo_.endo_simul(2,end)-oo_.endo_simul(2,1))/abs(oo_.endo_simul(2,1));
% v3d=100*(oo_.endo_simul(3,end)-oo_.endo_simul(3,1))/abs(oo_.endo_simul(3,1));
% v4d=100*(oo_.endo_simul(4,end)-oo_.endo_simul(4,1))/abs(oo_.endo_simul(4,1));
% v5d=100*(oo_.endo_simul(5,end)-oo_.endo_simul(5,1))/abs(oo_.endo_simul(5,1));
% v6d=100*(oo_.endo_simul(6,end)-oo_.endo_simul(6,1))/abs(oo_.endo_simul(6,1));
% v7d=100*(oo_.endo_simul(7,end)-oo_.endo_simul(7,1))/abs(oo_.endo_simul(7,1));
% v8d=100*(oo_.endo_simul(8,end)-oo_.endo_simul(8,1))/abs(oo_.endo_simul(8,1));
% v9d=100*(oo_.endo_simul(9,end)-oo_.endo_simul(9,1))/abs(oo_.endo_simul(9,1));
% v10d=100*(oo_.endo_simul(10,end)-oo_.endo_simul(10,1))/abs(oo_.endo_simul(10,1));
% v11d=100*(oo_.endo_simul(11,end)-oo_.endo_simul(11,1))/abs(oo_.endo_simul(11,1));
% 
% solcs=[v1d v2d v3d v4d v5d v6d v7d v8d v9d v10d v11d]';
% 
% 
% save('solcs3.mat','solcs');
% 
% 
% solcs1=load('solcs1.mat');
% solcs2=load('solcs2.mat');
% solcs3=load('solcs3.mat');
% solcs4=load('solcs4.mat');
% % Table: comparative statics
% tpath='/Users/bock/Google Drive/TDTE/Modeles/Codes/Simple OLG/Output/';
% TCSA = fopen([tpath,'TCSA_SS.tex'], 'w');
% fprintf(TCSA, '\\begin{tabular}{lcccc} \\hline \\hline \n');
% fprintf(TCSA, '  & \\multicolumn{4}{c}{Chocs} \\\\ \n\\cline{2-5}');
% fprintf(TCSA, '  & Productivit\\''e & Natalit\\''e & Sant\\''e & Immigration\\\\ \\hline \n ');
% fprintf(TCSA, '   & $A$ & $n$ & $m_{50}$ & $\\chi_{20}$\\\\ \\hline \n ');
% fprintf(TCSA, '  Production & %.4f & %.4f & %.4f & %.4f\\\\ \n', [solcs1.solcs(1) solcs2.solcs(1) solcs3.solcs(1) solcs4.solcs(1)]);
% fprintf(TCSA, '  Consommation & %.4f & %.4f & %.4f & %.4f\\\\ \n', [solcs1.solcs(2) solcs2.solcs(2) solcs3.solcs(2) solcs4.solcs(2)]);
% fprintf(TCSA, '  Investissement & %.4f & %.4f & %.4f & %.4f\\\\ \n', [solcs1.solcs(3) solcs2.solcs(3) solcs3.solcs(3) solcs4.solcs(3)]);
% fprintf(TCSA, '  D\\''epenses publiques & %.4f & %.4f & %.4f & %.4f\\\\ \n', [solcs1.solcs(4) solcs2.solcs(4) solcs3.solcs(4) solcs4.solcs(4)]);
% fprintf(TCSA, '  Taux d''int\\''er\\^et & %.4f & %.4f & %.4f & %.4f\\\\ \n', [solcs1.solcs(5) solcs2.solcs(5) solcs3.solcs(5) solcs4.solcs(5)]);
% fprintf(TCSA, '  PM du travail & %.4f & %.4f & %.4f & %.4f\\\\ \n', [solcs1.solcs(6) solcs2.solcs(6) solcs3.solcs(6) solcs4.solcs(6)]);
% fprintf(TCSA, '  Capital & %.4f & %.4f & %.4f & %.4f\\\\ \n', [solcs1.solcs(7) solcs2.solcs(7) solcs3.solcs(7) solcs4.solcs(7)]);
% fprintf(TCSA, '  Emploi & %.4f & %.4f & %.4f & %.4f\\\\ \n', [solcs1.solcs(8) solcs2.solcs(8) solcs3.solcs(8) solcs4.solcs(8)]);
% fprintf(TCSA, '  $\\tau_w$ & %.4f & %.4f & %.4f & %.4f\\\\ \n', [solcs1.solcs(9) solcs2.solcs(9) solcs3.solcs(9) solcs4.solcs(9)]);
% fprintf(TCSA, '  Ratio D\\''eficit/PIB & %.4f & %.4f & %.4f & %.4f\\\\ \n', [solcs1.solcs(10) solcs2.solcs(10) solcs3.solcs(10) solcs4.solcs(10)]);
% fprintf(TCSA, '  Ratio Dette/PIB & %.4f & %.4f & %.4f & %.4f\\\\ \n', [solcs1.solcs(11) solcs2.solcs(11) solcs3.solcs(11) solcs4.solcs(11)]);
% fprintf(TCSA, '\\hline \\hline \n \\end{tabular} \n');
% fclose(TCSA);
% 
% 
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% hFig = figure(1);
% % set(hFig, 'Position', [0 20 750 750])
% subplot(4,3,1)
% plot(tt,v1,'LineWidth',2,'Color',[0 0 1])
% axis([0 t -inf inf])
% % hline(0,'-k')
% title(name1)
% box on
% subplot(4,3,2)
% plot(tt,v2,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name2)
% axis([0 t -inf inf])
% box on
% subplot(4,3,3)
% plot(tt,v3,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name3)
% axis([0 t -inf inf])
% box on
% subplot(4,3,4)
% plot(tt,v4,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name4)
% axis([0 t -inf inf])
% box on
% subplot(4,3,5)
% plot(tt,v5,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name5)
% axis([0 t -inf inf])
% box on
% subplot(4,3,6)
% plot(tt,v6,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name6)
% axis([0 t -inf inf])
% box on
% subplot(4,3,7)
% plot(tt,v7,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name7)
% axis([0 t -inf inf])
% box on
% subplot(4,3,8)
% plot(tt,v8,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name8)
% axis([0 t -inf inf])
% box on
% subplot(4,3,9)
% plot(tt,v9,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name9)
% axis([0 t -inf inf])
% box on
% subplot(4,3,10)
% plot(tt,v10,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name10)
% axis([0 t -inf inf])
% box on
% subplot(4,3,11)
% plot(tt,v11,'LineWidth',2,'Color',[0 0 1])
% % hline(0,'-k')
% title(name11)
% axis([0 t -inf inf])
% box on
% subplot(4,3,12)
% plot(tt,v12,'LineWidth',2,'Color',[0 0 1])
% %hline(0,'-k')
% title(name12)
% axis([0 t -inf inf])
% box on




%
% figure;
% subplot(2,2,1)
% plot([oo_.endo_simul(215:274,1)],'LineWidth',2,'Color',[0 0 1])
% title('Epargne')
% subplot(2,2,2)
% plot([oo_.endo_simul(275:334,1)],'LineWidth',2,'Color',[0 0 1])
% title('Consommation')
% subplot(2,2,3)
% plot([oo_.endo_simul(135:174,1)],'LineWidth',2,'Color',[0 0 1])
% title('Offre de travail')
% subplot(2,2,4)
% plot([oo_.endo_simul(175:214,1)],'LineWidth',2,'Color',[0 0 1])
% title('Salaires')



%
%
% figure;
% subplot(2,2,1)
% plot([oo_.endo_simul(215:274,:)])
% title('Epargne')
% subplot(2,2,2)
% plot([oo_.endo_simul(275:334,:)])
% title('Consommation')
% subplot(2,2,3)
% plot([oo_.endo_simul(135:174,:)])
% title('Offre de travail')
% subplot(2,2,4)
% plot([oo_.endo_simul(175:214,:)])
% title('Salaires')
%






%
%
% Ha=theta_a*na;
% Hr=theta_r*nr;
% Hm=theta_m*nm;
% H=Ha+Hr+Hm;
% prod=Y./H;
% warp=wa./wr;
% wrmp=wr./wm;
% relar=Ha./Hr;
% relrm=Hr./Hm;
% khr=K./Hr;
% pk = 1./exp(zk);
% Ca=theta_a*ca;
% Cr=theta_r*cr;
% Cm=theta_m*cm;
% C=Ca+Cr+Cm;

% Y=(Y/Y(1,1)-ones(1,t))*100;
% wa=(wa/wa(1,1)-ones(1,t))*100;
% wr=(wr/wr(1,1)-ones(1,t))*100;
% wm=(wm/wm(1,1)-ones(1,t))*100;
% rk=(rk/rk(1,1)-ones(1,t))*100;
% K=(K/K(1,1)-ones(1,t))*100;
% In=(In/In(1,1)-ones(1,t))*100;
% H=(H/H(1,1)-ones(1,t))*100;
% Ha=(Ha/Ha(1,1)-ones(1,t))*100;
% Hr=(Hr/Hr(1,1)-ones(1,t))*100;
% Hm=(Hm/Hm(1,1)-ones(1,t))*100;
% warp=(warp/warp(1,1)-ones(1,t))*100;
% wrmp=(wrmp/wrmp(1,1)-ones(1,t))*100;
% relar=(relar/relar(1,1)-ones(1,t))*100;
% relrm=(relrm/relrm(1,1)-ones(1,t))*100;
% khr=(khr/khr(1,1)-ones(1,t))*100;
% C=(C/C(1,1)-ones(1,t))*100;
% Ca=(Ca/Ca(1,1)-ones(1,t))*100;
% Cr=(Cr/Cr(1,1)-ones(1,t))*100;
% Cm=(Cm/Cm(1,1)-ones(1,t))*100;
% prod=(prod/prod(1,1)-ones(1,t))*100;
% pk=(pk/pk(1,1)-ones(1,t))*100;

% % Graphs
%
% hFig = figure(1);
% set(hFig, 'Position', [0 20 750 750])
% subplot(4,3,1)
% plot(tt,Y,'LineWidth',2,'Color',[0 0 1])
% axis([0 t -inf inf])
% hline(0,'-k')
% title('Output')
% box on
% subplot(4,3,2)
% plot(tt,In,'LineWidth',2,'Color',[0 0 1])
% hline(0,'-k')
% title('Investment')
% axis([0 t -inf inf])
% box on
% subplot(4,3,3)
% plot(tt,C,'LineWidth',2,'Color',[0 0 1])
% hline(0,'-k')
% title('Consumption')
% axis([0 t -inf inf])
% box on
% subplot(4,3,4)
% plot(tt,H,'LineWidth',2,'Color',[0 0 1])
% hline(0,'-k')
% title('Total hours')
% axis([0 t -inf inf])
% box on
% subplot(4,3,5)
% plot(tt,K,'LineWidth',2,'Color',[0 0 1])
% hline(0,'-k')
% title('Capital')
% axis([0 t -inf inf])
% box on
% subplot(4,3,6)
% plot(tt,rk,'LineWidth',2,'Color',[0 0 1])
% hline(0,'-k')
% title('Rental rate')
% axis([0 t -inf inf])
% box on
% subplot(4,3,7)
% plot(tt,prod,'LineWidth',2,'Color',[0 0 1])
% hline(0,'-k')
% title('Productivity')
% axis([0 t -inf inf])
% box on
% subplot(4,3,8)
% plot(tt,warp,'LineWidth',2,'Color',[0 0 1])
% hold on
% plot(tt,wrmp,'--','LineWidth',2,'Color',[0 0 1])
% hold off
% hline(0,'-k')
% title('Task premiums')
% legend({'W_a/W_r','W_r/W_m'},'Location','best','fontsize',6);
% axis([0 t -inf inf])
% box on
% subplot(4,3,9)
% plot(tt,relar,'LineWidth',2,'Color',[0 0 1])
% hold on
% plot(tt,relrm,'--','LineWidth',2,'Color',[0 0 1])
% hold off
% hline(0,'-k')
% title('Relative hours')
% legend({'H_a/H_r','H_r/H_m'},'Location','best','fontsize',6);
% axis([0 t -inf inf])
% box on
% subplot(4,3,10)
% plot(tt,Ca,'-x','Color',[1 0 0])
% hold on
% plot(tt,Cr,'-d','Color',[0 0 1])
% hold on
% plot(tt,Cm,'-o','Color',[0 1 0])
% hold off
% hline(0,'-k')
% title('Consumption by task')
% legend({'A','R','M'},'Location','best','fontsize',6)
% axis([0 t -inf inf])
% box on
% subplot(4,3,11)
% plot(tt,wa,'-x','Color',[1 0 0])
% hold on
% plot(tt,wr,'-d','Color',[0 0 1])
% hold on
% plot(tt,wm,'-o','Color',[0 1 0])
% hold off
% hline(0,'-k')
% title('Wage by task')
% axis([0 t -inf inf])
% box on
% subplot(4,3,12)
% plot(tt,Ha,'-x','Color',[1 0 0])
% hold on
% plot(tt,Hr,'-d','Color',[0 0 1])
% hold on
% plot(tt,Hm,'-o','Color',[0 1 0])
% hold off
% hline(0,'-k')
% title('Hours by task')
% axis([0 t -inf inf])
% box on

% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% if sh==1
% print(gcf, '-dpdf', [fpath 'trans_rbc_za.pdf']);
% elseif sh==2
% print(gcf, '-dpdf', [fpath 'trans_rbc_zk.pdf']);
% elseif sh==3
% print(gcf, '-dpdf', [fpath 'trans_rbc_ba.pdf']);
% elseif sh==4
% print(gcf, '-dpdf', [fpath 'trans_rbc_br.pdf']);
% elseif sh==5
% print(gcf, '-dpdf', [fpath 'trans_rbc_bm.pdf']);
% end
