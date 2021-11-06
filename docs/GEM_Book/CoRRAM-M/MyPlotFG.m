function  MyPlotFG(Data,Names,Wo,fig);

    lfh=fopen('MyPlot_LogFile.txt','w');
    fprintf(lfh,strcat(char(datetime('now')),'\n'));
    save('MyPLotFG','Data','Names','Wo','fig');
    fprintf(lfh,'Data saved\n');

    [nobs,tmp]=size(Data);
    [nvar,tmp]=size(Names);
    maxp=max(Wo);
    maxl=zeros(maxp,1);
    for i=1:1:nvar;
        if Wo(i)>0; maxl(Wo(i))=maxl(Wo(i))+1; end;
    end;
    Pvar=struct('L','','I',[]);
    for i=1:1:maxp;
        Pvar(i).L=cell(1,maxl(i)); % create cell array
       % Pvar(i).I=zeros(nobs,1);
    end;
   % for i=1:1:maxp;
   %     Pvar(i).L=cell(1,maxl(i)+1); % always plot the zero line
   %     Pvar(i).L{1}='Trend Line';
   %     Pvar(i).I=zeros(nobs,1);
   % end;

    k=ones(maxp,1);

    j=Wo(nvar); % plot the shock
    %k(j)=k(j)+1;    
    Pvar(Wo(nvar)).L{k(j)}=Names(nvar,:);
    Pvar(Wo(nvar)).I=[Pvar(Wo(nvar)).I,Data(:,nvar)];
    k(j)=k(j)+1;
    
    for i=1:nvar-1;
        if Wo(i)>0;  
            j=Wo(i);
            Pvar(Wo(i)).L{k(j)}=Names(i,:);
            Pvar(Wo(i)).I=[Pvar(Wo(i)).I,Data(:,i)];
            k(j)=k(j)+1;
        end;
    end;
  fprintf(lfh,'Plot information written\n');      
% Plot
% Graphic settings
    myColor=    ['k','b','g','c','m','r', 'k', 'b', 'g', 'c', 'm'];
    myLineStyle={'-','-','-','-','-','-','--','--','--','--','--'};
    figname=strcat('Shock ',num2str(fig))
    figure('Name',figname);
    if maxp==1; % one panel
        hf=axes();
        hl=plot(1:nobs,Pvar(1).I);
        xlabel('Quarter','FontSize',8,'FontWeight','normal');
        ylabel('%Dev. from Trend','FontSize',8,'FontWeight','normal');
        for j=1:1:length(hl); hl(j).LineStyle=myLineStyle{j}; hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
        hf.FontSize=7;
        legend(Pvar(1).L,'Location','northeast');
        grid on;
    end;

    if maxp==2; % two panels    
        for i=1:1:2;
            ha=subplot(2,1,i);
            hl=plot(1:nobs,Pvar(i).I);
            xlabel('Quarter','FontSize',8,'FontWeight','normal');
            ylabel('%Dev. from Trend','FontSize',8,'FontWeight','normal');
            for j=1:1:length(hl); hl(j).LineStyle=myLineStyle{j};hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
            ha.FontSize=7;
            legend(Pvar(i).L,'Location','northeast');
        end;
    end;

    if maxp==3; % three panels
        for i=1:1:3;
            ha=subplot(2,2,i);
            hl=plot(1:nobs,Pvar(i).I);
            xlabel('Quarter','FontSize',8,'FontWeight','normal');
            ylabel('%Dev. from Trend','FontSize',8,'FontWeight','normal');
            for j=1:1:length(hl); hl(j).LineStyle=myLineStyle{j};hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
            ha.FontSize=7;
            legend(Pvar(i).L,'Location','northeast');
        end;
    end;

    if maxp==4; % four panels
        for i=1:1:4;
            ha=subplot(2,2,i);
            hl=plot(1:nobs,Pvar(i).I);
            xlabel('Quarter','FontSize',8,'FontWeight','normal');
            ylabel('%Dev. from Trend','FontSize',8,'FontWeight','normal');
            for j=1:1:length(hl); hl(j).LineStyle=myLineStyle{j};hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
            ha.FontSize=7;
            legend(Pvar(i).L,'Location','northeast');
        end;
    end;

    if maxp==5; % five panels
        for i=1:1:5;
            ha=subplot(3,2,i);
            hl=plot(1:nobs,Pvar(i).I);
            xlabel('Quarter','FontSize',8,'FontWeight','normal');
            ylabel('%Dev. from Trend','FontSize',8,'FontWeight','normal');
            for j=1:1:length(hl);hl(j).LineStyle=myLineStyle{j}; hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
            ha.FontSize=7;
            legend(Pvar(i).L,'Location','northeast');
        end;    
    end;

    if maxp==6; % six panels
        for i=1:1:6;
            ha=subplot(3,2,i);
            hl=plot(1:nobs,Pvar(i).I);
            xlabel('Quarter','FontSize',8,'FontWeight','normal');
            ylabel('%Dev. from Trend','FontSize',8,'FontWeight','normal');
            for j=1:1:length(hl);hl(j).LineStyle=myLineStyle{j}; hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
            ha.FontSize=7;
            legend(Pvar(i).L,'Location','northeast');
            legend('boxoff');
            ax=gca;
            ax.XLim=[1 nobs];
            grid on;
        end;    
    end;
%    filename=[pwd '\Pic\TestModel2_Fig',num2str(fig),'.fig'];
    saveas(gcf,strcat('Figure',num2str(fig),'.fig'));
    fprintf(lfh,'Finished\n');
    fclose(lfh);
    
    return;

end

