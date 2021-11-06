function  MyPlot(xt,Model,iz)
% Plot impulse responses

%{
    Copyright : Alfred Maußner
    First version:   02 February 2017 (from MyPlot_T)
    Revision:        15 February 2018, bug fixed in using 7 and 8 panels
    Revision:        21 February 2018, bug fixed (counter was not increased for plots of shocks)
    
    Purpose : Plot the impulse responses in xt according to the information given in Var.
              Different from previous versions, this version removes any non consequtive
              panels.

    Input   : . xt nobs by nx+ny+1 matrix which stores the impulse respones of the state and
                the control variables to the impulse, whose time profile is given in the
                rightmost column of xt.
              . Model, an instance of the DSGE class
              . iz, scalar, the index of the shock
             
%}

%Fontname=Defaults.Fontname;
Fontname='Source Sans Pro';

Var=Model.Var;
nx=Model.nx;
nz=Model.nz;
ny=Model.ny;
nvar=nx+nz+ny;
nobs=Model.inobs;

% Graphic settings
if Model.Trendline;
    myColor=    ['r','k','b','g','c','m','k','b','g','c','m'];
    myLineStyle={'--','-','-','-','-','-','--','--','--','--','--'};
else
    myColor=    ['k','b','g','c','r','m','k','b','g','c','m'];
    myLineStyle={'-','-','-','-','-','--','--','--','--','--','--'};
end;

% find the number of panels to plot to and the number of lines to be plotted in each panel.
maxl=GetMaxPlotNumber(Var,iz);
maxp=size(nonzeros(maxl),1);

k=ones(8,1); % counter for lines per panel

% Gather information about each panel in a structure
Pvar=struct('L','','I',[]);
for i=1:1:8;
    if maxl(i)>0;
        if Model.Trendline;
            Pvar(i).L=cell(1,maxl(i)+1); 
            Pvar(i).L{1}='Trend Line';
            Pvar(i).I=zeros(nobs,1);
            k(i)=k(i)+1;
        else
            Pvar(i).L=cell(1,maxl(i)); 
        end;
    end;
end;

% write information in Pvar

for i=1:1:nvar;
    if Var(i).Plotno>0;
        if ~strcmp(Var(i).Type,'z');  
            j=Var(i).Plotno;
            Pvar(Var(i).Plotno).L{k(j)}=Var(i).Name;
            if Var(i).Type=='x';
                Pvar(Var(i).Plotno).I=[Pvar(Var(i).Plotno).I,xt(:,Var(i).Pos)];
            else
                Pvar(Var(i).Plotno).I=[Pvar(Var(i).Plotno).I,xt(:,Var(i).Pos+nx)]; 
            end;          
            k(j)=k(j)+1;
        else
            if Var(i).Pos==iz;
                j=Var(nx+iz).Plotno;
                Pvar(Var(nx+iz).Plotno).L{k(j)}=Var(nx+iz).Name;
                Pvar(Var(nx+iz).Plotno).I=[Pvar(Var(nx+iz).Plotno).I,xt(:,nx+ny+1)];
                k(j)=k(j)+1;
            end;
        end;
    end;
end;

% now plot
figname=strcat('Shock ',num2str(iz));
figure('Name',figname);

ip=find(maxl); % the panels to plot to

if maxp==1; % one panel
   hl=plot(1:nobs,Pvar(ip.I));
   ax=gca;
   ax.FontName=Fontname;
   ax.FontSize=10;
   ax.LabelFontSizeMultiplier=1.2;
   xlabel('Quarter');
   ylabel('%Dev. from Trend');
   for j=1:1:length(hl); hl(j).LineStyle=myLineStyle{j}; hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
   legend(Pvar(ip).L,'Location','northeast','FontName','Source Sans Pro','FontSize',10);      
   if Model.LegendBox; legend('boxon'); else legend('boxoff'); end;
   if Model.Grid; grid on; else grid off; end;
end;

if maxp==2; % two panels   
    for i=1:1:2;
        ha=subplot(2,1,i);
        hl=plot(1:nobs,Pvar(ip(i)).I);
        ax=gca;
        ax.FontName='Source Sans Pro';
        ax.FontSize=9;
        ax.LabelFontSizeMultiplier=1.2;
        xlabel('Quarter');
        ylabel('%Dev. from Trend');        
        for j=1:1:length(hl); hl(j).LineStyle=myLineStyle{j};hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
        legend(Pvar(ip(i)).L,'Location','northeast','FontName','Source Sans Pro','FontSize',9);   
        if Model.LegendBox; legend('boxon'); else legend('boxoff'); end;
        if Model.Grid; grid on; else grid off; end;
    end;
end;

if maxp==3; % three panels
    for i=1:1:3;
        subplot(2,2,i);
        hl=plot(1:nobs,Pvar(ip(i)).I);
        ax=gca;
        ax.FontName='Source Sans Pro';
        ax.FontSize=9;
        ax.LabelFontSizeMultiplier=1.2;
        xlabel('Quarter');
        ylabel('%Dev. from Trend');  
        for j=1:1:length(hl); hl(j).LineStyle=myLineStyle{j};hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
        legend(Pvar(ip(i)).L,'Location','northeast','FontName','Source Sans Pro','FontSize',9);   
        if Model.LegendBox; legend('boxon'); else legend('boxoff'); end;
        if Model.Grid; grid on; else grid off; end;
    end;
end;

if maxp==4; % four panels
    for i=1:1:4;
        ha=subplot(2,2,i);
        hl=plot(1:nobs,Pvar(ip(i)).I);
        ax=gca;
        ax.FontName='Source Sans Pro';
        ax.FontSize=10;
        ax.LabelFontSizeMultiplier=1.2;
        xlabel('Quarter');
        ylabel('%Dev. from Trend');  
        for j=1:1:length(hl); hl(j).LineStyle=myLineStyle{j};hl(j).LineWidth=1; hl(j).Color=myColor(j); end;        
        legend(Pvar(ip(i)).L,'Location','northeast','FontName','Source Sans Pro','FontSize',10);   
        if Model.LegendBox; legend('boxon'); else legend('boxoff'); end;
        if Model.Grid; grid on; else grid off; end;
    end;
end;

if maxp==5; % five panels
    for i=1:1:5;
        ha=subplot(3,2,i);
        hl=plot(1:nobs,Pvar(i).I);
        ax=gca;
        ax.FontName='Source Sans Pro';
        ax.FontSize=9;
        ax.LabelFontSizeMultiplier=1.2;
        xlabel('Quarter');
        ylabel('%Dev. from Trend');  
        for j=1:1:length(hl);hl(j).LineStyle=myLineStyle{j}; hl(j).LineWidth=1; hl(j).Color=myColor(j); end;        
        legend(Pvar(i).L,'Location','northeast','FontName','Source Sans Pro','FontSize',10);   
        if Model.LegendBox; legend('boxon'); else legend('boxoff'); end;
        if Model.Grid; grid on; else grid off; end;
    end;    
end;

if maxp==6; % six panels
    for i=1:1:6;
        ha=subplot(3,2,i);
        hl=plot(1:nobs,Pvar(i).I);        
        ax=gca;
       % ax.FontName='Source Sans Pro';
        ax.FontName=Fontname;
        %ax.FontSize=9;
        ax.FontSize=8;
        %ax.LabelFontSizeMultiplier=1.2;
        ax.LabelFontSizeMultiplier=1.1;
        xlabel('Quarter');
        ylabel('%Dev. from Trend');
        for j=1:1:length(hl);hl(j).LineStyle=myLineStyle{j}; hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
        %legend(Pvar(i).L,'Location','northeast','FontName','Helvetica','FontSize',10);   
        legend(Pvar(i).L,'Location','northeast','FontName',Fontname,'FontSize',9);   
        if Model.LegendBox; legend('boxon'); else legend('boxoff'); end;
        if Model.Grid; grid on; else grid off; end;
    end;    
end;
if maxp==7; % seven panels
    for i=1:1:7;
        ha=subplot(4,2,i);
        hl=plot(1:nobs,Pvar(i).I);        
        ax=gca;
        ax.FontName=Fontname;
        ax.FontSize=7;
        ax.LabelFontSizeMultiplier=1;
        xlabel('Quarter');
        ylabel('%Dev. from Trend');
        for j=1:1:length(hl);hl(j).LineStyle=myLineStyle{j}; hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
        legend(Pvar(i).L,'Location','northeast','FontName',Fontname,'FontSize',8);   
        if Model.LegendBox; legend('boxon'); else legend('boxoff'); end;
        if Model.Grid; grid on; else grid off; end;
    end;
end;
    
if maxp==8; % eight panels
    for i=1:1:8;
        ha=subplot(4,2,i);
        hl=plot(1:nobs,Pvar(i).I);        
        ax=gca;
        ax.FontName=Fontname;
        ax.FontSize=7;
        ax.LabelFontSizeMultiplier=1;
        xlabel('Quarter');
        ylabel('%Dev. from Trend');
        for j=1:1:length(hl);hl(j).LineStyle=myLineStyle{j}; hl(j).LineWidth=1; hl(j).Color=myColor(j); end;
        legend(Pvar(i).L,'Location','northeast','FontName',Fontname,'FontSize',8);   
        if Model.LegendBox; legend('boxon'); else legend('boxoff'); end;
        if Model.Grid; grid on; else grid off; end;
    end;    
end;

end

