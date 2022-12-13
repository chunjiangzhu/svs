function createfigure_line(ymatrix1, std, xvalue, yflag, XTEXT,file_name)

% ymatrix1 is for cost and #edges, 4*6 ( 4 different options, 6 lines)
% XTEXT = "epsilon", "number of sites", "sampling rate"
% file_name: save as pdf


figure1 = figure('Name',file_name);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
axes1.FontSize = 16

% barname = ["Comm. Cost (LocalSS)", "Comm. Cost (LocalSC)", "# Edges (LocalSS)", "# Edges (LocalSC)", "\epsilon (LocalSS)", "\epsilon (LocalSC)"]
% colors = ['#0000FF', '#00FFFF', '#A2142F', '#0072BD', '#4DBEEE', '#EDB120']
% colors = [[0 0 1]; [0 1 1]; [0.6350 0.0780 0.1840];  [0.9290 0.6940 0.1250] ; [0 0.4470 0.7410]; [0.3010 0.7450 0.9330];]
colors = [[0 0.4470 0.7410]; [0 0.4470 0.7410];
%           [0.4660 0.6740 0.1880]; [0.4660 0.6740 0.1880];
          [0.6350 0.0780 0.1840];[0.6350 0.0780 0.1840];
          [0.4660 0.6740 0.1880]; [0.4660 0.6740 0.1880];
%           [0.9290 0.6940 0.1250];
%           [0.4940 0.1840 0.5560];[0.4940 0.1840 0.5560];
%           [0.3010 0.7450 0.9330]; [0.3010 0.7450 0.9330]
]

styles = [ "-.","-", "-.","-", "-.", "-"
    
         ]
markers = ["o", "o", "+", "+", "*", "*"]

% Error bar: https://www.mathworks.com/matlabcentral/answers/438514-adding-error-bars-to-a-grouped-bar-plot
% ngroups = size(ymatrix1, 1);
% nbars = size(ymatrix1, 2);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));


% 2 y axis: https://www.mathworks.com/help/matlab/creating_plots/plotting-with-two-y-axes.html
% https://www.mathworks.com/matlabcentral/answers/239236-bar-plot-with-2-y-axes-and-same-x-axis

yyaxis left
% bar1 = bar(ymatrix1,'Parent',axes1);
ylabel("Comm. Cost");
x = xvalue(2:5)
ymax = max(max(ymatrix1));
if ymax<16000
    ymax = fix(ymax/2000) + 1;
    ytick = 0 :2: 2*ymax;
elseif ymax < 30000
    ymax = fix(ymax/4000) + 1;
    ytick = 0 :4: 4*ymax;
else
    ymax = fix(ymax/10000) + 1;
    ytick = 0 :10: 10*ymax;
end
ytick = ytick*1000
set(axes1,'YTick',ytick);

yticklabel = cell(1,size(ytick,2))
for i =1:size(ytick,2)
    yticklabel{1,i} = [num2str(ytick(i)/1000) 'K'];
end
yticklabel{1,1} = '0';
% yticklabel{1,size(ytick,2)} = '';
yticklabels(yticklabel)



plots = cell(6,1);
for i =[1,2,5,6]
    data = squeeze(ymatrix1(:,i)');
    error = squeeze(std(:,i)');
    plots{i} = line(x, data, 'Color', [colors(i,:)], 'LineStyle', styles(i), 'LineWidth', 1.5, 'Marker', markers(i))
    errorbar(x, data, error, 'Color', [colors(i,:)], 'LineStyle', styles(i), 'LineWidth', 1.5, 'Marker', markers(i))
end

yyaxis right
ylabel("Approx. Quality");
% if yflag == 1
%     ylim([0 0.05])
%     ytick = 0:0.01:0.05;
%     set(axes1,'YTick',ytick);
% else
%     ylim([0 0.06])
%     ytick = 0:0.02:0.06;
%     set(axes1,'YTick',ytick);
% end
ylim([0 0.1])
    
for i =3:4
    data = squeeze(ymatrix1(:,i)');
    error = squeeze(std(:,i)');
    plots{i} = line(x,data, 'Color', [colors(i,:)], 'LineStyle', styles(i), 'LineWidth', 1.5, 'Marker', markers(i))
    errorbar(x, data, error, 'Color', [colors(i,:)], 'LineStyle', styles(i), 'LineWidth', 1.5, 'Marker', markers(i))
end



%xlabel('Time step');
% xlabel(XTEXT);
box(axes1,'on');



if XTEXT == "Number of Sites"
    set(axes1,'XTick',xvalue);
    xlabel([XTEXT '{\it s}'])
    xlim([xvalue(1) xvalue(6)]);
    xticklabels({'', '3','5','10','15', ''})
elseif XTEXT == "Sampling Rate"
    set(axes1,'XTick',xvalue);
    xlabel([XTEXT '{\it r}'])
    xlim([xvalue(1) xvalue(6)]);
    xticklabels({'', '0.01','0.05','0.10','0.15', ''})
elseif XTEXT == "Number of Vertices"
    set(axes1,'XTick',xvalue);
    xlabel([XTEXT '{\it n}'])
    xlim([xvalue(1) xvalue(6)]);
    xticklabels({'', '0.5K','1K','2K','3K', ''})
end
    
% legend1 = legend(axes1,'show', "Comm. Cost (LocalSS)", "Comm. Cost (LocalSC)", "Approx. Quality (LocalSS)", "Approx. Quality (LocalSC)", "SC Size (LocalSS)", "SC Size (LocalSC)");
% legend1 = legend([bar1(1:4), bar2(5:6)], "Comm. Cost (LocalSS)", "Comm. Cost (LocalSC)", "Approx. Quality (LocalSS)", "Approx. Quality (LocalSC)", "# Edges (LocalSS)", "# Edges (LocalSC)")
legend1 = legend([plots{1}, plots{2},plots{3},plots{4},plots{5},plots{6}], ['Comm. Cost (' '{\itLocalSS}' ')'], ['Comm. Cost (' '{\itLocalSC}' ')'], ['Approx. Quality (' '{\itLocalSS}' ')'], ['Approx. Quality (' '{\itLocalSC}' ')'], ['SC Size (' '{\itLocalSS}' ')'], ['SC Size (' '{\itLocalSC}' ')'])

set(legend1,'Orientation','horizontal','NumColumns',2,...
    'Location','northoutside',...
    'FontSize',15);
saveas(figure1,file_name,"pdf")
hold off
