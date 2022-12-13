function createfigure_bar(ymatrix1, ymatrix2, std1, std2, y_limit, XTEXT,file_name)

% ymatrix1 is for cost and #edges, 4*6 and last 2 columns are nan.
% ymatrix2 is for eps, 4*6 and first 4 columns are nan.
% XTEXT = "epsilon", "number of sites", "sampling rate"
% file_name: save as pdf


figure1 = figure('Name',file_name);
axes1 = axes('Parent',figure1);
axes1.FontSize = 16
hold(axes1,'on');

barname = ["Comm. Cost (LocalSS)", "Comm. Cost (LocalSC)", "# Edges (LocalSS)", "# Edges (LocalSC)", "\epsilon (LocalSS)", "\epsilon (LocalSC)"]
% colors = ['#0000FF', '#00FFFF', '#A2142F', '#0072BD', '#4DBEEE', '#EDB120']
% colors = [[0 0 1]; [0 1 1]; [0.6350 0.0780 0.1840];  [0.9290 0.6940 0.1250] ; [0 0.4470 0.7410]; [0.3010 0.7450 0.9330];]
colors = [[0 0.4470 0.7410];
          [0.4660 0.6740 0.1880];
          [0.6350 0.0780 0.1840];
          [0.9290 0.6940 0.1250];
          
          [0.4940 0.1840 0.5560];
          [0.3010 0.7450 0.9330]
]


% Error bar: https://www.mathworks.com/matlabcentral/answers/438514-adding-error-bars-to-a-grouped-bar-plot
ngroups = size(ymatrix1, 1);
nbars = size(ymatrix1, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));


% 2 y axis: https://www.mathworks.com/help/matlab/creating_plots/plotting-with-two-y-axes.html
% https://www.mathworks.com/matlabcentral/answers/239236-bar-plot-with-2-y-axes-and-same-x-axis
yyaxis left
bar1 = bar(ymatrix1,'Parent',axes1);
ylabel("Comm. Cost");
% ymax = fix(max(max(ymatrix1))/2000) + 1;
% ytick = 0 :2: 2*ymax;
% ytick = ytick*1000
% set(axes1,'YTick',ytick);
% 
% yticklabel = cell(1,size(ytick,2))
% for i =1:size(ytick,2)
%     yticklabel{1,i} = [num2str(ytick(i)/1000) 'K'];
% end
% yticklabel{1,1} = '0';
% % yticklabel{1,size(ytick,2)} = '';
% yticklabels(yticklabel)


ymax = max(max(ymatrix1));
if ymax<14000
    ymax = fix(ymax/2000) + 1;
    ytick = 0 :2: 2*ymax;
elseif ymax<30000
    ymax = fix(ymax/4000) + 2;
    ytick = 0 :4: 4*ymax;
else
    ymax = fix(ymax/20000) + 1;
    ytick = 0 :20: 20*ymax;
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


for j=1:size(barname,2)
    i=size(barname,2)-j+1
    if (i ~= 3) & (i~=4)
        set(bar1(i),'DisplayName',barname(i));
%         hp = findobj(bar1(i),'type','patch'); 
%         hatch(hp(1),45,'k','-',4,2);
    else
        set(bar1(i),'DisplayName','off');
    end
    bar1(i).FaceColor = colors(i,:);
    
    
end
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    if i~=3 & i~=4
        errorbar(x, ymatrix1(:,i), std1(:,i), '.');
    end
end



yyaxis right
bar2 = bar(ymatrix2,'Parent',axes1);
ylabel("Approx. Quality");
ylim(y_limit)

for j=1:size(barname,2)
    i=size(barname,2)-j+1
    if i == 3 || i==4
        set(bar2(i),'DisplayName',barname(i));
    end
    bar2(i).FaceColor = colors(i,:);
end  

for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    if i==3 || i==4
        errorbar(x, ymatrix2(:,i), std2(:,i), '.');
    end
end


ymax = max(max(ymatrix2));
if ymax> 0.25
    ymax = fix(ymax/0.1) + 2;
    ytick = 0 :0.1: 0.1*ymax;
elseif ymax> 0.15 & ymax <= 0.25
    ymax = fix(ymax/0.1) + 1;
    ytick = 0 :0.1: 0.1*ymax;
    yticklabels({'0' '0.1', '0.2', '0.3'})
    % here need manually set up
else
    ymax = 0.1
    ytick = 0 :0.02: 0.1;
end
% ytick = ytick*1000
set(axes1,'YTick',ytick);

% yticklabel = cell(1,size(ytick,2))
% for i =1:size(ytick,2)
%     yticklabel{1,i} = [num2str(ytick(i)/1000) 'K'];
% end
% yticklabel{1,1} = '0';
% yticklabel{1,size(ytick,2)} = '';
% yticklabels(yticklabel)




%xlabel('Time step');
xlabel("Approx. Parameter \epsilon");
box(axes1,'on');
x_ticks = 1:size(ymatrix1)+3
% set(axes1,'XTick',[0 1 2 3 4 5]);
set(axes1,'XTick',x_ticks);
axes1.XColor='black';
axes1.XAxis.Color = 'black';
% xlim([0 5]);

xticklabel = cell(1,size(x_ticks,2))
for i =1:size(x_ticks,2)
    xticklabel{1,i} = [num2str((x_ticks(i)+1)/10)];
end
% xticklabel{1,1} = '';
xticklabel{1,size(x_ticks,2)} = '';

if XTEXT == "epsilon"
    xticklabels(xticklabel)
end


    
% legend1 = legend(axes1,'show');
legend1 = legend([bar1(1:4), bar2(5:6)],['Comm. Cost (' '{\itLocalSS}' ')'], ['Comm. Cost (' '{\itLocalSC}' ')'], ['Approx. Quality (' '{\itLocalSS}' ')'], ['Approx. Quality (' '{\itLocalSC}' ')'], ['SC Size (' '{\itLocalSS}' ')'], ['SC Size (' '{\itLocalSC}' ')'])

set(legend1,'Orientation','horizontal','NumColumns',2,...
    'Location','northoutside',...
    'FontSize',15);

saveas(figure1,file_name,"pdf")
hold off
