function [output] = friday13boxplot(csv_file_name,color_index,jitter)

%________________________________________________________________________________________________________________________
% Written by W. Davis Haselden
% M.D./Ph.D. Candidate, Department of Neuroscience
% The Pennsylvania State University & Herhsey College of Medicine
%________________________________________________________________________________________________________________________
%
%   Purpose: create a box plot showing individual data points for number of
%   trauma activations on friday the 13th and its neighboring fridays.
%________________________________________________________________________________________________________________________
%
%   Inputs: Friday 13th Graph - Sheet1.csv file containing the number of
%   trauma activations on friday the 13th and the fridays before and
%   following the week of.
%
%   Outputs: Box plot of trauma activations between the three fridays
%________________________________________________________________________________________________________________________

%% initialize variables
raw_data = importdata('Friday 13th Graph - Sheet1.csv');
color_index = {[242,240,247]./256,[158,154,200]./256,[74,20,134]./256};
jitter = 0.1;

%% Generate figure for illustrator: only plot data and refine later in the illustrator
monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
w = 600; ht = 400; % Define window size: width(w) and height(ht)
mon = 1; % Only one monitor

% Determine positioning for main figure4
pos = [((monitorPos(mon,3)-w)/2)+monitorPos(mon,1)/2,...
    ((monitorPos(mon,4)-ht)*16/20)+monitorPos(mon,2)/2,...
    w ht];

% turn off warning to prevent MATLAB throwing warning
warning('off','MATLAB:legend:PlotEmpty');
warning('off','MATLAB:legend:IgnoringExtraEntries');

figure('name','trauma activations around friday the 13th','NumberTitle','off',...
    'position',pos,'color','w');
hold on

%%

for ii = 1:length(raw_data.textdata)
DateString = raw_data.textdata{ii};
formatIn = 'dd-mmm-yy';
hold_date = datenum(DateString,formatIn);
date_vec_hold =  datevec(hold_date);
data_all(ii,:) = [date_vec_hold(1:3) raw_data.data(ii)];
end

prior_week_index = find(data_all(:,3)==6);
f13_week_index = find(data_all(:,3)==13);
next_week_index = find(data_all(:,3)==20);

prior = data_all(prior_week_index,:);
f13 = data_all(f13_week_index,:);
next = data_all(next_week_index,:);

activation_prior = sort(prior(:,4));
activation_f13 = sort(f13(:,4));
activation_next = sort(next(:,4));

x = [ones(length(activation_prior),1); 2.*ones(length(activation_f13),1); 3.*ones(length(activation_next),1)];
y = [activation_prior; activation_f13; activation_next];
% figure
% notBoxPlot(y,x,'markMedian',true)
% ylim([0 10])

output = boxplot(hAxes,y,x)
set(gca,'XTickLabels',{'prior week','friday the 13th','following week'})
%set(hBoxPlot,'Color','k')
ylabel('trauma activations')
ylim([0 10]);

jitter = gaussmf(linspace(0.75,1.25,length(activation_prior)),[0.1 1]);
s = scatter(ones(length(activation_f13),1),activation_prior,'fill','jitter','on','jitterAmount', jitter)
set(s,'MarkerFaceColor',color_index{1},'MarkerEdgeColor',[0 0 0])

s = scatter(2.*ones(length(activation_f13),1),activation_f13,'fill','jitter','on','jitterAmount', jitter)
set(s,'MarkerFaceColor',color_index{2},'MarkerEdgeColor',[0 0 0])

s = scatter(3.*ones(length(activation_next),1),activation_next,'fill','jitter','on','jitterAmount', jitter)
set(s,'MarkerFaceColor',color_index{3},'MarkerEdgeColor',[0 0 0])

end