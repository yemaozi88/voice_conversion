%
% 2010-04-26
% draws bar chart with standard deviation
%


%% bar graph
% data
X = 0.5:4.5;
Y = [50.0, 50.0; ... % nani1
    91.7, 100.0; ... % nani2
    58.3, 50; ... % iina1
    37.5, 37.5; ... % iina2
    11.1, 22.2];    % iie
 
hold on
    % bar chart
    h = bar(X, Y, 1.0);
    
    %xlabel('The number of mixtures', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Inteligibility per phoneme [%]', 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');
    tcksX={'nani1', 'nani2', 'iina1', 'iina2', 'iie'};
    set(gca, 'XTickLabel',tcksX ,'XTick',0.5:1:length(tcksX)+0.5);
    set(gca, 'FontName', 'Arial', 'FontSize', 18);
hold off