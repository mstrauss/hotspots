function psy_distribution_analyse( sortable_data, groupable_data, data_name )
[data_sorted, data_order] = sort(sortable_data);
v_sorted = groupable_data(data_order);
numGroups = 12;
grpSize = floor(numel(groupable_data)/numGroups);
v_group = repmat(1:numGroups, grpSize, 1 );
Vsorted = reshape( v_sorted(1:numel(v_group)), grpSize, numGroups );

% Kruskal-Wallis test
pKV = kruskalwallis(Vsorted, [], 'off')
fprintf('Kruskal-Wallis test (H0: all groups from same distr.): p=%.2g\n', pKV )

figure
set(gcf,'units','normalized','outerposition',[0.2 0.2 0.8 .5])
set(gcf, 'PaperPositionMode', 'auto')
subplot(131)
hist(data_sorted, numGroups)
title(sprintf('Histogramm of %s distribution', data_name))
xlabel(data_name);  ylabel('Number of Simulations')
subplot(132)
boxplot( reshape( data_sorted(1:numel(v_group)), grpSize, numGroups ) )
ylabel(data_name)
subplot(133)
boxplot( Vsorted )
ylabel('hospPsy regr. coeff.')
xlabel(sprintf('Kruskal-Wallis: p=%.2g', pKV))

%%
saveas(gcf, sprintf('boxplots hospPsy(%s).png', data_name))
end
