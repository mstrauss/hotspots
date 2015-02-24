function plot_simulation_envelopes( X, Y )
hold on
[x,y] = select_curve( X, Y, 'lower' );
plot(x,y, 'Color', [0.8 0.8 0.8], 'LineWidth', 1)
[x,y] = select_curve( X, Y, 'upper' );
plot(x,y,'Color', [0.8 0.8 0.8], 'LineWidth', 1)
[x,y] = select_curve( X, Y, @median );
plot(x,y,'Color', [.5 .5 .5], 'LineWidth', 2)
hold off

% formatting
set(gca,'FontSize',16);
set(gca,'LineWidth',2);

end
