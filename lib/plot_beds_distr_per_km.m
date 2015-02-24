function plot_beds_distr_per_km( n, Ph, Pctr, Ppkmh, Ppkmctr, text )

if ~exist('text', 'var')
  text = '';
end

clf

Pmax = max(Pctr);
Ppkmmax = max(Ppkmctr);

subplot(121) % psy beds per km
loglog( Ppkmctr, Ppkmh, 'x' );
xlim([0,Ppkmmax]);
xlabel('(Beds per Population) per (Railroad km per Pixel)');
ylabel('Pixel count')
title(sprintf('Distr. of Psy. Beds along the RR Network, i=%d, %s', n, text) );

subplot(122) % psy beds
semilogy( Pctr, Ph, 'x' );
xlim([0,Pmax]);
xlabel('(Beds per Population) per Pixel');
ylabel('Pixel count')
title( ['Distribution of Psychiatric Beds, i=', num2str(n)] );

saveas(gcf,sprintf('beds_distr_per_km_%03d.png',n));
end
