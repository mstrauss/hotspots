%% setup
clear all
close all
project_configuration
doPlots = true;
grid_width = 0.25;
kernel_function = 'gauss_convolute';
pop_bw = 1.6;
hospital_bw = pop_bw/2;  % hospital bandwidth

N = 999;  % number of shuffled configurations to create
acceptance_threshold = 5e3;

if mfilename
  diary(sprintf('%s_diary.txt',mfilename));
%   dos(sprintf('cp "%s.m" "%s"', mfilename, outdir));
end

%% initialize
assert( hospital_bw < pop_bw, 'hospital bw must be smaller than pop bw so that we have no strange division effects' );
load cases_mn
rr = cases_mn.net;
[~, G] = load_dataset( 'pop', grid_width );
pop = load_dataset( 'pop', G, 'bw', pop_bw, 'kernel_function', kernel_function );
L = G.line_density(rr.edges);

hosp = load_dataset( 'hosp', G, 'bw', hospital_bw, 'kernel_function', kernel_function );
hosp.totalbeds = hosp.psybeds + hosp.nonpsybeds;

%% initialize plot windows
if doPlots
  % smaller window
  figure(1)
  set(gcf,'units','normalized','outerposition',[0.2 0.2 0.8 .5])
  set(gcf, 'PaperPositionMode', 'auto');
end
%%
nBins = 15;  % number of bins chosen, such that the b coefficient of the fit flattens out

%%
P0 = hosp.denspsybeds;% ./pop.dens;
NP0 = hosp.densnonpsybeds;% ./pop.dens;

T = cell(N,1);  P = cell(N,1);  NP = cell(N,1);
X2total = zeros(N,1);  % chi2
KLtotal = zeros(N,1);
KLpsy = zeros(N,1);  % divergence
X2psy = zeros(N,1);  % chi2
S = cell(N,1);  % save shuffled configurations
not_acc = 0;

%%
Ppkm0 = P0./L;  Ppkm0(isinf(Ppkm0))=0;  Ppkm0(isnan(Ppkm0))=0;
[Ph0, Pctr] = hist( (P0(:)), nBins );
Ppkmctr = logspace(0, log10(max(Ppkm0(:))), nBins);
[Ppkmh0, Ppkmctr] = hist( (Ppkm0(:)), Ppkmctr );

if doPlots
  G.save_image( L + P0, 'beds_distr_000' );
  figure(1); plot_beds_distr_per_km( 0, Ph0, Pctr, Ppkmh0, Ppkmctr ); %, Th0, Tctr, Tpkmh0, Tpkmctr )
end

%%
n = 1;   % number of current configuration
while n <= N
  shuffled_hosp = G.shuffle_positions( hosp, pop.denstot );
  S{n} = shuffled_hosp;
  Z = G.project( shuffled_hosp.x, shuffled_hosp.y, hosp.psybeds );
  shuffled_hosp.psydens = G.convolution_smoothing( hospital_bw*eye(2), Z );
  
  p = shuffled_hosp.psydens;% ./pop.dens;
  
  % clean data
  p(isnan(p) | isinf(p))=0;
  
  P{n} = p;
  
  Ppkm = P{n}./L;
  Ppkm(L==0)=0;
  Ph = hist( (P{n}(:)), Pctr );
  Ppkmh = hist( (Ppkm(:)), Ppkmctr );
  
  % calculate distributional divergences on network
  KLpsy(n) = sum( kl( Ppkmh, Ppkmh0 ) );
  sel = ~( Ppkmh==0 & Ppkmh0==0 );
  % weighted chi-square
  X2psy(n) = sum( (Ppkmh(sel)-Ppkmh0(sel)).^2 ./ Ppkmh0(sel) .* Ppkmctr(sel) );
  
  % accept
  if X2psy(n) < acceptance_threshold
    % plot
    if doPlots
      G.save_image( L + p, sprintf('beds_distr_%03d.png',n) );
      fprintf( '%d accepted (chi2=%e, ratio=%e)\n', n, X2psy(n), n/(n+not_acc) );

      figure(1); plot_beds_distr_per_km( n, Ph, Pctr, Ppkmh, Ppkmctr, sprintf('chi2 = %.2e', X2psy(n) ));
    end
    n = n+1;
  else
    not_acc = not_acc + 1;
    %warning( '%d not accepted (ratio: %f, chi2=%e)', n, not_acc/(not_acc+n-1), X2psy(n) );
  end
  
end


%%
diary off
save('matlab.mat', '-v7.3');
