clearvars
close all

% load object space output from 'script_psych_distribution.m'
load matlab

% setup
varnames = {'pop' 'hospPsy' 'mainline'};
hospPsyIndex = 3;  % add 1 for the intercept
repetitions = N;
N_dummies = 3000;
wiggle_amount = 0.01;

if mfilename
  diary(sprintf('%s_diary.txt',mfilename));
end

%% initalize data
cases = load_dataset( 'cases', grid_width );
load cases_mn
pop = load_dataset('pop',G,'transform',@log,'bw',pop_bw);
gridPop = pop.denslogtot;

psy = load_dataset('hosp',G,'offset',1,'transform',@log,'bw',hospital_bw);
gridPsy = psy.denslogpsybeds;

[xx,yy]=ndgrid(pop.gx,pop.gy);

% space to hold results
results = cell(repetitions, 1);
v = [];

% loop counter
rep = 0;

%% loop
while rep <= repetitions
  % wiggle cases for numerical stability
  wiggle_ed = cases_mn.net.rmove( cases_mn.expanded_edge_network, wiggle_amount );
  
  % generate dummy points
  edgelist_dummies = rr.rpois(N_dummies);
  
  % fix coordinate ordering
  Y = [ones(cases.N,1); zeros(N_dummies,1)];
  
  % build marked network with data & dummies and update relevant vectors
  mn = MarkedNetwork( rr, [wiggle_ed; edgelist_dummies] );
  X = [mn.x mn.y];
  Y = Y(mn.edge_coordinates_ordering,:);

  % calculate weights
  W = mn.weights.weights;
  
  % validate weights
  if min(W)>0 && max(W)~= Inf
  else
    warning('Retrying because of invalid tesselation...');
    continue  % retry
  end
  
  % set outcome variable for Poisson regression
  YW = Y ./ W;

  % set psych parameters
  if rep == 0
    gridPsy = hosp.denspsybeds;
  else
    gridPsy = P{rep};
  end
  
  % calculate parameters
  interpolation_method = 'linear';  % no negative values allowed
  interPop = interpn( xx, yy, gridPop, X(:,1), X(:,2), interpolation_method );
  interPsy = interpn( xx, yy, gridPsy, X(:,1), X(:,2), interpolation_method );

  % calculate more parameters
  [~,Xdisp,Xedge] = rr.project( X );
  mainlineX = ( rr.edgedata.FCO(Xedge(:,1))==2 );
  
  % build regression input table
  T = array2table( [ interPop interPsy mainlineX YW ], ...
    'VariableNames', [varnames {'resp'}] );
  
  % Poisson regression
  model = fitglm( T, 'ResponseVar', 'resp', ...
    'PredictorVars', varnames, ...
    'Distribution', 'poisson', ...
    'DispersionFlag', true, 'Weights', W, ...
    'Intercept', true, 'Offset', 0, 'Link', 'log' )
 
  if rep == 0
    result0 = model;
    v0 = model.Coefficients.Estimate(hospPsyIndex);
  else
    results{rep} = model;
    vc = model.Coefficients.Estimate(hospPsyIndex);
    v(rep) = vc;
  end
  fprintf('%3d: psy = %.2e\n', rep, model.Coefficients.Estimate(hospPsyIndex) );
  if doPlots && numel(v) > 1
    figure(1);
    boxplot(v)
    hold on
    scatter( 1, v0, 100, 'fill' )
    scatter( 1, vc, 100, 'fill', 'MarkerFaceColor',[1 0 0] )
    hold off
    dat = [v(:); v0];
    R = range(dat);
    ylim([min(dat)-0.05*R, max(dat)+0.05*R]);
    title({'Reshuffled Location Test Statistic','and original loc. (blue dot)'})
    set(gca,'FontSize',8)
    drawnow
  end
  rep = rep + 1;
end

%% collect results
coeffs = zeros( length(results), sum(results{1}.VariableInfo.InModel) + 1 );
dispersions = zeros( length(results), 1 );
for i = 1:length(results)
  r = results{i};
  coeffs(i,:) = r.Coefficients.Estimate;
  dispersions(i) = r.Dispersion;
end
resultsTable = table( ...
  [mean(coeffs)'; mean(dispersions)], ...
  [std(coeffs)'; std(dispersions)], ...
  'VariableNames', {'Mean'; 'SD'}, 'RowNames', ...
  ['Intercept'; r.VariableNames(r.VariableInfo.InModel); 'Dispersion'])

%% plot results
clf
range = 1:nnz(coeffs(:,hospPsyIndex));
v = coeffs(range,hospPsyIndex);
boxplot( v )
v0 = result0.Coefficients.Estimate(hospPsyIndex);
hold on
scatter( 1, v0, 100, 'fill' )
hold off
title('Boxplot from Surrogate Data (blue dot marks coeff. for non-surrogate)')
saveas(gcf, 'boxplot_surrogate.png')

%% find p-value
vNormal = fitdist( v, 'normal' );
pValue = cdf(vNormal,v0,'upper')

%% Group Comparison Plot
if doPlots
  psy_distribution_analyse( X2psy, v, 'chi2' )
end

%% save results
save('results.mat', '-v7.3');
diary off
