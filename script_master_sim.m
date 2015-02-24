%% init
clearvars
close all
project_configuration

%% config

% Save intermediate stuff?
save_intermediates = true;

% Gaussian noise model bandwidth (km)
wiggle_amount = 1;

% number of dummy points
N_dummies = 2e4;

% number of initial fitting runs
N_wiggle_fit = 39;

%%
% (quadratic) pixel width (km)
grid_width = 0.25;

% psychiatric beds smoothing bandwidth (km)
hospital_bw = 0.8;
pop_bw = 1.6;

% number of simulation runs of fitted model
N_sim_model = 199;

%% setup
if mfilename
  diary(sprintf('%s_diary.txt',mfilename));
end
timer = tic;

%% define variables
disp('**************************')
disp('*** Building Variables ***')
disp('**************************')

stdopt = {'transform', 'log', 'extrapolate', true, 'domain'};

%%% gemeinde variables
SocioDemoDomain = Domain.Polygon(fullfile(datadir,'Gemeinden2011_mit_SozioOekonomDaten-simplified.shp'));
gemedu1    = Variable('gemedu1',    stdopt{:}, SocioDemoDomain);
gemedu2    = Variable('gemedu2',    stdopt{:}, SocioDemoDomain);
gemedu3    = Variable('gemedu3',    stdopt{:}, SocioDemoDomain);

%%% district variables
DistrictDataDomain = Domain.Polygon(fullfile(datadir,'BezirkSuizidraten2001-2008-simplified.shp'));
suirate = Variable('suirate', stdopt{:}, DistrictDataDomain );
income = Variable('wkoinc03', stdopt{:}, DistrictDataDomain );

%%% network variables
load cases_mn
NetworkDomain = Domain.NetworkEdge(cases_mn.net);
fco = Variable( 'FCO', 'categorical', true, 'domain', NetworkDomain, 'transform', @(x) x==2 );

%%% grid variables
PopDomain = Domain.Scattered(fullfile(datadir,'population-epsg3035.mat'));
PopDomain.no_smoothing_method = 'natural';  % to help having no NaNs
pop = Variable( 'pop', stdopt{:}, PopDomain, 'smoothing', pop_bw );

psy = Variable( 'psy', 'transform', 'log', ...
  'smoothing', hospital_bw, 'offset', 1, ...
  'domain', Domain.Scattered(fullfile(datadir,'hospitals.mat')) );

ErwerbsstatDomain = Domain.Scattered(fullfile(datadir,'AERWST_10km_2008.mat'));
ErwerbsstatDomain.no_smoothing_method = 'natural';  % to help having no NaNs
erwopt = [stdopt {ErwerbsstatDomain, 'offset', 1, 'fixNan', 1}];
m = Variable( 'm', erwopt{:} );
w = Variable( 'w', erwopt{:} );
age1 = Variable( 'age1', erwopt{:} );
age2 = Variable( 'age2', erwopt{:} );
age3 = Variable( 'age3', erwopt{:} );
popausl = Variable( 'popausl', erwopt{:} );
jobless = Variable( 'jobless', erwopt{:} );

%% PCA
pcavars = {income gemedu1 gemedu2 gemedu3 m w age1 age2 age3 suirate popausl jobless};
PCA = Domain.PCA( pcavars, cases_mn.net, N_dummies, pop );
PCA.coefftab
comp1 = Variable('comp01', 'domain', PCA );
comp2 = Variable('comp02', 'domain', PCA );

%% build model
disp('***********************')
disp('*** Building Models ***')
disp('***********************')
variables = {'intercept' psy fco pop comp1 comp2};
M = Grid.Model( cases_mn.net, 'loglinear', variables, grid_width );
if save_intermediates; save('model','M'); end
if ismethod(M,'save_variable_images'); M.save_variable_images; end

%% fit model
disp('**********************')
disp('*** Fitting Models ***')
disp('**********************')
fitOptions = { ...
  'name', sprintf('Wiggled Cases by %d km, pop_bw=%.2g', wiggle_amount, pop_bw), ...
  'SaveImages', true, 'SaveRealizationImages', false, ...
  'SaveState', true, ...
  'wiggle_amount', wiggle_amount, ...
  'N_simulations', N_wiggle_fit, ...
  'N_dummies', N_dummies };
fit = M.fit(cases_mn, fitOptions{:} );
fit.resultsTable
if save_intermediates; save('fit', 'fit', '-v7.3'); end
if ismethod(fit,'saveplots'); fit.saveplots; end

%% simulate model
disp('*************************')
disp('*** Simulating Models ***')
disp('*************************')
simOptions = {'name', 'Simulated Cases from Model', ...
  'SaveImages', true, 'SaveRealizationImages', false, ...
  'SaveState', true};
N_sim_cases = cases_mn.nummarks;
[simfit, model_intensity] = fit.simulate(N_sim_cases, N_sim_model, N_dummies, simOptions{:});
simfit.resultsTable
if save_intermediates; save('simfit', 'simfit', 'model_intensity'); end
if ismethod(simfit,'saveplots'); simfit.saveplots; end

%% plot and save pair distributions
disp('*****************************')
disp('*** Analyzing Simulations ***')
disp('*****************************')
simfit.plot_pair_distribution( cases_mn, 'histogram_plot_max', 400 );

%% done
disp('***************** COMPLETE *****************')
toc( timer )
disp('********************************************')
diary off
