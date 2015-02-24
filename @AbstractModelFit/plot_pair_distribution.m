function [] = plot_pair_distribution( simfit, marked_network, varargin )
% the marked_network is required!

nummarks = simfit.configvars.N_cases;
assert( nummarks == marked_network.nummarks );
net = simfit.model.network;

% parse options
argi = 1;
while argi <= numel(varargin)
  switch varargin{argi}
    case 'model_reps'
      argi = argi + 1; model_reps = varargin{argi};
    case 'baseline_reps'
      argi = argi + 1; baseline_reps = varargin{argi};
    case 'wiggle_reps'
      argi = argi + 1; wiggle_reps = varargin{argi};
    case 'wiggle_amount'
      argi = argi + 1; wiggle_amount = varargin{argi};
    case 'histogram_plot_max'
      argi = argi + 1; histogram_plot_max = varargin{argi};
    case 'histogram_calc_max'
      argi = argi + 1; histogram_calc_max = varargin{argi};
    case 'kl_bandwidth'
      argi = argi + 1; kl_bandwidth = varargin{argi};
    otherwise
      error('Unkown option %s',varargin{argi});
  end
  argi = argi + 1;
end

% set defaults
if ~exist('model_reps','var'); model_reps = numel( simfit.realizations ); end
if ~exist('baseline_reps','var'); baseline_reps = model_reps; end
if ~exist('wiggle_reps','var'); wiggle_reps = model_reps; end
if ~exist('wiggle_amount','var'); wiggle_amount = simfit.parent.configvars.wiggle_amount; end
if ~exist('histogram_calc_max','var'); histogram_calc_max = 800; end
if ~exist('histogram_plot_max','var'); histogram_plot_max = histogram_calc_max; end
if ~exist('kl_bandwidth','var'); kl_bandwidth = 1; end

%% do stuff
[ctr, poisson_fields, gHPP, gIPP, gM, nnIPP, nnM] = init_fields;
close all  % reset figures
plot_save_hpp
plot_save_ipp

%% functions

  function [ctr, poisson_fields, gHPP, gIPP, gM, nnIPP, nnM] = init_fields
    % initialize histogram centers
    ctr = histogram_bins( 1, histogram_calc_max );

    % Simulate Baseline
    fprintf('Generating random poisson fields (N=%d) %d times.\n', nummarks, baseline_reps);
    [gHPP,dHPP,nnHPP,poisson_fields] = generate_random_field( baseline_reps, @(i) poisson_random_field(net,nummarks), ctr );
    
    % Simulate Model
    fprintf('Analyzing %d existing model fields (N=%d).\n', model_reps, nummarks);
    ippfun = @(i) MarkedNetwork(net,simfit.realizations{i});
    [gIPP,dIPP,nnIPP] = generate_random_field( model_reps, ippfun, ctr );
    
    % wiggle original data
    fprintf('Wiggling original data %d times.\n', wiggle_reps);
    [gM,dM,nnM] = generate_random_field( wiggle_reps, @(i) wiggled_marks_random_field(marked_network,wiggle_amount), ctr );
  end

  function plot_save_hpp
    % Pair distribution function
    clf
    s = ctr<histogram_plot_max;    
    plot_simulation_envelopes( ctr(s), gHPP(:,s) );
    savefig(filename(simfit.name,'HPP histogram.fig'));
    clf
    
    plot_simulation_envelopes( ctr(s), centered(gHPP(:,s), gHPP(:,s)) );
    savefig(filename(simfit.name,'HPP histogram normalized.fig'));
    clf
  end

  function plot_save_ipp
    s = ctr<histogram_plot_max;
    plot_simulation_envelopes( ctr(s), centered(gIPP(:,s), gHPP(:,s)) ); hold on
    hline = refline(0,1); set(hline,'LineStyle','--', 'Color',[.8 .8 .8], 'LineWidth', 1.5)
    plot( ctr(s), centered(mean(gM(:,s),1), gHPP(:,s)), 'LineWidth', 2 )
    %grid on
    xlabel('Distance r (km)')
    ylabel('Relative Frequency');
    savefig(filename(simfit.name,'g(r).fig'));
    xlim([0,150])
    savefig(filename(simfit.name,'g(r)_closeup.fig'));
    xlim([0,15])
    savefig(filename(simfit.name,'g(r)_closeup_2.fig'));
  end

  function plot_save_knn
    %% mean distance to k-th neighbor
    clf
    k = 2:histogram_plot_max;
    %plot_simulation_envelopes( k-1, centered(nnHPP(:,k))); hold on
    plot_simulation_envelopes( k-1, nnIPP(:,k) ); hold on
    plot(k-1, mean(nnM(:,k))); hold on
    xlabel('Neighborhood rank k')
    ylabel('Mean distance (km)');
    title({'Mean distance to k-th neighbor, nn(k)',...
      'gray: fitted inhom. Poisson process with 95% sim. env.',...
      'blue: actual cases'});
    saveas(gcf, filename(simfit.name,'nn(k).png'));
    xlim([1,100])
    saveas(gcf, filename(simfit.name,'nn(k)_closeup.png'));
    xlim([1,15])
    saveas(gcf, filename(simfit.name,'nn(k)_closeup_2.png'));
  end

  function plot_save_contact
    %% spherical contact distribution
    ctr = histogram_bins( 1, 350 );
    m = numel(ctr);
    D = zeros(numel(poisson_fields),m);
    DN = zeros(numel(poisson_fields),m);
    rnd_fld_size = 1e3;
    for i = 1:numel(poisson_fields);
      rnd_fld = poisson_random_field(net,rnd_fld_size);
      % distances btw. random field and cases
      Ds=net.pdist2(rnd_fld.edge_network, marked_network.edge_network);
      Ds=sort(Ds');
      D(i,:) = hist( Ds(1,:), ctr );
      % distances btw. two random fields (null model)
      Ds=net.pdist2( rnd_fld.edge_network, ...
        poisson_fields{i}.edge_network);
      Ds=sort(Ds');
      DN(i,:) = hist( Ds(1,:), ctr );
    end
    
    clf
    cumD = cumsum(D,2);
    cumDN = cumsum(DN,2);
    plot_simulation_envelopes( ctr, cumD/mean(cumD(:,end)) ); hold on
    plot( ctr, mean(cumDN)/mean(cumDN(:,end)), '-' );
    ylim([0.1,1]);
    %xlim([0,40]);
    set( gca, 'XScale', 'log' );
    set( gca, 'YScale', 'log' );
    title({'Spherical Contact Distribution, H_s(r)', ...
      'gray: actual cases with 95% sim. env.',...
      'blue: null model (homog. Poisson)'});
    xlabel('Distance r to nearest test location (km)')
    ylabel('H_s(r)');
    grid on
    save_figure(filename(simfit.name,'H_s(r)'), 'rnd_fld_size', rnd_fld_size);
  end

  function plot_save_kldiv
    %% KL Divergence Cases / Model
    t = 'Cases/Model';
    %z = kl(cases.netdens', simfit.intensity');
    z = kl(simfit.intensity', simfit.parent.intensity');
    KL = sum(z(:));
    fprintf('KL(%s)=%.2g\n', t, KL);
    simfit.model.G.save_image( z, filename(simfit.name,'kl_div_cases_vs_model_intensities'), 'bw', kl_bandwidth, 'KL', KL );
  end

end
