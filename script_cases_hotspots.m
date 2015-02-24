function script_cases_hotspots

%% configuration
alpha = 0.05;
save_single_images = false;

grid_width = 1;
Ngrids = 16;

wiggle_amount = 1;
wiggle_reps = 1000;

simfile = 'cases_hotspots.mat';
mainVars = {'GC','Z','mn', 'Ngrids','wiggle_reps', 'Z0','L0', 'coords'};

%% load network
global rr
load cases_mn
rr = cases_mn.net;
mn = cases_mn;

%% load or create other data
if exist(simfile,'file')
  fprintf('Loading variables "%s" from file...\n', strjoin(mainVars,', '));
  load(simfile, mainVars{:});
else

  %% administrative
  if mfilename
    diary(sprintf('%s_diary.txt',mfilename));
  end
  
  %% construct marked network
  [~,G] = load_dataset( 'pop', grid_width );
  
  %% Grids (will be combined)
  GC = GridCollection( G, Ngrids );
  L0 = cell(1,GC.N);  % railroad length per pixel
  Z0 = cell(1,GC.N);
  Z = cell(1, wiggle_reps);
  coords = cell(1, wiggle_reps);
  
  %% rep == 0
  % generate line densities for each grid
  fprintf('Generate line densities...\n');
  parfor Gi = 1:GC.N
    G = GC.Grids{Gi};
    L0{Gi} = sparse(GC.Grids{Gi}.line_density(rr.edges));
    Z0{Gi} = inside_loop(G, mn.coords, L0{Gi});
  end
  
  %% rep > 0
  fprintf('Wiggling...\n');
  parfor rep = 1:wiggle_reps
    
    Z{rep} = cell(1,GC.N);
    
    wiggle_ed = mn.net.rmove( mn.edge_network, wiggle_amount, 'normal' );
    coords{rep} = mn.net.xy_from_edge( wiggle_ed(:,1), wiggle_ed(:,2) );
    
    %% Suicides per Pixel
    for Gi = 1:GC.N
      G = GC.Grids{Gi};
      Z{rep}{Gi} = inside_loop(G, coords{rep});
    end
    fprintf('%2d/%2d √\n', rep, wiggle_reps);
  end
  
  %% save results
  fprintf('Saving results...\n');
  save(simfile, '-v7.3');
end

%%
if exist('save_single_images','var') && save_single_images
  % export single-images, but with a consistent color coding
  fprintf('Saving images...\n');
  basename = 'cases';
  
  for Gi = 1:GC.N
    GC.Grids{Gi}.save_image( L0{Gi}, 'rr_len_per_pixel',...
      'partial', Gi, 'value-range', value_range(L0) );
    for rep = 1:wiggle_reps
      GC.Grids{Gi}.save_image( Z{rep,Gi}, basename, 'rep', rep, ...
        'partial', Gi, 'value-range', value_range(Z) );
    end
  end
end

%% load quantile data, if possible
quantileVars = {'ZqLo', 'ZqLo2', 'ZqLo3', 'Zmedian', 'Zmean', 'ZqHi3', 'ZqHi2', 'ZqHi', 'combLmedian', 'combLmean'};
basename = 'cases_wiggled';
if exist(simfile,'file')
  simfile_mat = matfile(simfile);
  if all( ismember( quantileVars, who(simfile_mat) ) )
    fprintf('Loading variables "%s" from file...\n', strjoin(quantileVars,', '));
    load(simfile, quantileVars{:} );
  else
    %% split grid into smaller pieces and calculate quantile
    fprintf('Calculating quantiles...\n');
    blockfun = @(bbi) GC.quantile( Z, [alpha/2 alpha 2*alpha 0.5 0 1-2*alpha 1-alpha 1-alpha/2], bbi );
    [ZqLo, ZqLo2, ZqLo3, Zmedian, Zmean, ZqHi3, ZqHi2, ZqHi] = GC.G.splitapply( numel(Z), blockfun );
    combL = GC.quantile({L0}, [0.5 0]);
    combLmedian = squeeze( combL(1,:,:) );
    combLmean = squeeze( combL(2,:,:) );
    
    %% save results
    fprintf('Saving results...\n');
    save(simfile, quantileVars{:}, '-append');
    
    %% save result images
    fprintf('Saving result images...\n');
    GC.G.save_image(combLmedian, 'rr_len_per_pixel_median');
    GC.G.save_image(combLmean, 'rr_len_per_pixel_mean');
    GC.G.save_image(ZqLo, [basename,'_qLo'], 'q', alpha/2);
    GC.G.save_image(ZqLo2, [basename,'_qLo'], 'q', alpha);
    GC.G.save_image(ZqLo3, [basename,'_qLo'], 'q', 2*alpha);
    GC.G.save_image(Zmedian, [basename,'_median']);
    GC.G.save_image(Zmean, [basename,'_mean']);
    GC.G.save_image(ZqHi3, [basename,'_qHi'], 'q', 2*alpha);
    GC.G.save_image(ZqHi2, [basename,'_qHi'], 'q', alpha);
    GC.G.save_image(ZqHi, [basename,'_qHi'], 'q', alpha/2);
  end
end


%% analyze lower quantile data for clusters
% load fit to include in analyzationsize
fprintf('Loading fit...\n');
load('fit.mat','fit');

%%
fprintf('Analyzing fit...\n');
T=analyse_pixel_clusters(struct( ...
  'lo',ZqLo3, 'med',Zmedian, 'hi',ZqHi3, 'coords',{coords}, ...
  'grid',GC.G, 'netlen',combLmedian, ...
  'name',[basename,'-qLo-all'], 'grow_cluster',3), ...
  Inf, mn, fit)

%% finito
diary off
return
end


function [CM, cases_per_railkilometer] = inside_loop(G, coords, L)

[~,~,~,CM] = G.boxcount( coords(:,1), coords(:,2) );

if exist('L','var')
  cases_per_railkilometer = sparse( CM ./ L );
  cases_per_railkilometer(L==0) = 0;
  nCasesOutsideNetwork = nnz( CM > 0 & L == 0 );
  assert( nCasesOutsideNetwork == 0, 'there are %d cases outside network', nCasesOutsideNetwork );
end

CM = sparse(CM);

end
