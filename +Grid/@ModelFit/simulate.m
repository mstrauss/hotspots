function [simfit, model_intensity] = simulate( obj, N_sim_cases, repetitions, N_dummies, varargin )
% simulate the given model repetitions times

assert( isnumeric(N_sim_cases) );
assert( isnumeric(repetitions) );
assert( isnumeric(N_dummies) );

% parse options
argi = 1;
while argi <= numel(varargin)
  switch varargin{argi}
    case 'SaveImages'
      argi = argi + 1; SaveImages = varargin{argi};
  end
  argi = argi + 1;
end

% default options
if ~exist('SaveImages','var'); SaveImages = false; end

  
if ~exist('N_dummies','var')
  %N_dummies = obj.results{1}.fitglm_output.NumObservations - size(obj.realizations{1},1);
  N_dummies = obj.configvars.N_dummies;
end

model_intensity = obj.intensity( SaveImages );
assert( ~any(isinf(model_intensity(:))) );
assert( ~any(isnan(model_intensity(:))) );
assert( all(model_intensity(:) >= 0 ), 'require non-negative model intensities' );

grid = obj.model.grid;
gridded_network = GriddedNetwork( obj.model.network, grid );
simOptions = {'N_dummies', N_dummies, 'N_simulations', repetitions};
simOptions = [simOptions, varargin];
simfit = obj.model.fit( ...
  @() gridded_network.rpoisinhom( N_sim_cases, model_intensity ), simOptions{:} );
simfit.parent = obj;

% if SaveImages
%   simfit.saveplots;
% end
