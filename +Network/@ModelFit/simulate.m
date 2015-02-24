function [simulated_point_patterns, model_intensity] = simulate( obj, varargin )
% simulate the given model repetitions times

assert( ~isempty(obj.marked_network), 'Invalid marked_network. Cannot simulate.' );

% parse options
argi = 1;
while argi <= numel(varargin)
  switch varargin{argi}
    case 'SaveImages'
      argi = argi + 1; SaveImages = varargin{argi};
    case 'repetitions'
      argi = argi + 1; repetitions = varargin{argi};
    case 'N_dummies'
      argi = argi + 1; N_dummies = varargin{argi};
  end
  argi = argi + 1;
end

% default options
if ~exist('SaveImages','var'); SaveImages = false; end
if ~exist('repetitions','var'); repetitions = obj.configvars.N_simulations; end
if ~exist('N_dummies','var'); N_dummies = obj.configvars.N_dummies; end
% if ~exist('N_dummies','var')
%   N_dummies = obj.results{1}.fitglm_output.NumObservations - size(obj.realizations{1},1);
% end

model_intensity = obj.intensity;
assert( ~any(isinf(model_intensity(:))) );
assert( ~any(isnan(model_intensity(:))) );
assert( all(model_intensity(:) >= 0 ), 'require non-negative model intensities' );

mn = obj.marked_network;
assert( mn.net == obj.model.network );

N_sim_cases = mn.N;  % number of marks

% for now, calculate ALL realizations of point patterns from THE SAME model
% intensity object.
for rep = 1:repetitions
  simulated_point_patterns(rep) = mn.net.rpoisinhom( N_sim_cases, model_intensity );
end
