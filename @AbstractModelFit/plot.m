function plot( obj, theta, theta_text )
if ~exist('theta','var')
  theta = obj.theta;
  theta_text = 'fitting/simulation mean';
end
m = obj.number_of_meta_parameters();
plot_scatter( obj.resultsMatrix(:,1:end-m), obj.varnames, theta, theta_text );
end
