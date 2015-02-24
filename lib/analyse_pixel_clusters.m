function [T, threshold] = analyse_pixel_clusters( datastruct, nClusters, cases_network, fit )

assert( isa(cases_network,'MarkedNetwork') );
assert( isa(fit,'Grid.ModelFit') );

%%

if isstruct(datastruct)
  imlo = datastruct.lo;   % lower bound data
  immed = datastruct.med; % median data
  imhi = datastruct.hi;   % upper bound data
  G = datastruct.grid;
  if isfield(datastruct,'netlen')
    L = datastruct.netlen;
    assert( all( size(L) == size(G) ) );
  end
  file_basename = datastruct.name;
else
  error('need struct input');
end

assert( all( size(G) == size(imlo) ) );
assert( all( size(G) == size(imhi) ) );
assert( all( size(G) == size(immed) ) );

if nClusters < Inf
  % find optimal threshold so we find exactly nClusters clusters
  threshold = fzero( @objective, max(imlo(:)) );
else
  threshold = 0;
end

imlo_thresholded = imlo;
imlo_thresholded(imlo_thresholded<threshold) = 0;

%% identify connected components (= clusters) from convoluted image
% (grown by border of width 1)
if isfield(datastruct,'grow_cluster') && datastruct.grow_cluster
  imlo_conv = conv2( imlo_thresholded, ones(datastruct.grow_cluster+1), 'same' );
else
  imlo_conv = imlo_thresholded;
end
cc = bwconncomp( imlo_conv, 8);
Ncc = cc.NumObjects;

% keep only pixels which exist on the original image
for cci = 1:Ncc
  pixels = cc.PixelIdxList{cci};
  keepPixels = ismember( pixels, find(imlo_thresholded) );
  cc.PixelIdxList{cci} = pixels(keepPixels);
end

%% find quantiles of cases for identified clusters
numReal = numel(datastruct.coords);
cluster_case_counts = cell(numReal,Ncc);
for real=1:numReal  % loop over realizations
  coords = datastruct.coords{real};
  for cci=1:Ncc  % loop over clusters
    [~,~,~,CountsMatrix] = G.boxcount( coords(:,1), coords(:,2) );
    cluster_case_counts{real,cci} = sum( CountsMatrix( cc.PixelIdxList{cci} ) );
  end
end

sorted_cluster_case_counts = sort( cell2mat( cluster_case_counts ), 1 );

% cluster cases
alpha = 0.05;
clcLo = sorted_cluster_case_counts( numReal*alpha, : );
clcMed = sorted_cluster_case_counts( numReal/2, : );
clcHi = sorted_cluster_case_counts( numReal*(1-alpha), : );

%% project cases onto pixel grid
cases_coords = cases_network.coords;
[~,~,~,CountsMatrix] = G.boxcount( cases_coords(:,1), cases_coords(:,2) );

% save results to table T
Covariates = zeros( Ncc, numel(fit.theta) );
C = cell(Ncc,1);
M = cell(Ncc,1);
S = struct.empty;
for i = 1:Ncc
  ci = cc.PixelIdxList{i};
  ni = length(ci);
  centroid = sqrt(mean( G.ind2coord(ci).^2, 1 ));
  maxoidI = ci( find( imlo(ci)==max(imlo(ci)) ) );
  maxoid = sqrt(mean( G.ind2coord( maxoidI ).^2, 1 ));
  total_weight_lo = sum( imlo(ci) );
  total_weight_hi = sum( imhi(ci) );
  S(i).id = i;
  S(i).median_cluster_cases = [clcLo(i) clcMed(i) clcHi(i)];
  S(i).actual_cluster_cases = sum( CountsMatrix(ci) );
  S(i).cluster_pixels = ni;
  if exist('L','var')
    netlen = sum( L(ci) );
    S(i).cluster_netlen = netlen;
  end
  
  S(i).no_maxoid_pixels = numel(maxoidI);

  for j = 1:numel(fit.theta)
    covariate = G.interpolate( maxoidI, fit.model.G, fit.model.GX(j) );
    Covariates(i,j) = mean( exp( fit.theta(j) * covariate ) );
  end
  C{i} = sprintf('%d,%d', int32(1000*centroid(1)), int32(1000*centroid(2)));
  M{i} = sprintf('%d,%d', int32(1000*maxoid(1)), int32(1000*maxoid(2)));
end
T = struct2table(S);
T.centroid_epsg3035 = C;
T.maxoid_epsg3035 = M;

% excess risk
for j = 1:numel(fit.theta)
  var = fit.model.varnames(j);
  T.(['er_',var{:}]) = Covariates(:,j);
end
[~,order] = sort(T.median_cluster_cases(:,3),'descend');
T = T(order,:);

%% save to textfile
writetable( T, [file_basename,'.csv'], 'Delimiter', ';' );

%% plot
G.save_image(imlo_thresholded, file_basename);

  function y = objective( threshold )
    imlo_thresholded = imlo;
    imlo_thresholded(imlo_thresholded<threshold) = 0;
    cc = bwconncomp(imlo_thresholded, 8);
    y = cc.NumObjects - nClusters;
  end

end
