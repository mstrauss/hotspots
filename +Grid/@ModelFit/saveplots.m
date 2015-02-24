function saveplots( obj )

% save model intensity
thisModelInt = obj.intensity(true);  % intensity saves the corresponding figures itself
% Scatter plot
f=figure('units','normalized','outerposition',[0 0 1 1]);
if isempty(obj.parent)
  plot( obj, obj.theta, 'theta = mean parameters' );
else
  % we have a parent
  plot( obj, obj.parent.theta, 'original theta = sim parameters' );
  % save model intensities on the same color scale for comparison
  parentModelInt = obj.parent.intensity;
  vr = value_range({thisModelInt parentModelInt});
  obj.model.G.save_image( thisModelInt, filename(obj.name,'model intensity'), 'value-range', vr );
  obj.model.G.save_image( parentModelInt, filename(obj.name,'parent model intensity'), 'value-range', vr );
end
save_figure(filename(obj.name, 'scattermatrix'));
close(f)
end
