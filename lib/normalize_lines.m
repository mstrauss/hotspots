function [ lines, sortindex ] = normalize_lines( lines )
% NORMALIZE_LINES reorganizes the data such that dx is non-negative
% and sorts by x, y, dx, dy.
%   % lines: rows of [x y dx dy]

% % pre-diagnostic
% x = lines(:,1); y = lines(:,2);
% dx = lines(:,3); dy = lines(:,4);
% figure; plot( [x x+dx]', [y y+dy]', 'x-b' ); hold on

with_negative = lines(:,3) < 0;
lines( with_negative, 1) = lines( with_negative, 1) + lines( with_negative, 3);
lines( with_negative, 3) = -lines( with_negative, 3);
lines( with_negative, 2) = lines( with_negative, 2) + lines( with_negative, 4);
lines( with_negative, 4) = -lines( with_negative, 4);

% sorting
[lines, sortindex] = sortrows( lines );

% % post-diagnostic
% x = lines(:,1); y = lines(:,2);
% dx = lines(:,3); dy = lines(:,4);
% plot( [x x+dx]', [y y+dy]', 'x-r' ); hold off
% figure
end
