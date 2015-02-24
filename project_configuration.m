%% Project Variables
global projectdir datadir pop_default_file pop_epsg3035_default_file ...
    cases_epsg3035_default_file hospitals_default_file cachedir
projectdir = fileparts( mfilename('fullpath') );

datadir = fullfile(projectdir,'data');
pop_default_file = fullfile(datadir, 'population.mat');
pop_epsg3035_default_file = fullfile(datadir, 'population-epsg3035.mat');
cases_epsg3035_default_file = fullfile(datadir, 'suizide-epsg3035.mat');
hospitals_default_file = fullfile(datadir, 'hospitals.mat');

cachedir = fullfile(projectdir,'cache');

%% add code and data paths
addpath(genpath(fullfile(projectdir,'lib')));
addpath(datadir);

%% Mini functions
h = @(x) x(2)-x(1);
