classdef Regular2GridTests < matlab.unittest.TestCase
  
  properties
    gridSize
    grid
    lines
    linesLength  % total lenght of all lines
  end
  
  methods (TestMethodSetup)
    function createGrid( tc )
      tc.gridSize = 100;
      tc.grid = Regular2Grid(1,1,1,1,tc.gridSize,tc.gridSize);
    end
    
    function createRandomLines(tc)
      N = 100;
      lineLength = 10;  % length of a single line
      tc.lines = [0.5+tc.gridSize*rand(N,2) 2*lineLength*rand(N,2)-lineLength];
      tc.linesLength = sum( vnorm(tc.lines(:,3:4)) );
    end
  end
  
  methods (Test)
    
    function testLineDensity1(tc)
      import matlab.unittest.constraints.IsEqualTo;
      import matlab.unittest.constraints.RelativeTolerance;
      [L,O] = tc.grid.line_density(tc.lines);
      tc.assertThat( tc.linesLength, IsEqualTo( sum(L(:)) + O, 'Within', RelativeTolerance(100*eps) ) );
      %if tc.gridSize <= 10; flipud(L'); end
    end
        
    function testDensity(tc)
      L = tc.grid.line_density(tc.lines);
      S = tc.grid.convolution_smoothing([0.8 0; 0 0.8], L);
      figure
      subplot(121) ; imagesc(L); axis xy; colorbar
      subplot(122) ; imagesc(S); axis xy; colorbar
    end
    
  end
  
end

