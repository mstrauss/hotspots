classdef NetworkTests < matlab.unittest.TestCase
  
  properties
  end
  
  methods (TestMethodSetup)
  end
  
  methods (Test)
    function testMoveFix(testCase)
      import matlab.unittest.constraints.IsEqualTo;
      import matlab.unittest.constraints.AbsoluteTolerance;
      import tests.NetworkFactory.VerySimpleLoopfreeNetwork;
      % network without loops
      n = VerySimpleLoopfreeNetwork;
      figure; plot(n, true);
      testCase.assertThat( n.move_fix([1 3]), IsEqualTo([3 3-sqrt(2)-1]) );
      testCase.assertThat( n.move_fix([3 -2]), IsEqualTo([1 sqrt(2)-1]) );
      testCase.assertThat( n.move_fix([1 -0.2]), IsEqualTo([1 0.2]) );
      testCase.assertThat( n.move_fix([3 2]), IsEqualTo([3 sqrt(2)-(2-sqrt(2))]) );
      testCase.assertThat( n.move_fix([1 4]), IsEqualTo([3 sqrt(2)-(4-2*sqrt(2)-1)]) );
      testCase.assertThat( n.move_fix([1 6]), IsEqualTo([2 sqrt(32)-5], 'Within', AbsoluteTolerance(2*eps)) );
    end
    
    function testMoveFixLoopNetwork(testCase)
      import matlab.unittest.constraints.IsEqualTo;
      import matlab.unittest.constraints.AbsoluteTolerance;
      import tests.NetworkFactory.VerySimpleSingleloopNetwork;
      % network with single loop
      n2 = VerySimpleSingleloopNetwork;
      figure; plot(n2,true);
      n2.rmove([1 0; 1 1], 400);
    end
    
    function testShortestDistanceCalculation1(testCase)
      import matlab.unittest.constraints.IsEqualTo;
      n = Network([0 0; 1 0; 2 1; 3 0; 4 0], struct( 'adjacency', [0 1 0 0 0; 1 0 1 0 0; 0 1 0 0 0; 0 0 0 0 1; 0 0 0 1 0], ...
        'distance', 'dijkstra') );
      l1 = sqrt(2); l2 = l1+1;
      testCase.assertThat( n.D, IsEqualTo([0 1 l2 Inf Inf; 1 0 l1 Inf Inf; l2 l1 0 Inf Inf; Inf Inf Inf 0 1; Inf Inf Inf 1 0]) );
    end
    
    function testShortestDistanceCalculation2(testCase)
      import matlab.unittest.constraints.IsEqualTo;
      import matlab.unittest.constraints.AbsoluteTolerance;
      g = importdata('/Users/markus/Documents/MATLAB/Apps/gaimc/graphs/cores_example.mat');
      n = Network( g.xy, struct('adjacency', g.A, 'distance', 'dijkstra') );
      load +tests/testShortestDistanceCalculation2.mat
      testCase.assertThat( n.D, IsEqualTo(A, 'Within', AbsoluteTolerance(10*eps)) );
    end
    
  end
  
end

