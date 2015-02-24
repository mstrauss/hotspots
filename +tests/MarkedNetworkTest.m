classdef MarkedNetworkTest < matlab.unittest.TestCase
  
  properties
     mn
  end
  
  methods (TestMethodSetup)
    function createMarkedNetwork(tc)
      import tests.NetworkFactory.SimpleSyntheticMarkedNetwork;
      tc.mn = SimpleSyntheticMarkedNetwork;
    end
  end
  
  methods (Test)
    
    function testMarksDistances(testCase)
      import matlab.unittest.constraints.IsEqualTo;
      import matlab.unittest.constraints.AbsoluteTolerance;
      v = [ 0 Inf sqrt(2)+1.2 3 2.2 1.2;  Inf 0 Inf Inf Inf Inf;  sqrt(2)+1.2 Inf 0 sqrt(2)+2.2 sqrt(2)+1.4 sqrt(2)+0.4;  3 Inf sqrt(2)+2.2 0 0.8 1.8;  2.2 Inf sqrt(2)+1.4 0.8 0 1;  1.2 Inf sqrt(2)+0.4 1.8 1 0 ]';
      order = testCase.mn.edge_coordinates_ordering;
      v = v(order,order);
      testCase.assertThat( testCase.mn.MD, IsEqualTo(v, 'Within', AbsoluteTolerance(2*eps)) )
    end
    
    function testMarksDistancesWithDuplicateMarks(testCase)
      import matlab.unittest.constraints.IsEqualTo;
      import matlab.unittest.constraints.AbsoluteTolerance;
      c = [ 5.1 2.2; 5.1 3.3; 5.1 3.3 ];
      [~,~,edges]=testCase.mn.net.project(c);
      mn = MarkedNetwork( testCase.mn.net, edges );
      v=[sqrt(2)+1.2 sqrt(2)+1.2 0];
      testCase.assertThat( mn.MD, IsEqualTo(squareform(v), 'Within', AbsoluteTolerance(2*eps)) );
      testCase.assertThat( mn.nummarks, IsEqualTo(3) );
    end
    
    function testEdgeNetworkRepresentation(testCase)
      import matlab.unittest.constraints.IsEqualTo;
      c =[ 5.1 2.2; 7.1 3.3; 5.1 3.3; 2.1 2; 2.9 2; 3.9 2; 2 -.5; 6.1 1.9];
      [~,~,edges]=testCase.mn.net.project(c);
      mn = MarkedNetwork( testCase.mn.net, edges );
      testCase.assertThat( mn.coords_from_edge_network, IsEqualTo( mn.coords ) );
    end
    
    function testWeights0(testCase)
      import matlab.unittest.constraints.IsEqualTo;
      import matlab.unittest.constraints.AbsoluteTolerance;
      w = testCase.mn.weights;
      testCase.assertThat( sum(w), IsEqualTo( sum( testCase.mn.net.L ), 'Within', AbsoluteTolerance(100*eps) ) );
    end
    
    function testWeights1(testCase)
      import matlab.unittest.constraints.IsEqualTo;
      import matlab.unittest.constraints.RelativeTolerance;
      mn = MarkedNetwork( Network( [ 0 0; 1 0; 1.5 0.5 ], struct('adjacency', [ 0 1 0; 1 0 1; 0 1 0 ]) ), ...
        [ 1 0.6 ] );
      w = mn.weights;
      % the whole network length must weigh on the single mark
      testCase.assertThat( sum(w), IsEqualTo( 1/sqrt(2)+1, 'Within', RelativeTolerance(2*eps) ) );
    end
    
    function testWeights2(testCase)
      import matlab.unittest.constraints.IsEqualTo;
      import matlab.unittest.constraints.RelativeTolerance;
      mn = MarkedNetwork( ...
        Network( [ 0 0; 1 0; 1+1/sqrt(2) 1/sqrt(2); 1+1/sqrt(2) -1/sqrt(2) ], ...
        struct('adjacency',[ 0 1 0 0; 1 0 1 1; 0 1 0 0; 0 1 0 0 ]) ), ...
        [ 1 0.6; 2 0.8; 3 0.2 ] );
      w = mn.weights;
      testCase.assertThat( sum(w), IsEqualTo( 3, 'Within', RelativeTolerance(2*eps) ) );
    end
    
    function testWeightsWithDuplicateMarks(testCase)
      import matlab.unittest.constraints.IsEqualTo;
      import matlab.unittest.constraints.RelativeTolerance;
      c = [ 5.1 2.2; 5.1 3.3; 5.1 3.3 ];
      [~,~,edges]=testCase.mn.net.project(c);
      mn = MarkedNetwork( testCase.mn.net, edges );
      [w,unacc] = mn.weights;
      testCase.assertThat( sum(w), IsEqualTo( sum(mn.graph.Edges.Weight(:))-unacc, 'Within', RelativeTolerance(2*eps) ) );
    end
    
    function testPlot(tc)
      figure; tc.mn.plot
    end
    
%     function testBipartite(tc)
%       tc.assertTrue( tc.mn.bipartite.is_bipartite )
%     end
  end
  
end

