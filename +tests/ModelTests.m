classdef ModelTests < matlab.unittest.TestCase

  properties
  end
  
  methods (TestMethodSetup)
  end
  methods (Test)
    function createModel(tc)
      import tests.NetworkFactory.SimpleSyntheticNetwork;
      net = SimpleSyntheticNetwork;
      Grid.Model( net, 'loglinear', {}, 1 );
    end
  end
  
end
