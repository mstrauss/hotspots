classdef NetworkFactory
  
  properties
  end
  
  methods (Static)
    function net = VerySimpleLoopfreeNetwork
      net = Network([0 0; 1 0; 2 1; 3 0], struct('adjacency', [0 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 1 0]));
    end
    
    function net = VerySimpleSingleloopNetwork
      net = Network([0 0; 1 0; 2 1; 3 0], struct('adjacency', [0 1 0 0; 1 0 1 1; 0 1 0 1; 0 1 1 0]) );
    end
    
    function mn = VerySimpleSingleloopMarkedNetwork
      mn = MarkedNetwork( tests.NetworkFactory.VerySimpleSingleloopNetwork, ...
        [2 0.1; 3 0.4; 3 1.9] );
    end
    
    function net = SimpleSyntheticNetwork
      g = importdata('+tests/SimpleSyntheticNetwork.mat');
      net = Network( g.xy, struct('adjacency', g.A) );
    end
    
    function mn = SimpleSyntheticMarkedNetwork
      net = tests.NetworkFactory.SimpleSyntheticNetwork;
      c =[ 5.1 2.2; 7.1 3.3; 5.1 3.3; 2.1 2; 2.9 2; 3.9 2 ];
      [~,~,edges] = net.project(c);
      mn = MarkedNetwork( net, edges );
    end
    
    function mn = RealisticMarkedNetwork
      load cases_mn
      mn = cases_mn;
    end
    
    function net = ChicagoStreetNetwork
      % credits : extracted from R's SPATSTAT package, manually digitized
      % by Adrian Baddeley (see '? chicago' in R);  export command:
      % write.mat( chicago$domain$lines$ends, 'chicago_lines.mat')
      lines = load('+tests/data/chicago_lines.mat');
      assert( isstruct(lines) );
      lines = cell2mat(struct2cell(lines))';
      net = Network.createFromLines( lines );
    end
    
    function mn = ChicagoMarkedStreetNetwork
      % credits : extracted from R's SPATSTAT package, manually digitized
      % by Adrian Baddeley (see '? chicago' in R);  export command:
      % write.mat(as.data.frame(chicago$data),'chicago_marks.mat')
      net = tests.NetworkFactory.ChicagoStreetNetwork;
      marks = load('+tests/data/chicago_marks.mat');
      [~,~,edges] = net.project( [marks.x' marks.y'] );
      mn = MarkedNetwork( net, edges );
    end
  end
  
end
