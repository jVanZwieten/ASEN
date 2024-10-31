classdef Triangle3d
    properties
        Vertices
    end
    
    properties (Dependent)
        nV
        nHat
        A 
    end
    
    methods
        function obj = Triangle3d(x1, x2, x3)
            obj.Vertices = {x1 x2 x3};
        end

        function nV = get.nV(obj)
            nV = cross(obj.Vertices{2} - obj.Vertices{1}, obj.Vertices{3} - obj.Vertices{1});
        end
        
        function nHat = get.nHat(obj)
            nHat = Utilities.unitVector(obj.nV);
        end
        
        function A = get.A(obj)
            A = 0.5 * norm(obj.nV);
        end
    end
end
