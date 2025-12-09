classdef tests < matlab.unittest.TestCase
    % The 'properties' section is where you can store variables,
    % but for simple tests, it's often not needed.
    
    methods (Test)
        % This section contains the actual test methods
        
        function testBasicAddition(testCase)
            % This is a test method. It must accept 'testCase' as its first argument.
            
            dvec = DesignVector();
            wingDesign = WingDesign(dvec);

        end
        
        
    end
end