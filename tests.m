classdef tests < matlab.unittest.TestCase
    % The 'properties' section is where you can store variables,
    % but for simple tests, it's often not needed.
    
    methods (Test)
        % This section contains the actual test methods
        
        function testBasicAddition(testCase)
            % This is a test method. It must accept 'testCase' as its first argument.
            
            % 1. Arrange (set up inputs)
            inputA = 2;
            inputB = 3;
            expectedResult = 5;
            
            % 2. Act (call the function being tested)
            actualResult = myAddFunction(inputA, inputB); 
            
            % 3. Assert (check if the result is correct)
            testCase.verifyEqual(actualResult, expectedResult, ...
                'The addition result did not match the expected value.');
        end
        
        function testNegativeInput(testCase)
            inputA = -5;
            inputB = 10;
            expectedResult = 5;
            
            actualResult = myAddFunction(inputA, inputB);
            
            testCase.verifyEqual(actualResult, expectedResult, ...
                'Test failed for negative input.');
        end
    end
end