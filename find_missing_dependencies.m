function find_missing_dependencies(scriptName)
    % Read the content of the script
    scriptContent = fileread(scriptName);

    % Extract all function calls from the script
    functionCalls = extractFunctionCalls(scriptContent);

    % Initialize a cell array to hold missing files
    missingFiles = {};

    % Check each function call to see if its .m file exists
    for i = 1:length(functionCalls)
        funcName = functionCalls{i};
        % Try to find the function file
        funcFile = which(funcName);
        if isempty(funcFile)
            missingFiles{end+1} = [funcName '.m']; %#ok<AGROW>
        end
    end

    % Display the missing files
    if isempty(missingFiles)
        disp('All required files are present.');
    else
        disp('Missing files:');
        for i = 1:length(missingFiles)
            disp(missingFiles{i});
        end
    end
end

function functionCalls = extractFunctionCalls(scriptContent)
    % Regular expression to match function calls excluding array indexing
    funcPattern = '(?<=[^\.\w])\b\w+(?=\()';
    tokens = regexp(scriptContent, funcPattern, 'match');

    % Remove duplicates
    functionCalls = unique(tokens);
    
    % Exclude MATLAB built-in functions and keywords
    keywordsAndBuiltins = getMatlabKeywordsAndBuiltins();
    functionCalls = setdiff(functionCalls, keywordsAndBuiltins);
end

function keywordsAndBuiltins = getMatlabKeywordsAndBuiltins()
    % List of MATLAB keywords and common built-in functions
    keywords = {
   
    };

    keywordsAndBuiltins = unique(keywords);
end
