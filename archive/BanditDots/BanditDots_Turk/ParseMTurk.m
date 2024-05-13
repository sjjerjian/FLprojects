function data = ParseMTurk(inputFile, outputFile)
  % Take inputFile and parse it into a Matlab struct, then save that
  % struct with the name 'data' in outputFile, and return it. If no 
  % output file is given, just return it.
  
  % Parse this file and save
  inputMatrix = readtext(inputFile, ',', '', '"');
  inputMatrix = cellfun(@StripQuotes, inputMatrix, 'UniformOutput', false);
	data = ParseMatrix([], inputMatrix);
  data = ConvertToArrays(data);
  if nargin == 2
    save(outputFile, 'data');
  end
end

function data = ParseMatrix(data, inputMatrix)
  % For each column
	for i=1:size(inputMatrix,2)
    
    % Name of column
    columnName = inputMatrix{1,i};
    
    % Remove the header row
		curCells = inputMatrix(2:end,i);
    
    % If all cells are NaNs or numeric, convert to a vector rather than a
    % cell array
		if all(cellfun(@(x)(length(x)==1&&(isnan(x)||isnumeric(x))), curCells))
			curCells = cell2mat(curCells);
    end    
    
    % If the name of this column has a . in it, parse it into sub structs
    % Also split _ into sub structs.
    data = ParseColumn(columnName, curCells, data);
  end
end

function curStruct = ParseColumn(columnName, curCells, curStruct)
  if ~any(columnName == '.' | columnName == '_')
    % No substructure: just add this column to the struct, and return
    columnName = CleanToFieldName(columnName); 
    if isfield(curStruct, columnName)
      curStruct.(columnName) = ...
        [curStruct.(columnName); curCells];
    else
      curStruct.(columnName) = curCells;
    end    
    return;
  end
  
  % Split by '.' and '_' and recursively call again...
  splitLoc = min([find(columnName=='.',1), find(columnName=='_',1)]);
  beforeName = CleanToFieldName(columnName(1:splitLoc-1));
  afterSplit = columnName(splitLoc+1:end);
  if ~isfield(curStruct, beforeName)
    curStruct.(beforeName) = [];
  end
  curStruct.(beforeName) = ParseColumn(afterSplit, curCells, curStruct.(beforeName));
end

function data = ConvertToArrays(data)
  fieldN = fieldnames(data);
  for s = 1:length(fieldN)
    if isstruct(data.(fieldN{s}))
      % Recursively look for other struct arrays:
      data.(fieldN{s}) = ConvertToArrays(data.(fieldN{s}));
      
      % If this item is names abc12 or any similar, convert to an array:
      if isstrprop(fieldN{s}(end), 'digit')
        m = regexp(fieldN{s},'(.*?)(\d+)$','tokens');
        name = m{1}{1};
        num = str2double(m{1}{2});
        if num == 1
          data.(name) = data.(fieldN{s});
        else
          data.(name)(num) = data.(fieldN{s});
        end
        data = rmfield(data, fieldN{s});
      end
    end
  end
end

function c = CleanToFieldName(columnName)
  c = columnName(isstrprop(columnName,'alphanum'));
end

function x = StripQuotes(x)
  % Strip quotes from any cell with quotes around it
  if ischar(x)
    x = regexprep(x, '"', '');
  end
end
