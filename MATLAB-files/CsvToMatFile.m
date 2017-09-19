%This function takes a string as an argument, which contains the path to a
%csv file containing model parameters.
%Then it assigns the parameter values to variables and saves the work
%space.

function CsvToMatFile(sourceFile)
    
    csvHandle = fopen(sourceFile);
    csvText = textscan(csvHandle, '%s %f', 'Delimiter', ',');
    fclose(csvHandle)
    
    for i = 1:length(csvText{1})
        if strcmp(csvText{1}{i}, '_inputType')
            continue
        end
        
        eval([csvText{1}{i} ' = ' num2str(csvText{2}(i)) ';']);
          
        
    end
    
    %Workspace clean up
    clear csvHandle csvText i ans
    
    
    %This was done inline to keep the .mat file clean
    %matFileName = [csvFile(1:length(csvFile)-4) '.mat'];
    
    save([sourceFile(1:length(sourceFile)-4) '.mat'])
    