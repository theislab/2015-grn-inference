
function out = maxDataLength()

nfiles = 20;
T      = [1 10:10:100];

for temperature = T
    
    chain  = [];

    for i = 20:nfiles    
    
        filename = strcat('Risulta_',int2str(temperature),'_',int2str(i));
        estPara  = importdata(filename);
        chain    = [chain estPara'];
    
    end

    FinalChain = chain(:,1:100:end)';
    save(strcat('Risulta_',int2str(temperature)),'FinalChain','-ascii')
    
    fprintf('1/2 -- Iteration %d of %d\n',find(temperature==T),length(T))
    
end


val = [];

for j = 1:length(T)
    
    name = strcat('Risulta_',int2str(T(j)));
    estPara = importdata(name);
    
    val = [val length(estPara)];
    
    fprintf('2/2 -- Iteration %d of %d\n',j,length(T))
    
end

out = min(val);