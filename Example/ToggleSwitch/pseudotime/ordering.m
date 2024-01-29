
clear all
close all
clc

% ----------------------------------------------------------------------- %
%                              BRANCH CLUSTERING                          %
% ----------------------------------------------------------------------- %

startClustering = input('\n Start clustering: '); % e.g. 0 or 1
fprintf('\n')

if startClustering
    
    fprintf('\n   *** Clustering Algorithm ***\n')
    
    cd ..
    figure, load('diffusionInput'); data = DataDiff; clear DataDiff; load('psi');
    subplot(1,2,1), scatter(psi(:,1),psi(:,2),50,data(3,:)), xlabel('DC_1'), ylabel('DC_2')
    subplot(1,2,2), scatter(psi(:,1),psi(:,3),50,data(3,:)), xlabel('DC_1'), ylabel('DC_3')
    clear
    cd pseudotime
    
    % Select diffusion components to use for the clustering algorithm. The two
    % diffusion components where the branch is visually better defined, should
    % be chosen. Depending on the branch we want to cluster, we can choose DC1
    % and DC2, or DC1 and DC3. Theoretically we can include DC4, but we don't.
    
    chooseDim = input('\n Select two diffusion components: '); % e.g. [1 2] or [1 3] 
    fprintf('\n')
    
    chooseAxi = input('\n Select if the branch is on the positive or negative x-axis (1 for positive, -1 for negative): ');
    fprintf('\n')
    close
    
    nneigh = input('\n Select number of neighbours (default is 15): ');
    fprintf('\n')
    
    cd ..
    load('diffusionInput'); data = DataDiff'; clear DataDiff
    load('psi');
    cd pseudotime

    n = length(data);

    d1 = chooseDim(1); % xaxis
    d2 = chooseDim(2); % yaxis

    % Select starting and final cell. Enlarge the window size first, in order
    % to get a better approximate position for the starting and final cells.
    
    scatter(psi(:,d1),psi(:,d2),50,data(:,2)), title('Select starting cell C_S and then final cell C_F','FontSize',15)

    [xX,yY] = ginput(2);
    puntvia = [];

    for i = 1:length(xX)
        for j = 1:size(psi,1)
            NorP(j) = norm(psi(j,[d1 d2])-[xX(i) yY(i)]);
        end
        [aval aind] = min(NorP);
        puntvia = [puntvia aind];
    end    
    dataDM = psi(:,1:4);
    [~, idxs, ~] = nrsearch(dataDM', 1:n , 10, 0,struct('FastNNSearcher','best'));

    D = zeros(n,n);
    for i = 1:n
        nns = idxs{i};
        for j = 1:size(nns,2)
            D(i,nns(j))=dot(dataDM(i,:)-dataDM(nns(j),:),dataDM(i,:)-dataDM(nns(j),:));  
        end
    end
    Adj                   = max(max(D)) - D;
    Adj(Adj==max(max(D))) = 0;
    W2{1}                 = Adj + Adj';

    finalpath = puntvia(1);
    for i=1:length(puntvia)-1;
        [d,p] = dijk(W2{1},puntvia(i),puntvia(i+1));
        finalpath = [finalpath p(2:end)];
    end
    close all; figure
    scatter(psi(:,d1),psi(:,d2)), hold on, plot(psi(finalpath,d1),psi(finalpath,d2),'g.','MarkerSize',40)
    [~, idxsN, ~] = nrsearch(dataDM', finalpath , nneigh, 0,struct('FastNNSearcher','best'));

    br = cell2mat(idxsN);
    br = reshape(br,1,size(br,1)*size(br,2));

    Br = [];
    for i = 1:length(br)
        if chooseAxi==-1
            if psi(br(i),d1) < xX(1)
                Br = [Br br(i)];
            end
        end
        if chooseAxi==1
            if psi(br(i),d1) > xX(1)
                Br = [Br br(i)];
            end
        end
    end

    % Remove repeated cells (this happens because neighbors can be shared).
    sBr = sort(Br);
    BrN = sBr(1);
    diffBr = diff(sBr);
    for i = 1:length(diffBr)
        if diffBr(i)>0
            BrN = [BrN sBr(i+1)];
        end
    end

    scatter(psi(:,d1),psi(:,d2)), hold on, plot(psi(BrN,d1),psi(BrN,d2),'ro')
    figure, scatter3(psi(:,1),psi(:,2),psi(:,3)), hold on, scatter3(psi(BrN,1),psi(BrN,2),psi(BrN,3),90,'fill')

    clear Adj Br D NorP W2 aind aval br d d1 d2 diffBr finalpath NorP aind
    clear idxs idxsN j k n nneigh nns nway p sBr waypoints xX yY puntvia

    save WanderStart

end

% ----------------------------------------------------------------------- %
%                           WANDERLUST ALGORITHM                          %
% ----------------------------------------------------------------------- %

clear all
close all
clc

load('WanderStart')

fprintf('\n   *** Wanderlust Algorithm ***\n')

Nneigh = input('\n Select number of neighbours (e.g. 18): ');
fprintf('\n In case of error of the algorithm , increase the number of neighbours. \n\n')

d1 = chooseDim(1); % xaxis
d2 = chooseDim(2); % yaxis

psi0  = psi(BrN,:);
data0 = data(BrN,:)';

figure, scatter(psi0(:,d1),psi0(:,d2)), title('Select 10 consecutive waypoints along the branch, starting with initial cell and ending with final cell','FontSize',15)

[xX,yY]   = ginput(10);
waypoints = [];
for i = 1:10
    for j = 1:size(psi0,1)
        NorP(j) = norm(psi0(j,[d1 d2])-[xX(i) yY(i)]);
    end
    [aval aind] = min(NorP);
    waypoints = [waypoints aind];
end

init_cell_id = waypoints(1);
waypoints    = waypoints(2:end);
nl           = size(waypoints,2);
dataDM0      = data0';
%dataDM0      = psi0;

L        = max(dataDM0(:,1:size(dataDM0,2)))-min(dataDM0(:,1:size(dataDM0,2)));
psi      = dataDM0(:,1:size(dataDM0,2))./(repmat(L,size(dataDM0,1),1));
data     = psi;

n        = size(data,1);
k        = Nneigh;

[~, idxs, ~] = nrsearch(data', 1:n , k, 0,struct('FastNNSearcher','best'));

replica        = 5;
t1             = zeros(replica,n);
dist_waypoints = zeros(replica,n,nl);
W2             = cell(replica,1);

for replicate = 1:replica
    
    d2 = zeros(n,n);
    for i = 1:n
        choice = floor((rand( 1,ceil(k*1/3) )*k)+1);
        nns    = idxs{i}(choice);
        for j = 1:size(nns,2)
            d2(i,nns(j)) = dot(data(i,:)-data(nns(j),:),data(i,:)-data(nns(j),:));  
        end
    end


    Adj                    = max(max(d2))-d2;
    Adj(Adj==max(max(d2))) = 0;
    W2{replicate}          = Adj+Adj';


    for i = 1:n;
        for l = 1:nl
            [d,p] = dijk(W2{replicate},waypoints(l),i);
            dist_waypoints(replicate,i,l)=d;
        end
    end

end

counts = 10;
trj    = zeros(replica,counts,n);
t      = zeros(nl,n);

for replicate = 1:replica 
    trj(replicate,1,:) = dist_waypoints(replicate,:,1);
    for count = 2:counts
        for i = 1:n;
            for l = 1:nl
    
                if trj(replicate,count-1,i)<trj(replicate,count-1,waypoints(l))
                    t(l,i) = -dist_waypoints(replicate,i,l)+trj(replicate, count-1,waypoints(l) );
                else
                    t(l,i) = dist_waypoints(replicate,i,l)+trj (replicate, count-1,waypoints(l) );
                end
        
            end    
            trj(replicate,count,i) = mean(t(:,i));
        end
    end
end 

trjmean   = mean(trj(:,count,:),1);
trjvar    = var(trj(:,count,:),1);
[val,id]  = sort((trjmean(1,1,:)));
orderedID = reshape(id,1,n);

figure
plot(data0(:,orderedID(end:-1:1))','LineWidth',2)
set(gca,'FontSize',27,'FontName','Courier')
xlabel('pseudo-time'), ylabel('gene expression')
axis tight
legend('Gene A','Gene B','Gene C','Gene D','Gene E','Gene F')

replicate = [orderedID; data0];

% Change this every time to save 3 different replicates for each branch,
% Branch_1_rep_1, Branch_1_rep_2, Branch_1_rep_3
% Branch_2_rep_1, Branch_2_rep_2, Branch_2_rep_3
cd replicates
filename = 'Branch_1_rep_1';
%save(filename, 'replicate')
cd ..

