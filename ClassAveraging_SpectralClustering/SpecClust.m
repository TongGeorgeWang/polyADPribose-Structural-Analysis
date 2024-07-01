function C = SpecClust(AdjacencyMatrix, nKmeans)

warning('off','all') 

%% Perform spectral clustering on a graph, whose information is contained in an adjacency matrix 
%   AdjacencyMatrix - Input adjacency matrix (in this case of binary macromolecular pairings, needs to be square
%   nKmeans - Number of K-means clusters to compute
%
%   References: 
%       - Ng, Jordan, Weiss: On spectral clustering and an algorighm, NIPS 2001
%       - von Luxburg: A tutorial on spectral clustering, Statistics and Computing 2007
%       some code implementations adapted from Ingo Burk, 'areslp', Github, on Academic License
% 

%% Compute Laplacian matrix and normalize
degrees = sum(AdjacencyMatrix, 2);
D = sparse(1:size(AdjacencyMatrix, 1), 1:size(AdjacencyMatrix,2), degrees);
L = D - AdjacencyMatrix; 
D = spdiags(1./(degrees.^0.5), 0, size(D, 1), size(D, 2)); 
L = D * L * D;


%% Compute nKmeans smallest eigenvectors and perform K-means

diff  = eps;
if nKmeans ~= 0        
    [V, eigenvalues] = eigs(L, nKmeans, diff);
    V = bsxfun(@rdivide, V, sqrt(sum(V.^2, 2))); %normalize eigenvectors

    C = kmeans(V, nKmeans, 'start', 'cluster') - 1; %n-by-1 matrix containing the cluster number for each data point, minus 1 for the consistant with the true label
    %C = sparse(1:size(D, 1), C, 1); 


elseif nKmeans == 0 % Determine optimal number of clusters by maximizing the average silhouette
    
    kRange = 2:1:20; % nClusters search range

    for N = 1:10 % Do a few trials, as Kmeans result can vary between runs
        for nKmeans = kRange
            
            % Compute U as normal 
            V = eigs(L, nKmeans, diff);
            V = bsxfun(@rdivide, V, sqrt(sum(V.^2, 2)));

            C_vary{N,nKmeans-1} = kmeans(V, nKmeans, 'start', 'cluster') - 1;

            % Compute silhouette metric, which we aim to maximize 
            s{N,nKmeans-1} = silhouette(V,C_vary{N,nKmeans-1});
            s_avg(N,nKmeans-1) = mean(s{N,nKmeans-1});

        end
    end
    
    % Index the trial where average silhouette value is maximized and keep the clustering result corresponding to that trial 
    [maxSavg,whereMaxSavg] = max(s_avg(:)); 
    [whereMax_row, whereMax_col] = ind2sub(size(s_avg),whereMaxSavg);
    C = C_vary{whereMax_row,whereMax_col};

end


end
