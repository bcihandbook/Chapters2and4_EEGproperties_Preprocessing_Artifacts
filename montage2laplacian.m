function laplacian = montage2laplacian(montage, pattern)

% function laplacian = montage2laplacian(montage, pattern)
% Function computing the matrix for small laplacian spatial filtering
%
% Inputs:
% montage: electrode configuration in the form of a matrix containing the 
% electrode index where the electrodes are located  and zeros elsewhere. It
% should contain all electrode indices from 1 to K (electrode number).
%
% pattern: 3x3 matrix (first-order laplacina only) containing 1 where the
% electrode in question should be considered a neighbor of the central 
% electrode (the one in position (2,2)) Default is the standard NSWE 
% cross-neighbor pattern. pattern can also be an alphanumeric with values:
% 'cross', 'X' or 'all', for the three common configurations
%
% Outputs:
%
% laplacian: laplacian filter matrix (to be multiplied with the data matrix 
% to get spatially filtered data). laplacian is a KxK matrix where K the 
% number of electrodes/channels 
%
% Logic of the algorithm: after identifying the neighbors of each channel,
% in the laplacian matrix L, channel i will occupy column i, where element
% lii is always 1 and row element lji will have value -1/N if channel j is
% a neighbor of i and N is the total number of neighbors of channel i, and
% 0 otherwise.
%
% Algorithm pseudocode:
% 1) Check that montage and neighbor pattern are sane
% 2) Make a list of neighbors for each channel 
% 3) Create a KxK unit matrix as initialization of the laplacian matrix
% 4) For each column i, fill in positions lji (i~=j, since lii=1) with -1/N
% if j is a neighbor of i and N the total numberof neighbors and leave 0
% otherwise

if(nargin < 2)
    disp(['[montage2laplacian] No neighbor pattern provided, defaulting to cross Laplacian.'])
    pattern = [0 1 0; 1 0 1; 0 1 0];
end

if(nargin < 1)
    disp(['[montage2laplacian] No montage provided! Returning empty.']);
    laplacian = [];
    return;
end

if(ischar(pattern))
    switch(pattern)
        case 'cross'
            pattern = [0 1 0; 1 0 1; 0 1 0];
        case 'X'
            pattern = [1 0 1; 0 0 0;1 0 1];
        case 'all'            
            pattern = [1 1 1; 1 0 1;1 1 1];  
        otherwise
            disp(['[montage2laplacian] Unknown neighbor pattern, defaulting to cross Laplacian!']);            
            pattern = [0 1 0; 1 0 1; 0 1 0];            
    end
end

% Check if montage provided in 0/1 format
if(isequal(unique(montage(:)),[0 1]') || isequal(unique(montage(:)),[1]))
    disp(['[montage2laplacian] Warning: montage provided in 0/1 format. Converting assuming row-wise indices']);
    tmpmontage = montage';tmpmontage = tmpmontage(:);
    tmpmontage(tmpmontage==1) = [1:1:sum(tmpmontage(:))];
    montage = reshape(tmpmontage',[size(montage,2) size(montage,1)])'; 
end

% Find number of channels
N = max(montage(:));

% Check for montage sanity
if(~isequal(unique(montage(:)),[0:1:N]') && ~isequal(unique(montage(:)),[1:1:N]'))
    disp(['[montage2laplacian] Warning: channels are missing from the Laplacian montage!']);
end

% Check for pattern sanity
if(~isequal(size(pattern),[3 3]))
    disp(['[montage2laplacian] Neighbor pattern must be a 3x3 matrix!']);
    laplacian = [];
    return;
end

if(~isequal(unique(pattern(:)),[0 1]') && ~isequal(unique(pattern(:)),[1]) )
    disp(['[montage2laplacian] Invalid neighbor pattern! It must be a 3x3 matrix of 1s and 0s']);
    laplacian = [];
    return;
end

% Find set of neighbor coordinates wrt to a channel, given the pattern
NeighborLoc = [];
for x=-1:1
    for y=-1:1
        if(pattern(2+y,2+x)==1)
            NeighborLoc = [NeighborLoc ; [x y]];
        end
    end
end

% Make a list of neighbors for each channel
NoNeighbors = [];
for ch=1:N
    Neighbors{ch} = [];
    % Find location of channel
    [chy chx] = find(montage==ch);
    if(isempty(chx) || isempty(chy))
        continue;
    end
    % Check if possible neighbor locations exist, and have a channel
    for pnb=1:size(NeighborLoc,1)
        if(  (chx+NeighborLoc(pnb,1)>0) && (chx+NeighborLoc(pnb,1)<=size(montage,2)) && ... 
             (chy+NeighborLoc(pnb,2)>0) && (chy+NeighborLoc(pnb,2)<=size(montage,1)) ) 

            if(montage(chy+NeighborLoc(pnb,2),chx+NeighborLoc(pnb,1))~=0)
                    Neighbors{ch} = [Neighbors{ch} montage(chy+NeighborLoc(pnb,2),chx+NeighborLoc(pnb,1))];
            end
        end
    end
    
    if(isempty(Neighbors{ch}))
        NoNeighbors = [NoNeighbors; ch];
    end
end

if(~isempty(NoNeighbors))
    if(length(NoNeighbors)==1)
        disp(['[montage2laplacian] Warning: Channel ' num2str(NoNeighbors)...
        ' has no neighbors, your Laplacian derivation might be suboptimal.']);
    else
        Chnstr = sprintf('%d,',NoNeighbors);Chnstr=Chnstr(1:end-1);
        disp(['[montage2laplacian] Warning: Channels ' Chnstr...
        ' have no neighbors, your Laplacian derivation might be suboptimal.']);        
    end
    
end

% Initialize laplacian as a NxN unit matrix
laplacian = eye(N);
for col=1:N
    if(~isempty(Neighbors{col}))
        laplacian(sort(Neighbors{col},'ascend'),col) = -1.0/length(Neighbors{col});
    end
end