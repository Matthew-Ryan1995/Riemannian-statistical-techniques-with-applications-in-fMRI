%% Parameters
dimY = 15;
noise = 0.3;
npairs = 10;
udist = 3;

%% For figure
%X = [0:0.25:1 0:0.25:1; 0 0 0 0 0 1 1 1 1 1  ];
X = rand(2,10);
X = center(X);

%% Synthesized data
npivots = size(X,1); % Number of points except the base point
Yp = zeros(dimY, npivots);
while true
    Yp(:,1) = unitvec(randn(dimY,1));
    Yp(:,2) = addnoise_sphere(Yp(:,1),udist,2);
    Yp(:,3) = addnoise_sphere(Yp(:,1),udist,2);
    V = zeros(dimY,npivots);
    for j = 1:npivots
        V(:,j) = logmap_sphere(Yp(:,1),Yp(:,j+1));
    end
    
    %% Generate Ground Truth
    Y = zeros(dimY, size(X,2));
    issafe = true;
    for i = 1:length(X)
        Vtmp = V*X(:,i);
        if norm(Vtmp) > udist
            issafe = false;
            break
        end
        Y(:,i) = expmap_sphere(Yp(:,1),Vtmp);
    end

    if issafe
        break
    end
end

%% Add noise and make some samples.
Xsample = [];
Ysample = zeros(dimY,size(Y,2)*npairs);
isample = 1;
for i = 1:npairs
    for j = 1:size(Y,2)
        Ysample(:,isample) = addnoise_sphere(Y(:,j),noise);
        isample = isample + 1;
    end
    Xsample =[Xsample X];
end
X = Xsample;
Y = unitvec(Ysample);
