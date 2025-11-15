function [TsaiWuAnalysis] = TsaiWuv2(laminateAnalysis, principleStrength)
    
    % Extracting Principle Strengths from Input Structure
    Xc = principleStrength.Xc;
    Xt = principleStrength.Xt;
    Yc = principleStrength.Yc;
    Yt = principleStrength.Yt;
    S = principleStrength.S;

    % Computing Tsai-Wu Coefficients
    F1 = (1/Xt) - (1/Xc);
    F2 = (1/Yt) - (1/Yc);
    F11 = 1/(Xt*Xc);
    F22 = 1/(Yt*Yc);
    F66 = 1/(S^2);

    % Storing the Tsai-Wu Coefficients
    TsaiWuAnalysis.F1 = F1;
    TsaiWuAnalysis.F2 = F2;
    TsaiWuAnalysis.F11 = F11;
    TsaiWuAnalysis.F22 = F22;
    TsaiWuAnalysis.F66 = F66;

    % Displaying Tsai-Wu Coefficients
    disp("----------------- Tsai-Wu Results ----------------")
    disp(" ")
    disp("Tsai-Wu Coefficients:")
    disp(['F1 = ', num2str(F1)]);
    disp(['F2 = ', num2str(F2)]);
    disp(['F11 = ', num2str(F11)]);
    disp(['F22 = ', num2str(F22)]);
    disp(['F66 = ', num2str(F66)]);

    % Extracting Layer Stresses from Input Structure
    stressArray = table2array(laminateAnalysis.StressTable123);
    layerNum = stressArray(:,1);
    zCoordLayer = stressArray(:,2);
    sigma1 = stressArray(:,3);
    sigma2 = stressArray(:,4);
    tau12 = stressArray(:,5);
 
    posTsaiWu = zeros(length(layerNum), 1);
    negTsaiWu = zeros(length(layerNum), 1);
    a = NaN(length(layerNum), 1);
    b = NaN(length(layerNum), 1);

    % Loop to Compute the Limiting Loads using Tsai-Wu for each layer interface
    for i = 1:length(layerNum)
        % Computing the coefficients of the quadratic Tsai-Wu Criterion
        a(i) = F11.*sigma1(i).^2 - sqrt(F11.*F22).*sigma1(i).*sigma2(i) + F22.*sigma2(i).^2 + F66.*tau12(i).^2;
        b(i) = F1.*sigma1(i) + F2.*sigma2(i);
        c = - 1;
        % Computing the roots of the quadratic Tsai-Wu Criterion. Derived
        % using Mathematica.
        posTsaiWu(i) = ( -b(i) + sqrt( b(i).^2 - (4.*a(i).*c) ) ) ./ (2*a(i)) ;
        negTsaiWu(i) = ( -b(i) - sqrt( b(i).^2 - (4.*a(i).*c) ) ) ./ (2*a(i)) ;
    end
    
    % Storing the computed limiting loads in the TsaiWuAnalysis structure
    TsaiWuAnalysis.tensileLoads = posTsaiWu;
    TsaiWuAnalysis.compressiveLoads = negTsaiWu;
    
    % Tabular Output for Analysis Results
    TsaiWuTable = table(layerNum, zCoordLayer, negTsaiWu, posTsaiWu , 'VariableNames',{'Layer', 'z', '- N', '+ N'});
    disp("----------------- Tsai-Wu Table ----------------")
    disp(" ")
    disp(TsaiWuTable)
    TsaiWuAnalysis.TsaiWuTable = TsaiWuTable;

    % Determining the lowest positive and highest negative values of N
    tensileFailureLoad = min(posTsaiWu);
    compressiveFailureLoad = max(negTsaiWu);
    
    % Storing the failure loads in the TsaiWuAnalysis structure
    TsaiWuAnalysis.tensileFailureLoad = tensileFailureLoad;
    TsaiWuAnalysis.compressiveFailureLoad = compressiveFailureLoad;

    % Finding which layers drive failure.
    tensileRowIdx = find(posTsaiWu == tensileFailureLoad);
    compressiveRowIdx = find(negTsaiWu == compressiveFailureLoad);

    % This layer will be used for computations. 
    tensileRowIdxFirst = find(posTsaiWu == tensileFailureLoad, 1, 'first');
    compressiveRowIdxFirst = find(negTsaiWu == compressiveFailureLoad, 1, 'first');

    % Finding the layer(s) that failed
    tensileFailedLayer = layerNum(tensileRowIdx);
    compressiveFailedLayer = layerNum(compressiveRowIdx);
    
    % Storing the Failed Layer numbers in the output structure
    TsaiWuAnalysis.tensileFailedLayer = tensileFailedLayer;
    TsaiWuAnalysis.compressiveFailedLayer = compressiveFailedLayer;

    disp(" The layer(s) that control failure for a tensile load are: ")
    disp(tensileFailedLayer)

    disp(" The layer(s) that control failure for a compressive load are: ")
    disp(compressiveFailedLayer)
    
    %% Determining Failure Mode for Tsai-Wu
    
    % Tensile Load
    sigma1posN = sigma1(tensileRowIdxFirst).*tensileFailureLoad;
    sigma2posN = sigma2(tensileRowIdxFirst).*tensileFailureLoad;
    tau12posN = tau12(tensileRowIdxFirst).*tensileFailureLoad;

    % Compressive Load
    sigma1negN = sigma1(compressiveRowIdxFirst).*compressiveFailureLoad;
    sigma2negN = sigma2(compressiveRowIdxFirst).*compressiveFailureLoad;
    tau12negN = tau12(compressiveRowIdxFirst).*compressiveFailureLoad;

    % Calculating the sigma 1, sigma2, and tau12 contributions to failure
    sigma1posNcontr = F1.*sigma1posN + F11.*sigma1posN.^2 + 0.5.*(-sqrt(F11.*F22)).*sigma1posN.*sigma2posN; % the (0.5) splits the coupling term contribution
    sigma2posNcontr = F2.*sigma2posN + F22.*sigma2posN.^2 + 0.5.*(-sqrt(F11.*F22)).*sigma1posN.*sigma2posN;
    tau12posNcontr = F66.*tau12posN.^2;
    
    sigma1negNcontr = F1.*sigma1negN + F11.*sigma1negN.^2 + 0.5.*(-sqrt(F11.*F22)).*sigma1negN.*sigma2negN;
    sigma2negNcontr = F2.*sigma2negN + F22.*sigma2negN.^2 + 0.5.*(-sqrt(F11.*F22)).*sigma1negN.*sigma2negN;
    tau12negNcontr = F66.*tau12negN.^2;

    posNcontr = [sigma1posNcontr sigma2posNcontr tau12posNcontr];
    negNcontr = [sigma1negNcontr sigma2negNcontr tau12negNcontr];

    TsaiWuAnalysis.posNcontr = posNcontr;
    TsaiWuAnalysis.negNcontr = negNcontr; 

    TsaiWuAnalysis.sumposNcontr = sum(posNcontr);
    TsaiWuAnalysis.sumnegNcontr = sum(negNcontr);

    % Creating array to store failure data
    failureModeDataArray = [...
        compressiveFailureLoad, tensileFailureLoad;
        layerNum(compressiveRowIdxFirst), layerNum(tensileRowIdxFirst);
        sigma1negN, sigma1posN;
        sigma2negN, sigma2posN;
        tau12negN, tau12posN;
        F1.*sigma1negN, F1.*sigma1posN;
        F11.*sigma1negN.^2, F11.*sigma1posN.^2;
        F2.*sigma2negN, F2.*sigma2posN;
        F22.*sigma2negN.^2, F22.*sigma2posN.^2;
        F66.*tau12negN.^2, F66.*tau12posN.^2;
        (-sqrt(F11.*F22)).*sigma1negN.*sigma2negN, (-sqrt(F11.*F22)).*sigma1posN.*sigma2posN;
        sum(negNcontr),sum(posNcontr)];

    % Creating Table to Displayu Failure Contributions
    TsaiWuAnalysis.failureTable = array2table(failureModeDataArray,...
        'VariableNames',{'Compressive', 'Tensile'}, ... 
        'RowNames', {'P',... 
        'Layer',...
        'σ1',...
        'σ2',... 
        'τ12',...
        'F1σ1',...
        'F11σ1sq',...
        'F2σ2',...
        'F22σ2sq',...
        'F66τ12^2',...
        'Coupling',... 
        'Sum'});
 
    disp("------- Tsai-Wu Failure Mode Table -------")
    format short g 
    disp(TsaiWuAnalysis.failureTable)
    disp(' ')
    disp(['The First-Ply Failure Tensile Load = ', num2str(tensileFailureLoad), ' Newtons/meter. '])
    disp(' ')
    disp(['The First-Ply Failure Compressive Load = ', num2str(compressiveFailureLoad), ' Newtons/meter. '])

    %% Determining failure mode via code will be added in later. 


    
end