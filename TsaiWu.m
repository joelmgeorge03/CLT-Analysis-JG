function [TsaiWuAnalysis] = TsaiWu(laminateAnalysis, principleStrength)
    syms P

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
    sigma1 = stressArray(:,3).*P;
    sigma2 = stressArray(:,4).*P;
    tau12 = stressArray(:,5).*P;

    % Defining Tsai-Wu Failure Criterion
    TsaiWuCriterion = sym(zeros(length(layerNum), 1));
    TsaiWuLimitLoadsSym = sym(zeros(length(layerNum), 2));

    % Loop to Compute the Tsai-Wu Criterion for each layer interface
    for i =1:length(layerNum)
        % Plugging in the Tsai-Wu Coefficients and the Layer Stresses
        TsaiWuCriterion(i) = F1.*(sigma1(i)) + F2.*(sigma2(i)) + F11.*(sigma1(i).^2) + F22.*(sigma2(i).^2) + F66.*(tau12(i).^2) - (sqrt(F11.*F22)).*sigma1(i).*sigma2(i) - 1;
        % Solving the Tsai-Wu Quadratic Equation for the two values of P
        TsaiWuLimitLoadsSym(i,:) = solve(TsaiWuCriterion(i) == 0, P);
    end
    
    % Storing Tsai-Wu Analysis Results
    TsaiWuLimitLoadsNum = double(TsaiWuLimitLoadsSym);
    TsaiWuAnalysis.LimitLoads = TsaiWuLimitLoadsNum;
    compressiveLoads = TsaiWuAnalysis.LimitLoads(:,1);
    tensileLoads = TsaiWuAnalysis.LimitLoads(:,2);
    
    % Tabular Output for Analysis Results
    TsaiWuTable = table(layerNum, zCoordLayer, compressiveLoads, tensileLoads, 'VariableNames',{'Layer', 'z', '- P', '+ P'});
    disp("----------------- Tsai-Wu Results ----------------")
    disp(" ")
    disp(TsaiWuTable)
    TsaiWuAnalysis.TsaiWuTable = TsaiWuTable;

    % Calculate the limiting values of P
    compressiveFailureLoad = max(compressiveLoads);
    tensileFailureLoad = min(tensileLoads);
    TsaiWuAnalysis.CompressiveFailureLoad = compressiveFailureLoad;
    TsaiWuAnalysis.TensileFailureLoad = tensileFailureLoad;

    % Find the Layer(s) that are driving failure
    [tensileFailureIdx, ~] = find(TsaiWuLimitLoadsNum == tensileFailureLoad);
    [compressiveFailureIdx, ~] = find(TsaiWuLimitLoadsNum == compressiveFailureLoad);
    tensileFailedLayer = layerNum(tensileFailureIdx,1);
    compressiveFailedLayer = layerNum(compressiveFailureIdx,1);

    disp(" The layer(s) that control failure for a tensile load are: ")
    disp(tensileFailedLayer)

    disp(" The layer(s) that control failure for a compressive load are: ")
    disp(compressiveFailedLayer)
    disp("")
    % Storing Failed Layers
    TsaiWuAnalysis.tensileFailedLayer = tensileFailedLayer;
    TsaiWuAnalysis.compressiveFailedLayer = compressiveFailedLayer;

    % Storing the indices of the failed layers in the analysis results
    TsaiWuAnalysis.tensileFailureIdx = tensileFailureIdx;
    TsaiWuAnalysis.compressiveFailureIdx = compressiveFailureIdx;
    
    
end