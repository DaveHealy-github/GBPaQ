function GBPaQlength_new(traces, nPixelsPerMetre, northCorrection, xMin, yMin, xMax, yMax, ...
                                nBins, flag_histolength, flag_logloglength, flag_crossplot, ...
                                flag_mle, flag_censor, flag_revY, flag_revX, ...
                                sColour, nlc, nuc)
%   guiFracPaQ2Dlength.m 
%       calculates and plots statistics of line trace segment lengths  
%       
%   David Healy
%   July 2014 
%   ported from FracPaQ September 2021 
%   d.healy@abdn.ac.uk 

%% Copyright
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

maxPossTraceLength = ceil(sqrt(xMax^2 + yMax^2)) ; 

numTraces = length(traces) ; 

if flag_censor 
    
    numTraceSegmentsCensored = 0 ; 
    for i = 1:numTraces

        nSeg = traces(i).nSegments ; 

        %   if the first point of the trace is on the edge, don't use it in the
        %   length stats 
        if traces(i).Node(1).x == 0.0 || traces(i).Node(1).x == xMax 
            traces(i).segmentLength(1:nSeg) = NaN ;
            traces(i).totalLength = NaN ;
        elseif traces(i).Node(1).y == 0.0 || traces(i).Node(1).y == yMax 
            traces(i).segmentLength(1:nSeg) = NaN ;
            traces(i).totalLength = NaN ;
        end ; 

        %   if the last point of the trace is on the edge, don't use it in the
        %   length stats 
        if traces(i).Node(nSeg).x == 0.0 || traces(i).Node(nSeg).x == xMax 
            traces(i).segmentLength(1:nSeg) = NaN ;
            traces(i).totalLength = NaN ;
        elseif traces(i).Node(nSeg).y == 0.0 || traces(i).Node(nSeg).y == yMax 
            traces(i).segmentLength(1:nSeg) = NaN ;
            traces(i).totalLength = NaN ;
        end ; 

        if sum(isnan(traces(i).segmentLength(1:nSeg))) > 0 
            numTraceSegmentsCensored = numTraceSegmentsCensored + 1 ; 
        end ; 

    end ; 

else 
    numTraceSegmentsCensored = 0 ; 
end ; 

%   trace lengths
traceLengths = [ traces.totalLength ] ; 
traceLengths = traceLengths(~isnan(traceLengths)) ; 
minTraceLength = min(traceLengths) ; 
maxTraceLength = max(traceLengths) ; 

%   segment lengths 
traceSegmentLengths = [ traces.segmentLength ] ; 
len1 = length(traceSegmentLengths) ; 
traceSegmentLengths = traceSegmentLengths(~isnan(traceSegmentLengths)) ; 
len2 = length(traceSegmentLengths) ; 
minSegmentLength = min(traceSegmentLengths) ; 
maxSegmentLength = max(traceSegmentLengths) ; 

%   segment angles 
traceSegmentAngles = [ traces.segmentAngle ]' ; 
traceSegmentAngles = traceSegmentAngles(~isnan(traceSegmentLengths)) ; 

%   double the trace angle data over 360 degrees 
traceSegmentAngles2 = [ round(traceSegmentAngles - northCorrection); ...
                 round(traceSegmentAngles - northCorrection) + 180 ] ;
for i = 1:max(size(traceSegmentAngles2))
    if traceSegmentAngles2(i) < 0 
        traceSegmentAngles2(i) = traceSegmentAngles2(i) + 360 ; 
    end ; 
end ; 

disp(' ') ; 
disp('Length stats...') ; 
disp(['Total number of traces found: ', num2str(numTraces)]) ; 
disp(['Total number of segments found: ', num2str(length(traceSegmentLengths))]) ; 

disp(['Number of traces excluded (censored): ', num2str(numTraceSegmentsCensored)]) ; 
disp(['Number of segments excluded (censored): ', num2str(len1-len2)]) ; 

if nPixelsPerMetre > 0 
    disp(['Minimum trace length: ', num2str(minTraceLength, '%8.2f'), ' metres']) ; 
    disp(['Maximum trace length: ', num2str(maxTraceLength, '%8.2f'), ' metres']) ; 
    disp(['Average trace length: ', num2str(mean(traceLengths, 'omitnan'), '%8.2f'), ' metres']) ; 

    disp(['Minimum segment length: ', num2str(minSegmentLength, '%8.2f'), ' metres']) ; 
    disp(['Maximum segment length: ', num2str(maxSegmentLength, '%8.2f'), ' metres']) ; 
    disp(['Average segment length: ', num2str(mean(traceSegmentLengths, 'omitnan'), '%8.2f'), ' metres']) ; 
else 
    disp(['Minimum trace length: ', num2str(minTraceLength, '%8.2f'), ' pixels']) ; 
    disp(['Maximum trace length: ', num2str(maxTraceLength, '%8.2f'), ' pixels']) ; 
    disp(['Average trace length: ', num2str(mean(traceLengths, 'omitnan'), '%8.2f'), ' pixels']) ; 
    
    disp(['Minimum segment length: ', num2str(minSegmentLength, '%8.2f'), ' pixels']) ; 
    disp(['Maximum segment length: ', num2str(maxSegmentLength, '%8.2f'), ' pixels']) ; 
    disp(['Average segment length: ', num2str(mean(traceSegmentLengths, 'omitnan'), '%8.2f'), ' pixels']) ; 
end ; 

%   write the trace lengths to a text file
sorted_traceSegmentLengths = sort(traceSegmentLengths) ; 
fidLength = fopen('GBPaQsegmentlengths.txt', 'wt') ; 
for i = 1:max(size(traceSegmentLengths)) 
    fprintf(fidLength, '%6.2f\n', sorted_traceSegmentLengths(i)) ; 
end ; 
fclose(fidLength) ; 

%   write the trace lengths to a text file
sorted_traceLengths = sort(traceLengths) ; 
fidLength = fopen('GBPaQtracelengths.txt', 'wt') ; 
for i = 1:max(size(traceLengths)) 
    fprintf(fidLength, '%6.2f\n', sorted_traceLengths(i)) ; 
end ; 
fclose(fidLength) ; 

if flag_histolength 

    %   histogram of trace lengths 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    [ nTraceLengths, binTraceLengths ] = hist(traceLengths, ...
                                minTraceLength:(maxTraceLength-minTraceLength)/nBins:maxTraceLength) ; 

    yyaxis left ;                         
    bar(binTraceLengths, nTraceLengths, 1, 'FaceColor', sColour) ;
    hold on ; 
    plot([minTraceLength, minTraceLength], [0, max(nTraceLengths)*1.2], '-r', 'LineWidth', 1) ;  
    plot([maxTraceLength, maxTraceLength], [0, max(nTraceLengths)*1.2], '-r', 'LineWidth', 1) ;  
%     plot([maxPossTraceLength, maxPossTraceLength], [0, max(nTraceLengths)], '--r', 'LineWidth', 1) ;  
    ylim([0 max(nTraceLengths)*1.2]) ; 
    ylabel('Frequency') ; 
    yyaxis right ; 
    bar(binTraceLengths, (nTraceLengths/sum(nTraceLengths))*100, 1, 'FaceColor', sColour) ;
    hold off ; 
    xlim([0 maxTraceLength*1.1]) ; 
    ylim([0 max((nTraceLengths/sum(nTraceLengths))*100)*1.2]) ; 
    if nPixelsPerMetre > 0 
        xlabel('Trace length, metres') ; 
    else 
        xlabel('Trace length, pixels') ; 
    end ; 
    ylabel('Frequency, %') ; 
    axis on square ; 
    box on ; 
    grid on ; 
    title({['Trace lengths, n=', num2str(length(traceLengths))];''}) ; 

    %   save to file 
    guiPrint(f, 'GBPaQ_histotracelength') ; 

    %   histogram of segment lengths 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    [ nSegmentLengths, binSegmentLengths ] = hist(traceSegmentLengths, ...
                                minSegmentLength:(maxSegmentLength-minSegmentLength)/nBins:maxSegmentLength) ; 

    yyaxis left ; 
    bar(binSegmentLengths, nSegmentLengths, 1, 'FaceColor', sColour) ;
    hold on ; 
    plot([minSegmentLength, minSegmentLength], [0, max(nSegmentLengths)*1.2], '-r', 'LineWidth', 1) ;  
    plot([maxSegmentLength, maxSegmentLength], [0, max(nSegmentLengths)*1.2], '-r', 'LineWidth', 1) ;  
%     plot([maxPossTraceLength, maxPossTraceLength], [0, max(nSegmentLengths)], '--r', 'LineWidth', 1) ;  
    ylim([0 max(nSegmentLengths)*1.2]) ; 
    ylabel('Frequency') ; 
    yyaxis right ; 
    bar(binSegmentLengths, (nSegmentLengths/sum(nSegmentLengths))*100, 1, 'FaceColor', sColour) ;
    hold off ; 
    xlim([0 maxSegmentLength*1.1]) ; 
    ylim([0 max((nSegmentLengths/sum(nSegmentLengths))*100)*1.2]) ; 
    if nPixelsPerMetre > 0 
        xlabel('Segment length, metres') ; 
    else 
        xlabel('Segment length, pixels') ; 
    end ; 
    ylabel('Frequency, %') ; 
    axis on square ; 
    box on ; 
    grid on ; 
    title({['Segment lengths, n=', num2str(length(traceSegmentLengths))];''}) ; 

    %   save to file 
    guiPrint(f, 'GBPaQ_histosegmentlength') ; 

end ; 

%   cumulative log-log plot of trace lengths 
if flag_logloglength 

    %   trace lengths density distribution 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    [ nTraceLengths, binTraceLengths ] = hist(traceLengths, ... 
                                 minTraceLength:(maxTraceLength-minTraceLength)/nBins:maxTraceLength) ; 
    loglog(binTraceLengths, nTraceLengths, 's', ...
        'LineWidth', 1, 'Color', sColour, 'MarkerFaceColor', sColour) ;
    hold on ; 
    plot([minTraceLength, minTraceLength], [0.1, max(nTraceLengths)*1.2], '-r', 'LineWidth', 1) ;  
    plot([maxTraceLength, maxTraceLength], [0.1, max(nTraceLengths)*1.2], '-r', 'LineWidth', 1) ;  
    hold off ; 
    set(gca,'XTick', [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4]) ; 
    xlim([minTraceLength*0.9 maxTraceLength*1.1]) ; 
    ylim([0.1 max(nTraceLengths)*1.2]) ; 
    if nPixelsPerMetre > 0 
        xlabel('Trace length, metres') ; 
    else 
        xlabel('Trace length, pixels') ; 
    end ; 
    ylabel('Frequency, n(L)') ; 
    axis on square ; 
    box on ; 
    grid on ; 
    title({['Density distribution of trace lengths, n=', num2str(length(traceLengths))];''}) ; 
    
        %   save to file 
    guiPrint(f, 'GBPaQ_loglogplottracelengthdensity') ; 
    
    %   trace lengths cumulative distribution 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    if nPixelsPerMetre > 0 
        [ nTraceLengths, binTraceLengths ] = hist(traceLengths, minTraceLength:(maxTraceLength-minTraceLength)/1000:maxTraceLength) ; 
        cLengths = cumsum(nTraceLengths) ; 
        loglog((minTraceLength:(maxTraceLength-minTraceLength)/1000:maxTraceLength), ...
               (max(cLengths) - cLengths), 's', 'LineWidth', 1, 'Color', sColour, 'MarkerFaceColor', sColour) ; 
        hold on ; 
        plot([minTraceLength, minTraceLength], [0.1, max(cLengths)*1.2], '-r', 'LineWidth', 1) ; 
        plot([maxTraceLength, maxTraceLength], [0.1, max(cLengths)*1.2], '-r', 'LineWidth', 1) ; 
        plot([maxPossTraceLength, maxPossTraceLength], [0.1, max(cLengths)*1.2], '--r', 'LineWidth', 1) ; 
        hold off ; 
    else 
        [ nTraceLengths, binTraceLengths ] = hist(traceLengths, minTraceLength:maxTraceLength) ; 
        cLengths = cumsum(nTraceLengths) ; 
        loglog((minTraceLength:maxTraceLength), (max(cLengths) - cLengths), 's', ... 
               'LineWidth', 1, 'Color', sColour, 'MarkerFaceColor', sColour) ; 
        hold on ;   
        plot([minTraceLength, minTraceLength], [0.1, max(cLengths)*1.2], '-r', 'LineWidth', 1) ; 
        plot([maxTraceLength, maxTraceLength], [0.1, max(cLengths)*1.2], '-r', 'LineWidth', 1) ; 
        plot([maxPossTraceLength, maxPossTraceLength], [0.1, max(cLengths)*1.2], '--r', 'LineWidth', 1) ; 
        hold off ; 
    end ; 
    set(gca,'XTick', [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4]) ; 
    xlim([minTraceLength*0.9 maxPossTraceLength*1.1]) ; 
    ylim([0.1 max(cLengths)*1.2]) ; 
    if nPixelsPerMetre > 0 
        xlabel('Trace length, metres') ; 
    else 
        xlabel('Trace length, pixels') ; 
    end ; 
    ylabel('Cumulative frequency, C(L)') ; 
    axis on square ; 
    box on ; 
    grid on ; 
    title({['Cumulative distribution of trace lengths, n=', num2str(length(traceLengths))];''}) ; 

    %   save to file 
    guiPrint(f, 'GBPaQ_loglogplottracelength') ; 

    %   segment lengths density distribution 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    [ nSegmentLengths, binSegmentLengths ] = hist(traceSegmentLengths, ... 
                                 minSegmentLength:(maxSegmentLength-minSegmentLength)/nBins:maxSegmentLength) ; 
    loglog(binSegmentLengths, nSegmentLengths, 's', ...
        'LineWidth', 1, 'Color', sColour, 'MarkerFaceColor', sColour) ;
    hold on ; 
    plot([minSegmentLength, minSegmentLength], [0.1, max(nSegmentLengths)*1.2], '-r', 'LineWidth', 1) ;  
    plot([maxSegmentLength, maxSegmentLength], [0.1, max(nSegmentLengths)*1.2], '-r', 'LineWidth', 1) ;  
    hold off ; 
    set(gca,'XTick', [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4]) ; 
    xlim([minSegmentLength*0.9 maxSegmentLength*1.1]) ; 
    ylim([0.1 max(nSegmentLengths)*1.2]) ; 
    if nPixelsPerMetre > 0 
        xlabel('Segment length, metres') ; 
    else 
        xlabel('Segment length, pixels') ; 
    end ; 
    ylabel('Frequency, n(L)') ; 
    axis on square ; 
    box on ; 
    grid on ; 
    title({['Density distribution of segment lengths, n=', num2str(length(traceSegmentLengths))];''}) ; 
    
        %   save to file 
    guiPrint(f, 'GBPaQ_loglogplotsegmentlengthdensity') ; 

    %   segment lengths cumulative distribution 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    if nPixelsPerMetre > 0 
        [ nSegmentLengths, binSegmentLengths ] = hist(traceSegmentLengths, ...
                minSegmentLength:(maxSegmentLength-minSegmentLength)/1000:maxSegmentLength) ; 
        cLengths = cumsum(nSegmentLengths) ; 
        loglog((minSegmentLength:(maxSegmentLength-minSegmentLength)/1000:maxSegmentLength), ...
               (max(cLengths) - cLengths), 's', 'LineWidth', 1, 'Color', sColour, 'MarkerFaceColor', sColour) ; 
        hold on ; 
        plot([minSegmentLength, minSegmentLength], [0.1, max(cLengths)*1.2], '-r', 'LineWidth', 1) ; 
        plot([maxSegmentLength, maxSegmentLength], [0.1, max(cLengths)*1.2], '-r', 'LineWidth', 1) ; 
        plot([maxPossTraceLength, maxPossTraceLength], [0.1, max(cLengths)*1.2], '--r', 'LineWidth', 1) ; 
        hold off ; 
    else 
        [ nSegmentLengths, binSegmentLengths ] = hist(traceSegmentLengths, minSegmentLength:maxSegmentLength) ; 
        cLengths = cumsum(nSegmentLengths) ; 
        loglog((minSegmentLength:maxSegmentLength), (max(cLengths) - cLengths), 's', ... 
               'LineWidth', 1, 'Color', sColour, 'MarkerFaceColor', sColour) ;
        hold on ; 
        plot([minSegmentLength, minSegmentLength], [0.1, max(cLengths)*1.2], '-r', 'LineWidth', 1) ; 
        plot([maxSegmentLength, maxSegmentLength], [0.1, max(cLengths)*1.2], '-r', 'LineWidth', 1) ;
        plot([maxPossTraceLength, maxPossTraceLength], [0.1, max(cLengths)*1.2], '--r', 'LineWidth', 1) ;
        hold off ;        
    end ; 
    set(gca,'XTick', [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4]) ; 
    xlim([minSegmentLength*0.9 maxPossTraceLength*1.1]) ; 
    ylim([0.1 max(cLengths)*1.2]) ; 
    if nPixelsPerMetre > 0 
        xlabel('Segment length, metres') ; 
    else 
        xlabel('Segment length, pixels') ; 
    end ; 
    ylabel('Cumulative frequency, C(L)') ; 
    axis on square ; 
    box on ; 
    grid on ; 
    title({['Cumulative distribution of segment lengths, n=', num2str(length(traceSegmentLengths))];''}) ; 

    %   save to file 
    guiPrint(f, 'GBPaQ_loglogplotsegmentlength') ; 

end ; 

if flag_crossplot 
    
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    hold on ; 
    plot(traceSegmentAngles2, [traceSegmentLengths'; traceSegmentLengths'], 's', ...
                'Color', sColour, 'MarkerFaceColor', sColour, 'MarkerSize', 2) ; 
    plot([0, 360], [minSegmentLength, minSegmentLength], '-r', 'LineWidth', 1) ;  
    plot([0, 360], [maxSegmentLength, maxSegmentLength], '-r', 'LineWidth', 1) ;  
    plot([0, 360], [maxPossTraceLength, maxPossTraceLength], '--r', 'LineWidth', 1) ;  
    hold off ; 
    xlabel('Trace segment angle, degrees') ; 
    if nPixelsPerMetre > 0 
        ylabel('Trace segment length, metres') ; 
    else 
        ylabel('Trace segment length, pixels') ;
    end ; 
    xlim([-10 370]) ; 
    ylim([0 maxSegmentLength*1.1]) ; 
    set(gca,'XTick', 0:60:360) ; 
    axis on square ; 
    box on ; 
    grid on ; 
    title({['Trace segment lengths versus angles, n=', num2str(length(traceSegmentLengths))];''}) ; 

    %   save to file 
    guiPrint(f, 'GBPaQ_crossplotsegmentlengthangle') ; 

end ; 

% Call the functions for Maximum Likelihood (MLE) Statistical Analysis of lengths.
if flag_mle
    
    % Inputs
    % 'x'---> data set
    x = sorted_traceLengths;

    % 'uc'and 'lc'---> Cut off a certain percentage from the beginnig ('uc') 
    % or the end ('lc') of the data set. Default values are 0% cut-off.
    uc = nuc ;
    lc = nlc ;

    disp(' ') ; 
    disp('MLE for trace lengths...') ; 
    % Power Law Statistics
    [alpha,xmin,z,HpercentPL,PpercentPL]=fittingPL(x,uc,lc,sColour);
    if nPixelsPerMetre > 0 
        xlabel('Length, metres') ; 
    else 
        xlabel('Length, pixels') ;
    end ; 
    t = get(gca, 'Title') ; 
    title({[t.String, ', trace lengths'];''}) ; 
    %   save to file 
    guiPrint(gcf, 'GBPaQ_tracelengthsFittingPowerLaw') ;
    
    % Exponential Statistics
    [lambda,xmin2,z2,Hpercentexp,Ppercentexp]=fittingExponential(x,uc,lc,sColour);
    if nPixelsPerMetre > 0 
        xlabel('Length, metres') ; 
    else 
        xlabel('Length, pixels') ;
    end ; 
    t = get(gca, 'Title') ; 
    title({[t.String, ', trace lengths'];''}) ; 
    %   save to file 
    guiPrint(gcf, 'GBPaQ_tracelengthsFittingExponential') ; 
    
    % Log-Normal Statistics
    [mu,sigma,xmin3,z4,HpercentLN,PpercentLN]=fittingLognormal(x,uc,lc,sColour);
    if nPixelsPerMetre > 0 
        xlabel('Length, metres') ; 
    else 
        xlabel('Length, pixels') ;
    end ; 
    t = get(gca, 'Title') ; 
    title({[t.String, ', trace lengths'];''}) ; 
    %   save to file 
    guiPrint(gcf, 'GBPaQ_tracelengthsFittingLogNormal') ;
    
    % Inputs
    % 'x'---> data set
    x = sorted_traceSegmentLengths;

    % 'uc'and 'lc'---> Cut off a certain percentage from the beginnig ('uc') 
    % or the end ('lc') of the data set. Default values are 0% cut-off.
    uc = nuc ;
    lc = nlc ;

    disp(' ') ; 
    disp('MLE for segment lengths...') ; 
    % Power Law Statistics
    [alpha,xmin,z,HpercentPL,PpercentPL]=fittingPL(x,uc,lc,sColour);
    if nPixelsPerMetre > 0 
        xlabel('Length, metres') ; 
    else 
        xlabel('Length, pixels') ;
    end ; 
    t = get(gca, 'Title') ; 
    title({[t.String, ', segment lengths'];''}) ; 
    %   save to file 
    guiPrint(gcf, 'GBPaQ_segmentlengthsFittingPowerLaw') ;

    % Exponential Statistics
    [lambda,xmin2,z2,Hpercentexp,Ppercentexp]=fittingExponential(x,uc,lc,sColour);
    if nPixelsPerMetre > 0 
        xlabel('Length, metres') ; 
    else 
        xlabel('Length, pixels') ;
    end ; 
    t = get(gca, 'Title') ; 
    title({[t.String, ', segment lengths'];''}) ; 
    %   save to file 
    guiPrint(gcf, 'GBPaQ_segmentlengthsFittingExponential') ; 

    % Log-Normal Statistics
    [mu,sigma,xmin3,z4,HpercentLN,PpercentLN]=fittingLognormal(x,uc,lc,sColour);
    if nPixelsPerMetre > 0 
        xlabel('Length, metres') ; 
    else 
        xlabel('Length, pixels') ;
    end ; 
    t = get(gca, 'Title') ; 
    title({[t.String, ', segment lengths'];''}) ; 
    %   save to file 
    guiPrint(gcf, 'GBPaQ_segmentlengthsFittingLogNormal') ;

end ; 

end 