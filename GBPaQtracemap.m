function GBPaQtracemap(traces, nPixelsPerMetre, nNorth, xMin, yMin, xMax, yMax, ...
                flag_shownodes, flag_revY, flag_revX, sColour, ...
                fMulticolour, flag_tracemap, ...
                flag_tracesbylength, flag_segmentsbylength, flag_segmentsbystrike) 
%   guiFracPaQ2Dtracemap.m 
%       calculates and plots statistics of line trace segment lengths  
%       
%   David Healy
%   ported from FracPaQ September 2021 
%   July 2014 
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

nTraces = length(traces) ; 
nSegments = sum([traces(:).nSegments]) ; 
nNodes = sum([traces(:).nNodes]) ; 

if flag_tracemap

    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    hold on ; 

    for k = 1:nTraces

        plot( [ traces(k).Node.x ]', [ traces(k).Node.y ]', 'LineWidth', 0.75, 'Color', sColour) ;

    end ; 

    if flag_shownodes 

        for k = 1:nTraces

            for l = 1:traces(k).nSegments

                %   plot segment end points
                plot(traces(k).Segment(l).Point1(1), traces(k).Segment(l).Point1(2), 'ok', 'MarkerSize', 4) ; 
                plot(traces(k).Segment(l).Point2(1), traces(k).Segment(l).Point2(2), 'ok', 'MarkerSize', 4) ; 

                %   plot segment mid points 
                plot(traces(k).Segment(l).midpointX, traces(k).Segment(l).midpointY, 'sr') ;  

            end ; 

            %   plot trace mid points 
            plot(traces(k).midpointX, traces(k).midpointY, 'or') ;  

        end ; 

    end ; 

    hold off ; 
    axis on equal ; 
    box on ; 
    xlim([xMin xMax]) ; 
    ylim([yMin yMax]) ; 
    if flag_revX 
        set(gca, 'XDir', 'reverse') ; 
    end ; 
    if flag_revY 
        set(gca, 'YDir', 'reverse') ; 
    end ; 
    if nPixelsPerMetre > 0 
        xlabel('X, metres') ; 
        ylabel('Y, metres') ; 
    else
        xlabel('X, pixels') ; 
        ylabel('Y, pixels') ; 
    end ; 
    title({['Mapped traces (n=', num2str(nTraces), ...
           '), segments (n=', num2str(nSegments), ...
           ') & nodes (n=', num2str(nNodes), ')'];''}) ; 

    %   save to file 
    guiPrint(f, 'GBPaQ_tracemap') ; 

end ; 

%   multcolour tracemap, if necessary 
if fMulticolour
    
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    nTraces = length(traces) ; 
    nSegments = sum([traces(:).nSegments]) ; 
    nNodes = sum([traces(:).nNodes]) ; 

    hold on ; 
    for k = 1:nTraces
        plot( [ traces(k).Node.x ]', [ traces(k).Node.y ]', 'LineWidth', 0.75, 'Color', traces(k).Colour) ;
    end ; 
    hold off ; 
    axis on equal ; 
    box on ; 
    xlim([xMin xMax]) ; 
    ylim([yMin yMax]) ; 
    if flag_revX 
        set(gca, 'XDir', 'reverse') ; 
    end ; 
    if flag_revY 
        set(gca, 'YDir', 'reverse') ; 
    end ; 
    if nPixelsPerMetre > 0 
        xlabel('X, metres') ; 
        ylabel('Y, metres') ; 
    else
        xlabel('X, pixels') ; 
        ylabel('Y, pixels') ; 
    end ; 
    title({['Multicolour traces (n=', num2str(nTraces), ...
           '), segments (n=', num2str(nSegments), ...
           ') & nodes (n=', num2str(nNodes), ')'];''}) ; 

    %   save to file 
    guiPrint(f, 'GBPaQ_tracemap_multicolour') ; 

end ; 

if flag_tracesbylength 
    %   trace length by colour map 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    maxtraceLength = max([traces.totalLength]) ; 
    traceColourMap = cmocean('haline') ; 
    colormap(traceColourMap);

    hold on ; 
    for k = 1:nTraces
        iC = round( ( traces(k).totalLength / maxtraceLength ) * 256 ) ; 
        if iC < 1 || isnan(iC) 
           iC = 1 ; 
        end ; 
        plot( [ traces(k).Node.x ]', [ traces(k).Node.y ]', 'LineWidth', 0.75, 'color', traceColourMap(iC,:) ) ;
    end ; 
    hold off ;
    caxis([0 maxtraceLength]) ; 
    c = colorbar('southoutside') ; 
    axis on equal ; 
    box on ; 
    xlim([xMin xMax]) ; 
    ylim([yMin yMax]) ; 
    if flag_revX 
        set(gca, 'XDir', 'reverse') ; 
    end ; 
    if flag_revY 
        set(gca, 'YDir', 'reverse') ; 
    end ; 
    if nPixelsPerMetre > 0 
        xlabel('X, metres') ; 
        ylabel('Y, metres') ; 
        c.Label.String = 'Trace length, metres' ;  
    else
        xlabel('X, pixels') ; 
        ylabel('Y, pixels') ; 
        c.Label.String = 'Trace length, pixels' ;  
    end ; 
    title({['Trace length map, n=', num2str(nTraces)];''}) ; 

    %   save to file 
    guiPrint(f, 'GBPaQ_tracelengthmap') ; 
end ; 

if flag_segmentsbylength 
    %   segment length by colour map 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    maxsegmentLength = max([traces.segmentLength]) ; 
    segmentColourMap = cmocean('haline') ; 
    colormap(segmentColourMap);

    hold on ; 
    for k = 1:nTraces
        for l = 1:traces(k).nSegments
            iC = round( ( traces(k).segmentLength(l) / maxsegmentLength ) * 256 ) ; 
            if iC < 1 || isnan(iC)
                iC = 1 ; 
            end ; 
            plot( [ traces(k).Segment(l).Point1(1), traces(k).Segment(l).Point2(1) ]', ...
                  [ traces(k).Segment(l).Point1(2), traces(k).Segment(l).Point2(2) ]', ...
                        'LineWidth', 0.75, 'color', segmentColourMap(iC,:) ) ;
        end ; 
    end ; 
    hold off ;
    caxis([0 maxsegmentLength]) ; 
    c = colorbar('southoutside') ; 
    axis on equal ; 
    box on ; 
    xlim([xMin xMax]) ; 
    ylim([yMin yMax]) ; 
    if flag_revX 
        set(gca, 'XDir', 'reverse') ; 
    end ; 
    if flag_revY 
        set(gca, 'YDir', 'reverse') ; 
    end ; 
    if nPixelsPerMetre > 0 
        xlabel('X, metres') ; 
        ylabel('Y, metres') ; 
        c.Label.String = 'Segment length, metres' ;  
    else
        xlabel('X, pixels') ; 
        ylabel('Y, pixels') ; 
        c.Label.String = 'Segment length, pixels' ;  
    end ; 
    title({['Segment length map, n=', num2str(nSegments)];''}) ; 

    %   save to file 
    guiPrint(f, 'GBPaQ_segmentlengthmap') ; 
end ; 

if flag_segmentsbystrike
    %   segment strike by colour map 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

    segmentColours = cmocean('phase', 180) ; 

    hold on ; 
    for k = 1:nTraces

        for l = 1:traces(k).nSegments

            if traces(k).segmentLength(l) > 0

                iAngle = round(traces(k).segmentAngle(l) - nNorth) ; 

                if flag_revY 
                    iAngle = 180 - iAngle ; 
                end ; 

                if flag_revX
                    iAngle = 180 - iAngle ; 
                end ; 

                if iAngle < 0  
                    iAngle = iAngle + 360 ; 
                end ; 

                if iAngle > 180  
                    iAngle = iAngle - 180 ; 
                    if iAngle > 180  
                        iAngle = iAngle - 180 ; 
                    end ; 
                end ; 

                if iAngle < 1  
                    iAngle = 180 ; 
                end ; 

                plot( [ traces(k).Segment(l).Point1(1), traces(k).Segment(l).Point2(1) ]', ...
                      [ traces(k).Segment(l).Point1(2), traces(k).Segment(l).Point2(2) ]', ...
                            'LineWidth', 0.75, 'color', segmentColours(iAngle,:) ) ;

            end ; 

        end ;

    end ; 
    hold off ;
    colormap(segmentColours) ; 
    caxis([0 180]) ; 
    c = colorbar('southoutside') ; 
    c.Label.String = 'Segment strike, degrees' ;  
    axis on equal ; 
    box on ; 
    xlim([xMin xMax]) ; 
    ylim([yMin yMax]) ; 
    if flag_revX 
        set(gca, 'XDir', 'reverse') ; 
    end ; 
    if flag_revY 
        set(gca, 'YDir', 'reverse') ; 
    end ; 
    if nPixelsPerMetre > 0 
        xlabel('X, metres') ; 
        ylabel('Y, metres') ; 
    else
        xlabel('X, pixels') ; 
        ylabel('Y, pixels') ; 
    end ; 
    title({['Segment strike map, n=', num2str(nSegments)];''}) ; 

    %   save to file 
    guiPrint(f, 'GBPaQ_segmentstrikemap') ; 
end ; 

end 