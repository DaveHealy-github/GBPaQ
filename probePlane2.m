function probePlane2(x, xMin, yMin, xMax, yMax, traces, inc, angles, lengths)
%   probePlane.m
%
%   Uses stereology methods on 3 orthogonal 2D images to define 
%   3D fabric and crack tensors  
%
%   references:
%       Suzuki et al., 1998 IJRMMS v35 
%       Oda et al., 2002 JGR v107
%
%   David Healy
%   August 2012 
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

nq = 180 ; 

%   run scan traces over the image, centred on the image mid point
%   at nq different angles 
xMid = (xMax-xMin)/2 ; 
yMid = (yMax-yMin)/2 ; 
nLines = 1 ; 
if nLines == 1 
    disp(' ') ; 
    disp('*** Using radial scanlines...') ; 
    lenScanline = 0.95 * min(xMax-xMin, yMax-yMin) ; 
else       
%     disp(' ') ; 
%     disp('*** Using parallel scanlines...') ; 
    lenScanline = 0.95 * sqrt((xMax-xMin)^2 + (yMax-yMin)^2) ; 
end ; 
d = lenScanline / nLines ; 
id = -10:1:10 ; 
id = id' ; 
% mq = zeros(nq, 1) ; 
mq_new = zeros(nq, 1) ; 

traces = traces(~any(isnan(traces),2),:) ; 

i = 0 ; 
for a = 0:inc:nq 

    i = i + 1 ; 
    angleScan = a * pi / 180 ; 
    sina = sin(angleScan) ; 
    cosa = cos(angleScan) ; 
%         disp(['Angle of this scan line = ', num2str(angleScan*180/pi), ' deg']) ; 

    if nLines > 1  
        %   define parallel scan lines 
        deltaX = cosa * d ; 
        deltaY = sina * d ; 

        x1Scan = xMid - ( id * deltaX ) + ( cosa * lenScanline ) / 2 ; 
        y1Scan = yMid + ( id * deltaY ) + ( sina * lenScanline ) / 2 ; 

        x2Scan = xMid - ( id * deltaX ) - ( cosa * lenScanline ) / 2 ; 
        y2Scan = yMid + ( id * deltaY ) - ( sina * lenScanline ) / 2 ; 
    else 
        %   define radial scan line: length, angle and end points 
        x1Scan = xMid + ( cosa * lenScanline ) / 2 ;  
        y1Scan = yMid + ( sina * lenScanline ) / 2 ; 
        x2Scan = xMid - ( cosa * lenScanline ) / 2 ; 
        y2Scan = yMid - ( sina * lenScanline ) / 2 ; 
    end ; 

    %   find which traces intersected by scan line - their angles and lengths  
    scanXY = [ x1Scan, y1Scan, x2Scan, y2Scan ] ; 

    intersectionScan = lineSegmentIntersect(scanXY, traces) ; 

    %   count intersections of scan line with traces 
    numIntersections = max(size(find(intersectionScan.intMatrixX > 0))) ; 
    ni(i) = numIntersections ; 
%     if ni(i) > 10 
%         disp(traces(intersectionScan.intAdjacencyMatrix, :)) ; 
%     end 

    %   normalise count by total scan line length 
    mq(i) = numIntersections ./ ( nLines * lenScanline ) ; 
%         disp(['Count of intersections = ', num2str(numIntersections)]) ;
%         disp(['Normalised count of intersections = ', num2str(mq(a))]) ;

end ; 
    
thetaScan = deg2rad([0:inc:360]) ; 
ni = [ni(1:end-1), ni(1:end-1), ni(1)] ; 
mq = [mq(1:end-1), mq(1:end-1), mq(1)] ; 
ni_mean = zeros(1, max(size(ni))) ; 
ni_mean(:) = mean(ni) ; 
mq_mean = zeros(1, max(size(mq))) ; 
mq_mean(:) = mean(mq) ; 

[ ni_min, ni_min_ind ] = min(ni) ; 
[ ni_max, ni_max_ind ] = max(ni) ; 
[ mq_min, mq_min_ind ] = min(mq) ; 
[ mq_max, mq_max_ind ] = max(mq) ; 
ni_min_theta = thetaScan(ni_min_ind) ; 
ni_max_theta = thetaScan(ni_max_ind) ; 
mq_min_theta = thetaScan(mq_min_ind) ; 
mq_max_theta = thetaScan(mq_max_ind) ; 

ni_min_thetad = 90 - rad2deg(ni_min_theta) ;
if ni_min_thetad < 0 
    ni_min_thetad = 360 + ni_min_thetad ;
end 
ni_max_thetad = 90 - rad2deg(ni_max_theta) ;
if ni_max_thetad < 0
    ni_max_thetad = 360 + ni_max_thetad ;
end 
mq_min_thetad = 90 - rad2deg(mq_min_theta) ;
if mq_min_thetad < 0
    mq_min_thetad = 360 + mq_min_thetad ;
end 
mq_max_thetad = 90 - rad2deg(mq_max_theta) ; 
if mq_max_thetad < 0 
    mq_max_thetad = 360 + mq_max_thetad ; 
end 

ia = 0 ; 
for a = [ni_min_theta, ni_max_theta]  

    ia = ia + 1 ; 
    sina = sin(a) ; 
    cosa = cos(a) ; 

    %   define radial scan line: length, angle and end points 
    x1Scan = xMid + ( cosa * lenScanline ) / 2 ;  
    y1Scan = yMid + ( sina * lenScanline ) / 2 ; 
    x2Scan = xMid - ( cosa * lenScanline ) / 2 ; 
    y2Scan = yMid - ( sina * lenScanline ) / 2 ; 

    %   find which traces intersected by scan line - their angles and lengths  
    scanXY = [ x1Scan, y1Scan, x2Scan, y2Scan ] ; 
    intersectionScan = lineSegmentIntersect(scanXY, traces) ; 

    %   get angles of intersected segments
    numIntersections = max(size(find(intersectionScan.intMatrixX > 0))) ; 
    if numIntersections > 0 
        
        segAngles = angles(intersectionScan.intAdjacencyMatrix, :) ; 
        segLengths = lengths(intersectionScan.intAdjacencyMatrix, :) ; 
        
        segmentAngles = [ segAngles ] ; 
        segmentLengths = [ segLengths ] ; 
        
        %   double the trace angle data over 360 degrees 
        segmentAngles2 = [ round(segmentAngles) ; round(segmentAngles) + 180 ] ;
        for i = 1:max(size(segmentAngles2))
            if segmentAngles2(i) < 0 
                segmentAngles2(i) = segmentAngles2(i) + 360 ; 
            end
        end
        %   double the length data too, for length weighting of rose plot 
        segmentLengths2 = [ segmentLengths ; segmentLengths ] ; 

%         if flag_revX 
%             segmentAngles2 = 180 - segmentAngles2 ; 
%             for i = 1:max(size(segmentAngles2))
%                 if segmentAngles2(i) < 0 
%                     segmentAngles2(i) = segmentAngles2(i) + 360 ; 
%                 end
%             end
%         end 
% 
%         if flag_revY 
%             segmentAngles2 = 180 - segmentAngles2 ; 
%             for i = 1:max(size(segmentAngles2))
%                 if segmentAngles2(i) < 0 
%                     segmentAngles2(i) = segmentAngles2(i) + 360 ; 
%                 end 
%             end 
%         end

        f = figure ;
        set(gcf, 'PaperPositionMode', 'manual') ; 
        set(gcf, 'PaperUnits', 'inches') ; 
        set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

        roseEqualAreaColour(segmentAngles2, 10, 0, segmentLengths2, false, 'b') ; 

        if ia == 1 
            plot([ni_min*cos(a), -ni_min*cos(a)], ...
                 [ni_min*sin(a), -ni_min*sin(a)], '-g', 'LineWidth', 2) ; 
            title({['Segment angles intersected along minimum direction'] ; ...
                   ['(equal area, length weighted), n=', num2str(length(segLengths))]}) ; 
            %   save to file 
            guiPrint(f, 'GBPaQ_radialintersections_rose_minimum') ; 
        else
            plot([ni_max*cos(a), -ni_max*cos(a)], ...
                 [ni_max*sin(a), -ni_max*sin(a)], '-y', 'LineWidth', 2) ; 
            title({['Segment angles intersected along maximum direction'] ; ...
                   ['(equal area, length weighted), n=', num2str(length(segLengths))]}) ; 
            %   save to file 
            guiPrint(f, 'GBPaQ_radialintersections_rose_maximum') ; 
        end 
        
    end 

end

%   labels need adjusting to 'strike' values (0 at North)
a = [0:30:360] ; 
for i = 1:max(size(a)) 
    if 90 - a(i) >= 0 
        thetaLabels(i) = { num2str(90 - a(i)) } ; 
    else 
        thetaLabels(i) = { num2str(90 - a(i) + 360) } ; 
    end 
end

f = figure ;
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

polarplot(0,0, 'HandleVisibility', 'off') ; 
hold on ; 
polarplot(thetaScan, ni, '-r', 'LineWidth', 2) ; 
polarplot(thetaScan, ni_mean, '-b') ; 
polarplot([ni_min_theta, ni_min_theta+pi ; ni_min_theta, ni_min_theta+pi ], ...
          [0, 0; ni_min, ni_min], '-g', 'LineWidth', 2) ; 
polarplot([ni_max_theta, ni_max_theta+pi ; ni_max_theta, ni_max_theta+pi ], ...
          [0, 0; ni_max, ni_max], '-y', 'LineWidth', 2) ; 
hold off ; 
ax = gca ; 
ax.ThetaTickLabel = thetaLabels ;  
title('Number of intersections') ; 
legend(['Intersections (min=', num2str(ni_min), ' at ', num2str(ni_min_thetad, '%03d'), '\circ, max=', ...
        num2str(ni_max), ' at ', num2str(ni_max_thetad, '%03d'), '\circ)'], ...
       ['Average intersections (mean=', num2str(ni_mean(1)), ')'], ... 
        'Location', 'southoutside') ; 
%   save to file 
guiPrint(f, 'GBPaQ_radialintersections') ; 

f = figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

polarplot(0,0, 'HandleVisibility', 'off') ; 
hold on ; 
polarplot(thetaScan, mq, '-r', 'LineWidth', 2) ; 
polarplot(thetaScan, mq_mean, '-b') ; 
polarplot([mq_min_theta, mq_min_theta+pi ; mq_min_theta, mq_min_theta+pi ], ...
          [0, 0; mq_min, mq_min], '-g', 'LineWidth', 2) ; 
polarplot([mq_max_theta, mq_max_theta+pi ; mq_max_theta, mq_max_theta+pi ], ...
          [0, 0; mq_max, mq_max], '-y', 'LineWidth', 2) ; 
hold off ; 
ax = gca ; 
ax.ThetaTickLabel = thetaLabels ;  
title('Number of intersections per pixel') ; 
legend(['Intersections/pixel (min=', num2str(mq_min), ' at ', num2str(mq_min_thetad, '%03d'), '\circ, max=', ...
        num2str(mq_max), ' at ', num2str(mq_max_thetad, '%03d'), '\circ)'], ...
       ['Average intersections/pixel (mean=', num2str(mq_mean(1)), ')'], ... 
        'Location', 'southoutside') ; 

%   save to file 
guiPrint(f, 'GBPaQ_radialintersections_norm') ; 

%   write data to file 
fidInt = fopen('GBPaQintersections.txt', 'wt') ; 
for i = 1:max(size(ni)) 
    fprintf(fidInt, '%6.2f\t%8.2f\n', rad2deg(thetaScan(i)), ni(i)) ; 
end 
fclose(fidInt) ; 

fidIntper = fopen('GBPaQintersectionsperpixel.txt', 'wt') ; 
for i = 1:max(size(mq)) 
    fprintf(fidIntper, '%6.2f\t%10.4f\n', rad2deg(thetaScan(i)), mq(i)) ; 
end 
fclose(fidIntper) ; 

end 
