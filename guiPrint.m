function guiPrint(fig, fn) 
%   guiPrint.m 
%       prints figures in standard format, and adds tag to metadata  
%       
%   David Healy
%   July 2016 
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

print(fig, '-dpng', '-r300', [fn, '.png']) ;
% FolderName = 'TEST';                                    %first create a folder in the code's directory, put here the folder name
% x = evalin('base', 'user_file_name');                   %call the saved variable from the file guiFracPaQ2D.m line 600
% newStr = erase(x,".");                                  %delete the points in that character
% finalStr = strcat(newStr,fn);                           %build up the string
% saveas(fig, fullfile(FolderName, [finalStr,'.png']));   %save figures in the folder with the specified filename
% savefig(fig, fullfile(FolderName, [finalStr,'.fig'])); 
% %   add some FracPaQ metadata to the file 
% t = Tiff([fn, '.tif'], 'r+') ; 
% t.setTag('Artist', 'GBPaQ version 1.0') ; 
% t.rewriteDirectory() ; 
% t.close() ; 

end 