
clc
clear 

FolderName = 'TEST';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
FigNameUser = ["Segan_tst";"2";"3";"4";"5";"6"];
for iFig = 2:length(FigList)
  FigHandle = FigList(iFig);
  Name= FigNameUser(iFig);
  savefig(FigHandle, fullfile(FolderName, Name+'.fig'));
  saveas(FigHandle, fullfile(FolderName, Name+'.png'));
end

