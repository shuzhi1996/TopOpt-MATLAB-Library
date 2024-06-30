%%%%%%%% RENDER MULTI-MATERIAL TOPOLOGY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Render(x,D)
Name   = {'Void' 'A' 'B' 'C'};      TColor = ['w','b','r','k'];
[nely,nelx] = size(x);
% MEDIAN MATERIAL DENSITY
MedianD = D;
for i=1 : length(MedianD)-1
    MedianD(i)=0.5*(D(i)+D(i+1));
end
% DEFINITION OF IMAGE SIZE AND DATA WRITING
ES=3;   % ELEMENT SIZE ES
LS=1;   % LINE SIZE LS
m1 = nely*ES + (nely + 1)*LS;
n1 = nelx*ES + (nelx + 1)*LS;
Image=zeros(m1,n1,3);
for i=1:nely
    for j=1:nelx
        for k=1:length(MedianD)
            if x(i,j)<=MedianD(k)
                RGB=TextToRGB(TColor(k));
                break;
            end
        end
        for k=i*LS+1+(i-1)*ES:i*LS+1+i*ES
            for l=j*LS+1+(j-1)*ES:j*LS+1+j*ES
               Image(k,l,1)=RGB(1);
               Image(k,l,2)=RGB(2);
               Image(k,l,3)=RGB(3);
            end
        end
    end
end
% CREATE LEGENDS
% h=40;rs=10;dis=80;len=80; % LEGEND SIZE CONTROL
% Title=201*ones(h,n1,3);     % GRAY LEGEND BACKGROUND
MergeImage=[Image]./255.*255;
imshow(MergeImage);
% for i=1:length(MedianD)
%     text(len+rs+6+(i-1)*dis,0.5*h,Name(i),'color',TColor(i),'Fontsize',rs*0.8);
%     rectangle('Position',[len+(i-1)*dis,(h-rs)/2,rs,rs],'FaceColor',TColor(i));
% end
end
%%%%%%%% TRANSFER TEXT INTO RGB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RGB]=TextToRGB(TColor)
if TColor=='w'      % WHITE
    RGB(1)=255;   RGB(2)=255;   RGB(3)=255;
elseif TColor=='b'  % BLUE
    RGB(1)=0;   RGB(2)=255;   RGB(3)=255;
elseif TColor=='k'  % Green
    RGB(1)=0;   RGB(2)=204;   RGB(3)=0;
elseif TColor=='r'  % RED
    RGB(1)=255;   RGB(2)=0;   RGB(3)=0;
elseif TColor=='g'  % RED
    RGB(1)=100;   RGB(2)=100;   RGB(3)=100;
elseif TColor=='c'  % RED
    RGB(1)=255;   RGB(2)=255;   RGB(3)=0;
end
end