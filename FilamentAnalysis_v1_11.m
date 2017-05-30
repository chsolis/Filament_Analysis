%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
%%%FILAMENT LENGTH ANALYSIS%%%
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% FILENAME: FilamentAnalysis_v1_9.m, First Created on July 02, 2012
%
% DESCRIPTION: This script allows to determine filament leght distributions
% of regulated actins (Actin + Tn + Tm) and non-regulated actins.
%
% CHANGES TO VERSION: INcluding curvefitting of binding data base on
% intensity trajectories
%
%
% INPUT FILES: .tif images. It is necessary to have images for red, 
% green, and blue channels.3
% 
%
% OUTPUT FILES
%
% INSTRUCTIONS:
%   Set USER VARIABLES (filename, pathname) using the GUI. 
%   Upload red, green, and blue channel images in .tif format.
%   Its necesasary to have the same number of images for each color.
%
% Tested with Matlab R2011b
% Christopher Solis, South Dakota State University, updated January 2012
%
% BEFORE WE BEGIN

% Set constants
clear all;
close all;
clc;
CURR_DIRECT = pwd; %Scripts directory

% BEGIN
prompt = {'Enter pixels per microns:'};
dlg_title = 'Constants';
num_lines = 1;
def = {'15'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
ppum = str2double(answer);

if isequal(ppum,0)
 errordlg({'Input must an interger'},'modal');
end

ImageDataCell = {'N','rAcLength','nrAcLength','AcLength','Y1_TnTm','Y2_TnTm','Y3_TnTm','CorrTnTm','CorrTnAc','CorrTmAc','Y_IntRed','Y_IntGreen','rAc_CurvIdx','nrAc_CurvIdx','AcT_CurvIdx','X_Length','Ac_Length','rAc_Length','nrAc_Length','rAcLength','nrAcLength','AcLength'};

%% 1. Upload images (R G B)
% 
% Get names of image files to be imported

% Select Red images
DiagR = helpdlg('Load Red (Tm) channel images');
uiwait(DiagR);
[FileNameR,PathNameR] = uigetfile({'*.tif;*.tiff','TIFF(*.tif,*.tiff)'},...
    'SELECT RED CHANNEL IMAGES','MultiSelect','on'); 

cd(PathNameR); PathNameR
% Select Green images
DiagG = helpdlg('Load Green (Tn) channel images');
uiwait(DiagG);
[FileNameG,PathNameG] = uigetfile({'*.tif;*.tiff','TIFF(*.tif,*.tiff)'},...
    'SELECT GREEN CHANNEL IMAGES','MultiSelect','on'); 
%
% Select Blue images
DiagB = helpdlg('Load Blue (Actin) channel images');
uiwait(DiagB);
[FileNameB,PathNameB] = uigetfile({'*.tif;*.tiff','TIFF(*.tif,*.tiff)'},...
    'SELECT BLUE CHANNEL IMAGES','MultiSelect','on'); 
%
if isequal(FileNameB,0)
   errordlg({'Cannot continue without specifiying files',...
       'Files must be in TIFFF format'},'Software Error','modal');      
end

if not(length(FileNameG) == length(FileNameR))
   errordlg({'Image color series must have same length'},'modal');
end

%% 2. Catalog all the Images
% Import images into data structure:
cd(PathNameR)
for k = 1:length(FileNameR);
Red(:,:,k) = imread(char(FileNameR(k)));
end
%figure, imshow(Red(:,:,2));
cd(PathNameR)
for k = 1:length(FileNameG);
Green(:,:,k) = imread(char(FileNameG(k)));
end
%figure, imshow(Green(:,:,2));
cd(PathNameR)
for k = 1:length(FileNameB);
Blue(:,:,k) = imread(char(FileNameB(k)));
end
%figure, imshow(Blue(:,:,2));

%% 3. Median filter
% (change to Redu8, Greenu8 or Redu8 to use uint8 images)
MeaFilt = fspecial('average', [3 3]);

for k = 1:length(FileNameR);
Red(:,:,k) = medfilt2(Red(:,:,k),[3 3]);
Red(:,:,k) = imfilter(Red(:,:,k),MeaFilt);
end
%figure, imshow(Redb(:,:,6));
for k = 1:length(FileNameG);
Green(:,:,k) = medfilt2(Green(:,:,k),[3 3]);
Green(:,:,k) = imfilter(Green(:,:,k),MeaFilt);
end
%figure, imshow(Greenb(:,:,6));
for k = 1:length(FileNameB);
Blue(:,:,k) = medfilt2(Blue(:,:,k),[3 3]);
Blue(:,:,k) = imfilter(Blue(:,:,k),MeaFilt);
end
%figure, imshow(Blueb(:,:,6));

%% 4. Pixel Index
cd(CURR_DIRECT)
SLOPE = 0.2;

 for k = 1:length(FileNameR); % TmA equals regulated actins
rAcMask(:,:,k) = mask1plusmask2_2(Red(:,:,k),Green(:,:,k),SLOPE,SLOPE); 
 end

  for k = 1:length(FileNameB); % TmA equals regulated actins
AcTMask(:,:,k) = getMask2(Blue(:,:,k),SLOPE); 
  end

for k = 1:length(FileNameB) % this is to find the non regualted actins 
nrAcMask(:,:,k) = mask1minusmask2and3(Blue(:,:,k),Red(:,:,k),Green(:,:,k),SLOPE,SLOPE,SLOPE);
end

%% 5. Pixel Listing
for k = 1:length(FileNameR) % for regulated actins
rAcMask_k = rAcMask(:,:,k);
for p=1:rAcMask_k.NumObjects
PixListrAc{p,k} = rAcMask_k.PixelIdxList{p};
end 
end

for k = 1:length(FileNameB) % for total actins
AcTMask_k = AcTMask(:,:,k);
for p=1:AcTMask_k.NumObjects
PixListAcTMask{p,k} = AcTMask_k.PixelIdxList{p};
end 
end

for k = 1:length(FileNameB)% for non-regulated actins
nrAcMask_k = nrAcMask(:,:,k);
for p=1:nrAcMask_k.NumObjects
PixListnrAc{p,k} = nrAcMask_k.PixelIdxList{p};
end 
end

%% 6. Reduced Masks
MIN_VAL = 2; % filaments less or equal are not considered

ReducedListAcT = PixListAcTMask; 
for q = 1:length(ReducedListAcT(:,1))
        for w = 1:length(ReducedListAcT(1,:))
      if  length(ReducedListAcT{q,w}) <= MIN_VAL     
     ReducedListAcT{q,w} = [];
      else
      end
        end
end

ReducedListrAc = PixListrAc; 
for q = 1:length(ReducedListrAc(:,1))
        for w = 1:length(ReducedListrAc(1,:))
      if  length(ReducedListrAc{q,w}) <= MIN_VAL     
     ReducedListrAc{q,w} = [];
      else
      end
        end
end

ReducedListnrAc = PixListnrAc; 
for q = 1:length(ReducedListnrAc(:,1))
        for w = 1:length(ReducedListnrAc(1,:))
      if  length(ReducedListnrAc{q,w}) <= MIN_VAL     
     ReducedListnrAc{q,w} = [];
      else
      end
        end
end

%% 7. Check Images

IMAGE = 1; % select image to display
MAX_NUM = uint16(2^16-1);
color = logical([1 0 0]); % red!

% Adjust image
R = Red(:,:,IMAGE);
G = Green(:,:,IMAGE);
B = Blue(:,:,IMAGE);

LIM_R = stretchlim(R);
LIM_G = stretchlim(G);
LIM_B = stretchlim(B);
GammaR = 8;
GammaG = 6;
GammaB = 2;

R = round(imadjust(R,LIM_R,[],GammaR));
G = round(imadjust(G,LIM_G,[],GammaG));
B = round(imadjust(B,LIM_B,[],GammaB));
myRGB = cat(3,R,G,B); 
%TotalFilPix = cat(1,ReducedListAcT{:,IMAGE});
RegulatedPix = cat(1,ReducedListrAc{:,IMAGE});
nonRegulatedPix = cat(1,ReducedListnrAc{:,IMAGE});
pr = imposeMask2(myRGB,RegulatedPix,MAX_NUM,color); % regulated filaments
nrpr = imposeMask2(myRGB,nonRegulatedPix,MAX_NUM,color); %non-regulated filaments
%Acpr = imposeMask2(myRGB,TotalFilPix,MAX_NUM,color); %non-regulated filaments
%figure, imshow(myRGB); % shows image without mask!
figure, imshow(pr);
figure, imshow(nrpr);
%figure, imshow(Acpr);

%% 8. Check Image Approach No3
YES = 0;
IMAGE = 1;
if YES
MAX_NUM2 = uint8(255/1);
color2 = logical([1 0 0]);
mask1 = edge(Blue(:,:,IMAGE),'canny',0.15,2);
TotalPix = cat(1,ReducedListrAc{:,IMAGE});
myRGB3 = double(cat(3,mask1,mask1,mask1));
haha = imposeMask2(myRGB3,TotalPix,MAX_NUM2,color2);
figure, imshow(haha);
end

%% 9. Curvefitting Length and area
BINSIZE = 0.2;
NoBINS = BINSIZE:BINSIZE:10;
HrAcy_tot = [];
HAcTy_tot = [];
HnrAcy_tot = [];
LrAc_um3_tot = [];
LAcT_um3_tot = [];
LnrAc_um3_tot = [];

startingVals = [1,100];

for u = 1:length(FileNameR)

ImrAc = ReducedListrAc(:,u); % obtain number of pixels
for n = 1:length(ImrAc);
   LrAc(:,n) =length(ImrAc{n});
end
   LrAc_um = (LrAc./ppum)';
   LrAc_um2 = LrAc_um > 0;% remove fillaments with zero value length
   LrAc_um3 = LrAc_um(LrAc_um2);
   
ImAcT = ReducedListAcT(:,u);
for n = 1:length(ImAcT);
   LAcT(:,n) =length(ImAcT{n});
end
   LAcT_um = (LAcT./ppum)';
   LAcT_um2 = LAcT_um > 0; % remove fillaments with zero value length
   LAcT_um3 = LAcT_um(LAcT_um2);
   
ImnrAc = ReducedListnrAc(:,u);
for n = 1:length(ImnrAc);
   LnrAc(:,n) =length(ImnrAc{n});
end
   LnrAc_um = (LnrAc./ppum)';
   LnrAc_um2 = LnrAc_um > 0; % remove fillaments with zero value length
   LnrAc_um3 = LnrAc_um(LnrAc_um2);   
   
% accumulation of filament's length
LrAc_um3_tot(u) = sum(LrAc_um3);
LAcT_um3_tot(u) = sum(LAcT_um3);
LnrAc_um3_tot(u) = sum(LnrAc_um3);
   
  
% histograms
[HrAcy,HrAcx]=hist(LrAc_um3,NoBINS); %figure, hist(LrAc_um3,NoBINS);
[HAcTy,HAcTx]=hist(LAcT_um3,NoBINS); %figure, hist(LAcT_um3,NoBINS);
[HnrAcy,HnrAcx]=hist(LnrAc_um3,NoBINS); %figure, hist(LnrAc_um3,NoBINS);

HrAcy_tot(u,:) = HrAcy;
HAcTy_tot(u,:) = HAcTy;
HnrAcy_tot(u,:) = HnrAcy;

ImageDataCell{u+1,1} = u; % Image Number
end
 
%% Combiend histograms
HrAcy_sum = sum(HrAcy_tot);
HAcTy_sum = sum(HAcTy_tot);
HnrAcy_sum = sum(HnrAcy_tot);

% Combined lengths 
LrAc_um3_sum = sum(LrAc_um3_tot);
LAcT_um3_sum = sum(LAcT_um3_tot);
LnrAc_um3_sum =sum(LnrAc_um3_tot);

% plotting distributions
[FR_rAc,GOF_rAc,OUT_rAc] = fit(HrAcx',HrAcy_sum',fittype( 'exp1' ));
[FR_AcT,GOF_AcT,OUT_AcT] = fit(HAcTx',HAcTy_sum',fittype( 'exp1' ));
[FR_nrAc,GOF_nrAc,OUT_nrAc] = fit(HnrAcx',HnrAcy_sum',fittype( 'exp1' ));

% getting the fit datapoints
coefrAc_vals = FR_rAc(HrAcx);
coefrAcT_vals = FR_AcT(HAcTx);
coefnrAc_vals = FR_nrAc(HnrAcx);

coefrAc_resid = OUT_rAc.residuals;
coefAcT_resid = OUT_AcT.residuals;
coefnrAc_resid = OUT_nrAc.residuals;

% Coefficients
    rAc_Ca = FR_rAc.a;
    rAc_Cb = FR_rAc.b;
    AcT_Ca = FR_AcT.a;
    AcT_Cb = FR_AcT.b;  
    nrAc_Ca = FR_nrAc.a;
    nrAc_Cb = FR_nrAc.b;

% Error of the coefficients
FErr_rAc = confint(FR_rAc,0.95); 
FErr_AcT = confint(FR_AcT,0.95); 
FErr_nrAc = confint(FR_nrAc,0.95); 

    % Area calculations   
    FunrAc =  @(x) rAc_Ca.*exp(rAc_Cb.*x);%+ rAc_Cc.*exp(rAc_Cd.*x);
    FunAcT =  @(x) AcT_Ca.*exp(AcT_Cb.*x);% + AcT_Cc.*exp(AcT_Cd.*x);
    FunnrAc =  @(x) nrAc_Ca.*exp(nrAc_Cb.*x);% + nrAc_Cc.*exp(nrAc_Cd.*x); 
    ArAc = quad(FunrAc,BINSIZE,max(HrAcx));
    AAcT = quad(FunAcT,BINSIZE,max(HAcTx));
    AnrAc = quad(FunnrAc,BINSIZE,max(HnrAcx));


% Calculation of TnTm saturation or rAc ratio

% Approach 1: Using the sum of the filament lengths
YrAc_L = LrAc_um3_tot./LAcT_um3_tot;
YnrAc_L = LnrAc_um3_tot./LAcT_um3_tot;

% Approach 2: The areas of rAc relative to the total filaments
YrAc_A1 = ArAc/(AAcT);
YnrAc_A1 = AnrAc/(AAcT);

% Approach 3: The areas of rAc relative to the sume of rAc and nrAc
YrAc_A2 = ArAc/(ArAc+AnrAc);
YnrAc_A2 = AnrAc/(ArAc+AnrAc);

% Combination of 3 approachces
YrAc = geomean([YrAc_L,YrAc_A1,YrAc_A2]); % geometric mean
YnrAc = geomean([YnrAc_L,YnrAc_A1,YnrAc_A2]); % geometric mean

% Length calculations
    Lfit_rAc = abs(1/rAc_Cb);
    Lfit_AcT = abs(1/AcT_Cb);
    Lfit_nrAc = abs(1/nrAc_Cb);

% Error Length Calculations
LErrfit_rAc = abs(1/(rAc_Cb-abs(FErr_rAc(1,2))));
LErrfit_AcT = abs(1/(AcT_Cb-abs(FErr_AcT(1,2))));
LErrfit_nrAc = abs(1/(nrAc_Cb-abs(FErr_nrAc(1,2))));

for u = 1:length(FileNameR)
ImageDataCell{u+1,2} = Lfit_rAc;
ImageDataCell{u+1,3} = Lfit_nrAc;
ImageDataCell{u+1,4} = Lfit_AcT;
ImageDataCell{u+1,5} = YrAc_L;
ImageDataCell{u+1,6} = YrAc_A1;
ImageDataCell{u+1,7} = YrAc_A2;
end
% 10. Plot curves
X_MAXLIM = 5;
%Y_MAXLIM = 70;
Y_MAXLIM = max([0 max([HrAcy_sum HnrAcy_sum HAcTy_sum])]);  
RES_Y_LIM = 28;
figure,
subplot(3,2,1:4);  
    stairs(HrAcx,HrAcy_sum,'r');
    hold on
  stairs(HnrAcx,HnrAcy_sum,'b'); 
  stairs(HAcTx,HAcTy_sum,'m'); 
  plot(HrAcx,coefrAc_vals, '--r','Linewidth',2);
  plot(HnrAcx,coefnrAc_vals, '--b','Linewidth',2);
  plot(HAcTx,coefrAcT_vals, '--m','Linewidth',2);
    set(gca,'XTick',[]);    
    hold off
% set(gca, 'FontName', 'Helvetica');
% set(gca, 'FontSize', 18);
ylabel('Counts')
axis([0 X_MAXLIM 0 Y_MAXLIM])
subplot(3,2, 5:6)
plot(HrAcx,coefrAc_resid,'-r');
hold on
plot(HAcTx,coefAcT_resid,'-m');
plot(HnrAcx,coefnrAc_resid,'-b');
plot([0 X_MAXLIM],[0 0],'-k');
hold off

axis([0 X_MAXLIM -RES_Y_LIM RES_Y_LIM])
ylabel('Residuals');
xlabel('Filament length (\mum)');
daspect([0.023 1 1]) % Adjust the first value to fit both graphs: this is critical
text(2.5,155,[sprintf('<L>_A_c= %0.01f',Lfit_AcT),'\pm',num2str(LErrfit_AcT,'%0.01f'),' \muM'], 'FontName', 'Helvetica','FontSize',14,'BackgroundColor',[1 0.6 1]);
text(2.5,135,[sprintf('<L>_r_A_c= %0.01f',Lfit_rAc),'\pm',num2str(LErrfit_rAc,'%0.01f'),' \muM'], 'FontName', 'Helvetica','FontSize',14,'BackgroundColor',[1 0.6 0.6]);
text(2.5,115,[sprintf('<L>_n_r_A_c= %0.01f',Lfit_nrAc),'\pm',num2str(LErrfit_nrAc,'%0.01f'),' \muM'], 'FontName', 'Helvetica','FontSize',14,'BackgroundColor',[0.6 0.6 1]);
text(2.5,95,sprintf('Y_T_n_T_m= %0.2f',YrAc), 'FontName', 'Helvetica','FontSize',14,'BackgroundColor',[0.5 0.5 0.6]);

ImageDataCell{2,16} = HAcTx;
ImageDataCell{2,17} = HAcTy_sum;
ImageDataCell{2,18} = HrAcy_sum;
ImageDataCell{2,19} = HnrAcy_sum;

%% 11. Correlation data 

for u = 1:length(FileNameR);

OrderedPix = cat(1,ReducedListAcT{:,u});

RedI = spreadpointint(Red(:,:,u),OrderedPix);
GreenI = spreadpointint(Green(:,:,u),OrderedPix);
BlueI = spreadpointint(Blue(:,:,u),OrderedPix);

GreenI = GreenI - min(GreenI);  % yeah weird if I don't divide by the min 
RedI = RedI - min(RedI);        %there in not much relative change of one 
BlueI = BlueI - min(BlueI);     % intensity to the other an the correlation 
                                %is higher.                               

GreenI_T{u} = GreenI'; % Collection of all green intensities                                
RedI_T{u} = RedI';     % Collection of all red intensities 
BlueI_T{u} = BlueI';   % Collection of all red intensities 

GreenxRed = max(abs(xcorr(GreenI,RedI)))/(norm(GreenI,2)*norm(RedI,2));
GreenxBlue = max(abs(xcorr(GreenI,BlueI)))/(norm(GreenI,2)*norm(BlueI,2));
RedxBlue = max(abs(xcorr(RedI,BlueI)))/(norm(RedI,2)*norm(BlueI,2));

ImageDataCell{u+1,8} = GreenxRed;
ImageDataCell{u+1,9} = GreenxBlue;
ImageDataCell{u+1,10} = RedxBlue;
end

%% 12. Binding amount from intensity trajectories
No_BINS = 50;
%IMAGE = 3;
GAUSS = 1; %choose 1 for 2-gaussian curvefit or 0 for 3 gaussian
HisTot_RedCounts = [];
HisTot_GreenCounts = [];
HisTot_BlueCounts = [];
Red_bin_N = [];
Green_bin_N = [];
Blue_bin_N = [];

for IMAGE = 1:length(FileNameR);

Red_Ii = RedI_T{IMAGE};
Green_Ii = GreenI_T{IMAGE};
Blue_Ii = BlueI_T{IMAGE};

Red_Iinorm = Red_Ii/max(Red_Ii);
Green_Iinorm = Green_Ii/max(Green_Ii);
Blue_Iinorm = Blue_Ii/max(Blue_Ii);

[Red_counts,Red_bin] = hist(Red_Iinorm,No_BINS);
[Green_counts,Green_bin] = hist(Green_Iinorm,No_BINS);
[Blue_counts,Blue_bin] = hist(Blue_Iinorm,No_BINS);

Iinorm_List{IMAGE,1} = Red_Iinorm;
Iinorm_List{IMAGE,2} = Green_Iinorm;
Iinorm_List{IMAGE,3} = Blue_Iinorm;

HisTot_RedCounts(:,IMAGE) = Red_counts';
HisTot_GreenCounts(:,IMAGE) = Green_counts';
HisTot_BlueCounts(:,IMAGE) = Blue_counts';
%end


OptionsG1 = fitoptions('gauss1');
OptionsG1.Lower = [0 0 0]; 
OptionsG1.Upper = [10e5 1 1]; 
OptionsG1.Robust = 'off'; 
OptionsG2 = fitoptions('gauss2');
OptionsG2.Lower = [0 0 0 0 0.15 0]; 
OptionsG2.Upper = [10e5 0.3 0.15 10e5 1 1]; 
OptionsG2.Robust = 'off'; 
OptionsG3 = fitoptions('gauss3');
OptionsG3.Lower = [0 0 0 0 0.15 0 0 0 0 ];
OptionsG3.Upper = [10e5 0.5 1 10e5 0.15 1 10e5 1 1 ];
OptionsG3.Robust = 'off'; 

if GAUSS
FIT_RED = fit(Red_bin',HisTot_RedCounts(:,IMAGE),'gauss2',OptionsG2);
FIT_GREEN = fit(Green_bin',HisTot_GreenCounts(:,IMAGE),'gauss2',OptionsG2);
FIT_BLUE = fit(Blue_bin',HisTot_BlueCounts(:,IMAGE),'gauss1',OptionsG1);
else
FIT_RED = fit(Red_bin',HisTot_RedCounts(:,IMAGE),'gauss3',OptionsG3);
FIT_GREEN = fit(Green_bin',HisTot_GreenCounts(:,IMAGE),'gauss3',OptionsG3);
FIT_BLUE = fit(Blue_bin',HisTot_BlueCounts(:,IMAGE),'gauss1',OptionsG1);
end

Coef_RED = coeffvalues(FIT_RED);
Coef_GREEN = coeffvalues(FIT_GREEN);
Coef_BLUE = coeffvalues(FIT_BLUE);

if GAUSS
% separate the 2 gussian equations
EQN1_RED = @(x) Coef_RED(1)*exp(-((x-Coef_RED(2))/Coef_RED(3)).^2); 
EQN2_RED = @(x) Coef_RED(4)*exp(-((x-Coef_RED(5))/Coef_RED(6)).^2); 

% get the area ratios of the 2nd gaussian divided by the 2-gaussian area
A_EQN1_RED = quad(EQN1_RED,0,max(Red_bin));
A_EQN2_RED = quad(EQN2_RED,0,max(Red_bin));

% Condition to divide the area of bigger mean relative to the total
if Coef_RED(6) > Coef_RED(3)
A_ratio_RED = A_EQN2_RED./(A_EQN1_RED+A_EQN2_RED);
else
A_ratio_RED = A_EQN1_RED./(A_EQN1_RED+A_EQN2_RED);  
end

% separate the 2 gussian equations
EQN1_GREEN = @(x) Coef_GREEN(1)*exp(-((x-Coef_GREEN(2))/Coef_GREEN(3)).^2); 
EQN2_GREEN = @(x) Coef_GREEN(4)*exp(-((x-Coef_GREEN(5))/Coef_GREEN(6)).^2); 

% get the area ratios of the 2nd gaussian divided by the 2-gaussian area
A_EQN1_GREEN = quad(EQN1_GREEN,0,max(Red_bin));
A_EQN2_GREEN = quad(EQN2_GREEN,0,max(Red_bin));

% Condition to divide the area of bigger mean relative to the total
% use the std(sigma) since one peak is wider
if Coef_GREEN(6) > Coef_GREEN(3)
A_ratio_GREEN = A_EQN2_GREEN./(A_EQN1_GREEN+A_EQN2_GREEN);
else
A_ratio_GREEN = A_EQN1_GREEN./(A_EQN1_GREEN+A_EQN2_GREEN);   
end

else
% separate the 3 gussian equations
EQN1_RED = @(x) Coef_RED(1)*exp(-0.5*((x-Coef_RED(2))/Coef_RED(3)).^2); 
EQN2_RED = @(x) Coef_RED(4)*exp(-0.5*((x-Coef_RED(5))/Coef_RED(6)).^2); 
EQN3_RED = @(x) Coef_RED(7)*exp(-0.5*((x-Coef_RED(8))/Coef_RED(9)).^2);

% separate the 3 gussian equations
EQN1_GREEN = @(x) Coef_GREEN(1)*exp(-0.5*((x-Coef_GREEN(2))/Coef_GREEN(3)).^2); 
EQN2_GREEN = @(x) Coef_GREEN(4)*exp(-0.5*((x-Coef_GREEN(5))/Coef_GREEN(6)).^2); 
EQN3_GREEN = @(x) Coef_GREEN(7)*exp(-0.5*((x-Coef_GREEN(8))/Coef_GREEN(9)).^2);

% get the area ratios of the 3rd gaussian divided by the 3-gaussian area
A_EQN1_RED = quad(EQN1_RED,0,max(Red_bin));
A_EQN2_RED = quad(EQN2_RED,0,max(Red_bin));
A_EQN3_RED = quad(EQN3_RED,0,max(Red_bin));

% Condition to divide the area of bigger mean relative to the total
if Coef_RED(2) > ( Coef_RED(5) & Coef_RED(8) ) 
A_ratio_RED = A_EQN1_RED./(A_EQN1_RED+A_EQN2_RED+A_EQN3_RED);
elseif Coef_RED(5) > ( Coef_RED(8) & Coef_RED(2) )
A_ratio_RED = A_EQN2_RED./(A_EQN1_RED+A_EQN2_RED+A_EQN3_RED);  
else
A_ratio_RED = A_EQN3_RED./(A_EQN1_RED+A_EQN2_RED+A_EQN3_RED);    
end

% get the area ratios of the 3rd gaussian divided by the 3-gaussian area
A_EQN1_GREEN = quad(EQN1_GREEN,0,max(Red_bin));
A_EQN2_GREEN = quad(EQN2_GREEN,0,max(Red_bin));
A_EQN3_GREEN = quad(EQN3_GREEN,0,max(Red_bin));

% Condition to divide the area of bigger mean relative to the total
if Coef_GREEN(2) > ( Coef_GREEN(5) & Coef_GREEN(8) )
A_ratio_GREEN = A_EQN1_GREEN./(A_EQN1_GREEN+A_EQN2_GREEN+A_EQN3_GREEN);
elseif Coef_GREEN(5) > ( Coef_GREEN(8) & Coef_GREEN(2) )
A_ratio_GREEN = A_EQN2_GREEN./(A_EQN1_GREEN+A_EQN2_GREEN+A_EQN3_GREEN);  
else
A_ratio_GREEN = A_EQN3_GREEN./(A_EQN1_GREEN+A_EQN2_GREEN+A_EQN3_GREEN);    
end

end % from the gauss...

% get the datapoints of the curvefitting
DatapointsFIT_Red(:,IMAGE) = feval(FIT_RED,Red_bin'); % Red is 1
DatapointsFIT_Green(:,IMAGE) = feval(FIT_GREEN,Green_bin'); % Green is 2
DatapointsFIT_Blue(:,IMAGE) = feval(FIT_BLUE,Blue_bin'); % Blue is 3

% Normalized histograms

DatapointsFIT_Red_N(:,IMAGE) = DatapointsFIT_Red(:,IMAGE)/max(DatapointsFIT_Red(:,IMAGE));
DatapointsFIT_Green_N(:,IMAGE) = DatapointsFIT_Green(:,IMAGE)/max(DatapointsFIT_Green(:,IMAGE));
DatapointsFIT_Blue_N(:,IMAGE) = DatapointsFIT_Blue(:,IMAGE)/max(DatapointsFIT_Blue(:,IMAGE));
Hist_Counts_Red_N(:,IMAGE) = HisTot_RedCounts(:,IMAGE)./max(DatapointsFIT_Red(:,IMAGE));
Hist_Counts_Green_N(:,IMAGE) = HisTot_GreenCounts(:,IMAGE)./max(DatapointsFIT_Green(:,IMAGE));
Hist_Counts_Blue_N(:,IMAGE) = HisTot_BlueCounts(:,IMAGE)./max(DatapointsFIT_Blue(:,IMAGE));

Red_bin_N(:,IMAGE) = Red_bin/max(Red_bin);
Green_bin_N(:,IMAGE) = Green_bin/max(Green_bin);
Blue_bin_N(:,IMAGE) = Blue_bin/max(Blue_bin);
end

% Save in matrix
%ImageDataCell{1+1,11} = A_ratio_RED; % RED
%mageDataCell{1+1,12} = A_ratio_GREEN; % GREEN
%end

IMAGE = 2;
%     %%  to check each gaussian curve  %
%plot(Red_bin,EQN1_GREEN(Red_bin),'r');
%hold on
%plot(Red_bin,EQN2_GREEN(Red_bin),'g');
%plot(EQN3_GREEN(Red_bin),'b');

%     %%  to check curvefitting  %
figure,stairs(Red_bin_N(:,IMAGE),Hist_Counts_Red_N(:,IMAGE),'r');
hold on
plot(Red_bin_N(:,IMAGE),DatapointsFIT_Red_N(:,IMAGE),'k'); %ylim([0 100])
legend off
figure,stairs(Green_bin_N(:,IMAGE),Hist_Counts_Green_N(:,IMAGE),'g');
hold on
plot(Green_bin_N(:,IMAGE),DatapointsFIT_Green_N(:,IMAGE),'k'); %ylim([0 100]);
legend off

%figure,stairs(Blue_bin,Hist_Counts_Blue,'b');
%hold on
%plot(FIT_BLUE,'k'); %ylim([0 100])
%legend off
%%
%Gridss = ReducedListAcT(:,3);

%%
%for i = iFIL:fFIL
   
%    step = length(Gridss{3,:})
%    Lines = [i,1:0.1:10];
%end

%% 13. Intensity plots 
YES = 1;
if YES
IMAGE = 1;
iFIL = 1; % initial filament
fFIL = 600; % final filament

x_axis = (1:length(Iinorm_List{IMAGE,1}))/ppum; % all colors measure =
RedIN = Iinorm_List{IMAGE,1};
GreenIN = Iinorm_List{IMAGE,2};
BlueIN = Iinorm_List{IMAGE,3};
Hmax_plot = max([HisTot_RedCounts(:,IMAGE);HisTot_GreenCounts(:,IMAGE);HisTot_BlueCounts(:,IMAGE)]);

figure, 
% Intensity trarjectories 
subplot(3,4,[1,2,3]);
plot(x_axis(iFIL:fFIL),RedIN(iFIL:fFIL),'r'); 
set(gca,'XTickLabel',{''}); ylim([0 1]); ylabel({'Tm norm.';'intensity'}); 
axis fill;

subplot(3,4,[5,6,7]);
plot(x_axis(iFIL:fFIL),GreenIN(iFIL:fFIL),'g');
set(gca,'XTickLabel',{''}); ylim([0 1]); ylabel({'Tn norm.';'intensity'}); 
axis fill;

subplot(3,4,[9,10,11]);
plot(x_axis(iFIL:fFIL),BlueIN(iFIL:fFIL),'b');
xlabel('Filament coordinate (\mum)'); ylabel({'Ac norm.';'intensity'});
ylim([0 1]); axis fill;

% histograms
subplot(3,4,[4]);
barh(Red_bin,Hist_Counts_Red_N(:,IMAGE),'r','EdgeColor','r');
hold on
plot(DatapointsFIT_Red_N(:,IMAGE),Red_bin_N,'-k'); %xlim([0 Hmax_plot]);
set(gca,'XTickLabel',{''}); set(gca,'YTickLabel',{''}); axis fill;
xlim([0 1]);
text(0.1*Hmax_plot,0.85,['Y=',num2str(ImageDataCell{1+1,11},'%1.2f')],'FontSize',18);
hold off

subplot(3,4,[8]);
barh(Green_bin,Hist_Counts_Green_N(:,IMAGE),'g','EdgeColor','g');
hold on
plot(DatapointsFIT_Green_N(:,IMAGE),Green_bin_N,'-k'); %xlim([0 Hmax_plot]);
set(gca,'XTickLabel',{''}); set(gca,'YTickLabel',{''}); axis fill;
xlim([0 1]);
text(0.1*Hmax_plot,0.85,['Y=',num2str(ImageDataCell{1+1,12},'%1.2f')],'FontSize',18);
hold off

subplot(3,4,[12]);
barh(Blue_bin,Hist_Counts_Blue_N(:,IMAGE),'b','EdgeColor','b');
hold on
plot(DatapointsFIT_Blue_N(:,IMAGE),Blue_bin,'-k'); %xlim([0 Hmax_plot]);
xlim([0 1]);
xlabel('Norm. counts'); set(gca,'YTickLabel',{''}); axis fill;
hold off
end

%% 14. Curvature

rAc_curvature = [];
mean_rAc_curvature = [];
for u = 1:length(FileNameR)
for f = 1:length(ReducedListrAc(:,u));
FILrAc = [];
rAc_curvature_u = [];
rAc_curvature_u2 = [];

FILrAc = ReducedListrAc{f,u};
if numel(FILrAc) > 0;% removes small filaments
IND = [FILrAc(1);FILrAc(end)];
[I,J] = ind2sub(size(Red(:,:,1)),IND);
Distance = (((I(2)-I(1)))^2 + ((J(2)-J(1)))^2)^0.5;
Length = length(FILrAc);
rAc_curvature(f,u) = Length./Distance;
else
end
end
rAc_curvature_u = rAc_curvature(:,u) >0;
rAc_curvature_u2 = rAc_curvature(rAc_curvature_u,u);
mean_rAc_curvature = mean(rAc_curvature_u2);
ImageDataCell{u+1,13} = mean_rAc_curvature';
end

nrAc_curvature = [];
mean_nrAc_curvature = [];
for u = 1:length(FileNameR)
for f = 1:length(ReducedListnrAc(:,u));
FILnrAc = [];
nrAc_curvature_u = [];
nrAc_curvature_u2 = [];

FILnrAc = ReducedListnrAc{f,u};
if numel(FILnrAc) > 0;% removes small filaments
IND = [FILnrAc(1);FILnrAc(end)];
[I,J] = ind2sub(size(Red(:,:,1)),IND);
Distance = (((I(2)-I(1)))^2 + ((J(2)-J(1)))^2)^0.5;
Length = length(FILnrAc);
nrAc_curvature(f,u) = Length./Distance;
else
end
end
nrAc_curvature_u = nrAc_curvature(:,u) >0;
nrAc_curvature_u2 = nrAc_curvature(nrAc_curvature_u,u);
mean_nrAc_curvature = mean(nrAc_curvature_u2);
ImageDataCell{u+1,14} = mean_nrAc_curvature';
end

AcT_curvature = [];
mean_AcT_curvature = [];
for u = 1:length(FileNameR)
for f = 1:length(ReducedListAcT(:,u));
FILAcT = [];
AcT_curvature_u = [];
AcT_curvature_u2 = [];

FILAcT = ReducedListAcT{f,u};
if numel(FILAcT) > 0;% removes small filaments
IND = [FILAcT(1);FILAcT(end)];
[I,J] = ind2sub(size(Red(:,:,1)),IND);
Distance = (((I(2)-I(1)))^2 + ((J(2)-J(1)))^2)^0.5;
Length = length(FILAcT);
AcT_curvature(f,u) = Length./Distance;
else
end
end
AcT_curvature_u = AcT_curvature(:,u) >0;
AcT_curvature_u2 = AcT_curvature(AcT_curvature_u,u);
mean_AcT_curvature = mean(AcT_curvature_u2);
ImageDataCell{u+1,15} = mean_AcT_curvature';
end


CURV_rAC = ImageDataCell(:,13);
CURV_nrAC = ImageDataCell(:,14);
CURV_AcT = ImageDataCell(:,15);

figure, boxplot([cell2mat(CURV_rAC(2:length(FileNameR)+1)),cell2mat(CURV_nrAC(2:length(FileNameR)+1)),cell2mat(CURV_AcT(2:length(FileNameR)+1))],'labels',{'rAc','nrAc','AcT'});
ylabel('Curvature Index');
%ylim([1 1.6]);
%% 15. Saving
uisave('ImageDataCell');
