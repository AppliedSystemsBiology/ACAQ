/*
Created on Wed April 18 12:16:00 2018

@author: Zoltan Cseresnyes
@affiliation: Research Group Applied Systems Biology, Leibniz Institute for 
Natural Product Research and Infection Biology – Hans Knöll Institute (HKI),
Beutenbergstrasse 11a, 07745 Jena, Germany.
@email: zoltan.cseresnyes@leibniz-hki.de or zcseresn@gmail.com
For advice and rights, please contact via
@email: thilo.figge@leibniz-hki.de

This is an implementation of the Hessian-based algorithm for host and pathogens segmentation.The 
full details of the work can be found in Cseresnyes et al. (2018), "ACAQ: a Fiji and R toolkit targeted for
automated confrontation assay quantification", currently under review in Bioinformatics. 
 If any part of the code is used for academic purposes or publications, please cite the 
above mentioned paper.

Copyright (c) 2016-2018, 
Leibniz Institute for Natural Product Research and Infection Biology – 
Hans Knöll Institute (HKI)

Licence: BSD-3-Clause, see  
https://opensource.org/licenses/BSD-3-Clause for full details

__version__ = "ACAQ Analyzer"

*/

appName = "ACAQ ";
versionNumber="Analyzer";

function dateAndTime(){
	MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
    DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
    getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	TimeString ="Date: "+DayNames[dayOfWeek]+" ";
	if (dayOfMonth<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
	if (hour<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+hour+":";
	if (minute<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+minute+":";
	if (second<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+second;
	return TimeString;
}

function getTag(tag) {
      info = getImageInfo();
      index1 = indexOf(info, tag);
      if (index1==-1) return "";
      index1 = indexOf(info, ":", index1);
      if (index1==-1) return "";
      index2 = indexOf(info, "\n", index1);
      value = substring(info, index1+2, index2);
      return value;
  }

function isElement(element, array) {
      value = false;
      for(i=0; i<array.length;i++){
      	if(element == array[i]){
      		value = true;
      	}
      }
      return value;
  }

function colorThresholding(){
	run("RGB Color");
	min=newArray(3);
	max=newArray(3);
	filter=newArray(3);
	a=getTitle();
	run("RGB Stack");
	run("Stack to Images");
	selectWindow("Red");
	rename("0");
	selectWindow("Green");
	rename("1");
	selectWindow("Blue");
	rename("2");
	min[0]=0;
	max[0]=0;
	filter[0]="pass";
	min[1]=0;
	max[1]=0;
	filter[1]="pass";
	min[2]=0;
	max[2]=0;
	filter[2]="pass";
	for (ii=0;ii<3;ii++){
	  selectWindow(""+ii);
	  setThreshold(min[ii], max[ii]);
	  run("Convert to Mask");
	  if (filter[ii]=="stop")  run("Invert");
	}
	imageCalculator("AND create", "0","1");
	imageCalculator("AND create", "Result of 0","2");
	for (i=0;ii<3;ii++){
	  selectWindow(""+ii);
	  close();
	}
	selectWindow("Result of 0");
	close();
	selectWindow("Result of Result of 0");
	rename(a + "_colorThreshold");
}

function colorThresholdingHSB(){
	run("RGB Color");
	run("HSB Stack");
	min=newArray(3);
	max=newArray(3);
	filter=newArray(3);
	a=getTitle();
	run("Stack to Images");
	selectWindow("Hue");
	rename("0");
	selectWindow("Saturation");
	rename("1");
	selectWindow("Brightness");
	rename("2");
	min[0]=0;
	max[0]=255;
	filter[0]="pass";
	min[1]=0;
	max[1]=255;
	filter[1]="pass";
	min[2]=1;
	max[2]=255;
	filter[2]="pass";
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  setThreshold(min[i], max[i]);
	  run("Convert to Mask");
	  if (filter[i]=="stop")  run("Invert");
	}
	imageCalculator("AND create", "0","1");
	imageCalculator("AND create", "Result of 0","2");
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  close();
	}
	selectWindow("Result of 0");
	close();
	selectWindow("Result of Result of 0");
	rename(a + "_colorThreshold");
}

print("Algorithm for Confrontation Assay Quantification (" + versionNumber + ")");
print("**********************************************");
print(dateAndTime());
print(" ");
run("Clear Results");

//Initialize counters:
CLAHEmessage = " "; 
nROI_greenSpores = 0;
nROI_spores = 0;
nROI_Hough_spores = 0;
nROI_macrophages = 0;
nROI_redMacrophages = 0;
nROI_blueSpores = 0;
nROI_greenSporesInside = 0;
nROI_greenSporesOutside = 0;
nROI_greenSpores_BrightSpots = 0;
numberOfMacrophages = 0;
numberOfInsideGreenSporesInsideMacrophages_excludingEdges = 0;
numberOfAdherentBlueSpores_excludingEdges = 0;
numberOfMacrophagesWithPhagocytosedGreenSpores = 0;
numberOfMacrophages_excludingEdges = 0;
phagocytosisRatio_excludingEdges = 0;
uptakeRatio_excludingEdges = 0;
phagocyticIndex_excludingEdges = 0;
symmetrizedPhagocyticIndex_excludingEdges = 0;
pathogenAreaMean = 0;
pathogenAreaSTD = 0;
pathogenPerimMean = 0;
pathogenPerimSTD = 0;
pathogenARMean = 0;
pathogenARSTD = 0;
pathogenSolidityMean = 0;
pathogenSoliditySTD = 0;
pathogenPerimeterFeretRatioMean = 0;
pathogenPerimeterFeretRatioSTD = 0;
hostAreaMean = 0;
hostAreaSTD = 0;
hostPerimMean = 0;
hostPerimSTD = 0;
hostARMean = 0;
hostARSTD = 0;
hostSolidityMean = 0;
hostSoliditySTD = 0;
hostPerimeterFeretRatioMean = 0;
hostPerimeterFeretRatioSTD = 0;
hostFlscAreaMean = 0;
hostFlscAreaSTD = 0;
hostFlscPerimMean = 0;
hostFlscPerimSTD = 0;
hostFlscARMean = 0;
hostFlscARSTD = 0;
hostFlscSolidityMean = 0;
hostFlscSoliditySTD = 0;
hostFlscPerimeterFeretRatioMean = 0;
hostFlscPerimeterFeretRatioSTD = 0;
//end of counter initialization

//start of Dialog box 0:
Dialog.create(appName + versionNumber + ": Dialog 0");
Dialog.addMessage(appName + versionNumber + ": Algorithm for Confrontation Assay Quantification; Dialog 0: basic settings and preprocessing parameters");
Dialog.addCheckbox("FAST mode? (CLAHE calculations)", true);
Dialog.addCheckbox("FAST mode? (NO plots)", true);
Dialog.addCheckbox("Re-stiching Zeiss tilescan?", false);
Dialog.addCheckbox("Read multi-channel TIFF?", false);
Dialog.addCheckbox("Use MAX Hessian instead? (default: MIN)", false);
Dialog.addCheckbox("Do you want to save your results?", true);
Dialog.addCheckbox("Use subfolders?", false);
Dialog.addCheckbox("Analyse all images?", true);
Dialog.addCheckbox("Close all windows at end?", true);
Dialog.addCheckbox("Run it in Batch Mode?", false);
Dialog.addCheckbox("Illumination correction for TL images?", false);
Dialog.addCheckbox("Illumination correction for fluorescence images?", false);
Dialog.addNumber("Line thickness for ROIs: ",0);
Dialog.addNumber("X tile number:",3);
Dialog.addNumber("Y tile number:",3);
Dialog.addNumber("Gaussian sigma for illumination correction:",21);
Dialog.addNumber("Gaussian sigma for fluorescence illumination correction:",11);
Dialog.addNumber("Hessian smoothing factor for host cells:",2);
Dialog.addNumber("Hessian smoothing factor for pathogens:",2);
Dialog.addNumber("Gaussian smoothing factor for Hessian host cells:",5);
Dialog.addNumber("Gaussian smoothing factor for labelled host cells:",3);
Dialog.addNumber("Gaussian smoothing factor for Hessian pathogens:",2);
Dialog.addNumber("Dilate/erode steps for Hessian host cells:",2);
Dialog.addNumber("Additional erosion steps for Hessian host cells:",3);
Dialog.addNumber("Additional dilation steps for labelled host cells:",3);
Dialog.addNumber("Additional erosion steps for labelled host cells:",5);
Dialog.addNumber("Dilate/erode steps for pathogens:",2);
Dialog.addNumber("Gaussian smoothing factor for fluorescence pathogens:",3);
Dialog.addNumber("Lower threshold multiplier for pathogens (only for Internal Gradient):",1.0);
Dialog.addNumber("Upper threshold multiplier for pathogens (only for Internal Gradient):",1.125);
Dialog.addNumber("Lower threshold multiplier for host cells:",0.6);
Dialog.addNumber("Upper threshold multiplier for host cells:",1.00);
Dialog.show();

fastModeCLAHE=Dialog.getCheckbox();
fastModeNoPlots=Dialog.getCheckbox();
rearrangeZeiss=Dialog.getCheckbox();
useMultichannelTiffFormat=Dialog.getCheckbox();
useMAXHessian=Dialog.getCheckbox();
saveResults=Dialog.getCheckbox();
useSubfolders=Dialog.getCheckbox();
analyseAllImages=Dialog.getCheckbox();
closeAllWindows=Dialog.getCheckbox();
runInBatchMode=Dialog.getCheckbox();
correctIlluminationForTL=Dialog.getCheckbox();
correctIlluminationForRed=Dialog.getCheckbox();
lineThicknessForObjects=Dialog.getNumber();
XTileNumber=Dialog.getNumber();
YTileNumber=Dialog.getNumber();
illuminationCorrectionSigmaForTL=Dialog.getNumber();
illuminationCorrectionSigmaForRed=Dialog.getNumber();
hessianSmoothingForMacrophages=Dialog.getNumber();
hessianSmoothingForSpores=Dialog.getNumber();
gaussianSmoothingForMacrophages=Dialog.getNumber();
gaussianSmoothingForRedMacrophages=Dialog.getNumber();
gaussianSmoothingForSpores=Dialog.getNumber();
dilateErodeStepsForMacrophages=Dialog.getNumber();
additionalErodeStepsForMacrophages=Dialog.getNumber();
enlargeRedMacrophages=Dialog.getNumber();
additionalErodeStepsForLabelledMacrophages=Dialog.getNumber();
dilateErodeStepsForSpores=Dialog.getNumber();
gaussianBlurForSpores=Dialog.getNumber();
lowerThresholdMultiplier=Dialog.getNumber();
upperThresholdMultiplier=Dialog.getNumber();
lowerThresholdMultiplierMacrophages=Dialog.getNumber();
upperThresholdMultiplierMacrophages=Dialog.getNumber();
//end of Dialog box 0

//start of Dialog box 1:
Dialog.create(appName + versionNumber + ": Dialog 1");
Dialog.addMessage(appName + versionNumber + ": Algorithm for Confrontation Assay Quantification ; Dialog 1: image types and naming rules");
Dialog.addCheckbox("Labelled hosts?", true);
Dialog.addCheckbox("Labelled pathogens?", true);
Dialog.addCheckbox("Exclude edges?", true);
Dialog.addCheckbox("Gather results for entire image set?", false);
Dialog.addCheckbox("Use image file name for results?", true);
Dialog.addCheckbox("Use specific image groups based on filename?", false);
Dialog.addCheckbox("Save full path of image name in output file?", true);//allow to have the full image name saved
Dialog.addString("Image type: ", ".czi");
Dialog.addString("Search term #1: ", "");
Dialog.addString("Search term #2: ", "");
Dialog.addString("Search term #3: ", "");
Dialog.addString("Exclude term #1: ", "   ");
Dialog.addString("Exclude term #2: ", "   ");
Dialog.addString("Exclude term #3: ", "   ");
Dialog.addString("Name tag of saved files: ", "");
Dialog.addNumber("Background image number:",1);
Dialog.addNumber("First image number:",1);
Dialog.addNumber("Last image number:",1);
Dialog.addNumber("Test image number:",1);
Dialog.addNumber("Edge width (for Exclude Edges option):",5);
Dialog.addNumber("Labelled pathogen channel multiplier (for dim image data):",1);
Dialog.addNumber("Outside pathogen channel multiplier (for dim image data):",1);
Dialog.addNumber("Labelled host cells channel multiplier (for dim image data):",1);

Dialog.show();

labelledMacrophages=Dialog.getCheckbox();
labelledSpores=Dialog.getCheckbox();
excludeEdges=Dialog.getCheckbox();
gatherResultsForEntireImageSet=Dialog.getCheckbox();
useImageFilename=Dialog.getCheckbox();
useSpecificImageGroups=Dialog.getCheckbox();
saveFullPathImagename=Dialog.getCheckbox();
imageType=Dialog.getString();
searchTerm_1=Dialog.getString();
searchTerm_2=Dialog.getString();
searchTerm_3=Dialog.getString();
excludeTerm_1=Dialog.getString();
excludeTerm_2=Dialog.getString();
excludeTerm_3=Dialog.getString();
nameTagSavedFiles=Dialog.getString();
backgroundImageNumber=Dialog.getNumber();
firstImageNumber=Dialog.getNumber();
lastImageNumber=Dialog.getNumber();
testImageNumber=Dialog.getNumber();
edgeWidth=Dialog.getNumber();
greenChannelMultiplierForDimData=Dialog.getNumber();
blueChannelMultiplierForDimData=Dialog.getNumber();
redChannelMultiplierForDimData=Dialog.getNumber();
//end of Dialog box 1

//Start of Dialog box, #2:
Dialog.create(appName + versionNumber + ": Dialog 2");
Dialog.addMessage(appName + versionNumber + ": Algorithm for Confrontation Assay Quantification ; Dialog 2: segmentation parameters"); 
Dialog.addCheckbox("Watershed on host cells?", true);
Dialog.addCheckbox("Watershed on pathogens?", true);
Dialog.addCheckbox("CLAHE on unlabelled images?", false);
Dialog.addCheckbox("CLAHE on labelled images?", true);
Dialog.addCheckbox("Combine Internal Gradient and Bright Spot results?", false);
Dialog.addNumber("Internal Gradient radius for labelled pathogens (try 2; use 0 for No):", 0);
Dialog.addNumber("Internal Gradient radius for outside pathogens (try 2; use 0 for No):", 0);
Dialog.addNumber("Internal Gradient radius for labelled host cells (try 25; use 0 for No):", 25);
Dialog.addNumber("Rolling ball radius for pathogens background (try 20):",20);
Dialog.addNumber("Rolling ball radius for labelled host cells background (try 20):",10);
Dialog.addNumber("Local threshold radius for labelled host cells image (try b/w 30 and 40): ",30);
Dialog.addNumber("Remove outliers for host cell ROI, step 1 (try 10) : ", 10);
Dialog.addNumber("Remove outliers for host cell ROI, step 2 (try 20) : ", 20);
Dialog.addNumber("CLAHE blocks:",127);
Dialog.addNumber("CLAHE bins:",256);
Dialog.addNumber("CLAHE max slope:",3.0);
Dialog.addNumber("1st CLAHE max slope for labelled host cells:",0.0);
Dialog.addNumber("2nd CLAHE max slope for labelled host cells:",3.0);
Dialog.addNumber("1st dilation steps for labelled host cells:",3);
Dialog.addNumber("2nd dilation steps for labelled host cells:",3);
Dialog.show();

watershedOnMacrophages=Dialog.getCheckbox();
watershedOnSpores=Dialog.getCheckbox();
applyCLAHE=Dialog.getCheckbox();
applyCLAHEonFLSC=Dialog.getCheckbox();
combineInternalGradientAndBrightSpotResults=Dialog.getCheckbox();
internalGradientRadiusGreen=Dialog.getNumber();
internalGradientRadiusBlue=Dialog.getNumber();
internalGradientRadiusRed=Dialog.getNumber();
rollingBallRadius=Dialog.getNumber();
rollingBallRadiusRedMacrophages=Dialog.getNumber();
localThresholdRadius = Dialog.getNumber();
removeOutliersStep1=Dialog.getNumber();
removeOutliersStep2=Dialog.getNumber();
CLAHEblocks=Dialog.getNumber();
CLAHEbins=Dialog.getNumber();
CLAHEslope=Dialog.getNumber();
CLAHEslope1FlscMacrophages=Dialog.getNumber();
CLAHEslope2FlscMacrophages=Dialog.getNumber();
dilationsSteps1FlscMacrophages=Dialog.getNumber();
dilationsSteps2FlscMacrophages=Dialog.getNumber();
//end of Dialog box 2

//Start of Dialog box, #3:
Dialog.create(appName + versionNumber + ": Dialog 3");
Dialog.addMessage(appName + versionNumber + ": Algorithm for Confrontation Assay Quantification ; Dialog 3: size limits and thresholding methods"); 
Dialog.addNumber("Min host cell size:",3000);
Dialog.addNumber("Max host cell size:",30000);
Dialog.addNumber("Min host cell circularity:",0.3);
Dialog.addNumber("Max host cell circularity:",1.00);
Dialog.addNumber("Min pathogen size:",400);
Dialog.addNumber("Max pathogen size:",3000);
Dialog.addNumber("Min pathogen circularity:",0.4);
Dialog.addNumber("Max pathogen circularity:",1.00);
Dialog.addNumber("Outside pathogens threshold for all pathogen classifier (separate inside from outside) :",40);
Dialog.addNumber("Threshold for Bright Spots of all pathogen classifier :",20); 
thresholdMethodListHessian = newArray("Percentile", "Otsu", "Huang", "RenyiEntropy", "Triangle");
Dialog.addRadioButtonGroup("Threshold method for Hessian", thresholdMethodListHessian, 1, 5, "Otsu");
thresholdMethodListGreenFluorescence = newArray("Otsu", "Huang", "RenyiEntropy", "Triangle", "Li");
Dialog.addRadioButtonGroup("Threshold method for all pathogen fluorescence (NOT VALID for BS when combining IG and BS, see next button)", thresholdMethodListGreenFluorescence, 1, 5, "Otsu");
thresholdMethodListGreenFluorescenceBrightSpots = newArray("Otsu", "Huang", "RenyiEntropy", "Li");
Dialog.addRadioButtonGroup("Threshold method for all pathogen fluorescence Bright Spots method (only when combining IG and BS, otherwise see previous button)", thresholdMethodListGreenFluorescenceBrightSpots, 1, 5, "Li");
thresholdMethodListBlueFluorescence = newArray("Percentile", "Otsu", "Huang", "RenyiEntropy", "Triangle");
Dialog.addRadioButtonGroup("Threshold method for outside pathogen fluorescence", thresholdMethodListBlueFluorescence, 1, 5, "Otsu");
localThresholdMethodListRedFluorescence = newArray("Median", "Otsu", "Mean", "Niblack", "Phansalkar");
Dialog.addRadioButtonGroup("Local threshold method for labelled host cells fluorescence", localThresholdMethodListRedFluorescence, 1, 5, "Median");
thresholdMethodListRedFluorescence = newArray("Li", "Otsu", "Huang", "RenyiEntropy", "Triangle");
Dialog.addRadioButtonGroup("Threshold method for labelled host cells fluorescence", thresholdMethodListRedFluorescence, 1, 5, "Li");
Dialog.show();

minMacrophageSize=Dialog.getNumber();
maxMacrophageSize=Dialog.getNumber();
minMacrophageCircularity=Dialog.getNumber();
maxMacrophageCircularity=Dialog.getNumber();
minSporeSize=Dialog.getNumber();
maxSporeSize=Dialog.getNumber();
minSporeCircularity=Dialog.getNumber();
maxSporeCircularity=Dialog.getNumber();
blueThresholdForGreenSporeClassifier=Dialog.getNumber();
greenThresholdForGreenSporeClassifier=Dialog.getNumber();
thresholdMethodHessian = Dialog.getRadioButton;
thresholdMethodGreenFluorescence = Dialog.getRadioButton;
thresholdMethodGreenFluorescenceBrightSpots = Dialog.getRadioButton;
thresholdMethodBlueFluorescence = Dialog.getRadioButton;
localThresholdMethodRedFluorescence = Dialog.getRadioButton;
thresholdMethodRedFluorescence = Dialog.getRadioButton;
//end of Dialog box 3

//Start of Dialog box, #4:
Dialog.create(appName + versionNumber + ": Dialog 4");
Dialog.addMessage(appName + versionNumber + ": Algorithm for Confrontation Assay Quantification ; Dialog 4: Channel assignment"); 
Dialog.addCheckbox("Use GUI for channel info?", true);
//channel0 = newArray("TL", "Blue", "Green", "Red");
channel0 = newArray("Unlabelled hosts/pathogens", "Outside pathogens", "All pathogens", "Labelled hosts");
Dialog.addRadioButtonGroup("Channel 0: ", channel0, 1, 4, "Labelled hosts");
//channel1 = newArray("TL", "Blue", "Green", "Red");
channel1 = newArray("Unlabelled hosts/pathogens", "Outside pathogens", "All pathogens", "Labelled hosts");
Dialog.addRadioButtonGroup("Channel 1: ", channel1, 1, 4, "All pathogens");
//channel2 = newArray("TL", "Blue", "Green", "Red");
channel2 = newArray("Unlabelled hosts/pathogens", "Outside pathogens", "All pathogens all", "Labelled hosts");
Dialog.addRadioButtonGroup("Channel 2: ", channel2, 1, 4, "Outside pathogens");
//channel3 = newArray("TL", "Blue", "Green", "Red");
channel3 = newArray("Unlabelled hosts/pathogens", "Outside pathogens", "All pathogens all", "Labelled hosts");
Dialog.addRadioButtonGroup("Channel 3: ", channel3, 1, 4, "Unlabelled hosts/pathogens");
Dialog.show();

useGUIforChannelInfo=Dialog.getCheckbox();
channel00_new = Dialog.getRadioButton;
channel01_new = Dialog.getRadioButton;
channel02_new = Dialog.getRadioButton;
channel03_new = Dialog.getRadioButton;
print("New channel names from GUI: " + channel00_new + "; " + channel01_new + "; " + channel02_new + "; " + channel03_new + "; " );

if(channel00_new == "Unlabelled hosts/pathogens") channel00 = "TL";
else if(channel00_new == "Outside pathogens") channel00 = "Blue";
else if(channel00_new == "All pathogens") channel00 = "Green";
else channel00 = "Red";
if(channel01_new == "Unlabelled hosts/pathogens") channel01 = "TL";
else if(channel01_new == "Outside pathogens") channel01 = "Blue";
else if(channel01_new == "All pathogens") channel01 = "Green";
else channel01 = "Red";
if(channel02_new == "Unlabelled hosts/pathogens") channel02 = "TL";
else if(channel02_new == "Outside pathogens") channel02 = "Blue";
else if(channel02_new == "All pathogens") channel02 = "Green";
else channel02 = "Red";
if(channel03_new == "Unlabelled hosts/pathogens") channel03 = "TL";
else if(channel03_new == "Outside pathogens") channel03 = "Blue";
else if(channel03_new == "All pathogens") channel03 = "Green";
else channel03 = "Red";
//end of Dialog box 4

//Start of Dialog box, #5:
Dialog.create(appName + versionNumber + ": Dialog 5");
Dialog.addMessage(appName + versionNumber + ": Algorithm for Confrontation Assay Quantification ; Dialog 5: Hessian-based segmentation of unlabeled pathogens"); 
segmentationMethodListForUnlabeledPathogens = newArray("Hough-filter", "Guided Watershed", "Hough and Watershed");
Dialog.addRadioButtonGroup("Choose segmentation method for unlabeled pathogens: ", segmentationMethodListForUnlabeledPathogens, 1, 3, "Hough-filter");
Dialog.addNumber("Hessian smoothing for pathogens:",1.0);
Dialog.addNumber("Internal gradient radius for pathogens:",3);
Dialog.addNumber("Hough circles minimum radius:",7);
Dialog.addNumber("Hough circles maximum radius:",25);
Dialog.addNumber("Hough circles increment radius:",1);
Dialog.addNumber("Hough minimum number of circles:",0);
Dialog.addNumber("Hough maximum number of circles:",700);
Dialog.addNumber("Hough circles threshold:",0.6);
Dialog.addNumber("Hough circles resolution:",113);
Dialog.addNumber("Hough circles ratio:",1.0);
Dialog.addNumber("Hough circles bandwidth:",10);
Dialog.addNumber("Hough circles local radius:",10);
thresholdMethodListHessianSpores = newArray("Default", "Otsu", "Huang", "RenyiEntropy", "Li");
Dialog.addRadioButtonGroup("Threshold method for Hessian", thresholdMethodListHessianSpores, 1, 5, "Li");
Dialog.show();

segmentationMethodForUnlabeledPathogens=Dialog.getRadioButton;
smoothingHessianSpores=Dialog.getNumber();
gradientRadiusHessianSpores=Dialog.getNumber();
minimumRadiusHoughHessianSpores=Dialog.getNumber();
maximumRadiusHoughHessianSpores=Dialog.getNumber();
incrementRadiusHoughHessianSpores=Dialog.getNumber();
minimumNumberHoughHessianSpores=Dialog.getNumber();
maximumNumberHoughHessianSpores=Dialog.getNumber();
thresholdHoughHessianSpores=Dialog.getNumber();
resolutionHoughHessianSpores=Dialog.getNumber();
ratioHoughHessianSpores=Dialog.getNumber();
bandwidthRadiusHoughHessianSpores=Dialog.getNumber();
localradiusHoughHessianSpores=Dialog.getNumber();
thresholdMethodHessianSpores=Dialog.getRadioButton;
//end of Dialog box 5

//Set up Results file where we'll save the main results:
saveAllFile = File.open("");//final results in a columnar format
print(saveAllFile, "Image" + "	\t" + "Nh(total)" + "	\t" + "Np(total)" + "	\t" + "Np(phag)" + "	\t" + "Np(adh)" + "	\t" + "Nh(phag)" + "	\t" + "Fi(pathogen)" + "	\t" + "Fi(host)" + "	\t" + "Fi(i)" + "	\t" + "Fi(i_sym)" + "	\t" 
 	+ "Area_Mean(pathogen)" + "	\t" + "Area_STD(pathogen)" + "	\t" + "Perimeter_Mean(pathogen)" + "	\t" + "Perimeter_STD(pathogen)" + "	\t" + "AR_Mean(pathogen)" + "	\t" + "AR_STD(pathogen)" + "	\t" + "Solidity_Mean(pathogen)" + "	\t" + "Solidity_STD(pathogen)" +  "	\t" 
	+ "Area_Mean(host)" + "	\t" + "Area_STD(host)" + "	\t" + "Perimeter_Mean(host)" + "	\t" + "Perimeter_STD(host)" + "	\t" + "AR_Mean(host)" + "	\t" + "AR_STD(host)" + "	\t" + "Solidity_Mean(host)" + "	\t" + "Solidity_STD(host)" + "	\n");	
//end of Results file set-up
		
dir = getDirectory("Choose a Directory ");
list = getFileList(dir);
print("Image folder: " + dir);
print("File list 0: " + list[0]);
print("List length: " + list.length);
numSubfolders = list.length;
subFolders = newArray(numSubfolders);//debugging

if(useSubfolders){
	for(i=0; i<list.length; i++) {
		if(endsWith(list[i], "/")) {
			//subFolders[i] = " "+dir+list[i];
			subFolders[i]=dir+list[i];
          	print("Subfolder found: " + subFolders[i]);
		}
	}
	numSubfolders=list.length;
	print("Number of subfolders = " + numSubfolders); //debugging
} else {
	numSubfolders=1;
	subFolders[0] = dir;
}

if(runInBatchMode==true){
   setBatchMode(true);
}


for(s=0; s<numSubfolders; s++){
	dir = subFolders[s];
	list = getFileList(dir);
	if(analyseAllImages==true){
		firstImageNumber = 0;
		lastImageNumber = list.length-1;
	}
	for(m=firstImageNumber;m<=lastImageNumber;m++){	
		testImageNumber = m;
		image_m = "" + dir + list[m];
		fullImagename = image_m;
		if(endsWith(image_m,imageType)==1){	
			if(useSpecificImageGroups==false || (useSpecificImageGroups==true && indexOf(image_m, searchTerm_1)>=0 && indexOf(image_m, searchTerm_2)>=0 && indexOf(image_m, searchTerm_3)>=0 && indexOf(image_m, excludeTerm_1)==-1 && indexOf(image_m, excludeTerm_2)==-1 && indexOf(image_m, excludeTerm_3)==-1)){
				if(useMultichannelTiffFormat==true){
					open(image_m);
					print("=================================================================================");
					print("Image name: " + image_m);
					run("Stack to Images");
					numOfChannels=nImages;//save original channel number to avoid confusion later
					currentImagename=list[m];//save image name for parameter file saving
					getDimensions(width, height, channels, slices, frames);
					imageWidth=width;
					imageHeight=height;
					imageBitDepth = bitDepth();
					print("Image size = " + width + " x " + height + " ; bit depth: " + imageBitDepth);
					width = getWidth();
					height = getHeight();
					channelNames = newArray(numOfChannels);
				} 
				else if (labelledMacrophages && imageType == ".czi"){
					run("Bio-Formats Importer", "open=[" + image_m + "] color_mode=Default rois_import=[ROI manager] split_channels view=[Standard ImageJ] stack_order=Default series_1");				
					print("=================================================================================");
					print("Image name: " + image_m);
					numOfChannels=nImages;//save original channel number to avoid confusion later
					currentImagename=list[m];//save image name for parameter file saving
					getDimensions(width, height, channels, slices, frames);
					imageWidth=width;
					imageHeight=height;
					imageBitDepth = bitDepth();
					print("Image size = " + width + " x " + height + " ; bit depth: " + imageBitDepth);
					width = getWidth();
					height = getHeight();
					channelNames = newArray(numOfChannels);
				}
				else if (!labelledMacrophages && imageType == ".czi"){
					run("Bio-Formats Importer", "open=[" + image_m + "] color_mode=Default rois_import=[ROI manager] split_channels view=[Standard ImageJ] stack_order=Default series_1");				
					print("=================================================================================");
					print("Image name: " + image_m);
					numOfChannels=nImages;//save original channel number to avoid confusion later
					currentImagename=list[m];//save image name for parameter file saving
					getDimensions(width, height, channels, slices, frames);
					imageWidth=width;
					imageHeight=height;
					imageBitDepth = bitDepth();
					print("Image size = " + width + " x " + height + " ; bit depth: " + imageBitDepth);
					width = getWidth();
					height = getHeight();
					channelNames = newArray(numOfChannels);
				}
				else if (!labelledMacrophages && imageType == ".czi" && !labelledSpores){
					run("Bio-Formats Importer", "open=[" + image_m + "] color_mode=Default rois_import=[ROI manager] split_channels view=[Standard ImageJ] stack_order=Default series_1");				
					print("=================================================================================");
					print("Image name: " + image_m);
					numOfChannels=nImages;//save original channel number to avoid confusion later#
					print("Number of channels = " + numOfChannels);
					currentImagename=list[m];//save image name for parameter file saving
					getDimensions(width, height, channels, slices, frames);
					imageWidth=width;
					imageHeight=height;
					imageBitDepth = bitDepth();
					print("Image size = " + width + " x " + height + " ; bit depth: " + imageBitDepth);
					width = getWidth();
					height = getHeight();
					channelNames = newArray(numOfChannels);
				}
				else {
					run("Bio-Formats Importer", "open=[" + image_m + "] color_mode=Default split_channels view=[Standard ImageJ] stack_order=Default");				
					print("=================================================================================");
					print("Image name: " + image_m);
					numOfChannels=nImages;//save original channel number to avoid confusion later
					currentImagename=list[m];//save image name for parameter file saving
					getDimensions(width, height, channels, slices, frames);
					imageWidth=width;
					imageHeight=height;
					imageBitDepth = bitDepth();
					print("Image size = " + width + " x " + height + " ; bit depth: " + imageBitDepth);
					width = getWidth();
					height = getHeight();
					channelNames = newArray(numOfChannels);
				}
				
				if(useMultichannelTiffFormat==true){
					channelNames = newArray(4);//here we use the number of channels (4) from the GUI, not from the actual image! 
					print("Channel names from GUI: " + channel00_new + "; " + channel01_new + "; " + channel02_new + "; " + channel03_new + "; " );
					print("Channel number = " + numOfChannels);
					if(useGUIforChannelInfo == false){
						print("Please use the GUI to assign the channels!!");
					}
					else {
						channelNames[0]=channel00;
						channelNames[1]=channel01;
						channelNames[2]=channel02;
						channelNames[3]=channel03;
					}
					listOfImages = getList("image.titles");
					print("Number of image windows = "+listOfImages.length);
				    if (listOfImages.length==0)
				       print("No image windows are open");
				    else {
				       print("Image windows:");
				       for (iCh=0; iCh<listOfImages.length; iCh++){
				          print("   "+listOfImages[iCh]);
				          print("   "+channelNames[iCh]);
				          selectWindow(listOfImages[iCh]);
				          rename(channelNames[iCh]);
				          }
				    }
				}
				else{
					channelNames = newArray(4);//here we use the number of channels (4) from the GUI, not from the actual image! 
					print("Channel names from GUI: " + channel00_new + "; " + channel01_new + "; " + channel02_new + "; " + channel03_new + "; " );
					print("Channel number = " + numOfChannels);
					if(useGUIforChannelInfo == false){
						print("Please use the GUI to assign the channels!!");
					}
					else {
						channelNames[0]=channel00;
						channelNames[1]=channel01;
						channelNames[2]=channel02;
						channelNames[3]=channel03;
					}
					listOfImages = getList("image.titles");
					print("Number of image windows = "+listOfImages.length);
				    if (listOfImages.length==0)
				       print("No image windows are open");
				    else {
				       print("Image windows:");
				       for (iCh=0; iCh<listOfImages.length; iCh++){
				          print("   "+listOfImages[iCh]);
				          print("   "+channelNames[iCh]);
				          selectWindow(listOfImages[iCh]);
				          rename(channelNames[iCh]);
				          }
				    }			
				}
				
				//*** now the re-stitching in the proper order: ***
				if(rearrangeZeiss==true){
					panelsizeX = width/XTileNumber;
					panelsizeY = height/YTileNumber;
					for(ch=0;ch<numOfChannels;ch++){ 
						selectWindow(channelNames[ch]);
						rename("RawData");
						for(mi=0;mi<XTileNumber;mi++){
							for(mj=0;mj<YTileNumber;mj++){
								selectWindow("RawData");
								makeRectangle(mi*panelsizeX + 1, mj*panelsizeY + 1, panelsizeX, panelsizeY);
								li=XTileNumber - mj - 1;
								lj=YTileNumber - mi - 1; 
								//tilename = "panel_" + li + "_" + lj;
								tilename = "panel_" + mi + "_" + mj;
								command1="title=" + tilename;
								run("Duplicate...", command1);
							}
						}
						if(XTileNumber == 3 && YTileNumber == 3){ //run for 3x3 matrix:
							run("Combine...", "stack1=panel_2_2 stack2=panel_2_1"); 
							run("Duplicate...", "title=combo01");
							run("Combine...", "stack1=combo01 stack2=panel_2_0");
							run("Duplicate...", "title=HorizontalLine1");
							
							run("Combine...", "stack1=panel_1_2 stack2=panel_1_1"); 
							run("Duplicate...", "title=combo02");
							run("Combine...", "stack1=combo02 stack2=panel_1_0");
							run("Duplicate...", "title=HorizontalLine2");
							
							run("Combine...", "stack1=panel_0_2 stack2=panel_0_1"); 
							run("Duplicate...", "title=combo03");
							run("Combine...", "stack1=combo03 stack2=panel_0_0");
							run("Duplicate...", "title=HorizontalLine3");
							
							//now the vertical stitching
							run("Combine...", "stack1=HorizontalLine1 stack2=HorizontalLine2 combine"); 
							run("Duplicate...", "title=VerticalCombo1");
							run("Combine...", "stack1=VerticalCombo1 stack2=HorizontalLine3 combine");
						}
						else if(XTileNumber == 5 && YTileNumber == 5){ //run for 5x5 matrix:
							run("Combine...", "stack1=panel_4_4 stack2=panel_4_3"); 
							run("Duplicate...", "title=combo01");
							run("Combine...", "stack1=combo01 stack2=panel_4_2");
							run("Duplicate...", "title=combo02");
							run("Combine...", "stack1=combo02 stack2=panel_4_1");
							run("Duplicate...", "title=combo03");
							run("Combine...", "stack1=combo03 stack2=panel_4_0");
							run("Duplicate...", "title=HorizontalLine0");
				
							run("Combine...", "stack1=panel_3_4 stack2=panel_3_3"); 
							run("Duplicate...", "title=combo11");
							run("Combine...", "stack1=combo11 stack2=panel_3_2");
							run("Duplicate...", "title=combo12");
							run("Combine...", "stack1=combo12 stack2=panel_3_1");
							run("Duplicate...", "title=combo13");
							run("Combine...", "stack1=combo13 stack2=panel_3_0");
							run("Duplicate...", "title=HorizontalLine1");
							
							run("Combine...", "stack1=panel_2_4 stack2=panel_2_3"); 
							run("Duplicate...", "title=combo21");
							run("Combine...", "stack1=combo21 stack2=panel_2_2");
							run("Duplicate...", "title=combo22");
							run("Combine...", "stack1=combo22 stack2=panel_2_1");
							run("Duplicate...", "title=combo23");
							run("Combine...", "stack1=combo23 stack2=panel_2_0");
							run("Duplicate...", "title=HorizontalLine2");
				
							run("Combine...", "stack1=panel_1_4 stack2=panel_1_3"); 
							run("Duplicate...", "title=combo31");
							run("Combine...", "stack1=combo31 stack2=panel_1_2");
							run("Duplicate...", "title=combo32");
							run("Combine...", "stack1=combo32 stack2=panel_1_1");
							run("Duplicate...", "title=combo33");
							run("Combine...", "stack1=combo33 stack2=panel_1_0");
							run("Duplicate...", "title=HorizontalLine3");
				
							run("Combine...", "stack1=panel_0_4 stack2=panel_0_3"); 
							run("Duplicate...", "title=combo41");
							run("Combine...", "stack1=combo41 stack2=panel_0_2");
							run("Duplicate...", "title=combo42");
							run("Combine...", "stack1=combo42 stack2=panel_0_1");
							run("Duplicate...", "title=combo43");
							run("Combine...", "stack1=combo43 stack2=panel_0_0");
							run("Duplicate...", "title=HorizontalLine4");
				
							//now the vertical stitching
							run("Combine...", "stack1=HorizontalLine0 stack2=HorizontalLine1 combine"); 
							run("Duplicate...", "title=VerticalCombo0");
							run("Combine...", "stack1=VerticalCombo0 stack2=HorizontalLine2 combine");
							run("Duplicate...", "title=VerticalCombo1");
							run("Combine...", "stack1=VerticalCombo1 stack2=HorizontalLine3 combine");
							run("Duplicate...", "title=VerticalCombo2");
							run("Combine...", "stack1=VerticalCombo2 stack2=HorizontalLine4 combine");
							run("Duplicate...", "title=VerticalCombo3");
						}
						else {
							print("Unhandled array size!  Must use 3x3 or 5x5 arrays!");
						}
		
						if(saveResults==true){
						   	if(useImageFilename==true){
								saveAs("png", dir + list[m] + "__" + nameTagSavedFiles + "__Channel_" + channelNames[ch] + "__Restitched.png");
							}
							else{
								saveAs("png", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Channel_" + channelNames[ch] + "__Restitched.png");
							}
						}
						close();
						close("RawData");
					}
				}

				//*** reopen the correctly stitched images: ***
				if(rearrangeZeiss == true){
					run("Close All");//debugging
					for(cho=0;cho<numOfChannels;cho++){
					   	if(useImageFilename==true){
							open(dir + list[m] + "__" + nameTagSavedFiles + "__Channel_" + channelNames[cho] + "__Restitched.png");
						}
						else{
							open("png", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Channel_" + channelNames[cho] + "__Restitched.png");
						}
						rename(channelNames[cho]);
						run("32-bit");
					}
				} else {
						for(cho=0;cho<numOfChannels;cho++){
							selectWindow(channelNames[cho]);
							run("Properties...", "channels=1 slices=1 frames=1 unit=pixel pixel_width=1 pixel_height=1 voxel_depth=1");
							run("32-bit");
						}
				}

				//*** Illumination correction if desired for TL: ***
				if(correctIlluminationForTL==true){
					selectWindow("TL");
					run("32-bit");
					run("Duplicate...", "title=Illumination_mask ");
					run("Gaussian Blur...", "sigma=" + illuminationCorrectionSigmaForTL);
					getStatistics(TL_area, TL_mean, TL_min, TL_max);
					run("Divide...", "value=" + TL_max);
					imageCalculator("Divide create 32-bit", "TL","Illumination_mask");
					rename("TL_IlluminationCorrected");
					run("Duplicate...", "title=Macrophages");
					close("Normalized_TL");
				}
				else{
					selectWindow("TL");
					run("Duplicate...", "title=TL_IlluminationCorrected");
					run("Duplicate...", "title=Macrophages");
				}

				//*** Illumination correction if desired for red channel: ***
				if(correctIlluminationForRed==true){
					selectWindow("Red");
					rename("Red_original");
					run("32-bit");
					run("Duplicate...", "title=RedIllumination_mask ");
					run("Gaussian Blur...", "sigma=" + illuminationCorrectionSigmaForRed);
					getStatistics(Red_area, Red_mean, Red_min, Red_max);
					run("Divide...", "value=" + Red_max);
					imageCalculator("Divide create 32-bit", "Red_original","RedIllumination_mask");
					rename("Red_IlluminationCorrected");
					run("Duplicate...", "title=Red");
					close("Normalized_Red");
				}
				

				if(applyCLAHE==true){
					selectWindow("Macrophages");
					if(fastModeCLAHE){
						CLAHEmessage = "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope + " mask=*None* fast_(less_accurate");
					} else {
						CLAHEmessage = "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope + " mask=*None*";
					}					
					run("Enhance Local Contrast (CLAHE)", CLAHEmessage);
				}

				selectWindow("Macrophages");
				run("Duplicate...", "title=Spores_TL");

				//*** Create ROIs for excludeEdges option: ***
				roiManager("reset");
				makeRectangle(0, 0, edgeWidth, imageHeight);
				roiManager("Add");
				makeRectangle(imageWidth-edgeWidth, 0, edgeWidth, imageHeight);
				roiManager("Add");
				makeRectangle(0, 0, imageWidth, edgeWidth);
				roiManager("Add");
				makeRectangle(0, imageHeight-edgeWidth, imageWidth, edgeWidth);
				roiManager("Add");
				run("Select All");
				roiManager("Combine");
				roiManager("Add");
				if(useImageFilename==true){
					roiManager("save", dir + list[m] + "__Edges_ROIs.zip");
				}
				else{
					roiManager("save",  dir + "/Image_" + testImageNumber + "__Edges_ROIs.zip");
				}
				
//***************************************************************************************************************
				
				//*** Pathogen  finder with Hessian starts: ***
				if(labelledSpores == false) {
					selectWindow("TL");//not using illum correction
					run("Duplicate...", "title=Spores_TL");
					selectWindow("Spores_TL");
					run("FeatureJ Hessian", "largest absolute smoothing=" + smoothingHessianSpores);
					run("Morphological Filters", "operation=[Internal Gradient] element=Octagon radius=" + gradientRadiusHessianSpores);
					setAutoThreshold(thresholdMethodHessianSpores + " dark");
					run("Convert to Mask");
					run("Despeckle");
					run("Despeckle");
					rename("Pathogens_Hessian_Hough");
					run("Duplicate...", "title=Pathogens_Hessian_BrightSpot");
					run("Duplicate...", "title=Pathogens_Hessian_GuidedWatershed");

					//*** Here we use the Hough-filter to find circles in the thresholded image: ***
					if(segmentationMethodForUnlabeledPathogens=="Hough-filter"|| segmentationMethodForUnlabeledPathogens=="Hough and Watershed"){
						selectWindow("Pathogens_Hessian_Hough");
						run("Hough Circle Transform","minRadius=" + minimumRadiusHoughHessianSpores + ", maxRadius=" + maximumRadiusHoughHessianSpores + ", inc=" + incrementRadiusHoughHessianSpores + ", minCircles=" + minimumNumberHoughHessianSpores + ", maxCircles=" + maximumNumberHoughHessianSpores + ", threshold=" + thresholdHoughHessianSpores + ", resolution=" + resolutionHoughHessianSpores + ", ratio=" + ratioHoughHessianSpores + ", bandwidth=" + bandwidthRadiusHoughHessianSpores + ", local_radius=" + localradiusHoughHessianSpores + ",  reduce show_mask show_centroids show_scores results_table");
						while(!isOpen("Centroid map")){};
						wait(200);
						selectWindow("Centroid map");
						colorThresholdingHSB();
						selectWindow("Centroid map_colorThreshold");
						//run("Invert");//NOT needed when using colorThresholdingHSB()
						run("Open");
						run("Watershed");
						roiManager("reset");
						if(excludeEdges==true){
							if(gatherResultsForEntireImageSet==true){
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude summarize add");	
							}
							else{
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude clear summarize add");
							}	
						}
						else{
							if(gatherResultsForEntireImageSet==true){
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display summarize add");				
							}
							else{
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display clear summarize add");
							}
						}
						
						run("Invert");
						run("Dilate");
						rename("Pathogens_Hessian_Hough_ROIs");
						run("Duplicate...", "title=Pathogens_Hessian_Hough_ROIs_copy");
						selectWindow("Pathogens_Hessian_Hough_ROIs");
					
						nROI_spores = roiManager("count");
						nROI_Hough_spores =  = roiManager("count");//save a copy for the ROI merging
						print("Number of Hessian-based pathogens found= " + nROI_spores);
						if(saveResults==true && nROI_spores>0){
								if(useImageFilename==true){
									roiManager("save", dir + list[m] + "__Pathogens_ROIs.zip");
								}
								else{
									roiManager("save",  dir + "/Image_" + testImageNumber + "__Pathogens_ROIs.zip");
								}
						}
					}

					//*** Here we apply the Bright Spot method to the Hessian pathogens image, then we merge the two ROIs: ***
					if(segmentationMethodForUnlabeledPathogens=="Hough and Watershed"){
						selectWindow("Pathogens_Hessian_BrightSpot");
						run("Invert");
						roiManager("reset");
						if(excludeEdges==true){
							if(gatherResultsForEntireImageSet==true){
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude summarize add");	
							}
							else{
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude clear summarize add");
							}	
						}
						else{
							if(gatherResultsForEntireImageSet==true){
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display summarize add");				
							}
							else{
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display clear summarize add");
							}
						}
						run("Invert");
						run("Dilate");
						
						rename("Pathogens_Hessian_BrightSpot_ROIs");
						run("Duplicate...", "title=Pathogens_Hessian_BrightSpot_ROIs_copy");
						selectWindow("Pathogens_Hessian_BrightSpot_ROIs");
						
						nROI_BrightSpots_spores = roiManager("count");
						print("Number of Hessian-based Bright Spots pathogens found= " + nROI_BrightSpots_spores);
						if(saveResults==true && nROI_BrightSpots_spores>0){
								if(useImageFilename==true){
									roiManager("save", dir + list[m] + "__Pathogens_BrightSpots_ROIs.zip");
								}
								else{
									roiManager("save",  dir + "/Image_" + testImageNumber + "__Pathogens_BrightSpots_ROIs.zip");
								}
						}
						//*** End of Bright Spot pathogen finder ***
	
						//*** Merging Hough and Bright Spot ROIs if such option selected: ***
						roiManager("reset");
						if(nROI_Hough_spores>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__Pathogens_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__Pathogens_ROIs.zip");
							}
						
							run("Select All");
							if(nROI_Hough_spores>=2){
								roiManager("Combine");
							}
							roiManager("Add");	//now the merged ROI is in position "nROI_spores" of the ROI array
						}
						
						if(nROI_BrightSpots_spores>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__Pathogens_ROIs_BrightSpots.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__Pathogens_ROIs_BrightSpots.zip");
							}
							GreenSporesAlreadyCounted = newArray(nROI_BrightSpots_spores);
							GreenSporesAlreadyCountedCounter = 0;
							for (ibrightspot=nROI_Hough_spores+1; ibrightspot<nROI_BrightSpots_spores+nROI_Hough_spores+1; ibrightspot++) {
								roiManager("Select", newArray(ibrightspot,nROI_Hough_spores));
								roiManager("AND");
								if(selectionType != (-1)){ //already counted
									GreenSporesAlreadyCounted[GreenSporesAlreadyCountedCounter]=ibrightspot;
									GreenSporesAlreadyCountedCounter++;	
								}
							}
		
							if(GreenSporesAlreadyCountedCounter>0){
								roiManager("Select", GreenSporesAlreadyCounted);
								if(GreenSporesAlreadyCountedCounter >= 2){ //Combine only works with >=2 ROIs
									roiManager("Combine");
								}
								roiManager("Delete");
							}
							roiManager("Select", nROI_Hough_spores);
							roiManager("Delete");	//delete the merged original Hough ROIs
							nROI_spores = roiManager("count");	//refresh and save the total ROI number, now including Bright Spot pathogens 
							print("New total pathogen count after Bright Spots added = " + nROI_spores);
							if(saveResults==true && nROI_spores>0){
								if(useImageFilename==true){
									roiManager("save", dir + list[m] + "__AllPathogens_ROIs.zip");
								}
								else{
									roiManager("save",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
								}
							}
						}
						//*** End of merging Hough-filter and Bright Spot ROIs ***
					}
					
					//*** Here we apply only the Bright Spot method, where we use the center of each pathogen as marker for watershed: ***
					if(segmentationMethodForUnlabeledPathogens=="Guided Watershed"){
						selectWindow("Pathogens_Hessian_GuidedWatershed");
						run("Duplicate...", "title=Pathogens_Hessian_GuidedWatershed_copy");
						selectWindow("Pathogens_Hessian_GuidedWatershed");
						/*run("Duplicate...", "title=Mask");
						selectWindow("Mask");
						run("Gaussian Blur...", "sigma=3");
						run("Invert");
						setAutoThreshold("Otsu");
						setOption("BlackBackground", true);
						run("Convert to Mask");
						run("Dilate");
						run("Dilate");
						run("Fill Holes");
						run("Divide...", "value=255");
						imageCalculator("Multiply create", "Pathogens_Hessian_GuidedWatershed","Mask");
						run("Erode");
						run("Remove Outliers...", "radius=3 threshold=50 which=Bright");
						run("Analyze Particles...", "size=3-300 circularity=0.10-1.00 show=[Bare Outlines] display exclude clear summarize add");
						selectWindow("Drawing of Result of Pathogens_Hessian_GuidedWatershed");
						run("Invert");
						run("Fill Holes");
						rename("Seeds");*/

						selectWindow("Pathogens_Hessian_GuidedWatershed_copy");
						run("Invert");
							//run("Analyze Particles...", "size=20-3000 circularity=0.10-1.00 show=[Bare Outlines] display exclude clear summarize add");
						if(excludeEdges==true){
							if(gatherResultsForEntireImageSet==true){
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude summarize add");	
							}
							else{
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude clear summarize add");
							}	
						}
						else{
							if(gatherResultsForEntireImageSet==true){
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display summarize add");				
							}
							else{
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display clear summarize add");
							}
						}
						/*selectWindow("Drawing of Pathogens_Hessian_GuidedWatershed_copy");
						run("Invert");
						run("Fill Holes");
						run("Gaussian Blur...", "sigma=3");
						setAutoThreshold("Otsu");
						//run("Threshold...");
						run("Convert to Mask");
						
						roiManager("reset");
						if(excludeEdges==true){
							if(gatherResultsForEntireImageSet==true){
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude summarize add");	
							}
							else{
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude clear summarize add");
							}	
						}
						else{
							if(gatherResultsForEntireImageSet==true){
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display summarize add");				
							}
							else{
								run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display clear summarize add");
							}
						}
						//selectWindow("Drawing of Drawing of Pathogens_Hessian_GuidedWatershed_copy");*/
						run("Invert");
						run("Dilate");
						rename("Pathogens_Hessian_GuidedWatershed_ROIs");
						run("Duplicate...", "title=Pathogens_Hessian_GuidedWatershed_ROIs_copy");
						selectWindow("Pathogens_Hessian_GuidedWatershed_ROIs");
						
						nROI_GuidedWatershed_spores = roiManager("count");
						print("Number of Hessian-based Guided Watershed pathogens found= " + nROI_GuidedWatershed_spores);
						nROI_spores = roiManager("count");	//refresh and save the total ROI number, now the Guided Watershed pathogens 
						if(saveResults==true && nROI_GuidedWatershed_spores>0){
							if(useImageFilename==true){
								roiManager("save", dir + list[m] + "__Pathogens_GuidedWatershed_ROIs.zip");
							}
							else{
								roiManager("save",  dir + "/Image_" + testImageNumber + "__Pathogens_GuidedWatershed_ROIs.zip");
							}
						}
					}

					//*** Now back to handling all possible ROIs together, so the rest of the code runs regardless: ***
					rename("SegmentedPathogens");
					run("Duplicate...", "title=SegmentedAllPathogens");
					selectWindow("SegmentedPathogens");
					
					if(saveResults==true){
						if(useImageFilename==true){
							saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__SegmentedPathogens.jpg");
						}
						else{
							saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__SegmentedPathogens.jpg");
						}
					}

					//*** Now create pathogen ROI as if with green labelling, for further analysis: ***
					nROI_greenSpores=nROI_spores;
					if(saveResults==true && nROI_greenSpores>0){
							if(useImageFilename==true){
								roiManager("save", dir + list[m] + "__AllPathogens_ROIs.zip");
							}
							else{
								roiManager("save",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
							}
					}
					
					selectWindow("SegmentedAllPathogens");
					run("Duplicate...", "title=Green");//need this to create blue-green image later
					selectWindow("SegmentedAllPathogens");
					if(saveResults==true){
						if(useImageFilename==true){
							saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "SegmentedAllPathogens.jpg");
						}
						else{
							saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "SegmentedAllPathogens.jpg");
						}
					} 
				}
				//*** End of pathogen finder with Hessian ***

//***************************************************************************************************************
				
				//*** Object analysis for pathogens with Hessian starts: ***
				if(labelledSpores == false) {
					nROI_spores = roiManager("count");
					if(nROI_spores>0){
						roiManager("Show All");
						updateResults();
						run("Summarize");
						updateResults();
						selectWindow("Results");
						for(ii=0;ii<nROI_spores;ii++){
							feretAngle = getResult("FeretAngle",ii);
							if(feretAngle > 90){
								normFeretAngle = 180-feretAngle;
							} 
							else{
								normFeretAngle = feretAngle;
							}
							setResult("NormalisedFeretAngle",ii,normFeretAngle);
							perimeter = getResult("Perim.",ii);
							feret = getResult("Feret",ii);
							perimFeretRatio=perimeter/feret;
							setResult("PerimeterFeretRatio",ii,perimFeretRatio);
							setResult("Label", ii, list[m]);//add image name to results table so that the results can be analysed image by image
						}
						
						updateResults();
						run("Summarize");
						updateResults();
	
						//Save pathogen morphology values for the output file:
						pathogenAreaMean = getResult("Area", nROI_greenSpores);
						pathogenAreaSTD = getResult("Area", nROI_greenSpores+1);
						pathogenPerimMean = getResult("Perim.", nROI_greenSpores);
						pathogenPerimSTD = getResult("Perim.", nROI_greenSpores+1);
						pathogenARMean = getResult("AR", nROI_greenSpores);
						pathogenARSTD = getResult("AR", nROI_greenSpores+1);
						pathogenSolidityMean = getResult("Solidity", nROI_greenSpores);
						pathogenSoliditySTD = getResult("Solidity", nROI_greenSpores+1);
						pathogenPerimeterFeretRatioMean = getResult("PerimeterFeretRatio", nROI_greenSpores);
						pathogenPerimeterFeretRatioSTD = getResult("PerimeterFeretRatio", nROI_greenSpores+1);
						//End of saving pathogen morphology values		
										
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("results", dir + list[m] + "__" + nameTagSavedFiles + "__Segmentation_Results_Pathogens.txt");
							}
							else{
								saveAs("results", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Segmentation_Results_Pathogens.txt");
							}
						}
	
						if(!fastModeNoPlots){
							run("Distribution...", "parameter=Circ. automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Circularity_Distribution_Spores.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Circularity_Distribution_Spores.jpg");
								}
							}
							run("Distribution...", "parameter=Round automatic");
							if(saveResults==true){
							   if(useImageFilename==true){
							   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Roundness_Distribution_Spores.jpg");
							   }
							   else{
							   		saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Roundness_Distribution_Spores.jpg");
							   }
							}
							run("Distribution...", "parameter=AR automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__AspectRatio_Distribution_Spores.jpg");	
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AspectRatio_Distribution_Spores.jpg");
								}
							}
							run("Distribution...", "parameter=FeretAngle automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__FeretAngle_Distribution_Spores.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__FeretAngle_Distribution_Spores.jpg");
								}
							}
							run("Distribution...", "parameter=NormalisedFeretAngle automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__NormalisedFeretAngle_Distribution_Spores.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__NormalisedFeretAngle_Distribution_Spores.jpg");
								}
							}
							run("Distribution...", "parameter=PerimeterFeretRatio automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__PerimeterFeretRatio_Distribution_Spores.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__PerimeterFeretRatio_Distribution_Spores.jpg");
								}
							}
							run("Distribution...", "parameter=Area automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Area_Distribution_Spores.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Area_Distribution_Spores.jpg");
								}
							}
							run("Distribution...", "parameter=Perim. automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Perimeter_Distribution_Spores.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Perimeter_Distribution_Spores.jpg");
								}
							}
						}
						selectWindow("TL_IlluminationCorrected");
						run("Duplicate...", "title=Unlabelled_image_with_final_objects");
						roiManager("Show All without labels");
						run("Colors...", "foreground=yellow background=black selection=yellow");
						run("Flatten");
						if(saveResults==true){
							if(useImageFilename==true){
						   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Unlabelled_image_with_pathogens.jpg");
						   	}
						   	else{
						   		saveAs("Jpeg", dir + "/Image_" + "__" + testImageNumber + "__" + nameTagSavedFiles + "__Unlabelled_image_with_pathogens.jpg");
						   	}
						}
					} //end of 'if(nROI_spores>0)' block
				} //end of "if(labelledSpores == false)" block
				//*** End of object analysis for pathogens with Hessian***


//***************************************************************************************************************
			
				//*** Fluorescence-based segmentation for host cells starts: ***
				if(labelledMacrophages == true) {
					selectWindow("Red");
					run("Duplicate...", "title=Macrophages_red");
					run("16-bit");
					// Offer choice between Internal Gradient (to find objects with a hole inside) and simple bright spot analysis:
					if(internalGradientRadiusRed >= 1){ 
						if(fastModeCLAHE){
							CLAHEmessage =  "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope2FlscMacrophages + " mask=*None*  fast_(less_accurate)";	
						} else {
							CLAHEmessage =  "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope2FlscMacrophages + " mask=*None*";
						}
						run("Enhance Local Contrast (CLAHE)",CLAHEmessage);
						run("Gaussian Blur...", "sigma=" + gaussianSmoothingForRedMacrophages);
						run("8-bit");
						run("Morphological Filters", "operation=[Internal Gradient] element=Octagon radius=" + internalGradientRadiusRed);
						run("Enhance Local Contrast (CLAHE)",CLAHEmessage);
						setAutoThreshold(thresholdMethodRedFluorescence + " dark");
						setOption("BlackBackground", true);
						run("Convert to Mask");
						run("Dilate");
						run("Fill Holes");
						run("Dilate");
						run("Dilate");
						run("Dilate");
						run("Fill Holes");
						run("Erode");
						run("Erode");
						if(watershedOnMacrophages==true){
							run("Watershed");
						}
						for(i=0;i<additionalErodeStepsForLabelledMacrophages;i++){
							run("Erode");
						}
					}
					else{
						setAutoThreshold(thresholdMethodRedFluorescence + " dark");
						setOption("BlackBackground", true);
						run("Convert to Mask");
						for(i=0;i<dilateErodeStepsForMacrophages;i++){
							run("Dilate");
						}
						run("Fill Holes");
						for(i=0;i<dilateErodeStepsForMacrophages;i++){
							run("Erode");
						}
						if(gaussianSmoothingForMacrophages>=1){
							run("Gaussian Blur...", "sigma=" + gaussianSmoothingForMacrophages);
						}
						setAutoThreshold(thresholdMethodRedFluorescence + " dark");
						setOption("BlackBackground", true);
						run("Convert to Mask");
						if(watershedOnMacrophages==true){
							run("Watershed");
						}
					}
					
					roiManager("reset");
					if(excludeEdges==true){
						if(gatherResultsForEntireImageSet==true){
							run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display exclude summarize add");	
						}
						else{
							run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display exclude clear summarize add");
						}	
					}
					else{
						if(gatherResultsForEntireImageSet==true){
							run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display summarize add");				
						}
						else{
							run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display clear summarize add");
						}
					}

					//Now rerun the particle finder on enlarged areas:
					run("Invert");
					run("Fill Holes");
					for(i=0;i<enlargeRedMacrophages;i++){
						run("Dilate");
					}
					if(watershedOnMacrophages==true){
						run("Watershed");
					}
					if(excludeEdges==true){
						if(gatherResultsForEntireImageSet==true){
							run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display exclude summarize add");	
						}
						else{
							run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display exclude clear summarize add");
						}	
					}
					else{
						if(gatherResultsForEntireImageSet==true){
							run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display summarize add");				
						}
						else{
							run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display clear summarize add");
						}
					}
					
					nROI_redMacrophages = roiManager("count");
					print("Total labelled host cell count  = " + nROI_redMacrophages);
					if(saveResults==true && nROI_redMacrophages>0){
							if(useImageFilename==true){
								roiManager("save", dir + list[m] + "__RedMacrophages_ROIs.zip");
							}
							else{
								roiManager("save",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs.zip");
							}
							
							if(internalGradientRadiusRed >= 1){ 
								if(useImageFilename==true){
									roiManager("save", dir + list[m] + "__RedMacrophages_ROIs_InternalGradient.zip");
								}
								else{
									roiManager("save",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs_InternalGradient.zip");
								}
							}						
					}
					
					rename("SegmentedRedMacrophages");
					if(saveResults==true){
						if(useImageFilename==true){
							saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "SegmentedRedMacrophages.jpg");
						}
						else{
							saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "SegmentedRedMacrophages.jpg");
						}
					}
				}
				//*** End of fluorescence-based segmentation for red host cells ***

//***************************************************************************************************************
			
				//*** Object analysis for red host cells with fluorescence starts: ***
				if(labelledMacrophages == true) {
					roiManager("reset");
					if(nROI_redMacrophages>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__RedMacrophages_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs.zip");
						}
						
						if(internalGradientRadiusRed >= 1){ 
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__RedMacrophages_ROIs_InternalGradient.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs_InternalGradient.zip");
							}
						}		
					}
	
					numberOfRedMacrophages = nROI_redMacrophages; //save for later use
					if(nROI_redMacrophages>0){
						roiManager("Show All");
						updateResults();
						run("Summarize");
						updateResults();
						selectWindow("Results");
						for(ii=0;ii<nROI_redMacrophages;ii++){
							feretAngle = getResult("FeretAngle",ii);
							if(feretAngle > 90){
								normFeretAngle = 180-feretAngle;
							} 
							else{
								normFeretAngle = feretAngle;
							}
							setResult("NormalisedFeretAngle",ii,normFeretAngle);
							perimeter = getResult("Perim.",ii);
							feret = getResult("Feret",ii);
							perimFeretRatio=perimeter/feret;
							setResult("PerimeterFeretRatio",ii,perimFeretRatio);
							setResult("Label", ii, list[m]);//add image name to results table so that the results can be analysed image by image
						}
						
						updateResults();
						run("Summarize");
						updateResults();
	
						//Save host morphology values for the output file:
						hostFlscAreaMean = getResult("Area", nROI_redMacrophages);
						hostFlscAreaSTD = getResult("Area", nROI_redMacrophages+1);
						hostFlscPerimMean = getResult("Perim.", nROI_redMacrophages);
						hostFlscPerimSTD = getResult("Perim.", nROI_redMacrophages+1);
						hostFlscARMean = getResult("AR", nROI_redMacrophages);
						hostFlscARSTD = getResult("AR", nROI_redMacrophages+1);
						hostFlscSolidityMean = getResult("Solidity", nROI_redMacrophages);
						hostFlscSoliditySTD = getResult("Solidity", nROI_redMacrophages+1);
						hostFlscPerimeterFeretRatioMean = getResult("PerimeterFeretRatio", nROI_redMacrophages);
						hostFlscPerimeterFeretRatioSTD = getResult("PerimeterFeretRatio", nROI_redMacrophages+1);
						//End of saving host morphology values
					
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("results", dir + list[m] + "__" + nameTagSavedFiles + "__Segmentation_Results_RedMacrophages.txt");
							}
							else{
								saveAs("results", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Segmentation_Results_RedMacrophages.txt");
							}
						}

						if(!fastModeNoPlots){
							run("Distribution...", "parameter=Circ. automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Circularity_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Circularity_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=Round automatic");
							if(saveResults==true){
							   if(useImageFilename==true){
							   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Roundness_Distribution.jpg");
							   }
							   else{
							   		saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Roundness_Distribution.jpg");
							   }
							}
							run("Distribution...", "parameter=AR automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__AspectRatio_Distribution.jpg");	
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AspectRatio_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=FeretAngle automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__FeretAngle_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__FeretAngle_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=NormalisedFeretAngle automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__NormalisedFeretAngle_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__NormalisedFeretAngle_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=PerimeterFeretRatio automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__PerimeterFeretRatio_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__PerimeterFeretRatio_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=Area automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Area_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Area_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=Perim. automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Perimeter_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Perimeter_Distribution.jpg");
								}
							}
						}
						
						selectWindow("TL_IlluminationCorrected");
						run("Duplicate...", "title=Unlabelled_image_with_redMacrophages");
						roiManager("Show All without labels");
						run("Colors...", "foreground=yellow background=black selection=yellow");
						run("Flatten");
						if(saveResults==true){
							if(useImageFilename==true){
						   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Unlabelled_image_with_redMacrophages.jpg");
						   	}
						   	else{
						   		saveAs("Jpeg", dir + "/Image_" + "__" + testImageNumber + "__" + nameTagSavedFiles + "__Unlabelled_image_with_redMacrophages.jpg");
						   	}
						}
	
						selectWindow("Red");
						run("Duplicate...", "title=Red_image_with_redMacrophages");
						roiManager("Show All without labels");
						run("Colors...", "foreground=yellow background=black selection=yellow");
						run("Flatten");
						if(saveResults==true){
							if(useImageFilename==true){
						   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Red_image_with_redMacrophages.jpg");
						   	}
						   	else{
						   		saveAs("Jpeg", dir + "/Image_" + "__" + testImageNumber + "__" + nameTagSavedFiles + "__Red_image_with_redMacrophages.jpg");
						   	}
						}
					} //end of 'if(nROI_macrophages>0)' block
				} //end of 'if(labelledMacrophages)' block
				//*** End of object analysis for labelled host cells with fluorescence ***

//***************************************************************************************************************
				
				//*** Macrophage finder with Hessian: ***
				selectWindow("Macrophages");
				run("Sharpen");
				run("FeatureJ Hessian", "largest smallest absolute smoothing=" + hessianSmoothingForMacrophages);
				if(useMAXHessian == true) {
					selectWindow("Macrophages largest Hessian eigenvalues");
				}
				else {
					selectWindow("Macrophages smallest Hessian eigenvalues");
				}
				if(gaussianSmoothingForMacrophages>=1){
					run("Gaussian Blur...", "sigma=" + gaussianSmoothingForMacrophages);
				}
				setAutoThreshold(thresholdMethodHessian + " dark");
				setOption("BlackBackground", true);
				getThreshold(lower, upper);
				setThreshold(lower*lowerThresholdMultiplierMacrophages, upper*upperThresholdMultiplierMacrophages);
				run("Convert to Mask");

				for(i=0;i<dilateErodeStepsForMacrophages;i++){
					run("Dilate");
				}
				run("Fill Holes");
				for(i=0;i<dilateErodeStepsForMacrophages;i++){
					run("Erode");
				}

				run("Remove Outliers...", "radius=" + removeOutliersStep1 + " threshold=50 which=Bright");
				run("Remove Outliers...", "radius=" + removeOutliersStep2 + " threshold=50 which=Bright");
				
				if(watershedOnMacrophages==true){
					run("Watershed");
				}
				
				for(i=0;i<lineThicknessForObjects;i++){
					run("Erode");
				}
				for(i=0;i<additionalErodeStepsForMacrophages;i++){ 
						run("Erode");
					}
				roiManager("reset");
				if(excludeEdges==true){
					if(gatherResultsForEntireImageSet==true){
						run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display exclude summarize add");	
					}
					else{
						run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display exclude clear summarize add");
					}	
				}
				else{
					if(gatherResultsForEntireImageSet==true){
						run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display summarize add");				
					}
					else{
						run("Analyze Particles...", "size=" + minMacrophageSize + "-" + maxMacrophageSize + " circularity=" + minMacrophageCircularity + "-" + maxMacrophageCircularity + " show=[Bare Outlines] display clear summarize add");
					}
				}
				nROI_macrophages = roiManager("count");
				if(saveResults==true && nROI_macrophages>0){
						if(useImageFilename==true){
							roiManager("save", dir + list[m] + "__HostCells_ROIs.zip");
						}
						else{
							roiManager("save",  dir + "/Image_" + testImageNumber + "__HostCells_ROIs.zip");
						}
				}
				for(i=0;i<lineThicknessForObjects;i++){
					run("Erode");
				}
				rename("SegmentedHostCells");
				if(saveResults==true){
					if(useImageFilename==true){
						saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__SegmentedHostCells.jpg");
					}
					else{
						saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__SegmentedHostCells.jpg");
					}
				}
				//*** End of host cell finder with Hessian ***

//***************************************************************************************************************
			
				//*** Object analysis for host cells with Hessian starts: ***
				nROI_macrophages = roiManager("count");
				numberOfMacrophages = nROI_macrophages; //save for later use
				if(nROI_macrophages>0){
					roiManager("Show All");
					updateResults();
					run("Summarize");
					updateResults();
					selectWindow("Results");
					for(ii=0;ii<nROI_macrophages;ii++){
						feretAngle = getResult("FeretAngle",ii);
						if(feretAngle > 90){
							normFeretAngle = 180-feretAngle;
						} 
						else{
							normFeretAngle = feretAngle;
						}
						setResult("NormalisedFeretAngle",ii,normFeretAngle);
						perimeter = getResult("Perim.",ii);
						feret = getResult("Feret",ii);
						perimFeretRatio=perimeter/feret;
						setResult("PerimeterFeretRatio",ii,perimFeretRatio);
						setResult("Label", ii, list[m]);//add image name to results table so that the results can be analysed image by image
					}
					
					updateResults();
					run("Summarize");
					updateResults();

					//Save host morphology values for the output file:
					hostAreaMean = getResult("Area", nROI_macrophages);
					hostAreaSTD = getResult("Area", nROI_macrophages+1);
					hostPerimMean = getResult("Perim.", nROI_macrophages);
					hostPerimSTD = getResult("Perim.", nROI_macrophages+1);
					hostARMean = getResult("AR", nROI_macrophages);
					hostARSTD = getResult("AR", nROI_macrophages+1);
					hostSolidityMean = getResult("Solidity", nROI_macrophages);
					hostSoliditySTD = getResult("Solidity", nROI_macrophages+1);
					hostPerimeterFeretRatioMean = getResult("PerimeterFeretRatio", nROI_macrophages);
					hostPerimeterFeretRatioSTD = getResult("PerimeterFeretRatio", nROI_macrophages+1);
					//End of saving host morphology values
					
					if(saveResults==true){
						if(useImageFilename==true){
							saveAs("results", dir + list[m] + "__" + nameTagSavedFiles + "__Segmentation_Results_HostCells.txt");
						}
						else{
							saveAs("results", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Segmentation_Results_HostCells.txt");
						}
					}

					if(!fastModeNoPlots){
						run("Distribution...", "parameter=Circ. automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Circularity_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Circularity_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=Round automatic");
						if(saveResults==true){
						   if(useImageFilename==true){
						   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Roundness_Distribution.jpg");
						   }
						   else{
						   		saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Roundness_Distribution.jpg");
						   }
						}
						run("Distribution...", "parameter=AR automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__AspectRatio_Distribution.jpg");	
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AspectRatio_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=FeretAngle automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__FeretAngle_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__FeretAngle_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=NormalisedFeretAngle automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__NormalisedFeretAngle_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__NormalisedFeretAngle_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=PerimeterFeretRatio automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__PerimeterFeretRatio_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__PerimeterFeretRatio_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=Area automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Area_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Area_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=Perim. automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Perimeter_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Perimeter_Distribution.jpg");
							}
						}
					}
		
					selectWindow("TL_IlluminationCorrected");
					run("Duplicate...", "title=Unlabelled_image_with_final_objects");
					roiManager("Show All without labels");
					run("Colors...", "foreground=yellow background=black selection=yellow");
					run("Flatten");
					if(saveResults==true){
						if(useImageFilename==true){
					   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Unlabelled_image_with_hostCells.jpg");
					   	}
					   	else{
					   		saveAs("Jpeg", dir + "/Image_" + "__" + testImageNumber + "__" + nameTagSavedFiles + "__Unlabelled_image_with_hostCells.jpg");
					   	}
					}
				} //end of 'if(nROI_macrophages>0)' block
				//*** End of object analysis for macrophages with Hessian ***

//***************************************************************************************************************
			
				//*** Fluorescence-based segmentation for green pathogens starts: ***
				if(labelledSpores == true) {
					selectWindow("Green");
					run("Duplicate...", "title=Spores_green");
					if(applyCLAHEonFLSC==true){
						selectWindow("Spores_green");
						if(fastModeCLAHE){
							CLAHEmessage = "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope + " mask=*None* fast_(less_accurate)";
						} else {
							CLAHEmessage = "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope + " mask=*None*";						
						}
						//run("Enhance Local Contrast (CLAHE)", "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope + " mask=*None*");
						run("Enhance Local Contrast (CLAHE)", CLAHEmessage);
					}
					run("Subtract Background...", "rolling=" + rollingBallRadius);
					run("Multiply...", "value="  + greenChannelMultiplierForDimData);
					
					// Offer choice between Internal Gradient (to find objects with a hole inside) and simple bright spot analysis:
					if(internalGradientRadiusGreen >= 1){ 
						run("Morphological Filters", "operation=[Internal Gradient] element=Square radius=" + internalGradientRadiusGreen);
						run("Invert");
						setAutoThreshold(thresholdMethodGreenFluorescence);
						getThreshold(lower, upper);
						setThreshold(lower*lowerThresholdMultiplier, upper*upperThresholdMultiplier);
						run("Convert to Mask");
						run("Fill Holes");
						if(watershedOnSpores==true){
							run("Watershed");
						}
					}
					else{
						setAutoThreshold(thresholdMethodGreenFluorescence + " dark");
						setOption("BlackBackground", true);
						run("Convert to Mask");
						for(i=0;i<dilateErodeStepsForSpores;i++){
							run("Dilate");
						}
						run("Fill Holes");
						for(i=0;i<dilateErodeStepsForSpores;i++){
							run("Erode");
						}
						
						if(gaussianBlurForSpores>=1){
							run("Gaussian Blur...", "sigma=" + gaussianBlurForSpores);
						}
						setAutoThreshold(thresholdMethodGreenFluorescence + " dark");
						setOption("BlackBackground", true);
						run("Convert to Mask");
						if(watershedOnSpores==true){
							run("Watershed");
						}
					}
					for(i=0;i<lineThicknessForObjects;i++){
						run("Erode");
					}
	
					roiManager("reset");
					if(excludeEdges==true){
						if(gatherResultsForEntireImageSet==true){
							run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude summarize add");	
						}
						else{
							run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude clear summarize add");
						}	
					}
					else{
						if(gatherResultsForEntireImageSet==true){
							run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display summarize add");				
						}
						else{
							run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display clear summarize add");
						}
					}
					nROI_greenSpores = roiManager("count");
					print("Total pathogen count  = " + nROI_greenSpores);
					if(saveResults==true && nROI_greenSpores>0){
							if(useImageFilename==true){
								roiManager("save", dir + list[m] + "__AllPathogens_ROIs.zip");
							}
							else{
								roiManager("save",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
							}
							
							if(internalGradientRadiusGreen >= 1){ 
								if(useImageFilename==true){
									roiManager("save", dir + list[m] + "__AllPathogens_ROIs_InternalGradient.zip");
								}
								else{
									roiManager("save",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs_InternalGradient.zip");
								}
							}						
					}
					
					for(i=0;i<lineThicknessForObjects;i++){
						run("Erode");
					}
					rename("SegmentedAllPathogens");
					if(saveResults==true){
						if(useImageFilename==true){
							saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "SegmentedAllPathogens.jpg");
						}
						else{
							saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "SegmentedAllPathogens.jpg");
						}
					}
				}
				//*** End of fluorescence-based segmentation for green pathogens ***

//***************************************************************************************************************
		
				//*** Object analysis for green pathogens starts: ***
				if(labelledSpores == true) {
					nROI_greenSpores = roiManager("count");
					if(nROI_greenSpores>0){
						roiManager("Show All");
						updateResults();
						run("Summarize");
						updateResults();
						selectWindow("Results");
						for(ii=0;ii<nROI_greenSpores;ii++){
							feretAngle = getResult("FeretAngle",ii);
							if(feretAngle > 90){
								normFeretAngle = 180-feretAngle;
							} 
							else{
								normFeretAngle = feretAngle;
							}
							setResult("NormalisedFeretAngle",ii,normFeretAngle);
							perimeter = getResult("Perim.",ii);
							feret = getResult("Feret",ii);
							perimFeretRatio=perimeter/feret;
							setResult("PerimeterFeretRatio",ii,perimFeretRatio);
							if(saveFullPathImagename == true){
								setResult("Label", ii,fullImagename);//add full image name to results table so that the results can be analysed image by image
							} else {
								setResult("Label", ii, list[m]);
							}
						}
						updateResults();
						run("Summarize");
						updateResults();
	
						//Save pathogen morphology values for the output file:
						pathogenAreaMean = getResult("Area", nROI_greenSpores);
						pathogenAreaSTD = getResult("Area", nROI_greenSpores+1);
						pathogenPerimMean = getResult("Perim.", nROI_greenSpores);
						pathogenPerimSTD = getResult("Perim.", nROI_greenSpores+1);
						pathogenARMean = getResult("AR", nROI_greenSpores);
						pathogenARSTD = getResult("AR", nROI_greenSpores+1);
						pathogenSolidityMean = getResult("Solidity", nROI_greenSpores);
						pathogenSoliditySTD = getResult("Solidity", nROI_greenSpores+1);
						pathogenPerimeterFeretRatioMean = getResult("PerimeterFeretRatio", nROI_greenSpores);
						pathogenPerimeterFeretRatioSTD = getResult("PerimeterFeretRatio", nROI_greenSpores+1);
						//End of saving pathogen morphology values
						
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("results", dir + list[m] + "__" + nameTagSavedFiles + "__Segmentation_Results_AllPathogens.txt");
							}
							else{
								saveAs("results", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Segmentation_Results_AllPathogens.txt");
							}
						}
	
						if(!fastModeNoPlots){
							run("Distribution...", "parameter=Circ. automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Circularity_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Circularity_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=Round automatic");
							if(saveResults==true){
							   if(useImageFilename==true){
							   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Roundness_Distribution.jpg");
							   }
							   else{
							   		saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Roundness_Distribution.jpg");
							   }
							}
							run("Distribution...", "parameter=AR automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__AspectRatio_Distribution.jpg");	
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AspectRatio_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=FeretAngle automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__FeretAngle_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__FeretAngle_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=NormalisedFeretAngle automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__NormalisedFeretAngle_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__NormalisedFeretAngle_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=PerimeterFeretRatio automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__PerimeterFeretRatio_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__PerimeterFeretRatio_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=Area automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Area_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Area_Distribution.jpg");
								}
							}
							run("Distribution...", "parameter=Perim. automatic");
							if(saveResults==true){
								if(useImageFilename==true){
									saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Perimeter_Distribution.jpg");
								}
								else{
									saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Perimeter_Distribution.jpg");
								}
							}
						}
			
						selectWindow("TL_IlluminationCorrected");
						run("Duplicate...", "title=Unlabelled_image_with_greenSpores");
						roiManager("Show All without labels");
						run("Colors...", "foreground=green background=black selection=green");
						run("Flatten");
						if(saveResults==true){
							if(useImageFilename==true){
						   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "Unlabelled_image_with_pathogens.jpg");
						   	}
						   	else{
						   		saveAs("Jpeg", dir + "/Image_" + "__" + testImageNumber + "__" + nameTagSavedFiles + "Unlabelled_image_with_pathogens.jpg");
						   	}
						}
	
						selectWindow("Green");
						run("Duplicate...", "title=AllPathogens_image_with_pathogens");
						run("Enhance Contrast", "saturated=0.35");
						roiManager("Show All without labels");
						run("Colors...", "foreground=green background=black selection=green");
						run("Flatten");
						if(saveResults==true){
							if(useImageFilename==true){
						   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "AllPathogens_image_with_pathogens.jpg");
						   	}
						   	else{
						   		saveAs("Jpeg", dir + "/Image_" + "__" + testImageNumber + "__" + nameTagSavedFiles + "AllPathogens_image_with_pathogens.jpg");
						   	}
						}
					} //end of 'if(nROI_greenSpores>0)' block
				} //end of "if(labelledSpores == true)" block
				//*** End of object analysis for green pathogens ***
				
//***************************************************************************************************************
		
				//*** Fluorescence-based segmentation with Bright Spots for green pathogens starts (when combining Internal Gradient and Bright Spots): ***

				if(combineInternalGradientAndBrightSpotResults==true && labelledSpores == true){
					selectWindow("Green");
					run("Duplicate...", "title=Spores_green_BrightSpots");
					if(applyCLAHEonFLSC==true){
						selectWindow("Spores_green_BrightSpots");
						if(fastModeCLAHE){
							CLAHEmessage = "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope + " mask=*None* fast_(less_accurate)";
						} else {
							CLAHEmessage = "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope + " mask=*None*";
						}
						run("Enhance Local Contrast (CLAHE)", CLAHEmessage);
					}
					run("Subtract Background...", "rolling=" + rollingBallRadius);
					run("Multiply...", "value="  + greenChannelMultiplierForDimData);

					setAutoThreshold(thresholdMethodGreenFluorescenceBrightSpots + " dark"); //option to use different background method for Bright Spots method when combining IG and BS-based ROIs for green
					setOption("BlackBackground", true);
					run("Convert to Mask"); 
					for(i=0;i<dilateErodeStepsForSpores;i++){
						run("Dilate");
					}
					run("Fill Holes");
					for(i=0;i<dilateErodeStepsForSpores;i++){
						run("Erode");
					}
					
					if(gaussianBlurForSpores>=1){
						run("Gaussian Blur...", "sigma=" + gaussianBlurForSpores);
					}
					setAutoThreshold(thresholdMethodGreenFluorescence + " dark");
					setOption("BlackBackground", true);
					run("Convert to Mask");
					if(watershedOnSpores==true){
						run("Watershed");
					}
					
					for(i=0;i<lineThicknessForObjects;i++){
						run("Erode");
					}

					for(i=0;i<2;i++){ //try this for Bright Spots 'cause the ROIs are too loose
						run("Erode");
					}
	
					roiManager("reset");
					if(excludeEdges==true){
						if(gatherResultsForEntireImageSet==true){
							run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude summarize add");	
						}
						else{
							run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude clear summarize add");
						}	
					}
					else{
						if(gatherResultsForEntireImageSet==true){
							run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display summarize add");				
						}
						else{
							run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display clear summarize add");
						}
					}
					nROI_greenSpores_BrightSpots = roiManager("count");
					print("Bright Spots pathogen count = " + nROI_greenSpores_BrightSpots);
					if(saveResults==true && nROI_greenSpores_BrightSpots>0){
							if(useImageFilename==true){
								roiManager("save", dir + list[m] + "__AllPathogens_ROIs_BrightSpots.zip");
							}
							else{
								roiManager("save",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs_BrightSpots.zip");
							}						
					}
					
					rename("SegmentedAllPathogens_BrightSpots");
					if(saveResults==true){
						if(useImageFilename==true){
							saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "SegmentedAllPathogens_BrightSpots.jpg");
						}
						else{
							saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "SegmentedAllPathogens_BrightSpots.jpg");
						}
					}
				}
				//*** End of Fluorescence-based segmentation with Bright Spots for green pathogens (when combining Internal Gradient and Bright Spots) ***

//***************************************************************************************************************

				//*** Fluorescence-based segmentation for blue pathogens starts: ***
				selectWindow("Blue");
				run("Duplicate...", "title=Spores_blue");
				if(applyCLAHEonFLSC==true){
					selectWindow("Spores_blue");
					if(fastModeCLAHE){
						CLAHEmessage = "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope + " mask=*None* fast_(less_accurate)";
					} else {
						CLAHEmessage = "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope + " mask=*None*";
					}
					//run("Enhance Local Contrast (CLAHE)", "blocksize=" + CLAHEblocks + " histogram=" + CLAHEbins + " maximum=" + CLAHEslope + " mask=*None*");
					run("Enhance Local Contrast (CLAHE)", CLAHEmessage);
				}
				run("Subtract Background...", "rolling=" + rollingBallRadius);
				run("Multiply...", "value="  + blueChannelMultiplierForDimData);
				
				// Offer choice between Internal Gradient (to find objects with a hole inside) and simple bright spot analysis:
				if(internalGradientRadiusBlue >= 1){ 
					run("Morphological Filters", "operation=[Internal Gradient] element=Square radius=" + internalGradientRadiusBlue);
					run("Invert");
					setAutoThreshold(thresholdMethodBlueFluorescence);
					getThreshold(lower, upper);
					setThreshold(lower*lowerThresholdMultiplier, upper*upperThresholdMultiplier);
					run("Convert to Mask");
					run("Fill Holes");
					if(watershedOnSpores==true){
						run("Watershed");
					}
				}
				else{
					setAutoThreshold(thresholdMethodBlueFluorescence + " dark");
					setOption("BlackBackground", true);
					run("Convert to Mask");
					for(i=0;i<dilateErodeStepsForSpores;i++){
						run("Dilate");
					}
					run("Fill Holes");
					for(i=0;i<dilateErodeStepsForSpores;i++){
						run("Erode");
					}
					if(gaussianBlurForSpores>=1){
						run("Gaussian Blur...", "sigma=" + gaussianBlurForSpores);
					}
					setAutoThreshold("Otsu dark");
					setOption("BlackBackground", true);
					run("Convert to Mask");
					
					if(watershedOnSpores==true){
						run("Watershed");
					}
				}
				for(i=0;i<lineThicknessForObjects;i++){
					run("Erode");
				}
				roiManager("reset");
				if(excludeEdges==true){
					if(gatherResultsForEntireImageSet==true){
						run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude summarize add");	
					}
					else{
						run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display exclude clear summarize add");
					}	
				}
				else{
					if(gatherResultsForEntireImageSet==true){
						run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display summarize add");				
					}
					else{
						run("Analyze Particles...", "size=" + minSporeSize + "-" + maxSporeSize + " circularity=" + minSporeCircularity + "-" + maxSporeCircularity + " show=[Bare Outlines] display clear summarize add");
					}
				}
				nROI_blueSpores = roiManager("count");
				if(saveResults==true && nROI_blueSpores>0){
						if(useImageFilename==true){
							roiManager("save", dir + list[m] + "__OutsidePathogens_ROIs.zip");
						}
						else{
							roiManager("save",  dir + "/Image_" + testImageNumber + "__OutsidePathogens_ROIs.zip");
						}
				}
				for(i=0;i<lineThicknessForObjects;i++){
					run("Erode");
				}
				rename("SegmentedOutsidePathogens");
				if(saveResults==true){
					if(useImageFilename==true){
						saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "SegmentedOutsidePathogens.jpg");
					}
					else{
						saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "SegmentedOutsidePathogens.jpg");
					}
				}
				
				//*** End of fluorescence-based segmentation for blue pathogens ***
				
//***************************************************************************************************************
			
				//*** Object analysis for blue pathogens starts: ***
				nROI_blueSpores = roiManager("count");
				if(nROI_blueSpores>0){
					roiManager("Show All");
					updateResults();
					run("Summarize");
					updateResults();
					selectWindow("Results");
					for(ii=0;ii<nROI_blueSpores;ii++){
						feretAngle = getResult("FeretAngle",ii);
						if(feretAngle > 90){
							normFeretAngle = 180-feretAngle;
						} 
						else{
							normFeretAngle = feretAngle;
						}
						setResult("NormalisedFeretAngle",ii,normFeretAngle);
						perimeter = getResult("Perim.",ii);
						feret = getResult("Feret",ii);
						perimFeretRatio=perimeter/feret;
						setResult("PerimeterFeretRatio",ii,perimFeretRatio);
						setResult("Label", ii, list[m]);//add image name to results table so that the results can be analysed image by image
					}
					
					updateResults();
					if(saveResults==true){
						if(useImageFilename==true){
							saveAs("results", dir + list[m] + "__" + nameTagSavedFiles + "__Segmentation_Results_OutsidePathogens.txt");
						}
						else{
							saveAs("results", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Segmentation_Results_OutsidePathogens.txt");
						}
					}

					if(!fastModeNoPlots){
						run("Distribution...", "parameter=Circ. automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Circularity_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Circularity_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=Round automatic");
						if(saveResults==true){
						   if(useImageFilename==true){
						   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Roundness_Distribution.jpg");
						   }
						   else{
						   		saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Roundness_Distribution.jpg");
						   }
						}
						run("Distribution...", "parameter=AR automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__AspectRatio_Distribution.jpg");	
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AspectRatio_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=FeretAngle automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__FeretAngle_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__FeretAngle_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=NormalisedFeretAngle automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__NormalisedFeretAngle_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__NormalisedFeretAngle_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=PerimeterFeretRatio automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__PerimeterFeretRatio_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__PerimeterFeretRatio_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=Area automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Area_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Area_Distribution.jpg");
							}
						}
						run("Distribution...", "parameter=Perim. automatic");
						if(saveResults==true){
							if(useImageFilename==true){
								saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Perimeter_Distribution.jpg");
							}
							else{
								saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Perimeter_Distribution.jpg");
							}
						}
					}
		
					selectWindow("TL_IlluminationCorrected");
					run("Duplicate...", "title=Unlabelled_image_with_outsidePathogens");
					roiManager("Show All without labels");
					run("Colors...", "foreground=blue background=black selection=blue");
					run("Flatten");
					if(saveResults==true){
						if(useImageFilename==true){
					   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "Unlabelled_image_with_outsidePathogens.jpg");
					   	}
					   	else{
					   		saveAs("Jpeg", dir + "/Image_" + "__" + testImageNumber + "__" + nameTagSavedFiles + "Unlabelled_image_with_outsidePathogens.jpg");
					   	}
					}

					selectWindow("Blue");
					run("Duplicate...", "title=OutsidePathogen_image_with_outsidePathogens");
					run("Enhance Contrast", "saturated=0.35");
					roiManager("Show All without labels");
					run("Colors...", "foreground=blue background=black selection=blue");
					run("Flatten");
					if(saveResults==true){
						if(useImageFilename==true){
					   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "OutsidePathogen_image_with_outsidePathogens.jpg");
					   	}
					   	else{
					   		saveAs("Jpeg", dir + "/Image_" + "__" + testImageNumber + "__" + nameTagSavedFiles + "OutsidePathogen_image_with_outsidePathogens.jpg");
					   	}
					}
				} //end of 'if(nROI_blueSpores>0)' block
				//*** End of object analysis for blue pathogens ***

//***************************************************************************************************************
			
				//*** Hard thresholding the Bright Spots-method green ROIs for green channel fluorescence to eliminate dim pathogens, used only when merging Internal Gradient and Bright Spots: ***
				if(nROI_greenSpores_BrightSpots>0){
					GreenSporesAboveGreenThreshold = newArray(nROI_greenSpores_BrightSpots);
					GreenSporesBelowGreenThreshold = newArray(nROI_greenSpores_BrightSpots);
					GreenSporesAboveGreenThresholdCounter = 0;
					GreenSporesBelowGreenThresholdCounter = 0;
					roiManager("reset");
					run("Clear Results");
					selectWindow("Green");
					run("Duplicate...", "title=AllPathogens_image_with_aboveThreshold_BrightSpots_pathogens");
					selectWindow("AllPathogens_image_with_aboveThreshold_BrightSpots_pathogens");
					run("8-bit");
					if(nROI_greenSpores_BrightSpots>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__AllPathogens_ROIs_BrightSpots.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs_BrightSpots.zip");
						}
					}
					for (green=0; green<nROI_greenSpores_BrightSpots; green++) {
						roiManager("Select", green);
						roiManager("Measure");
						updateResults();
						meanGreenAtGreenROI = getResult("Mean",green);
						if(meanGreenAtGreenROI >= greenThresholdForGreenSporeClassifier){
							GreenSporesAboveGreenThreshold[green] = green;
							GreenSporesBelowGreenThreshold[green] = 0;
							GreenSporesAboveGreenThresholdCounter++;
						}
						else{
							GreenSporesBelowGreenThreshold[green] = green;
							GreenSporesAboveGreenThreshold[green] = 0;
							GreenSporesBelowGreenThresholdCounter++;
						}
					}
					roiManager("Set Line Width", 2);
					roiManager("Set Color", "green");
					run("Colors...", "foreground=blue background=black selection=green");
					roiManager("Select", GreenSporesAboveGreenThreshold);
					roiManager("Draw");
					run("Flatten");
					run("Colors...", "foreground=yellow background=black selection=yellow");
					roiManager("Select", GreenSporesBelowGreenThreshold);
					roiManager("Draw");
					run("Flatten");
					if(saveResults==true){
						if(useImageFilename==true){
					   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "AllPathogens_image_with_aboveThreshold_BrightSpots_pathogens.jpg");
					   	}
					   	else{
					   		saveAs("Jpeg", dir + "/Image_" + "__" + testImageNumber + "__" + nameTagSavedFiles + "AllPathogens_image_with_aboveThreshold_BrightSpots_pathogens.jpg");
					   	}
					}
					print("Number of above-threshold and below-threshold all pathogens = " + GreenSporesAboveGreenThresholdCounter + " ; " + GreenSporesBelowGreenThresholdCounter);
					
					if(GreenSporesBelowGreenThresholdCounter>0){
						roiManager("Select", GreenSporesBelowGreenThreshold);
						if(GreenSporesBelowGreenThresholdCounter >= 2){ //Combine only works with >=2 ROIs
							roiManager("Combine");
						}
						roiManager("Delete");
					}
					nROI_greenSporesAboveGreenThreshold = roiManager("count");
					nROI_greenSpores_BrightSpots = nROI_greenSporesAboveGreenThreshold;
					if(saveResults==true && nROI_greenSporesAboveGreenThreshold>0){
							if(useImageFilename==true){
								roiManager("save", dir + list[m] + "__AllPathogens_ROIs_BrightSpots.zip");
							}
							else{
								roiManager("save",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs_BrightSpots.zip");
							}
					}
				}

				//*** End of hard thresholding the Bright Spots-method green ROIs for green channel fluorescence to eliminate dim pathogens, used only when merging Internal Gradient and Bright Spots ***

//***************************************************************************************************************
				
				//*** Merging Inside Gradient and Bright Spot ROIs if such option selected: ***
				if(combineInternalGradientAndBrightSpotResults==true && labelledSpores==true){
					roiManager("reset");
					if(nROI_greenSpores>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__AllPathogens_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
						}
					
						run("Select All");
						if(nROI_greenSpores>=2){
							roiManager("Combine");
						}
						roiManager("Add");	//now the merged ROI is in position "nROI_greenSpores" of the ROI array
					}
					
					if(nROI_greenSpores_BrightSpots>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__AllPathogens_ROIs_BrightSpots.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs_BrightSpots.zip");
						}
						GreenSporesAlreadyCounted = newArray(nROI_greenSpores_BrightSpots);
						GreenSporesAlreadyCountedCounter = 0;
						for (ibrightspot=nROI_greenSpores+1; ibrightspot<nROI_greenSpores_BrightSpots+nROI_greenSpores+1; ibrightspot++) {
							roiManager("Select", newArray(ibrightspot,nROI_greenSpores));
							roiManager("AND");
							if(selectionType != (-1)){ //already counted green spore
								GreenSporesAlreadyCounted[GreenSporesAlreadyCountedCounter]=ibrightspot;
								GreenSporesAlreadyCountedCounter++;	
							}
						}
	
						if(GreenSporesAlreadyCountedCounter>0){
							roiManager("Select", GreenSporesAlreadyCounted);
							if(GreenSporesAlreadyCountedCounter >= 2){ //Combine only works with >=2 ROIs
								roiManager("Combine");
							}
							roiManager("Delete");
						}
						roiManager("Select", nROI_greenSpores);
						roiManager("Delete");	//delete the merged original green ROIs
						nROI_greenSpores = roiManager("count");	//refresh and save the green ROI number, now including Bright Spot pathogens 
						print("New total pathogen count after Bright Spots added = " + nROI_greenSpores);
						if(saveResults==true && nROI_greenSpores>0){
							if(useImageFilename==true){
								roiManager("save", dir + list[m] + "__AllPathogens_ROIs.zip");
							}
							else{
								roiManager("save",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
							}
						}
					}
				}
				
				//*** End of merging Inside Gradient and Bright Spot ROIs ***
				
//***************************************************************************************************************
			
				//*** Analysis of green ROIs for blue pathogen fluorescence to classify in- and outside pathogens: ***
				GreenSporesOutside = newArray(nROI_greenSpores);
				GreenSporesInside = newArray(nROI_greenSpores);
				GreenSporesOutsideCounter = 0;
				GreenSporesInsideCounter = 0;
				roiManager("reset");
				run("Clear Results");
				selectWindow("Blue");
				run("Duplicate...", "title=OutsidePathogens_image_with_classfied_pathogens");
				selectWindow("OutsidePathogens_image_with_classfied_pathogens");
				run("8-bit");
				if(nROI_greenSpores>0){
					if(useImageFilename==true){
						roiManager("open", dir + list[m] + "__AllPathogens_ROIs.zip");
					}
					else{
						roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
					}
				}
				for (greenSp=0; greenSp<nROI_greenSpores; greenSp++) {
					roiManager("Select", greenSp);
					roiManager("Measure");
					updateResults();
					meanBlueAtGreenROI = getResult("Mean",greenSp);
					if(meanBlueAtGreenROI >= blueThresholdForGreenSporeClassifier){
						GreenSporesOutside[greenSp] = greenSp;
						GreenSporesInside[greenSp] = 0;
						GreenSporesOutsideCounter++;
					}
					else{
						GreenSporesInside[greenSp] = greenSp;
						GreenSporesOutside[greenSp] = 0;
						GreenSporesInsideCounter++;
					}
				}
				roiManager("Set Line Width", 2);
				roiManager("Set Color", "green");
				run("Colors...", "foreground=blue background=black selection=green");
				roiManager("Select", GreenSporesInside);
				roiManager("Draw");
				run("Flatten");
				run("Colors...", "foreground=yellow background=black selection=yellow");
				roiManager("Select", GreenSporesOutside);
				roiManager("Draw");
				run("Flatten");
				if(saveResults==true){
					if(useImageFilename==true){
				   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "OutsidePathogens_image_with_classfied_pathogens.jpg");
				   	}
				   	else{
				   		saveAs("Jpeg", dir + "/Image_" + "__" + testImageNumber + "__" + nameTagSavedFiles + "OutsidePathogens_image_with_classfied_pathogens.jpg");
				   	}
				}
				print("Number of outside and inside all pathogens = " + GreenSporesOutsideCounter + " ; " + GreenSporesInsideCounter);
				
				selectWindow("Green");
				run("Duplicate...", "title=AllPathogens_image_with_classfied_pathogens");
				selectWindow("AllPathogens_image_with_classfied_pathogens");
				roiManager("Set Line Width", 2);
				run("Colors...", "foreground=blue background=black selection=green");
				if(GreenSporesInsideCounter>0){
					roiManager("Select", GreenSporesInside);
					roiManager("Draw");
					run("Flatten");
				}
				run("Colors...", "foreground=yellow background=black selection=yellow");
				if(GreenSporesOutsideCounter>0){
					roiManager("Select", GreenSporesOutside);
					roiManager("Draw");
					run("Flatten");
				}
				if(saveResults==true){
					if(useImageFilename==true){
				   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "AllPathogens_image_with_classfied_pathogens.jpg");
				   	}
				   	else{
				   		saveAs("Jpeg", dir + "/Image_" + "__" + testImageNumber + "__" + nameTagSavedFiles + "AllPathogens_image_with_classfied_pathogens.jpg");
				   	}
				}

				roiManager("reset");//save ROIs for green pathogens that are NOT blue, i.e. they are inside MPs
				if(nROI_greenSpores>0){
					if(useImageFilename==true){
						roiManager("open", dir + list[m] + "__AllPathogens_ROIs.zip");
					}
					else{
						roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
					}
				}
				if(GreenSporesOutsideCounter>0){
					roiManager("Select", GreenSporesOutside);
					roiManager("Update");
					if(GreenSporesOutsideCounter >= 2){ //Combine only works with >=2 ROIs
						roiManager("Combine");
					}
					roiManager("Delete");
				}
				nROI_greenSporesInside = roiManager("count");
				if(saveResults==true && nROI_greenSporesInside>0){
						if(useImageFilename==true){
							roiManager("save", dir + list[m] + "__AllPathogensInside_ROIs.zip");
						}
						else{
							roiManager("save",  dir + "/Image_" + testImageNumber + "__AllPathogensInside_ROIs.zip");
						}
				}

				roiManager("reset");//save ROIs for green pathogens that are also blue, i.e. they are outside MPs
				if(nROI_greenSpores>0){
					if(useImageFilename==true){
						roiManager("open", dir + list[m] + "__AllPathogens_ROIs.zip");
					}
					else{
						roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
					}
				}
				if(GreenSporesInsideCounter>0){
					roiManager("Select", GreenSporesInside);
					roiManager("Update");
					if(GreenSporesInsideCounter >= 2){ //Combine only works with >=2 ROIs
						roiManager("Combine");
					}
					roiManager("Delete");
				}
				nROI_greenSporesOutside = roiManager("count");
				if(saveResults==true && nROI_greenSporesOutside>0){
						if(useImageFilename==true){
							roiManager("save", dir + list[m] + "__AllPathogensOutside_ROIs.zip");
						}
						else{
							roiManager("save",  dir + "/Image_" + testImageNumber + "__AllPathogensOutside_ROIs.zip");
						}
				}
				
				//*** End of analysis of green ROIs for blue pathogen fluorescence to classify in- and outside pathogens: ***

//***************************************************************************************************************

				//*** Saving images begins: ***
				roiManager("reset");
				selectWindow("TL_IlluminationCorrected");
				run("Duplicate...", "title=TL_IlluminationCorrected_saveForLater");
				selectWindow("TL_IlluminationCorrected");
				if(saveResults==true){
					if(useImageFilename==true){
				   		saveAs("png", dir + list[m] + "__" + nameTagSavedFiles + "__TL_image_Illumination_corrected.png");
				   	}
				   	else{
				   		saveAs("png", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__TL_image_Illumination_corrected.png");
				   	}
				}
				
				if(nROI_macrophages>0){
					if(useImageFilename==true){
						roiManager("open", dir + list[m] + "__HostCells_ROIs.zip");
					}
					else{
						roiManager("open",  dir + "/Image_" + testImageNumber + "__HostCells_ROIs.zip");
					}
					nr_macrophages=roiManager("count");
					Macrophages_01 = newArray(nr_macrophages);
					roiManager("Set Fill Color", "green");
					roiManager("Set Color", "yellow");
					roiManager("Set Line Width", 2);
					for (i=0; i<nr_macrophages; i++) {
						Macrophages_01[i] = i;
						roiManager("select", Macrophages_01);
						roiManager("Draw");
					}
				}
				run("Flatten");
				if(saveResults==true){
					if(useImageFilename==true){
				   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Raw_image_with_hostCells.jpg");
				   	}
				   	else{
				   		saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Raw_image_with_hostCells.jpg");
				   	}
				}
				//*** Saving images ends ***

				//*** Merge all components into one image and save it ***
				if(useImageFilename==true){
			   		macrophageWindowName=list[m] + "__" + nameTagSavedFiles + "__Raw_image_with_hostCells.jpg";
			   		blueSporesWindowName=list[m] + "__" + nameTagSavedFiles + "SegmentedOutsidePathogens.jpg";
			   		greenSporesWindowName=list[m] + "__" + nameTagSavedFiles + "SegmentedAllPathogens.jpg";
			   		selectWindow(macrophageWindowName);
			   		selectWindow(blueSporesWindowName);
			   		run("Invert");
					run("Blue");
					run("RGB Color");
					run("Merge Channels...", "c3=" + blueSporesWindowName + " c4=" + macrophageWindowName + " keep");
					selectWindow(greenSporesWindowName);
					run("Invert");
					run("Green");
					run("RGB Color");
					selectWindow("RGB");
					run("Merge Channels...", "c2=" + greenSporesWindowName + " c4=RGB keep");
			   	}
			   	else{
			   		macrophageWindowName="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Raw_image_with_hostCells.jpg";
			   		blueSporesWindowName="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "SegmentedOutsidePathogens.jpg";
			   		greenSporesWindowName="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "SegmentedAllPathogens.jpg";
			   		selectWindow(macrophageWindowName);
			   		selectWindow(blueSporesWindowName);
			   		run("Invert");
					run("Blue");
					run("RGB Color");
					run("Merge Channels...", "c3=" + blueSporesWindowName + " c4=" + macrophageWindowName + " keep");
					selectWindow(greenSporesWindowName);
					run("Invert");
					run("Green");
					run("RGB Color");
					selectWindow("RGB");
					run("Merge Channels...", "c2=" + greenSporesWindowName + " c4=RGB keep");
			   	}
			   	if(saveResults==true){
					if(useImageFilename==true){
				   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__AllComponents_Merged.jpg");
				   	}
				   	else{
				   		saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AllComponents_Merged.jpg");
				   	}
				}
				//*** End of image merging and saving ***

//***************************************************************************************************************
//Firstly, use Hessian-based host cell ROIs to calculate phagocytic measures:
//***************************************************************************************************************

				//*** Calculate ROI overlap for hosts and pathogens *** 
				roiManager("reset");
				if(nROI_macrophages>0){
					if(useImageFilename==true){
						roiManager("open", dir + list[m] + "__HostCells_ROIs.zip");
					}
					else{
						roiManager("open",  dir + "/Image_" + testImageNumber + "__HostCells_ROIs.zip");
					}
					if(nROI_greenSpores>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__AllPathogens_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
						}
					}
					if(useImageFilename==true){
				   		mergedWindowName=list[m] + "__" + nameTagSavedFiles + "__AllComponents_Merged.jpg";
				   	}
				   	else{
				   		mergedWindowName="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AllComponents_Merged.jpg";
				   	}
					numberOfGreenSporesInsideMacrophages = 0;
					numberOfNotBlueGreenSporesInsideMacrophages = 0;
					alreadyCounted = newArray(nROI_greenSpores);//in order to avoid counting the same pathogen twice if it overlaps two host cells
					roiManager("Set Color", "red");
					roiManager("Set Line Width", 2);
					for (imacro=0; imacro<nROI_macrophages; imacro++) {
						for (ispore=nROI_macrophages; ispore<nROI_greenSpores+nROI_macrophages; ispore++) {
							roiManager("Select", newArray(imacro,ispore));
							roiManager("AND");
							if(selectionType>-1){ //pathogen overlaps with macrophage
								if(!isElement(ispore - nROI_macrophages,alreadyCounted)){
									numberOfGreenSporesInsideMacrophages++;
									alreadyCounted[ispore - nROI_macrophages]= ispore - nROI_macrophages;
									roiManager("Draw");
									if(GreenSporesInside[ispore - nROI_macrophages] != 0){//== ispore - nROI_macrophages){
										 numberOfNotBlueGreenSporesInsideMacrophages++;
									}
								}
							}
						}
					}
					print("Number of unclassified all pathogens inside or touching host cells = " + numberOfGreenSporesInsideMacrophages);
					numberOfUnclassifiedGreenSporesPerMacrophage = numberOfGreenSporesInsideMacrophages / nROI_macrophages;
					print("Number of unclassified all pathogens per host cell (inside or adherent) = " + numberOfUnclassifiedGreenSporesPerMacrophage);
					numberOfClassifiedGreenSporesPerMacrophage = GreenSporesInsideCounter / nROI_macrophages;
					print("Number of classified phagocytosed all pathogens per host cell = " + numberOfClassifiedGreenSporesPerMacrophage);
					print("Total number of host cells = " + nROI_macrophages);
					numberOfAdherentGreenSpores = numberOfGreenSporesInsideMacrophages - numberOfNotBlueGreenSporesInsideMacrophages;//GreenSporesInsideCounter;
					numberOfAdherentGreenSporesPerMacrophage = numberOfAdherentGreenSpores / nROI_macrophages;
					print("Total number of adherent all pathogens = " + numberOfAdherentGreenSpores);
					print("Number of adherent classified pathogens per host cell = " + numberOfAdherentGreenSporesPerMacrophage);
				}
				//*** End of ROI overlap calculation for host cells and green pathogens ***


				//*** Calculate ROI overlap for host cells and blue pathogens *** 
				roiManager("reset");
				if(nROI_macrophages>0){
					if(useImageFilename==true){
						roiManager("open", dir + list[m] + "__HostCells_ROIs.zip");
					}
					else{
						roiManager("open",  dir + "/Image_" + testImageNumber + "__HostCells_ROIs.zip");
					}
					if(nROI_greenSporesOutside>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__AllPathogensOutside_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogensOutside_ROIs.zip");
						}
					}
					if(useImageFilename==true){
				   		mergedWindowName=list[m] + "__" + nameTagSavedFiles + "__AllComponents_Merged.jpg";
				   	}
				   	else{
				   		mergedWindowName="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AllComponents_Merged.jpg";
				   	}
					numberOfBlueSporesOverlappingMacrophages= 0;
					alreadyCounted = newArray(nROI_greenSporesOutside);//in order to avoid counting the same pathogen twice if it overlaps two host cells
					roiManager("Set Color", "red");
					roiManager("Set Line Width", 2);
					for (imacro=0; imacro<nROI_macrophages; imacro++) {
						for (ispore=nROI_macrophages; ispore<nROI_greenSporesOutside+nROI_macrophages; ispore++) {
							roiManager("Select", newArray(imacro,ispore));
							roiManager("AND");
							if(selectionType>-1){ //pathogen overlaps with macrophage
								if(!isElement(ispore - nROI_macrophages,alreadyCounted)){
									numberOfBlueSporesOverlappingMacrophages++;
									alreadyCounted[ispore - nROI_macrophages]= ispore - nROI_macrophages;
									roiManager("Draw");
								}
							}
						}
					}
					numberOfAdherentBlueSpores = numberOfBlueSporesOverlappingMacrophages;//numberOfNotBlueGreenSporesInsideMacrophages;//GreenSporesInsideCounter;
					numberOfAdherentGreenSporesPerMacrophage = numberOfAdherentBlueSpores / nROI_macrophages;
					print("Total number of adherent outside pathogens = " + numberOfAdherentBlueSpores);
					print("Number of classified adherent pathogens per host cell = " + numberOfAdherentGreenSporesPerMacrophage);
				}
				//*** End of ROI overlap calculation for host cells and blue pathogens ***

				//*** Calculate ROI overlap for host cells and pathogens, excluding edges *** 
				if(excludeEdges == true){
					//Count all green pathogens (w/ and w/o blue) and non-edge host cells:
					roiManager("reset");
					if(nROI_macrophages>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__HostCells_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__HostCells_ROIs.zip");
						}
						if(nROI_greenSpores>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__AllPathogens_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
							}
						}
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__Edges_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__Edges_ROIs.zip");
						}
						if(useImageFilename==true){
					   		mergedWindowNameNonEdge=list[m] + "__" + nameTagSavedFiles + "__AllComponents_NonEdge_Merged.jpg";
					   	}
					   	else{
					   		mergedWindowNameNonEdge="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AllComponents_NonEdge_Merged.jpg";
					   	}
						numberOfGreenSporesInsideMacrophages_excludingEdges = 0;
						numberOfMacrophages_excludingEdges = 0;
						alreadyCounted = newArray(nROI_greenSpores);//in order to avoid counting the same pathogen twice if it overlaps two host cells
						for (imacro=0; imacro<nROI_macrophages; imacro++) {
							roiManager("Select", newArray(imacro,nROI_greenSpores+nROI_macrophages+4));
							roiManager("AND");
							if(selectionType == (-1)){ //non-edge macrophage
								numberOfMacrophages_excludingEdges++;
								for (ispore=nROI_macrophages; ispore<nROI_greenSpores+nROI_macrophages; ispore++) {
									roiManager("Select", newArray(imacro,ispore));
									roiManager("AND");
									if(selectionType>-1){ //pathogen and host cell overlap
										if(!isElement(ispore - nROI_macrophages,alreadyCounted)){
											numberOfGreenSporesInsideMacrophages_excludingEdges++;
											alreadyCounted[ispore - nROI_macrophages]= ispore - nROI_macrophages;
											roiManager("Draw");
										}
									}	
								}
							}
						}
						print("Number of unclassified all pathogens inside or touching host cells excluding the edge = " + numberOfGreenSporesInsideMacrophages_excludingEdges);
						numberOfUnclassifiedGreenSporesPerMacrophage_excludingEdges = numberOfGreenSporesInsideMacrophages_excludingEdges / numberOfMacrophages_excludingEdges;
						print("Number of unclassified all pathogens per host cell excluding the edge (inside or adherent) = " + numberOfUnclassifiedGreenSporesPerMacrophage);
					}

					//Count green-only pathogens and non-edge host cells:
					roiManager("reset");
					if(nROI_macrophages>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__HostCells_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__HostCells_ROIs.zip");
						}
						if(nROI_greenSporesInside>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__AllPathogensInside_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogensInside_ROIs.zip");
							}
						}
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__Edges_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__Edges_ROIs.zip");
						}
						if(useImageFilename==true){
					   		mergedWindowNameNonEdge=list[m] + "__" + nameTagSavedFiles + "__AllComponents_NonEdge_Merged.jpg";
					   	}
					   	else{
					   		mergedWindowNameNonEdge="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AllComponents_NonEdge_Merged.jpg";
					   	}
						numberOfGreenOnlySporesInsideMacrophages_excludingEdges = 0;
						alreadyCounted = newArray(nROI_greenSporesInside);//in order to avoid counting the same pathogen twice if it overlaps two host cells
						for (imacro=0; imacro<nROI_macrophages; imacro++) {
							roiManager("Select", newArray(imacro,nROI_greenSporesInside+nROI_macrophages+4));
							roiManager("AND");
							if(selectionType == (-1)){ //non-edge macrophage
								for (ispore=nROI_macrophages; ispore<nROI_greenSporesInside+nROI_macrophages; ispore++) {
									roiManager("Select", newArray(imacro,ispore));
									roiManager("AND");
									if(selectionType>-1){ //pathogen and host cell overlap
										if(!isElement(ispore - nROI_macrophages,alreadyCounted)){
											numberOfGreenOnlySporesInsideMacrophages_excludingEdges++;
											alreadyCounted[ispore - nROI_macrophages]= ispore - nROI_macrophages;
											roiManager("Draw");
										}
									}	
								}
							}
						}
						numberOfClassifiedGreenOnlySporesPerMacrophage_excludingEdges = numberOfGreenOnlySporesInsideMacrophages_excludingEdges / numberOfMacrophages_excludingEdges;
						print("Number of classified phagocytosed pathogens excluding the edge = " + numberOfGreenOnlySporesInsideMacrophages_excludingEdges);
						print("Number of classified phagocytosed pathogens per host cell excluding the edge = " + numberOfClassifiedGreenOnlySporesPerMacrophage_excludingEdges);
					}

					//Count blue pathogens and non-edge host cells:
					roiManager("reset");
					if(nROI_macrophages>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__HostCells_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__HostCells_ROIs.zip");
						}
						if(nROI_greenSporesOutside>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__AllPathogensOutside_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogensOutside_ROIs.zip");
							}
						}
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__Edges_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__Edges_ROIs.zip");
						}
						if(useImageFilename==true){
					   		mergedWindowNameNonEdge=list[m] + "__" + nameTagSavedFiles + "__OutsidePathogens_NonEdge_Merged.jpg";
					   	}
					   	else{
					   		mergedWindowNameNonEdge="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__BlueSPores_NonEdge_Merged.jpg";
					   	}
						numberOfBlueSporesOverlappingMacrophages_excludingEdges = 0;
						alreadyCounted = newArray(nROI_greenSpores);//in order to avoid counting the same pathogen twice if it overlaps two host cells
						for (imacro=0; imacro<nROI_macrophages; imacro++) {
							roiManager("Select", newArray(imacro,nROI_greenSporesOutside+nROI_macrophages+4));
							roiManager("AND");
							if(selectionType == (-1)){ //non-edge macrophage
								for (ispore=nROI_macrophages; ispore<nROI_greenSporesOutside+nROI_macrophages; ispore++) {
									roiManager("Select", newArray(imacro,ispore));
									roiManager("AND");
									if(selectionType>-1){ //pathogen and host cell overlap
										if(!isElement(ispore - nROI_macrophages,alreadyCounted)){
											numberOfBlueSporesOverlappingMacrophages_excludingEdges++;
											alreadyCounted[ispore - nROI_macrophages]= ispore - nROI_macrophages;
											roiManager("Draw");
										}
									}	
								}
							}
						}
						numberOfBlueSporesPerMacrophage_excludingEdges = numberOfBlueSporesOverlappingMacrophages_excludingEdges / numberOfMacrophages_excludingEdges;
						print("Number of overlapping outside pathogens per host cell excluding the edge = " + numberOfBlueSporesPerMacrophage_excludingEdges);
						numberOfAdherentBlueSpores_excludingEdges = numberOfBlueSporesOverlappingMacrophages_excludingEdges;
						numberOfAdherentBlueSporesPerMacrophage_excludingEdges = numberOfAdherentBlueSpores_excludingEdges / numberOfMacrophages_excludingEdges;
						print("Total number of adherent outside pathogens excluding the edge = " + numberOfAdherentBlueSpores_excludingEdges);
						print("Number of adherent outside pathogens per host cell excluding the edge = " + numberOfAdherentBlueSporesPerMacrophage_excludingEdges);
					}

					//Calculate phagocytosis ratio:
					roiManager("reset");
					if(nROI_macrophages>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__HostCells_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__HostCells_ROIs.zip");
						}
						if(nROI_greenSporesInside>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__AllPathogensInside_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogensInside_ROIs.zip");
							}
						}
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__Edges_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__Edges_ROIs.zip");
						}
						
						numberOfInsideGreenSporesInsideMacrophages_excludingEdges = 0;
						alreadyCounted = newArray(nROI_greenSpores);//in order to avoid counting the same pathogen twice if it overlaps two host cells
						for (imacro=0; imacro<nROI_macrophages; imacro++) {
							roiManager("Select", newArray(imacro,nROI_greenSporesInside+nROI_macrophages+4));
							roiManager("AND");
							if(selectionType == (-1) && nROI_greenSporesInside>=1){ //non-edge macrophage
								for (ispore=nROI_macrophages; ispore<nROI_greenSporesInside+nROI_macrophages; ispore++) {
									roiManager("Select", newArray(imacro,ispore));
									roiManager("AND");
									if(selectionType>-1){
										if(!isElement(ispore - nROI_macrophages,alreadyCounted)){
											numberOfInsideGreenSporesInsideMacrophages_excludingEdges++;
											alreadyCounted[ispore - nROI_macrophages]= ispore - nROI_macrophages;
											roiManager("Draw");
										}
									}	
								}
							}
						}
						numberOfPhagocytosedGreenSporesPerMacrophage = nROI_greenSporesInside / numberOfMacrophages;
						numberOfPhagocytosedGreenSporesPerMacrophage_excludingEdges = numberOfInsideGreenSporesInsideMacrophages_excludingEdges / numberOfMacrophages_excludingEdges;
						print("Number of phagocytosed pathogens, total = " + nROI_greenSporesInside);
						print("Number of phagocytosed pathogens, excluding edges = " + numberOfInsideGreenSporesInsideMacrophages_excludingEdges);
						print("Number of phagocytosed pathogens per host cell, total = " + numberOfPhagocytosedGreenSporesPerMacrophage);
						print("Number of phagocytosed pathogens per host cell, excluding edges = " + numberOfPhagocytosedGreenSporesPerMacrophage_excludingEdges);
						phagocytosisRatio_excludingEdges = numberOfInsideGreenSporesInsideMacrophages_excludingEdges / (numberOfInsideGreenSporesInsideMacrophages_excludingEdges + numberOfAdherentBlueSpores_excludingEdges);
						print("Phagocytosis ratio, excluding edges = " + phagocytosisRatio_excludingEdges);
					}

					//Calculate uptake ratio:
					roiManager("reset");
					if(nROI_macrophages>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__HostCells_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__HostCells_ROIs.zip");
						}
						if(nROI_greenSporesInside>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__AllPathogensInside_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogensInside_ROIs.zip");
							}
						}
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__Edges_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__Edges_ROIs.zip");
						}
						
						numberOfMacrophagesWithPhagocytosedGreenSpores = 0;
						alreadyCounted = newArray(nROI_greenSpores);//in order to avoid counting the same pathogen twice if it overlaps two host cells
						for (imacro=0; imacro<nROI_macrophages; imacro++) {
							numberOfInsideGreenSpores_excludingEdges = 0;
							roiManager("Select", newArray(imacro,nROI_greenSporesInside+nROI_macrophages+4));
							roiManager("AND");
							if(selectionType == (-1) && nROI_greenSporesInside>=1){ //non-edge macrophage
								for (ispore=nROI_macrophages; ispore<nROI_greenSporesInside+nROI_macrophages; ispore++) {
									roiManager("Select", newArray(imacro,ispore));
									roiManager("AND");
									if(selectionType>-1){
										if(!isElement(ispore - nROI_macrophages,alreadyCounted)){
											numberOfInsideGreenSpores_excludingEdges++;
											alreadyCounted[ispore - nROI_macrophages]= ispore - nROI_macrophages;
											roiManager("Draw");
										}
									}	
								}
								if(numberOfInsideGreenSpores_excludingEdges >= 1){
									numberOfMacrophagesWithPhagocytosedGreenSpores++;
								}
							}
						}
						uptakeRatio_excludingEdges = numberOfMacrophagesWithPhagocytosedGreenSpores / numberOfMacrophages_excludingEdges; 
						print("M_phag and M_total, excl edgs = " + numberOfMacrophagesWithPhagocytosedGreenSpores + " ; " + numberOfMacrophages_excludingEdges);
						print("Uptake ratio, excluding edges = " + uptakeRatio_excludingEdges);
						phagocyticIndex_excludingEdges = uptakeRatio_excludingEdges * numberOfInsideGreenSporesInsideMacrophages_excludingEdges / numberOfMacrophages_excludingEdges;
						print("Phagocytic index, excluding edges = " + phagocyticIndex_excludingEdges);
						symmetrizedPhagocyticIndex_excludingEdges = uptakeRatio_excludingEdges * phagocytosisRatio_excludingEdges;
						print("Symmetrized phagocytic index, excluding edges = " + symmetrizedPhagocyticIndex_excludingEdges);

						//Save the results file in a spreadsheet-compatible format:
						if(saveFullPathImagename == true){
							saveImagename = fullImagename;
						}
						else {
							saveImagename = currentImagename;
						}
						print(saveAllFile, saveImagename + "	\t" + d2s(numberOfMacrophages,0) + "	\t" + d2s(nROI_greenSpores,0) + "	\t" + d2s(numberOfInsideGreenSporesInsideMacrophages_excludingEdges,0) + "	\t" + d2s(numberOfAdherentBlueSpores_excludingEdges,0) + "	\t" + d2s(numberOfMacrophagesWithPhagocytosedGreenSpores,0) + "	\t" + d2s(phagocytosisRatio_excludingEdges,3) + "	\t" + d2s(uptakeRatio_excludingEdges,3) + "	\t" + d2s(phagocyticIndex_excludingEdges,3) + "	\t" + d2s(symmetrizedPhagocyticIndex_excludingEdges,3) + "	\t"	
							 + d2s(pathogenAreaMean,3) + "	\t" + d2s(pathogenAreaSTD,3) + "	\t" + d2s(pathogenPerimMean,3) + "	\t" + d2s(pathogenPerimSTD,3) + "	\t" + d2s(pathogenARMean,3) + "	\t" + d2s(pathogenARSTD,3) + "	\t" + d2s(pathogenSolidityMean,3) + "	\t" + d2s(pathogenSoliditySTD,3) + "	\t" 
							 + d2s(hostAreaMean,3) + "	\t" + d2s(hostAreaSTD,3) + "	\t" + d2s(hostPerimMean,3) + "	\t" + d2s(hostPerimSTD,3) + "	\t" + d2s(hostARMean,3) + "	\t" + d2s(hostARSTD,3) + "	\t" + d2s(hostSolidityMean,3) + "	\t" + d2s(hostSoliditySTD,3) + "	\n");
					}
				} //end of if(excludeEdges == true)
				//*** End of ROI overlap calculation, excluding edges ***


//***************************************************************************************************************
//Secondly, use Ab-based red host cell ROIs to calculate phagocytic measures:
//***************************************************************************************************************


				//*** Saving images with labelled host cells begins: ***
				roiManager("reset");
				selectWindow("TL_IlluminationCorrected_saveForLater");
				if(nROI_redMacrophages>0){
					if(useImageFilename==true){
						roiManager("open", dir + list[m] + "__RedMacrophages_ROIs.zip");
					}
					else{
						roiManager("open",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs.zip");
					}
					nr_redmacrophages=roiManager("count");
					labeledHosts_01 = newArray(nr_redmacrophages);
					roiManager("Set Fill Color", "green");
					roiManager("Set Color", "yellow");
					roiManager("Set Line Width", 2);
					for (i=0; i<nr_redmacrophages; i++) {
						labeledHosts_01[i] = i;
						roiManager("select", labeledHosts_01);
						roiManager("Draw");
					}
				}
				run("Flatten");
				if(saveResults==true){
					if(useImageFilename==true){
				   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__Raw_image_with_labeledHosts.jpg");
				   	}
				   	else{
				   		saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Raw_image_with_labeledHosts.jpg");
				   	}
				}
				//*** Saving images with labelled host cells ends ***

				//*** Merge all components into one image using red host cells and save it ***
				if(useImageFilename==true){
			   		redMacrophageWindowName=list[m] + "__" + nameTagSavedFiles + "__Raw_image_with_labeledHosts.jpg";
			   		blueSporesWindowName=list[m] + "__" + nameTagSavedFiles + "SegmentedOutsidePathogens.jpg";
			   		greenSporesWindowName=list[m] + "__" + nameTagSavedFiles + "SegmentedAllPathogens.jpg";
			   		selectWindow(redMacrophageWindowName);
			   		selectWindow(blueSporesWindowName);
			   		run("Merge Channels...", "c3=" + blueSporesWindowName + " c4=" + redMacrophageWindowName + " keep");
					rename("RGB_labeledHosts");
					selectWindow(greenSporesWindowName);
					selectWindow("RGB_labeledHosts");
					run("Merge Channels...", "c2=" + greenSporesWindowName + " c4=RGB_labeledHosts keep");
			   	}
			   	else{
			   		redMacrophageWindowName="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Raw_image_with_labeledHosts.jpg";
			   		blueSporesWindowName="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "SegmentedOutsidePathogens.jpg";
			   		greenSporesWindowName="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "SegmentedAllPathogens.jpg";
			   		selectWindow(macrophageWindowName);
			   		selectWindow(blueSporesWindowName);
			   		run("Merge Channels...", "c3=" + blueSporesWindowName + " c4=" + redMacrophageWindowName + " keep");
					rename("RGB_labeledHosts");
					selectWindow(greenSporesWindowName);
					selectWindow("RGB_labeledHosts");
					run("Merge Channels...", "c2=" + greenSporesWindowName + " c4=RGB_labeledHosts keep");
			   	}
			   	if(saveResults==true){
					if(useImageFilename==true){
				   		saveAs("Jpeg", dir + list[m] + "__" + nameTagSavedFiles + "__AllComponents_labeledHosts_Merged.jpg");
				   	}
				   	else{
				   		saveAs("Jpeg", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AllComponents_labeledHosts_Merged.jpg");
				   	}
				}
				//*** End of image merging and saving for labelled host cells***


				if(labelledMacrophages){
					//*** Calculate ROI overlap for host cells and green pathogens *** 
					roiManager("reset");
					if(nROI_redMacrophages>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__RedMacrophages_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs.zip");
						}
						if(nROI_greenSpores>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__AllPathogens_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
							}
						}
						if(useImageFilename==true){
					   		mergedWindowName=list[m] + "__" + nameTagSavedFiles + "__AllComponents_labeledHosts_Merged.jpg";
					   	}
					   	else{
					   		mergedWindowName="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AllComponents_labeledHosts_Merged.jpg";
					   	}
						numberOfGreenSporesInsideMacrophages = 0;
						numberOfNotBlueGreenSporesInsideMacrophages = 0;
						alreadyCounted = newArray(nROI_greenSpores);//in order to avoid counting the same pathogen twice if it overlaps two host cells
						roiManager("Set Color", "red");
						roiManager("Set Line Width", 2);
						for (imacro=0; imacro<nROI_redMacrophages; imacro++) {
							for (ispore=nROI_redMacrophages; ispore<nROI_greenSpores+nROI_redMacrophages; ispore++) {
								roiManager("Select", newArray(imacro,ispore));
								roiManager("AND");
								if(selectionType>-1){ //pathogen overlaps with macrophage
									if(!isElement(ispore - nROI_redMacrophages,alreadyCounted)){
										numberOfGreenSporesInsideMacrophages++;
										alreadyCounted[ispore - nROI_redMacrophages]= ispore - nROI_redMacrophages;
										roiManager("Draw");
										if(GreenSporesInside[ispore - nROI_redMacrophages] != 0){//== ispore - nROI_redMacrophages){
											 numberOfNotBlueGreenSporesInsideMacrophages++;
										}
									}
								}
							}
						}
						print("Number of unclassified all pathogens inside or touching labelled host cells = " + numberOfGreenSporesInsideMacrophages);
						numberOfUnclassifiedGreenSporesPerMacrophage = numberOfGreenSporesInsideMacrophages / nROI_redMacrophages;
						print("Number of unclassified all pathogens per labelled host cell (inside or adherent) = " + numberOfUnclassifiedGreenSporesPerMacrophage);
						numberOfClassifiedGreenSporesPerMacrophage = GreenSporesInsideCounter / nROI_redMacrophages;
						print("Number of classified phagocytosed pathogens per labelled host cell = " + numberOfClassifiedGreenSporesPerMacrophage);
						print("Total number of labelled host cells = " + nROI_redMacrophages);
						numberOfAdherentGreenSpores = numberOfGreenSporesInsideMacrophages - numberOfNotBlueGreenSporesInsideMacrophages;//GreenSporesInsideCounter;
						numberOfAdherentGreenSporesPerMacrophage = numberOfAdherentGreenSpores / nROI_redMacrophages;
						print("Total number of adherent all pathogens = " + numberOfAdherentGreenSpores);
						print("Number of classified adherent pathogens per labelled host cell = " + numberOfAdherentGreenSporesPerMacrophage);
					}
					//*** End of ROI overlap calculation for host cells and green pathogens ***
	
	
					//*** Calculate ROI overlap for host cells and blue pathogens *** 
					roiManager("reset");
					if(nROI_redMacrophages>0){
						if(useImageFilename==true){
							roiManager("open", dir + list[m] + "__RedMacrophages_ROIs.zip");
						}
						else{
							roiManager("open",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs.zip");
						}
						if(nROI_greenSporesOutside>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__AllPathogensOutside_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogensOutside_ROIs.zip");
							}
						}
						if(useImageFilename==true){
					   		mergedWindowName=list[m] + "__" + nameTagSavedFiles + "__AllComponents_RedMacrophage_Merged.jpg";
					   	}
					   	else{
					   		mergedWindowName="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AllComponents_RedMacrophage_Merged.jpg";
					   	}
						numberOfBlueSporesOverlappingMacrophages= 0;
						alreadyCounted = newArray(nROI_greenSporesOutside);//in order to avoid counting the same pathogen twice if it overlaps two host cells
						roiManager("Set Color", "red");
						roiManager("Set Line Width", 2);
						for (imacro=0; imacro<nROI_redMacrophages; imacro++) {
							for (ispore=nROI_redMacrophages; ispore<nROI_greenSporesOutside+nROI_redMacrophages; ispore++) {
								roiManager("Select", newArray(imacro,ispore));
								roiManager("AND");
								if(selectionType>-1){ //pathogen overlaps with macrophage
									if(!isElement(ispore - nROI_redMacrophages,alreadyCounted)){
										numberOfBlueSporesOverlappingMacrophages++;
										alreadyCounted[ispore - nROI_redMacrophages]= ispore - nROI_redMacrophages;
										roiManager("Draw");
									}
								}
							}
						}
						numberOfAdherentBlueSpores = numberOfBlueSporesOverlappingMacrophages;//numberOfNotBlueGreenSporesInsideMacrophages;//GreenSporesInsideCounter;
						numberOfAdherentGreenSporesPerMacrophage = numberOfAdherentBlueSpores / nROI_redMacrophages;
						print("Total number of adherent outside pathogens = " + numberOfAdherentBlueSpores);
						print("Number of classified adherent pathogens per host cell = " + numberOfAdherentGreenSporesPerMacrophage);
					}
					//*** End of ROI overlap calculation for host cells and blue pathogens ***
	
					//*** Calculate ROI overlap for host cells and pathogens, excluding edges *** 
					if(excludeEdges == true){
						//Count all green pathogens (w/ and w/o blue) and non-edge host cells:
						roiManager("reset");
						if(nROI_redMacrophages>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__RedMacrophages_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs.zip");
							}
							if(nROI_greenSpores>0){
								if(useImageFilename==true){
									roiManager("open", dir + list[m] + "__AllPathogens_ROIs.zip");
								}
								else{
									roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogens_ROIs.zip");
								}
							}
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__Edges_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__Edges_ROIs.zip");
							}
							if(useImageFilename==true){
						   		mergedWindowNameNonEdge=list[m] + "__" + nameTagSavedFiles + "__AllComponents_RedMacrophage_NonEdge_Merged.jpg";
						   	}
						   	else{
						   		mergedWindowNameNonEdge="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AllComponents_RedMacrophage_NonEdge_Merged.jpg";
						   	}
							numberOfGreenSporesInsideMacrophages_excludingEdges = 0;
							numberOfMacrophages_excludingEdges = 0;
							alreadyCounted = newArray(nROI_greenSpores);//in order to avoid counting the same pathogen twice if it overlaps two host cells
							for (imacro=0; imacro<nROI_redMacrophages; imacro++) {
								roiManager("Select", newArray(imacro,nROI_greenSpores+nROI_redMacrophages+4));
								roiManager("AND");
								if(selectionType == (-1)){ //non-edge macrophage
									numberOfMacrophages_excludingEdges++;
									for (ispore=nROI_redMacrophages; ispore<nROI_greenSpores+nROI_redMacrophages; ispore++) {
										roiManager("Select", newArray(imacro,ispore));
										roiManager("AND");
										if(selectionType>-1){ //pathogen and host cell overlap
											if(!isElement(ispore - nROI_redMacrophages,alreadyCounted)){
												numberOfGreenSporesInsideMacrophages_excludingEdges++;
												alreadyCounted[ispore - nROI_redMacrophages]= ispore - nROI_redMacrophages;
												roiManager("Draw");
											}
										}	
									}
								}
							}
							print("Number of unclassified all pathogens inside or touching labelled host cells excluding the edge = " + numberOfGreenSporesInsideMacrophages_excludingEdges);
							numberOfUnclassifiedGreenSporesPerMacrophage_excludingEdges = numberOfGreenSporesInsideMacrophages_excludingEdges / numberOfMacrophages_excludingEdges;
							print("Number of unclassified all pathogens per labelled host cell excluding the edge (inside or adherent) = " + numberOfUnclassifiedGreenSporesPerMacrophage);
						}
	
						//Count green-only pathogens and non-edge host cells:
						roiManager("reset");
						if(nROI_redMacrophages>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__RedMacrophages_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs.zip");
							}
							if(nROI_greenSporesInside>0){
								if(useImageFilename==true){
									roiManager("open", dir + list[m] + "__AllPathogensInside_ROIs.zip");
								}
								else{
									roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogensInside_ROIs.zip");
								}
							}
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__Edges_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__Edges_ROIs.zip");
							}
							if(useImageFilename==true){
						   		mergedWindowNameNonEdge=list[m] + "__" + nameTagSavedFiles + "__AllComponents_RedMacrophage_NonEdge_Merged.jpg";
						   	}
						   	else{
						   		mergedWindowNameNonEdge="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__AllComponents_RedMacrophage_NonEdge_Merged.jpg";
						   	}
							numberOfGreenOnlySporesInsideMacrophages_excludingEdges = 0;
							alreadyCounted = newArray(nROI_greenSporesInside);//in order to avoid counting the same pathogen twice if it overlaps two host cells
							for (imacro=0; imacro<nROI_redMacrophages; imacro++) {
								roiManager("Select", newArray(imacro,nROI_greenSporesInside+nROI_redMacrophages+4));
								roiManager("AND");
								if(selectionType == (-1)){ //non-edge macrophage
									for (ispore=nROI_redMacrophages; ispore<nROI_greenSporesInside+nROI_redMacrophages; ispore++) {
										roiManager("Select", newArray(imacro,ispore));
										roiManager("AND");
										if(selectionType>-1){ //pathogen and host cell overlap
											if(!isElement(ispore - nROI_redMacrophages,alreadyCounted)){
												numberOfGreenOnlySporesInsideMacrophages_excludingEdges++;
												alreadyCounted[ispore - nROI_redMacrophages]= ispore - nROI_redMacrophages;
												roiManager("Draw");
											}
										}	
									}
								}
							}
							numberOfClassifiedGreenOnlySporesPerMacrophage_excludingEdges = numberOfGreenOnlySporesInsideMacrophages_excludingEdges / numberOfMacrophages_excludingEdges;
							print("Number of classified phagocytosed pathogens excluding the edge = " + numberOfGreenOnlySporesInsideMacrophages_excludingEdges);
							print("Number of classified phagocytosed pathogens per labelled host cell excluding the edge = " + numberOfClassifiedGreenOnlySporesPerMacrophage_excludingEdges);
						}
	
						//Count blue pathogens and non-edge host cells:
						roiManager("reset");
						if(nROI_redMacrophages>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__RedMacrophages_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs.zip");
							}
							if(nROI_greenSporesOutside>0){
								if(useImageFilename==true){
									roiManager("open", dir + list[m] + "__AllPathogensOutside_ROIs.zip");
								}
								else{
									roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogensOutside_ROIs.zip");
								}
							}
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__Edges_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__Edges_ROIs.zip");
							}
							if(useImageFilename==true){
						   		mergedWindowNameNonEdge=list[m] + "__" + nameTagSavedFiles + "__OutsidePathogens_NonEdge_Merged.jpg";
						   	}
						   	else{
						   		mergedWindowNameNonEdge="/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__OutsidePathogens_NonEdge_Merged.jpg";
						   	}
							numberOfBlueSporesOverlappingMacrophages_excludingEdges = 0;
							alreadyCounted = newArray(nROI_greenSpores);//in order to avoid counting the same pathogen twice if it overlaps two host cells
							for (imacro=0; imacro<nROI_redMacrophages; imacro++) {
								roiManager("Select", newArray(imacro,nROI_greenSporesOutside+nROI_redMacrophages+4));
								roiManager("AND");
								if(selectionType == (-1)){ //non-edge macrophage
									for (ispore=nROI_redMacrophages; ispore<nROI_greenSporesOutside+nROI_redMacrophages; ispore++) {
										roiManager("Select", newArray(imacro,ispore));
										roiManager("AND");
										if(selectionType>-1){ //pathogen and host cell overlap
											if(!isElement(ispore - nROI_redMacrophages,alreadyCounted)){
												numberOfBlueSporesOverlappingMacrophages_excludingEdges++;
												alreadyCounted[ispore - nROI_redMacrophages]= ispore - nROI_redMacrophages;
												roiManager("Draw");
											}
										}	
									}
								}
							}
							numberOfBlueSporesPerMacrophage_excludingEdges = numberOfBlueSporesOverlappingMacrophages_excludingEdges / numberOfMacrophages_excludingEdges;
							print("Number of overlapping outside pathogens per labelled host cell excluding the edge = " + numberOfBlueSporesPerMacrophage_excludingEdges);
							numberOfAdherentBlueSpores_excludingEdges = numberOfBlueSporesOverlappingMacrophages_excludingEdges;
							numberOfAdherentBlueSporesPerMacrophage_excludingEdges = numberOfAdherentBlueSpores_excludingEdges / numberOfMacrophages_excludingEdges;
							print("Total number of adherent outside pathogens excluding the edge = " + numberOfAdherentBlueSpores_excludingEdges);
							print("Number of adherent outside pathogens per red host cell excluding the edge = " + numberOfAdherentBlueSporesPerMacrophage_excludingEdges);
						}
	
						//Calculate phagocytosis ratio:
						roiManager("reset");
						if(nROI_redMacrophages>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__RedMacrophages_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs.zip");
							}
							if(nROI_greenSporesInside>0){
								if(useImageFilename==true){
									roiManager("open", dir + list[m] + "__AllPathogensInside_ROIs.zip");
								}
								else{
									roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogensInside_ROIs.zip");
								}
							}
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__Edges_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__Edges_ROIs.zip");
							}
							
							numberOfInsideGreenSporesInsideMacrophages_excludingEdges = 0;
							alreadyCounted = newArray(nROI_greenSpores);//in order to avoid counting the same pathogen twice if it overlaps two host cells
							for (imacro=0; imacro<nROI_redMacrophages; imacro++) {
								roiManager("Select", newArray(imacro,nROI_greenSporesInside+nROI_redMacrophages+4));
								roiManager("AND");
								if(selectionType == (-1) && nROI_greenSporesInside>=1){ //non-edge macrophage
									for (ispore=nROI_redMacrophages; ispore<nROI_greenSporesInside+nROI_redMacrophages; ispore++) {
										roiManager("Select", newArray(imacro,ispore));
										roiManager("AND");
										if(selectionType>-1){
											if(!isElement(ispore - nROI_redMacrophages,alreadyCounted)){
												numberOfInsideGreenSporesInsideMacrophages_excludingEdges++;
												alreadyCounted[ispore - nROI_redMacrophages]= ispore - nROI_redMacrophages;
												roiManager("Draw");
											}
										}	
									}
								}
							}
							numberOfPhagocytosedGreenSporesPerMacrophage = nROI_greenSporesInside / numberOfMacrophages;
							numberOfPhagocytosedGreenSporesPerMacrophage_excludingEdges = numberOfInsideGreenSporesInsideMacrophages_excludingEdges / numberOfMacrophages_excludingEdges;
							print("Number of phagocytosed pathogens with labelled host cells, total = " + nROI_greenSporesInside);
							print("Number of phagocytosed pathogens with labelled host cells, excluding edges = " + numberOfInsideGreenSporesInsideMacrophages_excludingEdges);
							print("Number of phagocytosed pathogens per labelled host, total = " + numberOfPhagocytosedGreenSporesPerMacrophage);
							print("Number of phagocytosed pathogens per labelled host, excluding edges = " + numberOfPhagocytosedGreenSporesPerMacrophage_excludingEdges);
							phagocytosisRatio_excludingEdges = numberOfInsideGreenSporesInsideMacrophages_excludingEdges / (numberOfInsideGreenSporesInsideMacrophages_excludingEdges + numberOfAdherentBlueSpores_excludingEdges);
							print("Phagocytosis ratio with labelled host cells, excluding edges = " + phagocytosisRatio_excludingEdges);
						}
	
						//Calculate uptake ratio:
						roiManager("reset");
						if(nROI_redMacrophages>0){
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__RedMacrophages_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__RedMacrophages_ROIs.zip");
							}
							if(nROI_greenSporesInside>0){
								if(useImageFilename==true){
									roiManager("open", dir + list[m] + "__AllPathogensInside_ROIs.zip");
								}
								else{
									roiManager("open",  dir + "/Image_" + testImageNumber + "__AllPathogensInside_ROIs.zip");
								}
							}
							if(useImageFilename==true){
								roiManager("open", dir + list[m] + "__Edges_ROIs.zip");
							}
							else{
								roiManager("open",  dir + "/Image_" + testImageNumber + "__Edges_ROIs.zip");
							}
							
							numberOfMacrophagesWithPhagocytosedGreenSpores = 0;
							alreadyCounted = newArray(nROI_greenSpores);//in order to avoid counting the same pathogen twice if it overlaps two host cells
							for (imacro=0; imacro<nROI_redMacrophages; imacro++) {
								numberOfInsideGreenSpores_excludingEdges = 0;
								roiManager("Select", newArray(imacro,nROI_greenSporesInside+nROI_redMacrophages+4));
								roiManager("AND");
								if(selectionType == (-1) && nROI_greenSporesInside>=1){ //non-edge macrophage
									for (ispore=nROI_redMacrophages; ispore<nROI_greenSporesInside+nROI_redMacrophages; ispore++) {
										roiManager("Select", newArray(imacro,ispore));
										roiManager("AND");
										if(selectionType>-1){
											if(!isElement(ispore - nROI_redMacrophages,alreadyCounted)){
												numberOfInsideGreenSpores_excludingEdges++;
												alreadyCounted[ispore - nROI_redMacrophages]= ispore - nROI_redMacrophages;
												roiManager("Draw");
											}
										}	
									}
									if(numberOfInsideGreenSpores_excludingEdges >= 1){
										numberOfMacrophagesWithPhagocytosedGreenSpores++;
									}
								}
							}
							uptakeRatio_excludingEdges = numberOfMacrophagesWithPhagocytosedGreenSpores / numberOfMacrophages_excludingEdges; 
							print("M_phag and M_total with labelled host cells, excl edgs = " + numberOfMacrophagesWithPhagocytosedGreenSpores + " ; " + numberOfMacrophages_excludingEdges);
							print("Uptake ratio with labelled host cells, excluding edges = " + uptakeRatio_excludingEdges);
							phagocyticIndex_excludingEdges = uptakeRatio_excludingEdges * numberOfInsideGreenSporesInsideMacrophages_excludingEdges / numberOfMacrophages_excludingEdges;
							print("Phagocytic index with labelled host cells, excluding edges = " + phagocyticIndex_excludingEdges);
							symmetrizedPhagocyticIndex_excludingEdges = uptakeRatio_excludingEdges * phagocytosisRatio_excludingEdges;
							print("Symmetrized phagocytic index with labelled host cells, excluding edges = " + symmetrizedPhagocyticIndex_excludingEdges);
	
							//Save the results file in a spreadsheet-compatible format:
							if(saveFullPathImagename == true){
								saveImagename =  "RedMacrophage__" + fullImagename;
							}
							else {
								saveImagename =  "RedMacrophage__" + currentImagename;
							}
							print(saveAllFile, saveImagename + "	\t" + d2s(numberOfRedMacrophages,0) + "	\t" + d2s(nROI_greenSpores,0) + "	\t" + d2s(numberOfInsideGreenSporesInsideMacrophages_excludingEdges,0) + "	\t" + d2s(numberOfAdherentBlueSpores_excludingEdges,0) + "	\t" + d2s(numberOfMacrophagesWithPhagocytosedGreenSpores,0) + "	\t" + d2s(phagocytosisRatio_excludingEdges,3) + "	\t" + d2s(uptakeRatio_excludingEdges,3) + "	\t" + d2s(phagocyticIndex_excludingEdges,3) + "	\t" + d2s(symmetrizedPhagocyticIndex_excludingEdges,3) + "	\t"	
							 + d2s(pathogenAreaMean,3) + "	\t" + d2s(pathogenAreaSTD,3) + "	\t" + d2s(pathogenPerimMean,3) + "	\t" + d2s(pathogenPerimSTD,3) + "	\t" + d2s(pathogenARMean,3) + "	\t" + d2s(pathogenARSTD,3) + "	\t" + d2s(pathogenSolidityMean,3) + "	\t" + d2s(pathogenSoliditySTD,3) + " \t"
							 + d2s(hostFlscAreaMean,3) + "	\t" + d2s(hostFlscAreaSTD,3) + "	\t" + d2s(hostFlscPerimMean,3) + "	\t" + d2s(hostFlscPerimSTD,3) + "	\t" + d2s(hostFlscARMean,3) + "	\t" + d2s(hostFlscARSTD,3) + "	\t" + d2s(hostFlscSolidityMean,3) + "	\t" + d2s(hostFlscSoliditySTD,3) + "	\n");
						}
					} //end of if(excludeEdges == true)
					//*** End of ROI overlap calculation, excluding edges ***
				}//*** End of "labelledMacrophages==true" section

				if(closeAllWindows==true){
					run("Close All");
				}		
			} //end of image type 'if' block
		} //end of 'if' block that finds specific image file names
	} //end of main 'for' loop that goes thru all images in a folder
} //end of loop for subfolders 5.10.16

File.close(saveAllFile); 

//**** Save parameters from GUI:
print("\n");
print("################################ PARAMETERS ################################");
print("Version: ACAQ_" + versionNumber);
print("FAST mode? (CLAHE calculations) " + fastModeCLAHE);
print("FAST mode? (NO plots) " + fastModeNoPlots);
print("Re-stiching Zeiss tilescan? " +  rearrangeZeiss);
print("Read multi-channel TIFF? " +  useMultichannelTiffFormat);
print("Use MAX Hessian instead? (default: MIN) " +  useMAXHessian);
print("Do you want to save your results? " +  saveResults);
print("Use subfolders? " +  useSubfolders);
print("Analyse all images? " +  analyseAllImages);
print("Close all windows at end? " +  closeAllWindows);
print("Run it in Batch Mode? " +  runInBatchMode);
print("Illumination correction for TL images? " +  correctIlluminationForTL);
print("Illumination correction for red images? " +  correctIlluminationForRed);
print("Line thickness for ROIs: " + lineThicknessForObjects);
print("X tile number: " + XTileNumber);
print("Y tile number: " + YTileNumber);
print("Gaussian sigma for illumination correction: " + illuminationCorrectionSigmaForTL);
print("Gaussian sigma for red illumination correction: " + illuminationCorrectionSigmaForRed);
print("Hessian smoothing factor for host cells: " + hessianSmoothingForMacrophages);
print("Hessian smoothing factor for pathogens: " + hessianSmoothingForSpores);
print("Gaussian smoothing factor for host cells: " + gaussianSmoothingForMacrophages);
print("Gaussian smoothing factor for red host cells: " + gaussianSmoothingForRedMacrophages);
print("Gaussian smoothing factor for TL pathogens: " + gaussianSmoothingForSpores);
print("Dilate/erode steps for host cells: " + dilateErodeStepsForMacrophages);
print("Additional erosion steps for host cells: " + additionalErodeStepsForMacrophages);
print("Dilate/erode steps for pathogens: " + dilateErodeStepsForSpores);
print("Gaussian smoothing factor for fluorescence pathogens: " + gaussianBlurForSpores);
print("Lower threshold multiplier for pathogens (only for Internal Gradient): " + lowerThresholdMultiplier);
print("Upper threshold multiplier for pathogens (only for Internal Gradient): " + upperThresholdMultiplier);
print("Lower threshold multiplier for host cells: " + lowerThresholdMultiplierMacrophages);
print("Upper threshold multiplier for host cells: " + upperThresholdMultiplierMacrophages);

print("Labelled host cells? " +  labelledMacrophages);
print("Labelled pathogens? " +  labelledSpores);
print("Exclude edges? " +  excludeEdges);
print("Gather results for entire image set? " +  gatherResultsForEntireImageSet);
print("Use image file name for results? " +  useImageFilename);
print("Use specific image groups based on filename? " +  useSpecificImageGroups);
print("Save full path of image name in output file? " +  saveFullPathImagename);
print("Image type: " +  imageType);
print("Search term #1: " +  searchTerm_1);
print("Search term #2: " +  searchTerm_2);
print("Search term #3: " +  searchTerm_3);
print("Exclude term #1: " +  excludeTerm_1);
print("Exclude term #2: " +  excludeTerm_2);
print("Exclude term #3: " +  excludeTerm_3);
print("Name tag of saved files: " +  nameTagSavedFiles);
print("Background image number: " + backgroundImageNumber);
print("First image number? " + firstImageNumber);
print("Last image number? " + lastImageNumber);
print("Test image number? " + testImageNumber);
print("Edge width (for Exclude Edges option)? " + edgeWidth);
print("Green channel multiplier for dim image data (try 2): " + greenChannelMultiplierForDimData);
print("Blue channel multiplier for dim image data (left at 1): " + blueChannelMultiplierForDimData);
print("Red channel multiplier for dim image data (left at 1): " + redChannelMultiplierForDimData);

print("Watershed on host cells? " +  watershedOnMacrophages);
print("Watershed on pathogens? " +  watershedOnSpores);
print("CLAHE on TL? " +  applyCLAHE);
print("CLAHE on FLSC? " +  applyCLAHEonFLSC);
print("Combine Internal Gradient and Bright Spot results? " +  combineInternalGradientAndBrightSpotResults);
print("Internal Gradient radius for green pathogens (try 2; use 0 for No): " +  internalGradientRadiusGreen);
print("Internal Gradient radius for blue pathogens (try 2; use 0 for No): " +  internalGradientRadiusBlue);
print("Internal Gradient radius for red host cells (try 5; use 0 for No): " +  internalGradientRadiusRed);
print("Rolling ball radius for pathogens background (try 20): " + rollingBallRadius);
print("Rolling ball radius for red host cells background (try 20): " + rollingBallRadiusRedMacrophages);
print("Local threshold radius for red host cells background (try b/w 30 and 40): " + localThresholdRadius);
print("Remove outliers for host cell ROI +  step 1 (try 10): " +  removeOutliersStep1);
print("Remove outliers for host cell ROI +  step 2 (try 20): " +  removeOutliersStep2);
print("CLAHE blocks: " + CLAHEblocks);
print("CLAHE bins: " + CLAHEbins);
print("CLAHE max slope:" + CLAHEslope);
print("1st CLAHE max slope for labelled host cells: " + CLAHEslope1FlscMacrophages);
print("2nd CLAHE max slope for labelled host cells: " + CLAHEslope2FlscMacrophages);
print("1st dilation steps for labelled host cells: " + dilationsSteps1FlscMacrophages);
print("2nd dilation steps for labelled host cells: " + dilationsSteps2FlscMacrophages);

print("Min host cell size: " + minMacrophageSize);
print("Max host cell size: " + maxMacrophageSize);
print("Min host cell circularity: " + minMacrophageCircularity);
print("Max host cell circularity: " + maxMacrophageCircularity);
print("Min pathogen size: " + minSporeSize);
print("Max pathogen size: " + maxSporeSize);
print("Min pathogen circularity: " + minSporeCircularity);
print("Max pathogen circularity: " + maxSporeCircularity);
print("Blue threshold for green pathogen classifier (inside vs. outside): " + blueThresholdForGreenSporeClassifier);
print("Green threshold for Bright Spots green pathogen classifier: " + greenThresholdForGreenSporeClassifier); 
print("Threshold method for Hessian: " +  thresholdMethodHessian);
print("Threshold method for green fluorescence (NOT VALID for BS when combining IG and BS +  see next button): " +  thresholdMethodGreenFluorescence);
print("Threshold method for green fluorescence Bright Spots method (only when combining IG and BS +  otherwise see previous button): " +  thresholdMethodGreenFluorescenceBrightSpots);
print("Threshold method for blue fluorescence: " +  thresholdMethodBlueFluorescence);
print("Local threshold method for red fluorescence: " +  localThresholdMethodRedFluorescence);
print("Threshold method for red fluorescence: " +  thresholdMethodRedFluorescence);

print("Use GUI for channel info? " +  useGUIforChannelInfo);
print("Channel 0: " +  channel00);
print("Channel 1: " +  channel01);
print("Channel 2: " +  channel02);
print("Channel 3: " +  channel03);

selectWindow("Log");

if(saveResults==true){
	if(useImageFilename==true){
   		saveAs("Text", dir + currentImagename + "__" + nameTagSavedFiles +  "__Notes_and_Parameter_settings.txt");
   	}
   	else{
   		saveAs("Text", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Notes_and_Parameter_settings.txt");
   	}
}
//**** End of notes and parameter saving

//*** Now save the parameters separately, for later re-use:
print("\\Clear");
print("################################ PARAMETERS ################################");
print("Version: ACAQ_" + versionNumber);
print("FAST mode? (CLAHE calculations) " + fastModeCLAHE);
print("FAST mode? (without plots) " + fastModeNoPlots);
print("Re-stiching Zeiss tilescan? " +  rearrangeZeiss);
print("Read multi-channel TIFF? " +  useMultichannelTiffFormat);
print("Use MAX Hessian instead? (default: MIN) " +  useMAXHessian);
print("Do you want to save your results? " +  saveResults);
print("Use subfolders? " +  useSubfolders);
print("Analyse all images? " +  analyseAllImages);
print("Close all windows at end? " +  closeAllWindows);
print("Run it in Batch Mode? " +  runInBatchMode);
print("Illumination correction for unlabelled images? " +  correctIlluminationForTL);
print("Illumination correction for labelled images? " +  correctIlluminationForRed);
print("Line thickness for ROIs: " + lineThicknessForObjects);
print("X tile number: " + XTileNumber);
print("Y tile number: " + YTileNumber);
print("Gaussian sigma for illumination correction: " + illuminationCorrectionSigmaForTL);
print("Gaussian sigma for labelled illumination correction: " + illuminationCorrectionSigmaForRed);
print("Hessian smoothing factor for host cells: " + hessianSmoothingForMacrophages);
print("Hessian smoothing factor for pathogens: " + hessianSmoothingForSpores);
print("Gaussian smoothing factor for host cells: " + gaussianSmoothingForMacrophages);
print("Gaussian smoothing factor for labelled host cells: " + gaussianSmoothingForRedMacrophages);
print("Gaussian smoothing factor for unlabelled pathogens: " + gaussianSmoothingForSpores);
print("Dilate/erode steps for host cells: " + dilateErodeStepsForMacrophages);
print("Additional erosion steps for host cells: " + additionalErodeStepsForMacrophages);
print("Dilate/erode steps for pathogens: " + dilateErodeStepsForSpores);
print("Gaussian smoothing factor for labelled pathogens: " + gaussianBlurForSpores);
print("Lower threshold multiplier for pathogens (only for Internal Gradient): " + lowerThresholdMultiplier);
print("Upper threshold multiplier for pathogens (only for Internal Gradient): " + upperThresholdMultiplier);
print("Lower threshold multiplier for host cells: " + lowerThresholdMultiplierMacrophages);
print("Upper threshold multiplier for host cells: " + upperThresholdMultiplierMacrophages);

print("Labelled host cells? " +  labelledMacrophages);
print("Labelled pathogens? " +  labelledSpores);
print("Exclude edges? " +  excludeEdges);
print("Gather results for entire image set? " +  gatherResultsForEntireImageSet);
print("Use image file name for results? " +  useImageFilename);
print("Use specific image groups based on filename? " +  useSpecificImageGroups);
print("Save full path of image name in output file? " +  saveFullPathImagename);
print("Image type: " +  imageType);
print("Search term #1: " +  searchTerm_1);
print("Search term #2: " +  searchTerm_2);
print("Search term #3: " +  searchTerm_3);
print("Exclude term #1: " +  excludeTerm_1);
print("Exclude term #2: " +  excludeTerm_2);
print("Exclude term #3: " +  excludeTerm_3);
print("Name tag of saved files: " +  nameTagSavedFiles);
print("Background image number: " + backgroundImageNumber);
print("First image number: " + firstImageNumber);
print("Last image number: " + lastImageNumber);
print("Test image number: " + testImageNumber);
print("Edge width (for Exclude Edges option): " + edgeWidth);
print("Labelled pathogen channel multiplier for dim image data: " + greenChannelMultiplierForDimData);
print("Outside pathogen channel multiplier for dim image data: " + blueChannelMultiplierForDimData);
print("Labelled host cells channel multiplier for dim image data: " + redChannelMultiplierForDimData);

print("Watershed on host cells? " +  watershedOnMacrophages);
print("Watershed on pathogens? " +  watershedOnSpores);
print("CLAHE on unlabelled images? " +  applyCLAHE);
print("CLAHE on labelled images? " +  applyCLAHEonFLSC);
print("Combine Internal Gradient and Bright Spot results for pathogens? " +  combineInternalGradientAndBrightSpotResults);
print("Internal Gradient radius for labelled pathogens (try 2; use 0 for No): " +  internalGradientRadiusGreen);
print("Internal Gradient radius for outside pathogens (try 2; use 0 for No): " +  internalGradientRadiusBlue);
print("Internal Gradient radius for labelled host cells (try 5; use 0 for No): " +  internalGradientRadiusRed);
print("Rolling ball radius for pathogens background (try 20): " + rollingBallRadius);
print("Rolling ball radius for labelled host cells background (try 20): " + rollingBallRadiusRedMacrophages);
print("Local threshold radius for labelled host cells (try b/w 30 and 40): " + localThresholdRadius);
print("Remove outliers for host cell ROI  step 1 (try 10): " +  removeOutliersStep1);
print("Remove outliers for host cell ROI  step 2 (try 20): " +  removeOutliersStep2);
print("CLAHE blocks: " + CLAHEblocks);
print("CLAHE bins: " + CLAHEbins);
print("CLAHE max slope:" + CLAHEslope);
print("1st CLAHE max slope for labelled host cells: " + CLAHEslope1FlscMacrophages);
print("2nd CLAHE max slope for labelled host cells: " + CLAHEslope2FlscMacrophages);
print("1st dilation steps for labelled host cells: " + dilationsSteps1FlscMacrophages);
print("2nd dilation steps for labelled host cells: " + dilationsSteps2FlscMacrophages);

print("Min host cell size: " + minMacrophageSize);
print("Max host cell size: " + maxMacrophageSize);
print("Min host cell circularity: " + minMacrophageCircularity);
print("Max host cell circularity: " + maxMacrophageCircularity);
print("Min pathogen size: " + minSporeSize);
print("Max pathogen size: " + maxSporeSize);
print("Min pathogen circularity: " + minSporeCircularity);
print("Max pathogen circularity: " + maxSporeCircularity);
print("Outside pathogens threshold for all pathogen classifier (separate inside from outside): " + blueThresholdForGreenSporeClassifier);
print("All pathogens threshold for Bright Spots pathogen classifier: " + greenThresholdForGreenSporeClassifier); 
print("Threshold method for Hessian: " +  thresholdMethodHessian);
print("Threshold method for all pathogens fluorescence (NOT VALID for BS when combining IG and BS +  see next button): " +  thresholdMethodGreenFluorescence);
print("Threshold method for all pathogens fluorescence Bright Spots method (only when combining IG and BS +  otherwise see previous button): " +  thresholdMethodGreenFluorescenceBrightSpots);
print("Threshold method for outside pathogens fluorescence: " +  thresholdMethodBlueFluorescence);
print("Local threshold method for labelled hosts fluorescence: " +  localThresholdMethodRedFluorescence);
print("Threshold method for labelled hosts fluorescence: " +  thresholdMethodRedFluorescence);

print("Use GUI for channel info? " +  useGUIforChannelInfo);
print("Channel 0: " +  channel00_new);
print("Channel 1: " +  channel01_new);
print("Channel 2: " +  channel02_new);
print("Channel 3: " +  channel03_new);

print("Segmentation method for unlabelled pathogens: " + segmentationMethodForUnlabeledPathogens);
print("Hessian smoothing for unlabelled pathogens: " + smoothingHessianSpores);
print("Internal gradient radius for unlabelled pathogens: " + gradientRadiusHessianSpores);
print("Hough circles minimum radius: " + minimumRadiusHoughHessianSpores);
print("Hough circles maximum radius: " + maximumRadiusHoughHessianSpores);
print("Hough circles increment radius: " + incrementRadiusHoughHessianSpores);
print("Hough minimum number of circles: " + minimumNumberHoughHessianSpores);
print("Hough maximum number of circles: " + maximumNumberHoughHessianSpores);
print("Hough circles threshold: " + thresholdHoughHessianSpores);
print("Hough circles resolution: " + resolutionHoughHessianSpores);
print("Hough circles ratio: " + ratioHoughHessianSpores);
print("Hough circles bandwidth: " + bandwidthRadiusHoughHessianSpores);
print("Hough circles local radius: " + localradiusHoughHessianSpores);
print("Threshold method for unlabelled pathogens with Hessian filter: " + thresholdMethodHessianSpores);

selectWindow("Log");
if(saveResults==true){
	if(useImageFilename==true){
   		saveAs("Text", dir + currentImagename + "__" + nameTagSavedFiles +  "__Parameter_settings.txt");
   	}
   	else{
   		saveAs("Text", dir + "/Image_" + testImageNumber + "__" + nameTagSavedFiles + "__Parameter_settings.txt");
   	}
}
//**** End of parameter saving

//**** Clear up all images and windows if so desired:
if(closeAllWindows==true){
	run("Close All");
	if (isOpen("Log")){ 
		selectWindow("Log");
		run("Close");
	}
	if (isOpen("Summary")){ 
		selectWindow("Summary");
		run("Close");
	}
    if (isOpen("Results")){ 
		selectWindow("Results");
		run("Close");
    }
    if (isOpen("ROI Manager")){ 
		selectWindow("ROI Manager");
		run("Close");
    }
}
//**** Clear-up done, macro complete  
