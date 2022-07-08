################################################################################
# Plot scheduling or zenith angle lines. Uses functions of PointingPlotting
# TO-Do: Needs to be more generalized
################################################################################

from gammagwfollowup.include.PointingTools import TableImportCTA, TableImportCTA_Time
from gammagwfollowup.include.PointingPlotting import PointingPlottingGWCTA,PointingPlottingGW_ZenithSteps

def PlotScheduling(ID,InputFileName,dirName,pointingsFile,FOV):
    
    j=ID
    
    # GW file
    #InputFileName = '../../dataset/BNS-GW_onAxis5deg.txt'
    InputList = TableImportCTA(InputFileName)
    GWFile = "../../dataset/skymaps/" + InputList['run'][j] + '_' + InputList['MergerID'][j] + "_skymap.fits.gz"
    name =  InputList['run'][j] + '_' + InputList['MergerID'][j]
    
    InjectTimeFile = '../../dataset/BNS-GW-Time_onAxis5deg_postRome.txt'
    InputTimeList = TableImportCTA_Time(InjectTimeFile)
    
    #pointingsFile= '../../gw-follow-up-simulations-side-results/TestPointings/run0017_MergerID000261_GWOptimisation.txt'
    #FOV= 2.0
    
    PointingPlottingGWCTA(GWFile,name,dirName,pointingsFile,FOV,InputTimeList['Observatory'][j])


def PlotZenithAngleLines(ID,InputFileName,dirName,FOV):

    j = ID
    
    # GW file
    InputList = TableImportCTA(InputFileName)
    
    GWFile = "../../dataset/skymaps/" + InputList['run'][j] + '_' + InputList['MergerID'][j] + "_skymap.fits.gz"
    name = InputList['run'][j] + '_' + InputList['MergerID'][j]
    
    InjectTimeFile = '../../dataset/BNS-GW-Time_onAxis5deg_postRome.txt'
    InputTimeList = TableImportCTA_Time(InjectTimeFile)
    
    PointingPlottingGW_ZenithSteps(GWFile, name, dirName, FOV, InputTimeList[j])
