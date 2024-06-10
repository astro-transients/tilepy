from tilepy.tools.PlotCumulativeCoverage import ObtainCumulativeProbabilityPlot, ObtainCumulativeProbabilityPlotLog

##ObtainCumulativeProbabilityPlot('S190814bv/PGalinFoV_NObs/SuggestedPointings_GWOptimisation.txt','S190814bv', 'Pgal')

##ObtainCumulativeProbabilityPlotLog('S190814bv_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation_noHA.txt','S190814bv', 'Pgal','2023-09-15 01:30:10')
##ObtainCumulativeProbabilityPlot('S190814bv_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation_noHA.txt','S190814bv', 'PGW','2023-09-15 01:30:10')
#ObtainCumulativeProbabilityPlot('S190814bv_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation.txt','S190814bv', 'Pgal')

ObtainCumulativeProbabilityPlot('../bn231012231_PAPER/PGWinFoV_NObs/SuggestedPointings_GWOptimisation.txt','GRB 231012A', 'PGW')
ObtainCumulativeProbabilityPlot('../927563_lalinference_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation.txt','927563', 'Pgal')
ObtainCumulativeProbabilityPlotLog('../S190814bv_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation.txt','S190814bv', 'PGW')