from tilepy.tools.PlotCumulativeCoverage import ObtainCumulativeProbabilityPlotLST,ObtainCumulativeProbabilityPlotMWL,ObtainCumulativeProbabilityPlotCTA, ObtainCumulativeProbabilityPlotLog

#ObtainCumulativeProbabilityPlot('S190814bv/PGalinFoV_NObs/SuggestedPointings_GWOptimisation.txt','S190814bv', 'Pgal')
ObtainCumulativeProbabilityPlotLST('bn231012231_PAPER.fit/PGWinFoV_NObs/SuggestedPointings_GWOptimisation.txt','GRB 231012A', 'PGW')
##ObtainCumulativeProbabilityPlotLog('S190814bv_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation_noHA.txt','S190814bv', 'Pgal','2023-09-15 01:30:10')
##ObtainCumulativeProbabilityPlotLog('S190814bv_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation_noHA.txt','S190814bv', 'PGW','2023-09-15 01:30:10')
ObtainCumulativeProbabilityPlotMWL('S190814bv_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation_noHA.txt','S190814bv', 'Pgal')
ObtainCumulativeProbabilityPlotMWL('S190814bv_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation_noHA.txt','S190814bv', 'PGW')

##ObtainCumulativeProbabilityPlot('MS230826n_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation.txt','S230826n', 'Pgal')
##ObtainCumulativeProbabilityPlot('MS230826n_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation.txt','S230826n', 'PGW')3
ObtainCumulativeProbabilityPlotCTA('927563_lalinference_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation.txt','927563', 'PGW')
ObtainCumulativeProbabilityPlotCTA('927563_lalinference_PAPER/PGalinFoV_NObs/SuggestedPointings_GWOptimisation.txt','927563', 'Pgal')
