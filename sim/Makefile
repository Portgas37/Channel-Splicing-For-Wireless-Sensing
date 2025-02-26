
PATH_TO_WISPR         = ~/gitlab_wisense/packages/wispr
PATH_TO_SIM           = $(PATH_TO_WISPR)/matlab/simulators/channel_splicing
PATH_TO_SPLICING      = $(PATH_TO_WISPR)/matlab/wispr_matlab/splicing
PATH_TO_SPLICERS      = $(PATH_TO_WISPR)/matlab/wispr_matlab/splicing/splicers
PATH_TO_VISUALIZATION = $(PATH_TO_WISPR)/matlab/wispr_matlab/visualization
PATH_TO_80211         = $(PATH_TO_WISPR)/matlab/wispr_matlab/80211

download:
	cp $(PATH_TO_SIM)/*.m .
	cp $(PATH_TO_SIM)/README.md .
	cp $(PATH_TO_SPLICING)/*.m .
	cp $(PATH_TO_SPLICERS)/*.m .
	cp $(PATH_TO_VISUALIZATION)/plotComplexResponse.m .
	cp $(PATH_TO_80211)/generateChannel.m .
	cp $(PATH_TO_80211)/generateFreqAxis.m .
	cp $(PATH_TO_80211)/getFftDcIdx.m .
	cp $(PATH_TO_80211)/getPrunedDFTMatrix.m .
	cp $(PATH_TO_80211)/visualizeChannel.m .

	rm -f testCIR.m testDFT.m phaseOffsetEst.m spliceInitInfoUniform testnufft.m csiSplicerChronos.m testnufft_solver.m
