READ ME

.../slepian_bgetraer/scripts/

This directory contains scripts which implement 2D wavelet decomposition to analyze image data created from NASA's GRACE spherical harmonic data products. 

These scripts are supported by the suite of functions contained within:
	/slepian_bgetraer/functions
And write and read data-files from:
	/slepian_bgetaer/datafiles
	
These scripts are designed to be implemented in a specific order:

1	BOXGREENLAND.m
		IMPLEMENTS:		squaregrid(), lat_lon2Im() 
		CREATES:		jp02fig5, ptsGL.mat, im_tools.mat
2	IMAGERYSEQ.m
		IMPLEMENTS:		grace2plmt(), plm2grid() 
		CREATES:		im_seqSH.mat
3	CHOOSEWAVELET.m
		IMPLEMENTS:		iminvar(), imbias()
		CREATES:		jp02fig6
4	ANALYZEWAVELET.m
		IMPLEMENTS:		prctileThold(), histWTHCOEF()
		CREATES:		threshpassindex.mat
5	WAVEINPOLY.m
		IMPLEMENTS:		
		CREATES:
6	
7	
