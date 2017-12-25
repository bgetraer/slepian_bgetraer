# slepian_bgetraer

This directory contains MATLAB functions and scripts complementing the slepian_alpha, bravo, etc suite fount at https://github.com/csdms-contrib modified and written by bgetraer@princeton.edu, FALL 2017.

Common problems I ran into with very very simple solutions:
	
	Function continues to fail, asking for some file or directory you do not have.
		Make the directory folder being requested, find the requested file from http://geoweb.princeton.edu/people/simons/software.html and put it where it belongs.
	Slepian bases are not plotting where you expect them to.
		Some functions use latitude, some use co-latitude, and you are mixing them up.
	
	