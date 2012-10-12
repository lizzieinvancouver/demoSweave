2012 October 12

Thanks to all those who attended the Sweave workshop! 

A couple notes here. 

***********************
** DEMO FILES TO TRY **
***********************

Here I have included some slightly updated Sweave files here:
sweavedemo_RelativePaths.Rnw should run if you toss the pinot.txt file in the same folder as your Rnw (I don't use relative paths that often since my data and such overlap between projects most often, and I personally find them less tidy, but I should use them some places I don't)

sweavedemo.Rnw you will need to change for the location of pinot.txt at two points in the code (where the data are read in, twice), then it should run.

In both files:
- I turned off the line in the preamble that controls making a nice little folder for the figures and formatting them. But I highly suggest you try turning it on (remove the % in front of it) as it keeps things much cleaner than the default.

- I also fixed the one error, my LaTeX file is now compiling error-free. Thanks Flo!

***********************
** INSTALLING SWEAVE **
***********************

If you did not get Sweave set up successfully try following the instructions here: http://www.iugo-cafe.org/chinook/view_node.php?id=2214

Remember to start by downloading the files under 'files to download' (upper menu of this webpage, look to the right) then there are a couple commands in R, *then* follow the online slides. 

For Mac users I enclosed my Sweave.engine file in case that is tripping you up. (Also note that the library path you want is probably in your personal User account, for example mine is Lizzie/Library. . . )

*** If you are still having issues do not hesitate to come by my office (Beaty 241). Try not to get too frustrated! Once you have it running, things are much easier! ***

We managed to get everyone who stayed up and running, save for Rebecca where we tracked the issue back to an RStudio bug they are working on, and figured out a temporary work-around.

**********************
** Things mentioned **
**********************

latex2rtf
. . .is my favorite way to convert TeX to rtf (which you can then open with formatting in .doc). If you have trouble installing this (make, make install etc.) come by office (Beaty 241).

Good editors to try:
for Macs: Aquamacs, TextMate
for PCs: Notepad++

Thanks to Ross who found than in TeXShop you can trash extra LaTeX files via the menu File -> Trash aux files

