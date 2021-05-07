## Test environments
* local R installation, R 4.0.3
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* win-builder (devel)

## R CMD check results
0 errors ✓ | 0 warnings ✓ | 3 notes x

> checking CRAN incoming feasibility ... NOTE<br />
  New submission<br />
  Version contains large components (0.0.0.9000)<br />
  Possibly mis-spelled words in DESCRIPTION:  <br />
  Croucher (11:117)<br />
  BRATNextGen (11:188)<br />
  Maintainer: 'Chrispin Chaguza <Chrispin.Chaguza@gmail.com>'<br />
    Gubbins (11:108)<br />
    Marttinen (11:201)<br />
    al (11:129, 11:214)<br />
    et (11:126, 11:211)<br />
  Size of tarball: 7546355 bytes<br />

> checking installed package size ... NOTE<br />
    installed size is 12.0Mb<br />
    sub-directories of 1Mb or more:<br />
      doc       7.4Mb<br />
      extdata   1.6Mb<br />
      help      2.3Mb<br />

> checking for future file timestamps ... NOTE<br />
  unable to verify current time<br />

## COMMENTS
* All words in the DESCRIPTION file are spelled correctly.

* The documentation contain high resolution images and example dataset with size >1Mb.

* This is a new release.


