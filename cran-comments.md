## Test environments
* ubuntu 16.04 LTS, R 3.3.1
* ubuntu 12.04 LTS (on Travis-CI), R 3.4.0
* win-builder (devel, release, and oldrelease)

## R CMD check results
There were no ERRORs or WARNINGs. 

On Ubuntu

There was 1 NOTE:

* checking installed package size ... NOTE
  installed size is 21.5Mb
  sub-directories of 1Mb or more:
    libs  21.3Mb
 
  The package contains mostly compiled code to improve computation speed. This compiled objects can not be reduced in size.

On Windows

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Cedric Landerer <cedric.landerer@gmail.com>'

  New submission

  Possibly mis-spelled words in DESCRIPTION:
    Genomic (3:20)
    Stationarity (3:39)

  Words are correclty spelled and Maintainer is Cedric Landerer

## Downstream dependencies
There are no downstream depencies as this is a first submission.
