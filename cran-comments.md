## Test environments
* ubuntu 16.04 LTS, R 3.4.4
* ubuntu 14.04 LTS (on Travis-CI 3.4.0), R 3.4.2
* win-builder (devel, release, and oldrelease)

## R CMD check results
There were no ERRORs or WARNINGs. 

NOTEs:

On Ubuntu

* checking installed package size ... NOTE
  installed size is 25.2Mb
  sub-directories of 1Mb or more:
    libs  23.7Mb
 
  Response: The package contains mostly compiled code to improve computation speed. The compiled objects can not be reduced in size.

* Author field differs from that derived from Authors@R

  Response: Denizhan Pak started contributing to AnaCoDa since the last update.

On Windows

* Possibly mis-spelled words in DESCRIPTION:
  Codon (3:20)
  Drummond (27:20)
  FONSE (29:13)
  NonSense (29:49, 32:16)
  PANSE (31:57)
  PAusing (31:39, 31:64)
  Ribosome (25:24)
  al (26:12)
  checkpointing (23:22)
  codon (21:64, 24:64, 30:45)
  et (26:9)
  footprinting (32:52, 34:43)
  ribosome (32:43, 33:26, 34:34)

  Response: All words are correctly spelled but not part of a standard dictionary.

* checking installed package size ... NOTE
  installed size is  6.1Mb
  sub-directories of 1Mb or more:
    libs   4.8Mb

  Response: The package contains mostly compiled code to improve computation speed. The compiled objects can not be reduced in size.


## Downstream dependencies
There are no downstream dependencies
