# PEARL-SCALLOP
Implementation of PEARL-SCALLOP, a variant of SCALLOP with faster parameters and better security!

## To run:
This folder contains two implementations of PEARL-SCALLOP, one simple SageMath implementation, and one slightly more optimised C++ implementation.

### SageMath implementation

The SageMath implementation is tested on SageMath 10.5. To run a simple example (default at the 512-security level), simply do

```bash
$ cd implementation-sage/
$ sage pearl-scallop.sage
...
ALICE PUBLIC KEY DONE
...
BOB PUBLIC KEY DONE
...
ALICE SHARED SECRET DONE: 144185074944435751013877123195204564013106392817211619102261885777481986164124319439219439625571334553879424364589842140062588483282128984735836982548209219200758771425930893021970049778097898599568878915620521276513398571458172433598494415065502496810365273127*z2 + 80529370329683865931822615881702220186940142216728230855533937485569251097480219979832426754389802376872662172509614827254923947468011229675201999589019839539753833732282463979759122831083419722499174156451015286208949068160668989985120064745947261272557749468
...
BOB SHARED SECRET DONE: 144185074944435751013877123195204564013106392817211619102261885777481986164124319439219439625571334553879424364589842140062588483282128984735836982548209219200758771425930893021970049778097898599568878915620521276513398571458172433598494415065502496810365273127*z2 + 80529370329683865931822615881702220186940142216728230855533937485569251097480219979832426754389802376872662172509614827254923947468011229675201999589019839539753833732282463979759122831083419722499174156451015286208949068160668989985120064745947261272557749468
```

### C++ Implementation

The C++ implementation has the following requirement:
* [NTL](https://libntl.org/) (tested on version 11.5.1)

To run a simple example (default at 1536 security level), do
```bash
$ cd cpp-optimised-implementation/
$ make
$ ./main
...
j(E_AB) = 4443239459191301787939269931222523747810876428912937035423048707435469265394782689827453142941409231322758051398730534877911627198793190874158072897947396599600457997637180270094181495234668472206323328104296314286865658754419300845766801821628299505212185111401552042472632249151216340330071606584524333623795999352085039128891444735009766519872846967956871206248378750375742655464265606561853672486861702768246181785363328 + i*22610707962898645482795718395364724066643839401978279128118770621648206409267796466336284864301855430157046172949487215936043920270705624584072596505128872852174812696651067818763460163184067818041364092404821595077392768479390947445622253730970727109424506276731223428382237626829601140481875322088897727491562456985631083707535340350737307549473054421374823475920040079243365001628421696685512541756378488445410107343357682
j(E_AB) = 4443239459191301787939269931222523747810876428912937035423048707435469265394782689827453142941409231322758051398730534877911627198793190874158072897947396599600457997637180270094181495234668472206323328104296314286865658754419300845766801821628299505212185111401552042472632249151216340330071606584524333623795999352085039128891444735009766519872846967956871206248378750375742655464265606561853672486861702768246181785363328 + i*22610707962898645482795718395364724066643839401978279128118770621648206409267796466336284864301855430157046172949487215936043920270705624584072596505128872852174812696651067818763460163184067818041364092404821595077392768479390947445622253730970727109424506276731223428382237626829601140481875322088897727491562456985631083707535340350737307549473054421374823475920040079243365001628421696685512541756378488445410107343357682
```
Note that the examples uses simple CSIDH-like vectors, and not reduced vectors, from randomly sampled ideal classes. In particular, the coefficients at the 1536-level is significantly smaller than the typical reduced vectors.


## Overview of Repository:
This repository contain the following folders:
* `cpp-optimised-implementation/`: Contains the C++/NTL implementation of PEARL-SCALLOP.
* `implementation-sage/`: Contains the SageMath implementation of PEARL-SCALLOP.
* `params/`: Contains some code to generate a suitable discrimant, and the corresponding factorisation of the conductor.
* `relations/`: Contains the (reduced) lattice of relations at different security levels.
* `starting-curve/`: Contains various code for generating the starting curve.


