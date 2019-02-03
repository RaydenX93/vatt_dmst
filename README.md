# VATT DMST Code
Questo documento illustrerà brevemente il funzionamento e la struttura del codice MATLAB per la valutazione di performance di una turbina marina ad asse verticale (VATT) mediante teoria DMST con aggiunta di sottomodelli per tener conto di fenomeni idrodinamici ulteriori.

Per la teoria del modello, si rimanda al lavoro di tesi di Stefano Deluca scaricabile sul portale ETD dell’Università di Pisa ([Link](https://etd.adm.unipi.it/theses/browse/by_type/LM.html)).

Tutte le unità adottate in questo lavoro sono sempre SI, _m_ per le lunghezze, _s_ per il tempo, _kg_ per la massa.

Si consiglia di leggere il manuale nella sua interezza prima di effettuare simulazioni, possibilmente visionando il codice sorgente allo stesso tempo.

## Indice

## Prerequisiti
Il programma è stato sviluppato e testato su sistema operativo **Windows 10 (64 bit)**. L’uso su altre piattaforme potrebbe richiedere modifiche al codice.

È necessario installare MATLAB R2018b con pacchetti:

 - MATLAB Coder 4.1 Parallel Computing Toolbox 6.13 Image Processing
   Toolbox 10.3 (solo per simulazioni MIT)

Successivamente bisogna installare:

<![if !supportLists]>· <![endif]>MATLAB Support for MinGW-w64 C/C++ Compiler ([Link](https://it.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler))

L’uso di versioni diverse non è stato testato.
<!--stackedit_data:
eyJoaXN0b3J5IjpbMjA1MjQwNzM2MF19
-->