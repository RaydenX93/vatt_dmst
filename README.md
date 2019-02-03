# VATT DMST Code
Questo documento illustrerà brevemente il funzionamento e la struttura del codice MATLAB per la valutazione di performance di una turbina marina ad asse verticale (VATT) mediante teoria DMST con aggiunta di sottomodelli per tener conto di fenomeni idrodinamici ulteriori.

Per la teoria del modello, si rimanda al lavoro di tesi di Stefano Deluca scaricabile sul portale ETD dell’Università di Pisa ([Link](https://etd.adm.unipi.it/theses/browse/by_type/LM.html)).

Tutte le unità adottate in questo lavoro sono sempre SI, _m_ per le lunghezze, _s_ per il tempo, _kg_ per la massa.

Si consiglia di leggere il manuale nella sua interezza prima di effettuare simulazioni, possibilmente visionando il codice sorgente allo stesso tempo.

## Indice

## Prerequisiti
Il programma è stato sviluppato e testato su sistema operativo **Windows 10 (64 bit)**. L’uso su altre piattaforme potrebbe richiedere modifiche al codice.

È necessario installare MATLAB R2018b con pacchetti:

 - MATLAB Coder 4.1 
 - Parallel Computing Toolbox 6.13 
 - Image Processing Toolbox 10.3 (solo per simulazioni MIT)

Successivamente bisogna installare:

 - MATLAB Support for MinGW-w64 C/C++ Compiler ([Link](https://it.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler))

L’uso di versioni diverse non è stato testato.
<![endif]-->

## Come lanciare una simulazione normale

Questo paragrafo è volutamente discorsivo. Si rimanda all’analisi dei singoli file per informazioni sulla sintassi degli script.

 1. Aprire il file _init_input.m_ e definire le opzioni della simulazione.
 2. Lanciare la funzione _dmst_vatt.m_ specificando il profilo di velocità da usare, per simulazioni 3D, oppure semplicemente il valore di velocità del flusso, per simulazioni 2D. Facoltativamente, definire il nome del file di output.
 3. Realizzare i plot desiderati usando i dati di output.

<![endif]-->

## Come lanciare una simulazione con i dati MIT

Al momento sono disponibili due set di dati posizionati in due cartelle.
- Hz600mN010mw (griglia con risoluzione 600 m)
- Hz200mN010mw (griglia con risoluzione 200 m)

Il file con i dati è _out.`[<span></span](http://www<span></span)>nc_.

<![if !supportLists]>1. <![endif]>Per preparare i dati alla simulazione DMST, avviare lo script _MIT_read_nc_data.m_ modificando la variabile pe_file al suo interno con il percorso al file .nc desiderato. Verrà prodotto, insieme ad altri file, il file _%nome_file%_processed.mat._

<![if !supportLists]>2. <![endif]>Per avviare la simulazione MIT, aprire lo script MIT_run_simulation.m e modificare la variabile nome_file con il file _%nome_file%_processed.mat_ creato precedentemente (senza estensione).

<![if !supportLists]>3. <![endif]>Modificare la variabile sim_step per impostare la risoluzione desiderata.

<![if !supportLists]>4. <![endif]>Avviare lo script.

<![if !supportLists]>5. <![endif]>I risultati complessivi sono salvati in _%nome_file%_sub_finished.mat._ Sono anche presenti i file delle singole simulazioni in una sottocartella.
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTQ5MDI2NTM5MF19
-->