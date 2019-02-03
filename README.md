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

Il file con i dati è *out.<span></span>nc*.
1. Per preparare i dati alla simulazione DMST, avviare lo script *MIT_read_nc_data.m* modificando la variabile pe_file al suo interno con il percorso al file .nc desiderato. Verrà prodotto, insieme ad altri file, il file *%nome_file%_processed.mat*.
2. Per avviare la simulazione MIT, aprire lo script MIT_run_simulation.m e modificare la variabile nome_file con il file *%nome_file%_processed.mat* creato precedentemente (senza estensione).
3. Modificare la variabile sim_step per impostare la risoluzione desiderata.
4. Avviare lo script.
5. I risultati complessivi sono salvati in *%nome_file%_sub_finished.mat*. Sono anche presenti i file delle singole simulazioni in una sottocartella.

## Struttura dei file
In questo paragrafo verranno investigati i singoli script, di cui verranno spiegati i rispettivi input e output. Non si entrerà nei meriti del codice in maniera troppo tecnica in quanto gli script sono sufficientemente commentati. Per informazioni sugli script nella cartella _func_, data la loro natura molto specifica, si consiglia di leggere direttamente il codice sorgente.

### `vatt_dmst.m`
Questo file è la function da cui lanciare le simulazioni.
`[data_post, data_geom, data_vel, data_out_geom, data_out, data_dyn, sim_input, sim_settings] = vatt_dmst(vel_input, [output_file], [tsr_override]);`


|Nome Output|Tipo  | Descrizione |
|-|-|-|
|`data_post`  | Struttura | Contiene vettori utili per l’output |
|`data_geom`  | Struttura | Contiene grandezze informazioni sulla geometria calcolata del rotore |
|`data_vel`  | Struttura | Contiene grandezze relative al flusso in ingresso |
|`data_out_geom`  | Matrice | Contiene grandezze relative alla geometria del rotore e del flusso calcolate su ogni cella del rotore. Dimensioni: nz **X** n_ring **X** 16 |
|`data_out`  | Matrice | Contiene grandezze di output calcolate su ogni cella del rotore. Dimensioni: nz **X** n_ring **X** 8 |
|`data_dyn`  | Matrice | Contiene grandezze relative alla routine di stallo dinamico per ogni cella del rotore. È una matrice nz **X** n_ring **X** 8 |
|`sim_input`  | Struttura | Contiene impostazioni della simulazione |
|`sim_settings`  | Struttura| Contiene informazioni sui sottomodelli attivi |

|Nome Input|Tipo  | Descrizione |
|-|-|-|
|`vel_input`  | Matrice/Scalare | Contiene informazioni sul flusso indisturbato.
 _Simulazioni 2D_
Scalare della velocità.

 _Simulazioni 3D_
Matrice di 2 colonne.
La prima indica le posizioni z a cui la velocità viene misurata da 0 (pelo libero del mare) a <![if !msEquation]>  <![endif]> (fondale). La seconda, il valore di velocità. |

|`output_file`  | Struttura | Contiene grandezze informazioni sulla geometria calcolata del rotore |
|`tsr_override`  | Struttura | Contiene grandezze relative al flusso in ingresso |


<!--stackedit_data:
eyJoaXN0b3J5IjpbMjA2NTQ3MTczMywtNDkwMjY1MzkwXX0=
-->