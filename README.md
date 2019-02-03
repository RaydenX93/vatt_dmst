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

## Come lanciare una simulazione normale

Questo paragrafo è volutamente discorsivo. Si rimanda all’analisi dei singoli file per informazioni sulla sintassi degli script.

 1. Aprire il file _init_input.m_ e definire le opzioni della simulazione.
 2. Lanciare la funzione _dmst_vatt.m_ specificando il profilo di velocità da usare, per simulazioni 3D, oppure semplicemente il valore di velocità del flusso, per simulazioni 2D. Facoltativamente, definire il nome del file di output.
 3. Realizzare i plot desiderati usando i dati di output.

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
___
|Nome Input|Tipo  | Descrizione |
|-|-|-|
|`vel_input`  | Matrice/Scalare | Contiene informazioni sul flusso indisturbato.<ul><li>_Simulazioni 2D_<br>Scalare della velocità.</li><li> _Simulazioni 3D_<br>Matrice di 2 colonne.<br>La prima indica le posizioni z a cui la velocità viene misurata da 0 (pelo libero del mare) a -∞ (fondale). La seconda, il valore di velocità. |</li></ul>
|`output_file`  | Struttura | *Opzionale*. Indica il file di output dove salvare le simulazioni. Se non fornito, viene generato un nome dal programma. |
|`tsr_override`  | Struttura | *Opzionale*. Imposta un TSR diverso da quello specificato in `init_input.m` |

Per un esempio su come avviare una simulazione, fare riferimento allo script *ESEMPIO_lanciasim_singola.m* e *ESEMPIO_dmst_optimizer.m*

I nomi delle variabili all’interno delle strutture sono piuttosto intuitivi. In caso di dubbio, risalire dal codice sorgente alla grandezza calcolata.

Ora verranno illustrate le grandezze calcolate nelle matrici, riferite alla cella in posizione azimutale _i_, variabile fra 1 (piano di turbina più in basso) e *n_ring*, e posizione verticale _k_, variabile fra 1 e *nz*.

<table>  
<tr>  
	<code>data_out(k,i,x)</code>
</tr>  
<tr>  
	<td>1</td>  
	<td>Fattore di induzione assiale non corretto</td>  
	<td>5</td>  
	<td>Fattore di perdita alle punte</td>
</tr>  
<tr>  
	<td>2</td>  
	<td>Fattore di induzione assiale a monte, se la cella è in downstream</td>  
	<td>6</td>  
	<td>CL senza perdita alle punte</td>
</tr>  
<tr>
	<td>3</td>  
	<td>CL con perdita alle punte</td>  
	<td>7</td>  
	<td>CD senza perdita alle punte</td>
</tr>
<tr>
	<td>4</td>  
	<td>CD con perdita alle punte</td>  
	<td>8</td>  
	<td>Fattore di induzione assiale dopo correzione espansione streamtubes</td>
</tr>
</table>

___

<table>  
<tr>  
	<code>data_out_geom(k,i,x)</code>
</tr>  
<tr>  
	<td>1</td>  
	<td><i>R</i></td>  
	<td>Raggio</td>  
	<td>9</td>
	<td><i>cos(α)</i></td>
	<td>Coseno angolo di attacco</td>
</tr>  
<tr>  
	<td>2</td>  
	<td><i>cos(θ)</i></td>  
	<td>Coseno pos. azimutale</td>  
	<td>10</td>
	<td><i>α</i> [deg]</td>
	<td>Angolo di attacco</td>
</tr>  
<tr>
	<td>3</td>  
	<td><i>sin(θ)</i></td>  
	<td>Seno pos. azimutale</td>  
	<td>11</td>
	<td><i>α<sub>virt</sub></i></td>
	<td>Angolo di attacco virtuale</td>
</tr>
<tr>
	<td>4</td>  
	<td><i>W<sub>0</sub></i></td>  
	<td> Componente orizzontale (parallela flusso) di vel. relativa</td>  
	<td>12</td>
	<td><i>Δz</i></td>
	<td>Altezza del piano di turbina</td>
</tr>
<tr>  
	<td>5</td>  
	<td><i>W<sub>1</sub></i></td>  
	<td>Componente verticale (perpendicolare flusso) di vel. relativa</td>  
	<td>13</td>
	<td><i>θ</i> [rad]</td>
	<td>Posizione azimutale</td>
</tr>  
<tr>  
	<td>6</td>  
	<td><i>ModW</i></td>  
	<td>Modulo vel. relativa</td>  
	<td>14</td>
	<td><i>TSR<sub>loc</sub></i></td>
	<td>TSR locale nella cella</td>
</tr>  
<tr>
	<td>7</td>  
	<td><i>Re<sub>c</sub></i></td>  
	<td>Numero di Reynolds basato su corda</td>  
	<td>15</td>
	<td><i>U</i></td>
	<td>Componente orizzontale (parallela flusso) di vel. assoluta</td>
</tr>
<tr>
	<td>8</td>  
	<td><i>sin(α)</i></td>  
	<td>Seno angolo di attacco</td>  
	<td>16</td>
	<td><i>V</i></td>
	<td>Componente verticale (perpendicolare flusso) di vel. assoluta</td>
</tr>
</table>

La definizione della matrice data_dyn(k,i,x) dipende dal tipo di modello di stallo dinamico adottato. Si invita ad ispezionare _dmst_calc.m_ per maggiori informazioni. Queste informazioni non dovrebbero essere particolarmente rilevanti ai fini dell’output del codice DMST.
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTUxNzE5ODg3OSwyMDY1NDcxNzMzLC00OT
AyNjUzOTBdfQ==
-->