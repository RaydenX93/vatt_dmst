# VATT DMST Code
Questo documento illustrerà brevemente il funzionamento e la struttura del codice MATLAB per la valutazione di performance di una turbina marina ad asse verticale (VATT) mediante teoria DMST con aggiunta di sottomodelli per tener conto di fenomeni idrodinamici ulteriori.

Per la teoria del modello, si rimanda al lavoro di tesi magistrale di Stefano Deluca dal titolo _**Development and validation of a tool for vertical axis tidal turbine performance evaluation in marine environments**_, scaricabile sul portale *ETD* dell’*Università di Pisa* ([Link](https://etd.adm.unipi.it/theses/browse/by_type/LM.html)).

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
**Sintassi**
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

Ora verranno illustrate le grandezze calcolate nelle matrici, riferite alla cella in posizione azimutale _i_, variabile fra 1 (piano di turbina più in basso, vicino al fondale) e *n_ring*, e posizione verticale _k_, variabile fra 1 e *nz*.

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

La definizione della matrice `data_dyn(k,i,x)` dipende dal tipo di modello di stallo dinamico adottato. Si invita ad ispezionare _dmst_calc.m_ per maggiori informazioni. Queste informazioni non dovrebbero essere particolarmente rilevanti ai fini dell’output del codice DMST.

### `init_input.m`
Questo file è lo script che definisce le impostazioni della simulazione a livello di geometria di turbina, tipo di flusso, discretizzazione del rotore, ecc. Seguire i commenti nel file per configurare la simulazione come si desidera.

**Sintassi**
`[sim_settings, sim_input] = init_input`
Per gli output si rimanda al paragrafo _vatt_dmst.m_.

### `init_geom.m`
Questo file si occupa di realizzare la discretizzazione spaziale del rotore in nz*n_ring celle.

**Sintassi**
`[data_geom] = init_geom(sim_input)`
Per gli input e output si rimanda al paragrafo _vatt_dmst.m_.

### `init_vel.m`

Questo file si occupa di applicare il flusso indisturbato ad ogni cella del rotore.

**Sintassi**
`[data_vel] = init_vel(sim_input,  sim_settings, data_geom, vel_input)`
Per gli input e output si rimanda al paragrafo _vatt_dmst.m_.

### `dmst_update.m`

Questo file gestisce la singola iterazione del solutore DMST. Nello specifico, separa le informazioni disponibili per ogni piano di turbina _k_ e successivamente risolve in parallelo i singoli piani di turbina (sezione parfor).

**Sintassi**
`[data_out_geom, data_out, data_dyn] = dmst_update(sim_settings, sim_input, data_geom, data_vel)`
Per gli input e output si rimanda al paragrafo _vatt_dmst.m_.

Si sottolinea che la variabile interna `dyn_input` è dichiarata come variabile di tipo *persistent* e pertanto continua ad esistere fra le varie chiamate di _dmst_update.m._

### `dmst_par_loop.m`

Questo file gestisce la soluzione del singolo piano di turbina. Viene risolta prima la parte upstream, trovando il fattore di induzione che fa tornare i bilanci per ogni cella azimutale della prima metà di turbina, poi c’è l’eventuale correzione di espansione dei tubi di flusso e successivamente viene risolto il downstream.

**Sintassi**
`[geom_out_data, out_data, dyn_data] = dmst_par_loop(sim_settings, sim_input, dmst_input)`

|Nome Input|Tipo  | Descrizione |
|-|-|-|
|`dmst_input`  | Struttura | Contiene informazioni sul singolo piano da risolvere|

Per i restanti input e output si rimanda al paragrafo _vatt_dmst.m_.

### `dmst_calc.m`

Questo file risolve il singolo streamtube per il valore del fattore di induzione fornito.

**Sintassi**
`[f_eff, data_geom, data_out, data_dyn] = dmst_calc(sim_settings, sim_input, dmst_input, pos_theta, a, a_upstream)`

|Nome Output|Tipo  | Descrizione |
|-|-|-|
|`f_eff`  | Scalare| <![endif]--> Differenza fra forza di spinta calcolata con teoria Actuator Disk e Blade Element |
___
|Nome Input|Tipo  | Descrizione |
|-|-|-|
|`pos_theta`  | Scalare | Indice della posizione azimutale |
|`a`  | Scalare | Fattore di induzione da adottare |
|`a_upstream`  | Scalare | Fattore di induzione della cella a monte, se applicabile |

Per i restanti input e output si rimanda al paragrafo _vatt_dmst.m_.

### `dmst_post.m`

Questo file crea vettori utili per l’output

**Sintassi**
`post_data = dmst_post(sim_input, out_geom_data, geom_data, out_data, data_vel)`
Per gli input e output si rimanda al paragrafo _vatt_dmst.m_.

### `dmst_plot_update.m`

Questo file crea e aggiorna i plot di interesse ad ogni iterazione.

**Sintassi**
`dmst_plot_update(fig_out, data_post)`

|Nome Input|Tipo  | Descrizione |
|-|-|-|
|`fig_out`  | Figura | Figura dove creare/aggiornare plots |

Per i restanti input e output si rimanda al paragrafo _vatt_dmst.m_.

## Future work
Segue una lista di possibili spunti per il miglioramento del codice, in ordine di crescente difficoltà di implementazione (circa), con relativi commenti.

1. **Cambiare il database di dati di riferimento per profili in accordo con ([Link](https://www.sciencedirect.com/science/article/pii/S1364032114005784))**
Bisogna valutare se effettivamente ne valga la pena, ma è molto semplice da fare, bisogna solo modificare _static_data.m_ e il sorgente della routine di stallo dinamico con i nuovi valori.

2. **Raffinare modello di curvatura di streamlines**
Si consiglia di effettuare nuove simulazioni CFD 2D con vario numero di pale, TSR, ecc. In ognuna bisognerà impostare un report per ogni timestep della quantità _Velocity Angle_ su una circonferenza centrata nell’origine (dov’è il centro della turbina) di raggio pari al raggio nominale della turbina. Ogni report dovrà essere mediato nel tempo (per usare valori del velocity angle privi di fluttuazioni). L’idea è ricavare un database di velocity angle al variare di posizione azimutale, TSR e numero di pale ed eventualmente altro. Successivamente il fattore <![if !msEquation]>  <![endif]> dovrà essere ricavato da una interpolazione lineare in più dimensioni (ispirarsi al modello di perdite alle punte nuovo).

3. **Raffinare modello di espansione degli streamtube**
Migliorare la giustificazione fisica del modello corrente o trovarne uno nuovo. Sicuramente migliorare il fitting effettuato sul modello attuale.

4. **Passare ad una discretizzazione verticale non più uniforme, con infittimento verso le punte.**
Potrebbe essere piuttosto ostico, soprattutto perché bisognerebbe cambiare molte formule in giro per il codice.

5. **Implementare effetti di curvatura del flusso.**
Bisogna prima capire perché il modello presente non funziona per codici DMST ma va bene per UDF su Fluent. Se non dovesse andar bene, va reimplementato in altro modo.

6. **Migliorare modello di perdite alle punte**
Servono dati da simulazioni CFD 3D, possibilmente dipendenti anche dalla posizione azimutale. Bisogna capire quanto è importante. Non particolarmente difficile di per sé, ma richiede simulazioni CFD molto lunghe.
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE1NDQxNDAxMjEsMjA2NTQ3MTczMywtND
kwMjY1MzkwXX0=
-->