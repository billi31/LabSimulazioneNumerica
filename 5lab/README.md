Per compilare utilizzare il comando 'make all'. Quando si esegue il programma si può specificare lo stato 'arg_state' dell'idrogeno da utilizzare ('100' o '210') e la distribuzione di probabilità 'arg_mode' utilizzata per il campionamento delle posizioni ('unif' o 'gauss'). 

Per generare file contenenti punti distribuiti secondo le funzioni d'onda basta scrivere da terminale

    ./bin/main --type graph

per specificare quale stato o quale distribuzione utilizzare, basta aggiungere gli argomenti **facoltativi**

    ./bin/main --type graph --state arg_state --mode arg_mode


Per generare file contenenti i valori in media a blocchi ed istantanei della distanza dall'origine dei punti generati si sfruttano i file '.sh' nella cartella 'run'. Se si vogliono generare i file relativi unicamente ad uno stato e ad una distribuzione di probabilità specifici, basta scrivere a terminale:

    sh run/run.arg_state.arg_mode.sh

se si desidera specificare solo lo stato della funzione d'onda:

    sh run/run.arg_state.sh

se si desidera specificare solo la distribuzione:

    sh run/run.arg_mode.sh

se non si desidera specificare niente:

    sh run/run.all.sh