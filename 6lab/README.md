Per compilare ed eseguire il programma sia per l'algoritmo Metropolis, sia per l'algoritmo di Gibbs, con e senza campo magnetico esterno, scrivere da terminale:

	sh run.sh

Se si desidera, invece, eseguire il programma utilizzando l'algoritmo 'algo' ('metro' o 'gibbs') e un campo esterno 'h' = '0' o '02' è sufficiente compilare con

	sh prepare/prepare_algo_h.sh

se si desidera partire da una vecchia configurazione 'mode' = 'repeat', altrimenti 'mode' = 'start' ripartirà dalla configurazione fcc; per eseguire:

	./bin/main --mode 'mode'
