Per coompilare il programma per il problema del commesso viaggiatore eseguire

	sh prepare/prepare_arg_area.sh
dove 'arg_area'='circle','square' a seconda che si voglia eseguire la ricerca su città disposte su una circonferenza o 
in un quadrato, rispettivamente.

Per eseguire il programma:

	./bin/main --area arg_area --fitness arg_fitness
dove 'arg_fitness'='L1','L2' indica la funzione di costo da utilizzare.

Se si vuole svolgere anche una breve grid search sul numero di generazioni e di individui:
	
	./bin/main --area arg_area --fitness arg_fitness
Nei file 'output/arg_area/arg_fitness/best.fitness.dat' e 'output/arg_area/arg_fitness/ave.fitness.dat' si trovano, rispettivamente,
il valore della funzione di costo del miglior individuo e il valore mediato sulla miglior metà di popolazione, per ogni generazione.
Nel file 'output/arg_area/arg_fitness/solution.dat' si trova il percorso trovato al termine del programma.
Nel file 'output/arg_area/arg_fitness/grid.dat' sono riportati i valori della funzione di costo per le diverse combinazioni
di numero di individui e numero di generazioni.
