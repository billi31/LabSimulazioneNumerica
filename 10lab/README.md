Per compilare il programma in modo che si trovi la soluzione del problema del commesso viaggiatore con città
distribuite su una circonferenza, basta eseguire da terminale:

	sh prepare/prepare_circle.sh
	
Per trovare la soluzione del problema nel caso di città distribuite casualmente in un quadrato, eseguire

	sh prepare/prepare_square.sh

Per eseguire il programma:

	./bin/main --area arg_area --fitness arg_fitness

dove 'arg_area' = 'circle', 'square' e 'arg_fitness' = 'L1','L2'.
Il programma genererà il percorso trovato nel file 'output/arg_area/arg_fitness/solution.dat'; 
Nel file 'output/arg_area/arg_fitness/best.fitness.dat' si trovano i valori della funzione di costo ad ogni step di temperatura
(ottenuti dopo MC_step step Montecarlo).

Nel file 'output/arg_area/arg_fitness/cities.dat' si trovano le coordinate delle città distribuite randomicamente.
