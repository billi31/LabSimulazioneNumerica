Per compilare il programma di ricerca parallela per il problema del commesso viaggiatore è necessario eseguire 
il comando 'make all'. Per eseguire il programma:

	mpiexec -np 4 main

Nei file 'output/square/arg_fitness/arg_rank/best.fitness.dat' e 'output/square/arg_fitness/arg_rank/ave.fitness.dat'
(dove 'arg_fitness' = 'L1','L2'; 'arg_rank' = '0','1','2','3') sono contenuti i valori della funzione di costo del miglior individuo
e della migliore metà di popolazione, rispettivamente, per il rank 'arg_rank'. Nel file 'output/square/arg_fitness/arg_rank/solution.dat'
è mostrata la soluzione trovata per il rank 'arg_rank'.
Nel file 'output/square/arg_fitness/best.fitness.rank.dat' è contenuto il valore migliore della funzione di costo di ogni rank, per ogni generazione;
nel file 'output/square/arg_fitness/arg_rank/cities.dat' sono contenute le coordinate delle città distribuite randomicamente in un quadrato.

