Per compilare il programma, eseguire:

	sh run/run.sh
	
In questo modo si eseguirà una grid search per trovare i valori ottimali di mu e sigma (mostrata nel file 'output/grid.etot.0'
e si svolgerà il calcolo delle medie a blocchi (in 'output/output.etot.**.0') di valori nei pressi dei valori ottimali appena trovati;
inoltre, se svolge la media a blocchi per i valori di mu e sigma contenuti in input.dat, che viene registrata nel file 'output/output.etot.0'.

Se si desidera solo svolgere la grid search o il calcolo delle medie a blocchi per i valori fissati di mu e sigma, eseguire

	sh run/run_grid.sh
oppure
	
	sh run/run_fixed.sh

Eseguendo 
	
	sh run/run_insta.sh
il programma genera 5e5 valutazioni della funzione |psi|^2 (che si trovano nel file 'output/instantaneus.psi2_etot.0')
con i valori di mu e sigma presenti nel file input.dat.
