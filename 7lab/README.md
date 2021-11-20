Per compilare ed eseguire il programma ed ottenere 500000 valori istantanei di energia interna e pressione (per il calcolo dell'autocorrelazione), per tutte e tre le fasi, scrivere

	sh run_instantaneus.sh
	
Per compilare ed eseguire unicamente il calcolo delle medie a blocchi per le tre fasi, scrivere

	sh run.sh

Per ottenere i risultati di un'unica fase 'phase' = 'solid', 'liquid' o 'gas', cominciando da una vecchia configuarzione 'mode' = 'repeat' o dal reticolo fcc 'mode' = 'start', scrivere

	sh prepare/prepare_phase.sh
	sh prepare/prepare_mode.sh

per compilare (l'opzione '--insta' è facoltativa e procede nel calcolo dei valori istantanei di pressione ed energia interna):

	./bin/main --state 'phase' --insta
