Za pokretanje programa potrebno je instalirati CodeBlocks (naša verzija je 13.12),
RStudio (verzija 3.5.1) i Qt Creator(verzija 5.12 s pripadnim kompajlerom MinGW).

U programima projekt_x.c se izraèunava fitness pomoæu mravljih algoritama,
za pokretanje treba otvoriti programe u CodeBlocksu i pritisnuti Build and Run 
(ili tipku F9 na tipkovnici), csv datoteke moraju biti sadržane u istoj mapi( C-kodovi).

Za primjenu K-means metode otvore se programi tipa x.R, postavi se direktorij R-kodovi kao 
"radni" (set as working directory) i pokrenu se pogrami pomoæu Ctrl+Enter. 
Datoteke tipa x.csv moraju biti sadržane u istoj mapi R-kodovi. Program æe 
ispisati u datoteku x.txt niz klastera, a ta datoteka se ruèno prilagodi za rad
u C-u tako da sadrži samo klastere odvojene razmakom bez dodataka koje je program u orginalu 
ispisao. Takve prilagoðene datoteke su tipa output_x.txt (veæ dodane u datoteku C-kodovi) 
i koriste se u C-kodovima fitness_x.c gdje se izraèunava fitness nakon metode K-means.
Napomena: programi mall.R i ant_mall.R ne sadrže dio ispisivanja u output datoteke, veæ
grafom prikazuju rezultate(mall.R rezultate nakon k_means, ant_mall.R nakon ACO iz datoteke output.txt
koja je dobivena pomoæu projekt_mall.c).

Za pokretanje grafièkog suèelja potrebno je otvoriti program test.pro iz datoteke test u Qt Creator-u, 
pritisnuti run qmake, build i run.

Grafovi rezultata su implementirani u RStudiju, potrebno je samo pokrenuti programe x.R
iz datoteke grafovi-rezultati u RStudiju.

U mapi csv datoteke su samo pobrojani svi skupovi podataka koji su korišteni u projektu.


