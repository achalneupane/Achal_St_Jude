libname trk 		"C:\Users\ghyun\OneDrive - St. Jude Children's Research Hospital\UDrive\Data\SJLIFE 2020_20200430\Tracking Data" ;
libname cl  		"C:\Users\ghyun\OneDrive - St. Jude Children's Research Hospital\UDrive\Data\SJLIFE 2020_20200430\Clinical Data" ;
libname survey 		"C:\Users\ghyun\OneDrive - St. Jude Children's Research Hospital\UDrive\Data\SJLIFE 2020_20200430\Survey Data" ;
libname fmt64 		"C:\Users\ghyun\OneDrive - St. Jude Children's Research Hospital\UDrive\Data\SJLIFE 2020_20200430" ;

option fmtsearch=(fmt64) fmterr ;


proc sort data=trk.tracking			 	out=track ; 		by mrn ; 			run ; 
proc sort data=trk.visits				out=visits ; 		by mrn visitnum ; 	run ; 
proc sort data=trk.lstcondt				out=lstcondt ; 		by mrn ; 			run ; 
proc sort data=cl.demographics 			out=demo ; 			by mrn ; 			run ; 
proc sort data=survey.adult_behavior			out=adult_behavior;				by mrn datecomp ;	run ; 
proc sort data=survey.adolescent_behavior 		out=adolescent_behavior ; 		by mrn datecomp ; 	run ; 
proc sort data=survey.adult_healthhabits 		out=adult_healthhabits ;		by mrn datecomp ;	run ; 
proc sort data=survey.adolescent_healthhabits 	out=adolescent_healthhabits ; 	by mrn datecomp ; 	run ; 
proc sort data=survey.adult_home 				out=adult_home ;				by mrn datecomp ;	run ; 
proc sort data=survey.adolescent_home		 	out=adolescent_home ; 			by mrn datecomp ; 	run ; 
proc sort data=survey.adult_mw		 			out=adult_mw ; 					by mrn datecomp ; 	run ; 
proc sort data=survey.abbreviated		 		out=abbreviated ; 				by mrn datecomp ; 	run ; 


data firstvisit(keep=mrn datefirstvisit visitnum) ; 
	set visits(where=(studypop="Survivor" and visitnum=1)) ; 
	datefirstvisit=DateVisitStart ; 
	format datefirstvisit date9. ;
run ; 




****************************Smoking status ;
data adult_healthhabits2 ;
	set adult_healthhabits(keep=mrn cigmo evsm smnow cigd smyr survey datecomp) ; 
	if evsm=. and smnow=. then delete ;
run ; 
data adolescent_healthhabits2 ;
	set adolescent_healthhabits(keep=mrn relation smnvr smnow survey datecomp) ;
	if smnvr=. and smnow=. then delete ;
run ; 
data adolescent_home2 ;
	set adolescent_home(keep=mrn relation evsm smnow survey datecomp) ;
	if evsm=. and smnow=. then delete ;
run ; 
data abbreviated2 ;
	set abbreviated(keep=mrn relation smnow cigd smyr survey datecomp) ;
	if smnow=. and cigd=. then delete ;
run ; 
data Smk_stat ; set adolescent_healthhabits2 adolescent_home2 abbreviated2 adult_healthhabits2 ; run ;
proc sort data=Smk_stat ; by mrn datecomp ; run ; 
data Smk_stat2 ;
	merge Smk_stat(in=a) demo(keep=mrn dob) ;  
	by mrn ;
	if a ;
	a_smksurvey=(datecomp-dob)/365.25 ;
	*Overall correction ;
	if evsm=2 then do ; smnow=2 ; cigmo=2 ; cigd=0 ; smyr=0 ; end ;
	if smnvr=1 then do ; evsm=2 ; smnow=2 ; cigmo=2 ; cigd=0 ; smyr=0 ; end ;
	if evsm=. and smnow=1 then do ; evsm=1 ; end ;
	if evsm=1 and smnow=. and cigmo=1 then smnow=1 ;
	if evsm=1 and smnow=. and cigmo=2 then smnow=2 ;
run ;
/*
proc print data=Smk_stat2 noobs width=min ; where evsm=. or smnow=. ; run ; 
proc freq data=Smk_stat2 ; table count ; run ;
*/
data Smk_Stat3 ; 
	merge dexa(in=a keep=mrn labdt) Smk_stat2(rename=(datecomp=smk_dt)) ; 
	by mrn ; 
	if a ;
	*Smoking status ;
	length SmkStat $12 ;
	if 		evsm=2 				then SmkStat="1: Never" ;
	else if evsm=1 and smnow=2 	then SmkStat="2: Former" ;
	else if evsm=1 and smnow=1 	then SmkStat="3: Current" ;

	diff_smkdt=labdt-smk_dt ;
	absdiff=abs(diff_smkdt) ;
run ; 
proc sort data=Smk_Stat3(where=(SmkStat ne "")) out=Smk_Stat3_s ; by mrn absdiff ; run ; 
*in this study, need to use most close smoke data to the date of your study data ;
data Smk_Stat4 ; set Smk_Stat3_s ; by mrn ; if first.mrn ; run ; 

data Smk_Stat5 ; 
	merge dexa(in=a keep=mrn) Smk_Stat4 ; 
	by mrn ; 
	if a ;
run ; 
