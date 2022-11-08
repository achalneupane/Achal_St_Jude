libname s "Z:\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\2 Final Data SJLIFE\20200430\Survey Data";
options fmtsearch=(f);
run;

/*Get pa and vitality from questionnaires*/
/*Physical activity data*/

data phys_a (keep=mrn datecomp mvpawk mets meetcdc sedentary adolescent);
	set s.adult_healthhabits;
/*Get those with mvpa minutes first*/
if vpa10=1 and vpadays=. then do vpadays=1; end;
if vpa10=1 and vpamin=. then do vpamin=10; end;
if vpamin>360 then do vpamin=360; end; /*Cap six hours per day*/
if vpa10=1 then do wvpa=vpadays*vpamin; end;
if vpa10=2 then do wvpa=0; end;
if wvpa=. then do;
	if vpa10 in (.,2) and nopa=1 and pa20 in (.,0) then do wvpa=0; end;
	if vpa10 in (.,2) and pa20 not in (.,0) then do wvpa=pa20*20; end;
	if vpa10 in (.,2) and nopa=2 and pa20 in (.,0) then do wvpa=0; end;
end;

if mpa10=. and (mpadays ne . or mpamin ne .) then do mpa10=1; end;
if mpa10=1 and mpadays=. then do mpadays=1; end;
if mpa10=1 and mpamin=. then do mpamin=10; end;
if mpamin>360 then do mpamin=360; end; /*Cap six hours per day*/
if mpa10=1 then do wmpa=mpadays*mpamin; end;
if mpa10=2 then do wmpa=0; end;
if wmpa=. then do;
	if mpa10 in (.,2) and nopa=1 then do wmpa=0; end;
	if wvpa ne . then do wmpa=0; end;
end;

mvpawk=(wmpa)+(wvpa*2); if mvpawk>2520 then mvpawk=2520; /*cap six hours per day*/
mets=mvpawk*3;
if . < mvpawk <150 then meetcdc=0;
else if mvpawk >=150 then meetcdc=1;
if nopa in (.,1) and ltpa in (.,2) and mvpawk=0 then sedentary=1; 
else if nopa=2 or ltpa=1 or mvpawk > 0 then sedentary=0;
adolescent=0;
run;


