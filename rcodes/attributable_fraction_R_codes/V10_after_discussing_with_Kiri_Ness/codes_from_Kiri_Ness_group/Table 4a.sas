libname c "Z:\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\2 Final Data SJLIFE\20170630\Clinical Data";
libname t "Z:\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\2 Final Data SJLIFE\20170630\Tracking Data";
libname e "Z:\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\2 Final Data SJLIFE\20170630\Event Data";
libname s "Z:\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\2 Final Data SJLIFE\20170630\Survey Data";
libname Mess "Z:\SJShare\SJCOMMON\ECC\Ness Research Team\SJLIFE Data Cleaning\20170630\Final Appended Files";
libname HEI "\\stjude\data\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\1 Source Data\FFQ\HEI";
libname smk "\\stjude\data\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\1 Source Data";
libname data "Z:\SJShare\SJCOMMON\ECC\Ness Research Team\Rik Projects\Wilson C\MetSyn\Data";

OPTIONS NOFMTERR;

data inclusion;
	set data.inclusion;
		if part = 1;
run;

/* ____________________________________ DIET */

data hei2 ;
	set hei.HEI_AHEI_20190612;	
	mrndigits= compress (scan(RESPONDENTID,1,'_'),,'kd'); 
	mrn=mrndigits*1;
run;

Proc sql;
     create table hei3 as
          select *, max(TODAYSDATE) as maxdate
               from hei2
               group by mrn;
quit;

data hei4 ;
	set hei3;
	if TODAYSDATE = maxdate ;
run;

proc sort data=hei4 ; by mrn; run;
proc sort data=inclusion; by mrn; run;

data heimerge (drop = maxdate);
	merge inclusion (in=a) hei4 (drop=age);
	by mrn;
	if a;
run;

proc contents data=hei.hei_ahei_20210308 varnum; run;

proc means data= heimerge stackodsoutput n nmiss mean std median p25 p75 min max maxdec=1 nolabels;
	var 
	heiy1_totalveg
	heiy2_green_and_bean
	heiy3_totalfruit
	heiy4_wholefruit
	heiy5_wholegrain
	heiy6_totaldairy
	heiy7_totprot
	heiy8_seaplant_prot
	heiy9_fattyacid
	heiy10_sodium
	heiy11_refinedgrain
	heiy12_addsug
	heiy13_sfa
	hei2015_total_score
;
	ods output summary = heiexcel;
run;

data heiexcel2;
	set heiexcel; 
		MEANnew = strip(put(mean,5.1))||" "||"±"||" "||strip(put(StdDev,5.1));
run;

ods excel file = " C:\Users\rdhaduk\Desktop\random.xlsx";
proc print data=heiexcel2;
run;
ods excel close;


proc rank data=heimerge groups=4 out=dietrank(keep= mrn hei2015_total_score dietrank);
	var hei2015_total_score; ranks dietrank;
run;

proc means data=dietrank n nmiss median min max;
	var hei2015_total_score;
		class dietrank;
run;


/* ____________________________________ SMOKERS */

proc sort data=s.adult_healthhabits ; by mrn datecomp; run; 

data smoke3;
	set s.adult_healthhabits;;
		if evsm ne .;
run;

data smoke2;
	set smoke3;
		by mrn;	
			if last.mrn;
run;

data smoke;
	merge inclusion (in=a) smoke2;
		by mrn;
		if a;

	length smoker $10;
	
	if evsm=. and smnow=. then smoker = ' ';
	else if evsm = 2 then smoker='Never';
	else if evsm = 1 and smnow=1 then smoker='Current';
	else if evsm = 1 and smnow=2 or smnow=. then smoker='Former';
	
	*evsm (Smoked at least 100 cigarettes in your entire life );
	*smnow  (Do you smoke cigarettes now);

run;


/*Subsetting the 21 missing to check */

data smokemissing;
	set smoke;
	if smoker=' ';
run;


/*Kyla gave new source data for smoking information*/

data smkk;
	set smk.smoking_20180214;
run;

proc sort data=smkk; by mrn; run;
proc sort data=smokemissing; by mrn; run;

proc format;
value smoker
0="Never"
1="Ever"
2="Current";
run;


data smokemissingmerge;
	merge smokemissing (in=a) smkk;
	format smk smoker.;
	if a;
	by mrn;
	length smoker $10;
	if smk = . then smoker =' ';
	else if smk = 0 then smoker='Never';
	else if smk = 1 then smoker='Former';
	else if smk = 2 then smoker='Current';

run;

proc freq data=smokemissingmerge ;
	tables smk smoker/norow nocol ; 
run;

/*Got all 21 misssing*/
/*Merging these with the original merged smoke data*/

proc sort data=smoke; by mrn; run;
proc sort data=smokemissingmerge; by mrn; run;

data smokefinal;
	merge smoke (in=a) smokemissingmerge;
		if a;
			by mrn;
run;

proc freq data=smokefinal ;
	tables  smoker/norow nocol out=freqtables; 
run;

data freqtables2;
	set freqtables;
		FREQnew = strip(put(count,5.0))||" "||"("||strip(put(percent,5.2))||"%"||")";
run;

ods excel file = " C:\Users\rdhaduk\Desktop\random.xlsx";
proc print data=freqtables2;
run;
ods excel close;



/* ____________________________________ PHYSICAL ACTIVITY */

Proc sql;
     create table hh as
          select *, max(datecomp) as maxdate
               from s.adult_healthhabits
               group by mrn;
quit;

data hh2 ;
	set hh;
	if datecomp = maxdate ;

	*Correction partially ;
	if vpa10=2 and vpadays=. then vpadays=0 ; 
	if vpa10=2 and vpamin=. then vpamin=0 ;
	if mpa10=2 and mpadays=. then mpadays=0 ; 
	if mpa10=2 and mpamin=. then mpamin=0 ;

	*Counting minutes per week;
	WkTotVigAct=vpadays*vpamin ;
	WkTotModAct=mpadays*mpamin ;

	*Categorizing by CDC std;
	if 		. < WkTotModAct < 150 	and . < WkTotVigAct < 75  then PhysiAct=0 ;
	else if . < WkTotModAct < 150 	and 	WkTotVigAct >= 75 then PhysiAct=1 ;
	else if 	WkTotModAct >= 150 	and . < WkTotVigAct < 75  then PhysiAct=1 ;
	else if 	WkTotModAct >= 150 	and 	WkTotVigAct >= 75 then PhysiAct=1 ;
	else if 	WkTotModAct =.	 	and 	WkTotVigAct >= 75 then PhysiAct=1 ;
	else if 	WkTotModAct =.	 	and . < WkTotVigAct < 75  then PhysiAct=. ;
	else if 	WkTotModAct >= 150 	and 	WkTotVigAct =.	  then PhysiAct=1 ;
	else if . <	WkTotModAct < 150 	and 	WkTotVigAct =.	  then PhysiAct=. ;
run;


proc sort data=hh2; by mrn; run;
proc sort data=inclusion; by mrn; run;

proc contents data=hh2; run;

data pa;
	merge inclusion (in=a) hh2;
	by mrn; length PA3cat $10.;
	if a;
	if first.mrn;
	if WkTotVigAct = . and WkTotModAct ne . then WkTotVigAct = 0;
	if WkTotModAct = . and WkTotVigAct ne . then WkTotModAct = 0;

	totalmin = WkTotVigAct + WkTotModAct;
	if totalmin =. then PA3cat = "";
		else if totalmin lt 450 then PA3cat = "<450";
		else if 450 <= totalmin < 900 then PA3cat = "450 - 899";
		else PA3cat = ">900";

run;



proc freq data=pa ;
	tables  PhysiAct/norow nocol out=freqtables; 
run;

data freqtables2;
	set freqtables;
		FREQnew = strip(put(count,5.0))||" "||"("||strip(put(percent,5.2))||"%"||")";
run;

ods excel file = " C:\Users\rdhaduk\Desktop\random.xlsx";
proc print data=freqtables2;
run;
ods excel close;

proc means data=pa n nmiss mean median maxdec=1;
	var vpadays mpadays vpamin mpamin WkTotVigAct WkTotModAct ;
	class PhysiAct;
run;


/* Adding DT_KCAL intake */

data dtcal (keep = mrn datecomp DT_KCAL);
	set c.Ffqraw; 
		by mrn; if last.mrn;
run;

data cal ;
	merge inclusion (in=a) dtcal;
		by mrn; if a;
			if DT_KCAL = . then dailycal = .; 
				else if DT_KCAL le 1805.05 then dailycal = 0;
				else dailycal = 1;
run;

proc rank data=cal groups=4 out=calrank(keep= mrn DT_KCAL calrank);
	var DT_KCAL; ranks calrank;
run;

proc means data=calrank n nmiss median min max;
	var DT_KCAL;
		class calrank;
run;

/*Saving all three datasets with imp var for following table*/

data data.dietsmokepa (keep = mrn  hei2015_total_score dietrank smoker PhysiAct PA3cat DT_KCAL calrank);
	merge dietrank smokefinal pa calrank;
	by mrn;
run;
