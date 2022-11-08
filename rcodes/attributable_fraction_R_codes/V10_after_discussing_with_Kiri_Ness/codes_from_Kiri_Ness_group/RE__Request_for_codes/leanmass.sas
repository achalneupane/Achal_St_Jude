libname f "Z:\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\2 Final Data SJLIFE\20181231";
libname c "Z:\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\2 Final Data SJLIFE\20181231\Clinical Data";
libname s "Z:\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\2 Final Data SJLIFE\20181231\Survey Data";
libname t "Z:\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\2 Final Data SJLIFE\20181231\Tracking Data";
libname e "Z:\SJShare\SJCOMMON\ECC\SJLife\SJLIFE Data Freeze\2 Final Data SJLIFE\20181231\Event Data";
libname k "Z:\SJShare\SJCOMMON\ECC\Ness Research Team\SJLIFE Data Cleaning\20181231\Final Appended Files";
options fmtsearch=(f);
run;

data part;
	set t.participants (where=(newstatus=3)); run;

data qct;
	set c.qct;
run;

data qct1 (where=(qctdt>=baselinevisitstart));
	merge qct part (in=a);
by mrn; 
if a;
run;

proc sql;
CREATE TABLE QCT_NUM AS
SELECT
	mrn, qctdt, zqct,
	COUNT(*)
FROM
	qct1
GROUP BY
	mrn;

proc sql;
create table qct_numtot as
	select unique mrn, _TEMG001
from qct_num;

data have_qct (where=(qctnum ne 0));
	merge part (in=a) qct_numtot (in=b rename=(_TEMG001=qctnum));
by mrn;
if a;
if qctnum=. then do qctnum=0; end;
run;

data function;
	set k.function_combo_all;
run;

proc sql;
CREATE TABLE FXN_NUM AS
SELECT
	mrn, assmntdate,
	COUNT(*)
FROM
	function
GROUP BY
	mrn;

proc sql;
create table fxn_numtot as
	select unique mrn, _TEMG001
from fxn_num;

data have_fxn (where=(fxnnum ne 0));
	merge part (in=a) fxn_numtot (in=b rename=(_TEMG001=fxnnum));
by mrn;
if a;
if fxnnum=. then do fxnnum=0; end;
run;

/*For lean mass*/
/*Need to adjust for amputation*/
data ampu; 
    set c.surgery;
		array surg[4] ICD9_1 ICD9_2 ICD9_3 ICD9_4;
		  do i = 1 to 4;
			 if surg[i] in ("84.01" "84.02" "84.03" "84.04" "84.05" "84.06" "84.07" "84.08" "84.09"
					"84.11" "84.12" "84.13" "84.14" "84.15" "84.16" "84.17" "84.18" "84.19") then amp=1;
			 end;
		drop i;
		if amp=1;
		/*drop ICD9_2 ICD9_3 ICD9_4; * because no amp in ICD9_2 and later ones; */
run;

/* ADJUSTMENT VALUES */
data ampu_ICD;
    set ampu(in=a);
	by mrn;
            retain count_ICD percent_ICD_1 percent_ICD_2 percent_ICD_3 percent_ICD_4;
            if first.mrn=1 then do count_ICD=0;
                  percent_ICD_1=0; percent_ICD_2=0; percent_ICD_3=0; percent_ICD_4=0; 
            end;
            count_ICD=count_ICD+1;
array surg[4] ICD9_1 ICD9_2 ICD9_3 ICD9_4;
array pct[4] percent_ICD_1 percent_ICD_2 percent_ICD_3 percent_ICD_4;
  do i = 1 to 4;
            if surg[i]=  "84.01" then  pct[i]=0; 
            if surg[i]=  "84.02" then  pct[i]=0; 
            if surg[i]=  "84.03" then  pct[i]=0; 
            if surg[i]=  "84.05" then  pct[i]=2.3; 
            if surg[i]=  "84.06" then  pct[i]=2.3; 
            if surg[i]=  "84.07" then  pct[i]=3.7; 
            if surg[i]=  "84.08" then  pct[i]=5; 
            if surg[i]=  "84.09" then  pct[i]=5;
            if surg[i]=  "84.11" then  pct[i]=0; 
            if surg[i]=  "84.12" then  pct[i]=0;  
            if surg[i]=  "84.13" then  pct[i]=1.5; 
            if surg[i]=  "84.15" then  pct[i]=3.7; 
            if surg[i]=  "84.16" then  pct[i]=5.9; 
            if surg[i]=  "84.17" then  pct[i]=11;
            if surg[i]=  "84.18" then  pct[i]=16; 
            if surg[i]=  "84.19" then  pct[i]=16; 
end;
percent_ICD=max(of percent_ICD_1 percent_ICD_2 percent_ICD_3 percent_ICD_4);
/*These I had to look at because I had to decide if it was a revision or an amputation
of a different limb*/
/*Delete the first of those that are a revision - take the bigger percentage*/
if mrn in (2497,10482,10763,11252,11918,12501,12685,15691) and count_ICD=1 then delete;
if mrn=15691 and count_ICD=2 then delete;
run;
/*Sum the rest*/
proc means data=ampu_ICD sum noprint;
	class mrn;
var percent_icd;
output out=ampu_adjust sum=ampu_adjust; run; 
/*Drop the total*/
data ampu_adjust (keep=mrn ampu_adjust);
	set ampu_adjust (where=(mrn ne .)); run;
/*Here is the code that I used to figure out whetehr or not the 
procedure was a revision or a different limb 
proc sql;
	create table twoamp
as select unique mrn from ampu_icd
where count_icd>1;

data ampcheck;
	merge twoamp (in=a) ampu_icd ;
by mrn;
if a;
run;

proc print data=ampcheck;
var mrn surgproc surgdt percent_icd count_icd; run;*/

data dexa (keep=mrn dexadt FFMadj Eweight compfile);
	merge c.dexa_wholebody ampu_adjust(in=c);
compfile=1;
if c then do;
	FFMadj=(wbtot_lean-head_lean)/(1-ampu_adjust/100)/1000; 
 	Eweight=WBTOT_MASS/(1-ampu_adjust/100)/1000; end;
else do; FFMadj=(wbtot_lean-head_lean)/1000; Eweight=WBTOT_MASS/1000; end;
run;

data basicd (keep=mrn dob gender race racegrp hispanic);
	set c.demographics;
run;

data skinf (keep=mrn assmntdate FFMadj Eweight hgt bmiadj race2 compfile age gender);
	merge basicd (keep=mrn gender dob racegrp hispanic) function (in=a);
by mrn;
if a;
compfile=2;
Eweight=wgt;
age=(assmntdate-dob)/365.25;
if gender="Male" then do;
		chest=mean(of sfchest1 sfchest2);
		abdomen=mean(of sfabdom1 sfabdom2);
		thigh=mean(of sfthigh1m sfthigh2m);
		sskin = sum(of chest abdomen thigh);
		bdensity = 1.10938- (0.0008267*sskin) + (0.0000016*(sskin**2)) - (0.0002574*age);
		Bfat= (457/bdensity)-414.2; *first formula from ACSM's guidelines for Exercise testing and prescription, 8th edition - same formula for male and female;
		FFMadj = Eweight-(Eweight*bfat/100)/1000;
	end;

	if gender="Female" then do;
		triceps=mean(of sftric1 sftric2);
		suprai=mean(of sfiliac1 sfiliac2);
		thigh=mean(of sfthigh1f sfthigh2f);
		sskin=sum(of triceps suprai thigh);
		bdensity=1.099421-(0.0009929*sskin)+ (0.0000023*(sskin**2))-(0.0001392*age);
		Bfat = (457/bdensity)-414.2; *first formula from ACSM's guidelines for Exercise testing and prescription, 8th edition - same formula for male and female;
		FFMadj = Eweight-(Eweight*bfat/100)/1000;
	end;

if hispanic="Hispanic/Latino" then do race2="Mexican"; end;
else if hispanic ne "Hispanic/Latino" then do;
	if racegrp in ("","White","Other","Unknown") then race2="White";
	else if racegrp="Black" then race2="Black"; end;
run;


proc sort data=dexa; by mrn dexadt; run;
proc sort data=skinf; by mrn assmntdate; run;

data concat;
	set dexa (rename=(dexadt=assmntdate) in=in1) skinf (in=in2);
	by mrn assmntdate; 
	retain assmntdate1 assmntdate2 ffmadj1 ffmadj2 eweight1 eweight2;
	if first.mrn then do assmntdate1=.; assmntdate2=.; ffmadj1=.; ffmadj2=.; eweight1=.; eweight2=.; end;
	if first.mrn and in1=0 then put "***NO DEXA FOR "mrn=;
	if first.mrn and in2=0 then put "***NO SKINF FOR "mrn=;
	file=2-in1;
	if in1 then assmntdate1=assmntdate; else assmntdate2=assmntdate;
	if in1 then ffmadj1=ffmadj; else ffmadj2=ffmadj;
	if in1 then eweight1=eweight; else eweight2=eweight;
datedif=abs(assmntdate1-assmntdate2);
 if file=2 then output; 
run;

data lean;
	set concat (where=(age ge 18)); /*Only include adults here*/
by mrn;
   * group age for calulcation of lean mass to use for deriving fmetab from normatitive data;
	if age ne . then do;
		if age lt 25 then ager=20;
		else if age lt 30 then ager=25;
		else if age lt 35 then ager=30;
		else if age lt 40 then ager=35;
		else if age lt 45 then ager=40;
		else if age lt 50 then ager=45;
		else if age lt 55 then ager=50;
		else if age lt 60 then ager=55;
		else if age lt 65 then ager=60;
		else if age lt 70 then ager=65;
		else if age lt 75 then ager=70;
		else if age lt 80 then ager=75;
		else if age lt 85 then ager=80;
		else if age ge 85 then ager=85;
	end;	

	if FFMadj1 ne .  and . < datedif <=365 then do;
		LMadj=FFMadj1/((HGT/100)**2); *---> lean mass formula dexa; end;
	if FFMadj1 =. or datedif >365 then do;
		LMadj=FFMadj2/((HGT/100)**2); *---> lean mass formula skin folds; end;
	

if gender = "Male" then sex=1;
else if gender = "Female" then sex=2;
run;

PROC IMPORT OUT= WORK.norms
   DATAFILE= "Z:\SJShare\SJCOMMON\ECC\Ness Research Team\CarrieUDrive\AnalyticalProjects\Frailty_2019\RecentFrailtyCode\LeanMass_NHANES.xlsx" 
   DBMS=xlsx; 		     
RUN;

data norms1 (drop = sex   rename=(sex1=sex age=ager race=race2));
set norms;
if sex = "Male" then sex1=1;else sex1=2;
run;

proc sort data=norms1;by sex race2 ager;run;
proc sort data=lean;by sex race2 ager;run;

data rellean;
	merge lean (in=a) norms1;
by sex race2 ager;
		if a = 1;
		zREE= M*((LMadj/M)**L-1)/(L*sd);
		if . lt zree lt -1.0 then fmetab=1;  * we determined cut-off by looking at the distribution (see code fr_norms.sas); *updated to -1.0 per Kiri 06.02.2016;
		else if zree ge -1.0 then fmetab=0;
		label fmetab='low lean mass';
/*keep mrn assmntdate fmetab sex;*/
drop ager m sd l;
if .<bmiadj<18.5 then fmetab=1; *per kiri email on 07.26.2018;
run;

data walk (where=(age ge 18));
	merge function (keep=mrn assmntdate hgt ppt7 comments sxmwd in=a) basicd;
by mrn;
if a;
age=(assmntdate-dob)/365.25;
	if gender="Male" then do;
			if . < HGT <=173 then do;
				if PPT7 >=23.3 /*or PPT7=.*/ then fwalk=1;
				else if . < PPT7 < 23.3 then fwalk=0;
			end;
	   		if HGT >173 then do;
				if PPT7>=20 /*or PPT7=.*/ then fwalk=1;
				else if . < PPT7 < 20.0 then fwalk=0;
			end;
	  	end;

	if gender="Female" then do;
			if . < HGT <=159 then do;
				if PPT7 >=23.3 /*or PPT7=.*/ then fwalk=1;
				else if . < PPT7 < 23.3 then fwalk=0;
			end;
	    if HGT >159 then do;
			if PPT7>=20 /* or PPT7=.*/ then fwalk=1;
			else if . < PPT7 < 20.0 then fwalk=0;
			end;
	  	end;

wtime=(PPT7/50)*15;
/*Medical reasons to walk slow so not done*/
if mrn in (5464,11570,11063,8668,3555,7904,9915,615,769,3649,3969,4095,4274,4616,4699,6326,7863,
7932,8292,8542,8678,8729,9389,9591,9790,10361,10438,10620,10764,11430,12133,12397,12477,12501,
12632,14681,14852,15266,15357,16090,19373,16001,17099,10116,9073,6528,9336,10004,10178,10957,12775,
14412,17886,21381,22157,2582,4034,12334,13186,1317,1337,5651,6660,7674,8058,8153,8915
9456,9494,9850,10321,10435,10628,10715,10945,11128,11460,11806,11949,12074,12451,12803,13004,13650
13817,13903,13984,14416,14434,14671,14983,15053,15616,15967,16140,18110,18197,19720,20939,21002,
21835,22887,29678,30819) 
then do fwalk=1; end; /*Medical reason to not do walk*/
if mrn in (7278,3601,14117,15026,15290) then do fwalk=0; end; /*Did not do PPT7 but 6mw fine - other had missing hgt but is fine*/
if mrn=14065 then do; fwalk=0; end; *functional data not entered - looked up in chart - CH 06.02.2016;
if fwalk=. and sxmwd>399 then do fwalk=0; end; /*use six min walk to replace some missing*/
if fwalk=. and .< sxmwd<=399 then do fwalk=1; end;
label fwalk='slowness';
run;

data grip (where=(age ge 18));
	merge function (keep=mrn assmntdate grip_lf1 grip_lf2 grip_rt1 grip_rt2 bmiadj comments in=a) basicd;
by mrn;
if a;
age=(assmntdate-dob)/365.25;
peakleft=max (of grip_lf1 grip_lf2); *taking peak of left hand grip;
peakright=max (of grip_rt1 grip_rt2); *taking peak of right hand grip;
grip=mean (of peakleft peakright); *now taking the mean of these;

if gender="Male" then do;
	if . < bmiadj <=24 then do;
		if . < grip <=29 then fhand=1;
		else if grip >29 then fhand=0;
				end;
	if 24 < bmiadj <=26 then do;
		if . < grip <=30 then fhand=1;
		else if grip >30 then fhand=0;
				end;
	if 26 < bmiadj <=28 then do;
		if . < grip <=30 then fhand=1;
		else if grip >30 then fhand=0;
				end;
	if 28 < bmiadj  then do;
		if . < grip <=32 then fhand=1;
		else if grip >32 then fhand=0;
				end;
		  end;	

	if gender="Female" then do;
		if . < bmiadj <=23 then do;
			if . < grip <=17 then fhand=1;
			else if grip >17 then fhand=0;
				end;
		if 23 < bmiadj <=26 then do;
			if . < grip <=17.3 then fhand=1;
			else if grip >17.3 then fhand=0;
				end;
		if 26 < bmiadj <=29 then do;
			if . < grip <=18 then fhand=1;
			else if grip >18 then fhand=0;
				end;
		if 29 < bmiadj  then do;
			if . < grip <=21 then fhand=1;
			else if grip >21 then fhand=0;
				end;
		  end;	
if mrn in (6803,11063,8668,769,3649,3765,3809,4274,7863,8678,9073,9373,9389,9494,9591,9613,10116,10361,
12267,12477,14852,16090, 615, 10957,5140,5651,6660,7278,8058,8914,9047) then fhand=1; /*Medical reason*/
if mrn in (8995, 3601,9963,11128,12451,12803,13903,13984,15967,19044,19720,20939,23461,25069,26100,
26495,29296,29604,29692,30311,31504/*615 - removed due to determined as 1 per Kendra previous code*/) then fhand=0; /*Missed hand grip but normal mTNS*/
if mrn in (14065,15290) then fhand=0;
label fhand='weakness';
run;

Data gripn (where=(age ge 18));
	merge function (keep=mrn assmntdate grip_lf1 grip_lf2 grip_rt1 grip_rt2 bmiadj comments in=a) basicd;
by mrn;
if a;
age=(assmntdate-dob)/365.25;
		* group age;
	if age ne . then do;
		if age lt 20 then ager=19;
		else if age lt 30 then ager=29;
		else if age lt 40 then ager=39;
		else if age lt 50 then ager=49;
		else if age lt 60 then ager=59;
		else if age lt 70 then ager=69;
		else if age ge 70 then ager=70;
	end;

pkgrip=max(of grip_lf1 grip_lf2)+max(of grip_rt1 grip_rt2); *do not want to add if both values not available - this was how norms were calculated;
if pkgrip=. then do;
	if grip_lf1=. and grip_lf2=. then pkgrip=sum(of grip_rt1 grip_rt2);
	if grip_rt1=. and grip_rt2=. then pkgrip=sum(of grip_lf1 grip_lf2);
end;
pkgriplb=pkgrip*2.2; *convert to pounds;
run;

* Read in normative data from excel spreadsheet - ;

PROC IMPORT OUT= WORK.gripnorms
DATAFILE= "Z:\SJShare\SJCOMMON\ECC\Ness Research Team\CarrieUDrive\AnalyticalProjects\Frailty_2019\Final\Grip_NHANES_2.xlsx" 
DBMS=xlsx;
RUN;
	
data gripnorms;set gripnorms;if missing(age) then delete;run;
	
data gripnorms1 (drop = sex   rename=(sex1=sex age=ager));
			set gripnorms;
			if sex = "Females" then gender="Female";
			else gender="Male";
SD=SE*sqrt(n); *get sd value into data;
run;

proc sort data=gripnorms1;by gender ager;run;

proc sort data=gripn;by gender ager;run;

data handnorms;
	merge gripn(in=a) gripnorms1;
by gender ager;
		if a = 1;
		zgrip=(pkgriplb-M)/sd ;
		if . lt zgrip lt -1.5 then fgripn=1;  
		else if zgrip ge -1.5 then fgripn=0;
		label fgripn='weak grip via norm';
/*adding in determinations from code above since these missed hand grip*/;
if mrn in (6803,11063,8668,769,3649,3765,3809,4274,7863,8678,9073,9373,9389,9494,9591,9613,10116,10361,
12267,12477,14852,16090, 615, 10957,5140,5651,6660,7278,8058,8914,9047) then fgripn=1; /*Medical reason*/
if mrn in (8995, 3601,9963,11128,12451,12803,13903,13984,15967,19044,19720,20939,23461,25069,26100,
26495,29296,29604,29692,30311,31504/*615 - removed due to determined as 1 per Kendra previous code*/) then fgripn=0; /*Missed hand grip but normal mTNS*/
if mrn in (14065,15290) then fgripn=0;
run;

/*Merge the data from the lab*/
proc sort data=rellean; by mrn assmntdate; run;
proc sort data=walk; by mrn assmntdate; run;
proc sort data=grip; by mrn assmntdate; run;
proc sort data=handnorms; by mrn assmntdate; run;

data frail_lab (keep=mrn assmntdate age gender sex race2 hgt wgt bmiadj fmetab fwalk fhand fgripn comments studypop where=(age ge 18));
	merge rellean walk grip handnorms function (keep=mrn assmntdate wgt studypop);
by mrn assmntdate;
run;

/*Put functional data by core visit*/

data visitf (keep=mrn datevisitstart corenum);
	set t.visits (where=(corevisit=1));
run;

data visit1f;
	merge visitf (where=(corenum=1) in=a) frail_lab (where=(studypop="Survivor") in=b);
by mrn;
if b;
datedif=abs(datevisitstart-assmntdate);
run;

proc sort data=visit1; by mrn datedif; run;

data lab1;
	set visit1f;
by mrn;
if first.mrn;
run;

data visit2f;
	merge visitf (where=(corenum=2) in=a) frail_lab (where=(studypop="Survivor") in=b);
by mrn;
if a and b;
datedif=abs(datevisitstart-assmntdate);
run;

proc sort data=visit2f; by mrn datedif; run;

data lab2;
	set visit2f;
by mrn;
if first.mrn;
run;

data visit3f;
	merge visitf (where=(corenum=3) in=a) frail_lab (where=(studypop="Survivor") in=b);
by mrn;
if a and b;
datedif=abs(datevisitstart-assmntdate);
run;

proc sort data=visit3f; by mrn datedif; run;

data lab3;
	set visit3f;
by mrn;
if first.mrn;
run;

data visit_cf;
	merge frail_lab part (where=(studypop="Community Control") in=b);
by mrn;
if first.mrn;
if b;
datevisitstart=baselinevisitstart;
run;

proc sort data=lab1; by mrn datevisitstart; run;
proc sort data=lab2; by mrn datevisitstart; run;
proc sort data=lab3; by mrn datevisitstart; run;
proc sort data=visit_cf; by mrn datevisitstart; run;

data lab_core (keep=mrn corenum datevisitstart assmntdate age gender sex 
race2 hgt wgt bmiadj fmetab fwalk fhand fgripn comments studypop datedif_lab);
	merge lab1(in=a) lab2(in=b) lab3(in=c) visit_cf (in=d);
by mrn datevisitstart;
datedif_lab=datedif;
run;

/*Get pa and vitality from questionnaires*/
/*Physical activity data*/

data partpa (keep=mrn baselinevisitstart studypop newstatus);
	set t.tracking (where=(newstatus in (3,5)));
run;

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

data phys_b (keep=mrn datecomp mvpawk mets meetcdc sedentary adolescent);
	set s.adolescent_healthhabits (where=(survey in 
("11-17 Years of Age Parent Report","5-10 Years of Age Parent Report")));
/*Get those with mvpa minutes first*/
if vpa10=1 and vpadays=. then do vpadays=1; end;
if vpa10=1 and vpamin=. then do vpamin=10; end;
if vpamin>360 then do vpamin=360; end; /*Cap six hours per day*/
if vpa10=1 then do wvpa=vpadays*vpamin; end;
if vpa10=2 then do wvpa=0; end;
if wvpa=. then do;
	if vpa10 in (.,2) and activity=2 then do wvpa=0; end;
	if vpa10 in (.,2) and activity=1 then do wvpa=0; end;
end;

if mpa10=. and (mpadays ne . or mpamin ne .) then do mpa10=1; end;
if mpa10=1 and mpadays=. then do mpadays=1; end;
if mpa10=1 and mpamin=. then do mpamin=10; end;
if mpamin>360 then do mpamin=360; end; /*Cap six hours per day*/
if mpa10=1 then do wmpa=mpadays*mpamin; end;
if mpa10=2 then do wmpa=0; end;
if wmpa=. then do;
	if mpa10 in (.,2) and activity=2 then do wmpa=0; end;
	if wvpa ne . then do wmpa=0; end;
end;

mvpawk=(wmpa)+(wvpa*2); if mvpawk>2520 then mvpawk=2520; /*cap six hours per day*/
mets=mvpawk*3;
if . < mvpawk <300 then meetcdc=0;
else if mvpawk >=300 then meetcdc=1;
if activity=2 and mvpawk=0 then sedentary=1; 
else if activity=1 or mvpawk > 0 then sedentary=0;
adolescent=1;
run;

proc append base=phys_b data=phys_a; run;

proc sort data=phys_b; by mrn datecomp; run;

data visit (keep=mrn datevisitstart corenum);
	set t.visits (where=(corevisit=1));
run;

data visit1;
	merge visit (where=(corenum=1) in=a) phys_b partpa (where=(studypop="Survivor") in=b);
by mrn;
if b;
datedif=abs(datevisitstart-datecomp);
run;

proc sort data=visit1; by mrn datedif; run;

data activity1;
	set visit1;
by mrn;
if first.mrn;
run;

data visit2;
	merge visit (where=(corenum=2) in=a) phys_b partpa (where=(studypop="Survivor") in=b);
by mrn;
if a and b;
datedif=abs(datevisitstart-datecomp);
run;

proc sort data=visit2; by mrn datedif; run;

data activity2;
	set visit2;
by mrn;
if first.mrn;
run;

data visit3;
	merge visit (where=(corenum=3) in=a) phys_b partpa (where=(studypop="Survivor") in=b);
by mrn;
if a and b;
datedif=abs(datevisitstart-datecomp);
run;

proc sort data=visit3; by mrn datedif; run;

data activity3;
	set visit3;
by mrn;
if first.mrn;
run;

data visit_c;
	merge phys_b part (where=(studypop="Community Control") in=b);
by mrn;
if first.mrn;
if b;
datevisitstart=baselinevisitstart;
run;

proc sort data=activity1; by mrn datevisitstart; run;
proc sort data=activity2; by mrn datevisitstart; run;
proc sort data=activity3; by mrn datevisitstart; run;
proc sort data=visit_c; by mrn datevisitstart; run;

data activity (keep=mrn corenum datevisitstart datecomp mvpawk mets meetcdc sedentary adolescent datedif_pa);
	merge activity1(in=a) activity2(in=b) activity3(in=c) visit_c (in=d);
by mrn datevisitstart;
datedif_pa=datedif;
run;

proc sort data=lab_core; by mrn datevisitstart; run;
proc sort data=activity; by mrn datevisitstart; run;

data lab_act;
	merge lab_core (in=a) activity (in=b);
by mrn datevisitstart;
if a;
	kcal_week=mets*(3.5*wgt)/200;
			
			if gender="Male" then do;
				if . < kcal_week < 383 then fkcal=1;
				else if kcal_week >=383 then fkcal=0;
			end;

			if gender="Female" then do;
				if . < kcal_week < 270 then fkcal=1;
				else if kcal_week >=270 then fkcal=0;
			end;
			label fkcal='low activity level';
run;

proc freq data=lab_act;
	table fmetab fkcal fhand fgripn fwalk/missing;
run;

/*Get qol data*/
data qol (keep=mrn datecomp survey fvital);
	set s.adult_sf36;
		if . < vtnorm <=40 then fvital=1;
		else if vtnorm > 40 /*or vtnorm=.*/ then fvital=0;
		label fvital='poor endurance/energy';
run;

data visit (keep=mrn datevisitstart corenum);
	set t.visits (where=(corevisit=1));
run;

data visit1q;
	merge visit (where=(corenum=1) in=a) qol (where=(survey ne "Adult Core Version - Control") in=b);
by mrn;
if b;
datedif=abs(datevisitstart-datecomp);
run;

proc sort data=visit1; by mrn datedif; run;

data qol1;
	set visit1q;
by mrn;
if first.mrn;
run;

data visit2q;
	merge visit (where=(corenum=2) in=a) qol (where=(survey ne "Adult Core Version - Control") in=b);
by mrn;
if a and b;
datedif=abs(datevisitstart-datecomp);
run;

proc sort data=visit2q; by mrn datedif; run;

data qol2;
	set visit2q;
by mrn;
if first.mrn;
run;

data visit3q;
	merge visit (where=(corenum=3) in=a) qol (where=(survey ne "Adult Core Version - Control") in=b);
by mrn;
if a and b;
datedif=abs(datevisitstart-datecomp);
run;

proc sort data=visit3q; by mrn datedif; run;

data qol3;
	set visit3q;
by mrn;
if first.mrn;
run;

data visit_cq;
	merge qol part (where=(studypop="Community Control") in=b);
by mrn;
if first.mrn;
if b;
datevisitstart=baselinevisitstart;
run;

proc sort data=qol1; by mrn datevisitstart; run;
proc sort data=qol2; by mrn datevisitstart; run;
proc sort data=qol3; by mrn datevisitstart; run;
proc sort data=visit_cq; by mrn datevisitstart; run;

data vitality (keep=mrn corenum datevisitstart datecomp fvital datedif_q);
	merge qol1(in=a) qol2(in=b) qol3(in=c) visit_cq (in=d);
by mrn datevisitstart;
datedif_q=datedif;
run;

data lab_vit;
	merge lab_act (in=a) vitality (in=b);
by mrn datevisitstart;
if a;
if studypop="Community Control" then do datedif_q=0; datedif_pa=0; datedif_lab=0; corenum=1; end;
run;

proc freq data=lab_vit;
	table fmetab fkcal fhand fgripn fwalk fkcal fvital/missing;
run;

proc format;
	value sex 1="male" 2="female";
	value ynzf 0="no" 1="yes";
run;

data frail_core (keep=mrn studypop datevisitstart corenum age sex wgt hgt bmiadj frail_num_fried frail_num_handn fmetab fwalk fhand fgripn fkcal fvital 
mvpawk meetcdc sedentary datedif_lab datedif_pa datedif_q comments);
	retain mrn studypop datevisitstart corenum age sex wgt hgt bmiadj frail_num_fried frail_num_handn fmetab fwalk fhand fgripn fkcal fvital 
mvpawk meetcdc sedentary datedif_lab datedif_pa datedif_q comments;
set lab_vit;
frail_num_fried=sum(of fmetab fwalk fhand fkcal fvital);
frail_num_handn=sum(of fmetab fwalk fgripn fkcal fvital);
label frail_num_fried="number of frailty components Fried criteria";
label frail_num_handn="number of frailty componsnts using hand grip normative data";
label fmetab="meets criteria for low lean mass";
label fwalk="meets criteria for slow walking";
label fhand="meets criteria for weakness Fried";
label fgripn="meets criteria for weakness CDC norms";
label fkcal="meets criteria for low energy expenditure";
label mvpawk="weekly minutes of moderate or vigorous activity";
label meetcdc="at least 150 minutes mvpa/wk";
label sedentary="no activity in past month";
label datedif_lab="proximity of functional value to core visit";
label datedif_pa="proximity of reported pa value to core visit";
label datedif_q="proximity of reported vitality to core visit";
label corenum="core visit number"; 
format sex sex.;
format fmetab fwalk fhand fgripn fkcal fvital 
meetcdc sedentary ynzf.;
run;

data "Z:\SJShare\SJCOMMON\ECC\Ness Research Team\Goodenough_Ness\P2- Low BMD\frail_core.sas7bdat";
	set frail_core;
run;
