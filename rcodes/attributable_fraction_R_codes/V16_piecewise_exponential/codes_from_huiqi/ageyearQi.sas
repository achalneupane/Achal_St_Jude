/*
Feb 2016, Qi works on multiple events. For the resulted dataset from ageyear macro, I need to chop the small segment whenever there is an event (time depend variable).
In order to do this, I need the start and end time of each segment in the resulted dataset of ageyear macro. So I will modify the code so it outputs the start and end 
date of each segment as well.
*/


/* ************************************************************************ */
/* AGECAL macro                                                             */
/*                                                                          */
/* Calculates exact ages using datetime functions.                          */

%macro agecal(agevar,        /* age variable */
              startdt,       /* start date   */
              finishdt);     /* finish date  */

  if month(&startdt)<month(&finishdt) or
     (month(&startdt)=month(&finishdt) and day(&startdt)<=day(&finishdt))
   then &agevar=year(&finishdt)-year(&startdt);
  else if month(&startdt)>month(&finishdt) or
     (month(&startdt)=month(&finishdt) and day(&startdt)>day(&finishdt))
   then &agevar=year(&finishdt)-year(&startdt)-1;

%mend agecal;


/* ************************************************************************ */
/* PYR macro                                                                */
/*                                                                          */
/* Parameters FIRSTDT and LASTDT must be date variables                     */
/* Creates variable                                                         */
/*   PYEAR = (#days between month,day of LASTDT and month,day of FIRSTDT)   */
/*           / (# days in current year: YEAR1231-YEAR0101)                  */
/* and outputs the observation if PYEAR is non-zero.                        */

%macro pyr(firstdt,         /* start date */
           lastdt);         /* end date   */

  pyear=(mdy(month(&lastdt), day(&lastdt), calyear)
         - mdy(month(&firstdt), day(&firstdt), calyear))
        / (year1231-year0101);
	/*Qi added this part so the py macro not only ouput pyear but also the start and end date of each segment.*/
	d_start=mdy(month(&firstdt), day(&firstdt), calyear); 
	d_ending=mdy(month(&lastdt), day(&lastdt), calyear); 
  if pyear>0 then output;

%mend pyr;


/* ************************************************************************ */
/* PYR1 macro                                                               */
/*                                                                          */
/* Parameters FIRSTDT and LASTDT must be date variables                     */
/* Creates variable                                                         */
/*   PYEAR = ( 1 + (#days between month,day of LASTDT and month,day of FIRSTDT))   */
/*           / (# days in current year: YEAR1231-YEAR0101)                  */
/* and outputs the observation if PYEAR is non-zero.                        */
/* Used for final year at risk, if different to initial year at risk.       */

%macro pyr1(firstdt,         /* start date */
            lastdt);         /* end date   */

  pyear=(1 + mdy(month(&lastdt), day(&lastdt), calyear)
         - mdy(month(&firstdt), day(&firstdt), calyear))
        / (year1231-year0101);
	/*Qi added this part so the py macro not only ouput pyear but also the start and end date of each segment.*/
	d_start=mdy(month(&firstdt), day(&firstdt), calyear); 
	d_ending=mdy(month(&lastdt), day(&lastdt), calyear);   
  if pyear>0 then output;

%mend pyr1;


/* ************************************************************************ */
/* AGEYEAR macro                                                            */
/*                                                                          */
/* Parameters:                                                              */
/*   INSET:   input dataset                                                 */
/*   BIRDT:   birth date variable (in INSET)                                */
/*   STARTDT: cohort entry date variable (in INSET)                         */
/*   ENDDT:   cohort exit date variable (in INSET)                          */
/*   YEARS:   years-in-cohort variable (new variable)                       */
/*   TIMEST:  starting value of YEARS                                       */
/*   OUTSET:  output dataset                                                */
/* Macros used:                                                             */
/*   AGECAL, PYR                                                            */
/* Creates one record in the output dataset for each change in calendar     */
/* year (variable CALYEAR), person's age (variable AGEY) or number of years */
/* in the cohort (macro variable YEARS) between dates STARTDT and ENDDT.    */
/* Amended 9/04: fixed problem with end date=01/01/yyyy                     */

%macro ageyear(inset, birdt, startdt, enddt, years, timest, outset);

/*			%let inset=indata;
            %let   birdt=d_birth;
            %let   startdt=d5year;
            %let  enddt=d_end;
            %let  years=ysincedx;
            %let  timest=5;
            %let  outset=py1;	*/

  /* Determine minimum start year and maximum end year in dataset, and      */
  /* create macro variables STYEAR and ENYEAR.                              */
  
  data _null_;
  set &inset end=last;
  retain minst 2099 maxen 1800;
  if year(&startdt)<minst then minst=year(&startdt);
  if year(&enddt)>maxen then maxen=year(&enddt);
  if last then do;
    call symput('styear', minst);
    call symput('enyear', maxen);
  end;
  run;

  /*	%put &styear &enyear;	*/
  
  /* Create new dataset OUTSET */
  
  data &outset;
  set &inset;
  
  /* Calculate initial value of AGEY (= age at start date) and assign       */
  /* starting value TIMEST to years-in-cohort variable YEARS.               */

  %agecal(agey, &birdt, &startdt);
  &years=&timest;

  /* Create code for all years in dataset */
  
  %do y=%eval(&styear) %to %eval(&enyear);
    
    /* Create variables CALYEAR, YEAR0101 and YEAR1231                     */
    
    calyear=%eval(&y);
    year0101=mdy(1, 1, calyear);
    year1231=mdy(12, 31, calyear);
    
    /* Execute if patient was in the cohort during current CALYEAR         */
    
    if year(&startdt)<=calyear and year(&enddt)>=calyear then do;

      if %eval(&y)<year(&enddt) then do;

        if %eval(&y)=year(&startdt) then do;

    /* year<end year, year=start year, entry d/m later than birth d/m */
          if month(&startdt)>month(&birdt) or
             (month(&startdt)=month(&birdt) and day(&startdt)>=day(&birdt))
           then do;
            %pyr(&startdt, year1231);
          end;

    /* year<end year, year=start year, entry d/m earlier than birth d/m */
          else if month(&startdt)<month(&birdt) or
             (month(&startdt)=month(&birdt) and day(&startdt)<day(&birdt))
           then do;
            %pyr(&startdt, &birdt);
            agey=agey+1;
            %pyr(&birdt, year1231);
          end;
        end;

        else if %eval(&y)>year(&startdt) then do;
      
    /* year<end year, year>start year, entry d/m later than birth d/m */
          if month(&startdt)>month(&birdt) or
             (month(&startdt)=month(&birdt) and day(&startdt)>=day(&birdt))
           then do;
            %pyr(year0101, &birdt);
            agey=agey+1;
            %pyr(&birdt, &startdt);
            &years=&years+1;
            %pyr(&startdt, year1231);
          end;
        
    /* year<end year, year>start year, entry d/m earlier than birth d/m */
          else if month(&startdt)<month(&birdt) or
             (month(&startdt)=month(&birdt) and day(&startdt)<day(&birdt))
           then do;
            %pyr(year0101, &startdt);
            &years=&years+1;
            %pyr(&startdt, &birdt);
            agey=agey+1;
            %pyr(&birdt, year1231);
          end;
        end;
      end;
      
      else if %eval(&y)=year(&enddt) then do;
        if %eval(&y)=year(&startdt) then do;

    /* year=start year, year=end year, entry d/m later than birth d/m */
          if month(&startdt)>month(&birdt) or
             (month(&startdt)=month(&birdt) and day(&startdt)>=day(&birdt))
           then do;
            %pyr(&startdt, &enddt);
          end;

    /* year=start year, year=end year, entry d/m earlier than birth d/m */
          else if month(&startdt)<month(&birdt) or
             (month(&startdt)=month(&birdt) and day(&startdt)<day(&birdt))
           then do;

    /* ... and exit d/m later than birth d/m */
            if month(&enddt)>month(&birdt) or
              (month(&enddt)=month(&birdt) and day(&enddt)>=day(&birdt))
             then do;
              %pyr(&startdt, &birdt);
              agey=agey+1;
              %pyr(&birdt, &enddt);
            end;
          
    /* ... and exit d/m earlier than birth d/m */
            else if month(&enddt)<month(&birdt) or
              (month(&enddt)=month(&birdt) and day(&enddt)<day(&birdt))
             then do;
              %pyr(&startdt, &enddt);
            end;
          end;
        end;

        else if %eval(&y)>year(&startdt) then do;
      
    /* year>start year, year=end year, entry d/m later than birth d/m */
          if month(&startdt)>month(&birdt) or
             (month(&startdt)=month(&birdt) and day(&startdt)>=day(&birdt))
           then do;

    /* ... and birth d/m later than exit d/m */
             if month(&birdt)>month(&enddt) or
                (month(&birdt)=month(&enddt) and day(&birdt)>=day(&enddt))
              then do;
               %pyr1(year0101, &enddt);         /* changed from pyr to pyr1 9/10/04 */
             end;

    /* ... and entry d/m later than exit d/m */
             else if month(&startdt)>month(&enddt) or
                (month(&startdt)=month(&enddt) and day(&startdt)>=day(&enddt))
              then do;
               %pyr(year0101, &birdt);
               agey=agey+1;
               %pyr1(&birdt, &enddt);         /* changed from pyr to pyr1 9/10/04 */
             end;

    /* ... and exit d/m later than entry d/m */
             else if month(&startdt)<month(&enddt) or
                     (month(&startdt)=month(&enddt) and day(&startdt)<day(&enddt))
              then do;
               %pyr(year0101, &birdt);
               agey=agey+1;
               %pyr(&birdt, &startdt);
               &years=&years+1;
               %pyr1(&startdt, &enddt);         /* changed from pyr to pyr1 9/10/04 */
             end;
          end;
      
    /* year>start year, year=end year, birth d/m later than entry d/m */
          else if month(&startdt)<month(&birdt) or
             (month(&startdt)=month(&birdt) and day(&startdt)<day(&birdt))
           then do;

    /* ... and entry d/m later than exit d/m */
            if month(&startdt)>month(&enddt) or
               (month(&startdt)=month(&enddt) and day(&startdt)>=day(&enddt))
             then do;
              %pyr1(year0101, &enddt);         /* changed from pyr to pyr1 9/10/04 */
            end;

    /* ... and birth d/m later than exit d/m */
            else if month(&birdt)>month(&enddt) or
               (month(&birdt)=month(&enddt) and day(&birdt)>=day(&enddt))
             then do;
              %pyr(year0101, &startdt);
              &years=&years+1;
              %pyr1(&startdt, &enddt);         /* changed from pyr to pyr1 9/10/04 */
            end;

    /* ... and exit d/m later than birth d/m */
            else if month(&startdt)<month(&enddt) or
                    (month(&startdt)=month(&enddt) and day(&startdt)<day(&enddt))
             then do;
              %pyr(year0101, &startdt);
              &years=&years+1;
              %pyr(&startdt, &birdt);
              agey=agey+1;
              %pyr1(&birdt, &enddt);         /* changed from pyr to pyr1 9/10/04 */
            end;
          end;
        end;
      end;
    end;
  %end;
  drop year0101 year1231;
  run;
  
%mend ageyear;
