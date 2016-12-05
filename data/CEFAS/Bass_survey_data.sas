/********************************************************************
* Name of program: Beems_FSS_data.sas                               *
*                                                                   *
* Description: Code to extract length and count data from the Bass  *
* surveys: SOLENT & INAKTHAMES, YFS and BEEMS                       *
*                                                                   *
* Library location:                                                 *
* \\Lowfile1\lowestoftuserdata2\lr05\contract C3278                 *
*                                                                   *
* Creation date: 03/11/2008                                         *
*                                                                   *
* Created by: L Readdy                                              *
*                                                                   *
* Edit history                                                      *
* +---------------------------------------------------------------+ *
* | Programmer      |Date    |Description                         | *
* |---------------------------------------------------------------| *
* | L Readdy        |05/09/16|Adapted from Bass survey code       | *
* |                 |        |                                    | *
* +_________________+________+____________________________________+ *
*                                                                   *
********************************************************************/
/********************************************************************
* SOLENT survey has specific primary stations                       *
********************************************************************/
data DataStationLogs;
 set tblDataStationLogs;
 fldPrimeStation = upcase(fldPrimeStation);
 run;

proc sql;
 create table Bass_surveyS as
 	select a.fldSeriesName, a.fldStartDate, b.*, 
		c.fldGearCode, c.fldGearAdditional1, c.fldGearAdditional3, c.fldGearAdditional6, c.fldValidityCode,
		d.fldPrimeStationComment
	from (select fldSeriesName, fldCruiseName, fldStartDate from tblReferenceCruises 
			where fldSeriesName in ('CARLHELMAR','NWGFS','Q1SWBEAM','Q4SWIBTS','WCGFS')) as a 
	inner join
	(select fldICESRectangle, fldShotLonDecimalDegrees, fldShotLatDecimalDegrees, fldHaulLonDecimalDegrees, 
			fldHaulLatDecimalDegrees, fldCruiseName, fldCruiseStationNumber, fldDateTimeShot, fldTowDuration,
			fldPrimeStation, fldShotDepth, fldHaulDepth, fldDateTimeHaul 
	from DataStationLogs where fldTowDuration > 0) as b
	on a.fldCruiseName=b.fldCruiseName
	inner join
	(select * from tblDataGearDeployments /*where fldValidityCode ^='I'*/) as c
	on a.fldCruiseName=c.fldCruiseName
	& b.fldCruiseStationNumber=c.fldCruiseStationNumber

	left join
	(select fldPrimeStation, fldPrimeStationComment, fldSeriesName from tblSeriesPrimeStations) as d
	on b.fldPrimeStation = d.fldPrimeStation
	& a.fldSeriesName = d.fldSeriesName
;
quit;

Data Bass_surveyG (drop=fldGearAdditional3 fldPrimeStationComment);
 set Bass_surveyS;
 format Time HHMM5.;
 *where fldSeriesName in ('INAKTHAMES','SOLENT')
		or (fldSeriesName = 'BEEMS' 
					and fldGearCode in (101010101,101010301,107120401,107120201,107120101,107120301))
		or (fldSeriesName = 'YFS' 
					and fldValidityCode = 'V');
 /** extra removal YFS**/
 *if fldSeriesName = 'YFS' and substr(fldCruiseName,1,4) = 'ECST' and year(datepart(fldStartDate)) ^= '1999' then delete;
 /**/
 Month=month(datepart(fldDateTimeShot));
 Year=year(datepart(fldDateTimeShot));
 Day=Day(Datepart(fldDateTimeShot));
 Time=timepart(fldDateTimeShot);
 if fldSeriesName = 'BEEMS' and fldGearCode in (101010101,101010301) then Coverage = fldGearAdditional1;
 else if fldSeriesName = 'YFS' then Coverage = fldGearAdditional6;
 else Coverage =.;

 if fldSeriesName = 'YFS' then Area_name = fldGearAdditional3;
 else if fldSeriesName in ('INAKTHAMES','SOLENT') then Area_name = fldPrimeStationComment;
 run;

proc sql;
 create table Bass_surveyL as
 	select a.*, b.fldMainSpeciesCode, b.fldProcessCode, b.fldSex,
			c.fldCatchNumber, c.fldCatchWeight, c.fldMeasuringInterval,
			d.fldCategoryRaisedNumberAtLength, d.fldLengthGroup, 
			e.fldGearDescription,
			floor(d.fldLengthGroup/10) as Length,
			case when b.fldMainSpeciesCode ='BSE' then 'ESB' 
				 when b.fldMainSpeciesCode in ('CRC','CRH') then 'CRE' 
				 else b.fldMainSpeciesCode end as Species_Code
	from (select * from Bass_surveyG) as a
 
	left join
	(select * from tblDataCatchComponents) as b
	on a.fldCruiseName=b.fldCruiseName
	& a.fldCruiseStationNumber=b.fldCruiseStationNumber
	& a.fldGearCode=b.fldGearCode

	left join
	(select * from tblDataCategories) as c
	on a.fldCruiseName=c.fldCruiseName
	& a.fldCruiseStationNumber=c.fldCruiseStationNumber
	& a.fldGearCode=c.fldGearCode
	& b.fldMainSpeciesCode=c.fldMainSpeciesCode
	& b.fldSex=c.fldSex

	left join
	(select * from tblDataLengthSamples where fldCategoryRaisedNumberAtLength ^in (.,0)) as d
	on a.fldCruiseName=d.fldCruiseName
	& a.fldCruiseStationNumber=d.fldCruiseStationNumber
	& a.fldGearCode=d.fldGearCode
	& b.fldMainSpeciesCode=d.fldMainSpeciesCode
	& b.fldSex=d.fldSex
	& c.fldCategoryNumber=d.fldCategoryNumber
	
	left join
	(select * from tblReferenceMainGearCodes) as e
	on a.fldGearCode=e.fldGearCode
	;
	create table Bass_surveyN as
 	select a.*,	f.fldCommonName, f.fldScientificName, f.fldAlternateSpeciesCode			
	from (select * from Bass_surveyL) as a
	left join
	(select * from tblReferenceMainSpecies) as f
	on a.Species_Code=f.fldMainSpeciesCode
; 
	quit;

data Bass_SurveyS;
	format Numbers Best12.;
 set Bass_surveyN;
		 if fldProcessCode in ('WM','MO') then do;
			Numbers = fldCategoryRaisedNumberatLength;
			Weighed_observed_only = 'N';
		end;
		else if fldProcessCode in ('WC','CO') then do;
			Numbers = fldCatchNumber;
			Weighed_observed_only = 'N';
		end;
		else if fldProcessCode in  ('WO','OB') then do;
			Numbers = .;
			Weighed_observed_only = 'Y';
		end;
		
		if fldSeriesName = 'YFS' and fldProcessCode ^in ('MO','WM','') then do;
			Numbers =.;
			Weighed_observed_only ='';
			species_code = '';
			fldCommonName ='';
			fldScientificName ='';
			fldAlternateSpeciesCode ='';
			fldProcessCode ='';
		end;

		if numbers = 0 then do;
			Numbers =.;
			Weighed_observed_only ='';
			species_code = '';
			fldCommonName ='';
			fldScientificName ='';
			fldAlternateSpeciesCode = '';
		end;
		Area_name=translate(Area_name,":",",");

		if coverage ^=. then fldTowDuration=.;
 run;

Proc summary data = Bass_SurveyS nway missing;
 class fldSeriesName fldCruiseName Year Month Day Time fldPrimeStation Area_Name fldCruiseStationNumber fldDateTimeShot fldDateTimeHaul
		fldShotLonDecimalDegrees fldShotLatDecimalDegrees fldHaulLonDecimalDegrees fldHaulLatDecimalDegrees fldShotDepth fldHaulDepth
		fldValidityCode fldICESRectangle fldGearCode fldGearDescription fldTowDuration Coverage fldProcessCode 
		Weighed_observed_only fldMainSpeciesCode Species_Code fldCommonName fldScientificName fldMeasuringInterval fldLengthGroup;
 var Numbers;
 output out = Bass_SurveyZ (drop=_type_ _freq_) sum=;
 run;

data BASS;
 set Bass_Surveyz;
 *where fldSeriesName in ('INAKTHAMES','SOLENT');
 run;


proc export data=BASS
			outfile="C:\Temp\ForPaulD\WesternSurveys_V20160905.dat"
			DBMS=CSV replace;
			run;
