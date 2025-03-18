% ------------------------------------------------------------------------------
%
%                           function invjday
%
%  this function finds the year, month, day, hour, minute and second
%    given the julian date. tu can be ut1, tdt, tdb, etc.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    jd          - julian date                    days from 4713 bc
%
%  outputs       :
%    year        - year                           1900 .. 2100
%    mon         - month                          1 .. 12
%    day         - day                            1 .. 28,29,30,31
%    hr          - hour                           0 .. 23
%    min         - minute                         0 .. 59
%    sec         - second                         0.0 .. 59.999
%
%  locals        :
%    days        - day of year plus fractional
%                  portion of a day               days
%    tu          - julian centuries from 0 h
%                  jan 0, 1900
%    temp        - temporary real values
%    leapyrs     - number of leap years from 1900
%
%  coupling      :
%    days2mdhms  - finds month, day, hour, minute and second given days and year
%
%  references    :
%    vallado       2007, 208, alg 22, ex 3-13
%
% [year,mon,day,hr,min,sec] = invjday ( jd );
% -----------------------------------------------------------------------------

function [year,mon,day,hr,min,sec] = invjday ( jd )

     % ----------------- find year and days of the year ---------------
     temp   = jd-2415019.5;
     tu     = temp / 365.25;
     year   = 1900 + floor( tu );
     leapyrs= floor( ( year-1901 )*0.25 );
%     days   = temp - ((year-1900)*365.0 + leapyrs ) + 0.00000000001; % nudge by 8.64x10-7 sec to get even outputs
     days   = temp - ((year-1900)*365.0 + leapyrs );

     % ------------ check for case of beginning of a year -------------
     if days < 1.0
         year   = year - 1;
         leapyrs= floor( ( year-1901 )*0.25 );
         days   = temp - ((year-1900)*365.0 + leapyrs );
     end

     % ------------------- find remaining data  -----------------------
     [mon,day,hr,min,sec] = days2mdh( year,days );
%     sec= sec - 0.00000086400;

% ------------------------------------------------------------------------------
%
%                           function days2mdh
%
%  this function converts the day of the year, days, to the equivalent month
%    day, hour, minute and second.
%
%  author        : david vallado                  719-573-2600   22 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    year        - year                           1900 .. 2100
%    days        - julian day of the year         0.0  .. 366.0
%
%  outputs       :
%    mon         - month                          1 .. 12
%    day         - day                            1 .. 28,29,30,31
%    hr          - hour                           0 .. 23
%    minute      - minute                         0 .. 59
%    sec         - second                         0.0 .. 59.999
%
%  locals        :
%    dayofyr     - day of year
%    temp        - temporary extended values
%    inttemp     - temporary integer value
%    i           - index
%    lmonth(12)  - integer array containing the number of days per month
%
%  coupling      :
%    none.
%
% [mon,day,hr,minute,sec] = days2mdh ( year,days);
% -----------------------------------------------------------------------------

function [mon,day,hr,minute,sec] = days2mdh ( year,days);

        % --------------- set up array of days in month  --------------
        for i= 1 : 12
            lmonth(i) = 31;
            if i == 2
                lmonth(i)= 28;
              end;
            if i == 4 | i == 6 | i == 9 | i == 11
                lmonth(i)= 30;
              end;
        end

        dayofyr= floor(days );

        % ----------------- find month and day of month ---------------
        if rem(year-1900,4) == 0
            lmonth(2)= 29;
          end

        i= 1;
        inttemp= 0;
        while ( dayofyr > inttemp + lmonth(i) ) & ( i < 12 )
            inttemp= inttemp + lmonth(i);
            i= i+1;
          end

        mon= i;
        day= dayofyr - inttemp;

        % ----------------- find hours minutes and seconds ------------
        temp= (days - dayofyr )*24.0;
        hr  = fix( temp );
        temp= (temp-hr) * 60.0;
        minute = fix( temp );
        sec = (temp-minute) * 60.0;

