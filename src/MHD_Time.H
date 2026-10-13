#pragma once

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

inline bool gregorianLeapYear(int a_year)
{
    return a_year%4 == 0 && (a_year%100 != 0 || a_year%400 == 0);
}

inline int calendarMonthDays(int a_year, int a_month)
{
    static const int days[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    return days[a_month-1] + (a_month == 2 && gregorianLeapYear(a_year) ? 1 : 0);
}

inline long double calendarYearSeconds(int a_year)
{
    return (gregorianLeapYear(a_year) ? 366 : 365)*86400.0L;
}

struct CalendarEpoch
{
    int year = 0;
    long double seconds = 0.0L; // Seconds since January 1, UTC.

    double decimalYear() const
    {
        return static_cast<double>(year + seconds/calendarYearSeconds(year));
    }
};

// Date strings keep whole seconds exactly. Decimal years also bypass the
// shared parser's float conversion. A zero epoch retains the startup sentinel.
inline CalendarEpoch parseCalendarEpoch(const std::string& a_raw, const char* a_name)
{
    std::string value = a_raw.substr(0, a_raw.find('#'));
    auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) value = "0";
    else value = value.substr(first, value.find_last_not_of(" \t\r\n")-first+1);
    auto invalid = [&]() -> std::invalid_argument
    {
        return std::invalid_argument(std::string(a_name)
            + " must be a finite decimal year in [1, 10000), 0, or a valid UTC date (DD-MM-YYYYTHH:MM:SS).");
    };
    if (value.find('T') == std::string::npos)
    {
        try
        {
            size_t end = 0;
            long double year = std::stold(value, &end);
            if (end != value.size() || !std::isfinite(year)
                || (year != 0 && (year < 1 || year >= 10000))) throw invalid();
            int wholeYear = static_cast<int>(std::floor(year));
            return {wholeYear, (year-wholeYear)*calendarYearSeconds(wholeYear)};
        }
        catch (const std::exception&) { throw invalid(); }
    }
    if (value.size() != 19 || value[2] != '-' || value[5] != '-'
        || value[10] != 'T' || value[13] != ':' || value[16] != ':') throw invalid();
    auto digits = [&](int start, int count)
    {
        int number = 0;
        for (int i = start; i < start+count; ++i)
        {
            if (value[i] < '0' || value[i] > '9') throw invalid();
            number = 10*number + value[i]-'0';
        }
        return number;
    };
    int day = digits(0, 2), month = digits(3, 2), year = digits(6, 4);
    int hour = digits(11, 2), minute = digits(14, 2), second = digits(17, 2);
    if (year < 1 || month < 1 || month > 12 || day < 1
        || day > calendarMonthDays(year, month) || hour > 23
        || minute > 59 || second > 59) throw invalid();
    int dayOfYear = day-1;
    for (int m = 1; m < month; ++m) dayOfYear += calendarMonthDays(year, m);
    return {year, dayOfYear*86400.0L + hour*3600 + minute*60 + second};
}

inline double secondsBetweenEpochs(const CalendarEpoch& a_start, const CalendarEpoch& a_end)
{
    auto daysBeforeYear = [](int year) -> long long
    {
        long long y = year-1;
        return 365*y + y/4 - y/100 + y/400;
    };
    return static_cast<double>(
        (daysBeforeYear(a_end.year)-daysBeforeYear(a_start.year))*86400.0L
        + a_end.seconds-a_start.seconds);
}

inline CalendarEpoch advanceCalendarEpoch(CalendarEpoch a_epoch, long double a_seconds)
{
    a_epoch.seconds += a_seconds;
    while (a_epoch.seconds >= calendarYearSeconds(a_epoch.year))
    {
        a_epoch.seconds -= calendarYearSeconds(a_epoch.year);
        ++a_epoch.year;
    }
    while (a_epoch.seconds < 0)
    {
        --a_epoch.year;
        a_epoch.seconds += calendarYearSeconds(a_epoch.year);
    }
    return a_epoch;
}

inline std::string formatUTCEpoch(CalendarEpoch a_epoch)
{
    // Round before splitting the date, so carries cross month/year boundaries.
    a_epoch = advanceCalendarEpoch({a_epoch.year, 0}, std::floor(a_epoch.seconds+0.5L));
    int day = static_cast<int>(a_epoch.seconds/86400);
    int clock = static_cast<int>(a_epoch.seconds-day*86400.0L);
    int month = 1;
    while (day >= calendarMonthDays(a_epoch.year, month))
    {
        day -= calendarMonthDays(a_epoch.year, month);
        ++month;
    }
    char output[64];
    std::snprintf(output, sizeof(output), "%04d-%02d-%02dT%02d:%02d:%02d UTC",
        a_epoch.year, month, day+1, clock/3600, (clock%3600)/60, clock%60);
    return output;
}
