#pragma once

#include <limits>

class CMEInsertionSchedule
{
public:
    CMEInsertionSchedule(bool a_enabled, double a_insertionTime)
        : m_enabled(a_enabled),
          m_time(a_enabled ? a_insertionTime : std::numeric_limits<double>::infinity())
    {}

    double time() const { return m_time; }
    bool inserted() const { return m_inserted; }
    void setInserted(bool a_inserted) { m_inserted = a_inserted; }
    bool due(double a_time) const
    {
        return m_enabled && !m_inserted && a_time >= m_time;
    }

private:
    bool m_enabled;
    double m_time;
    bool m_inserted = false;
};
