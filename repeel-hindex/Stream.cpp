
#include "common.h"
#include "Stream.h"

/**
    find edges in [ts, te]
**/
void Stream::getEdges(uint32_t uiTS, uint32_t uiTE, vector<pair<uint32_t, uint32_t> > &vDesEdges)
{
    /*map<int, vector<pair<uint32_t, uint32_t> > >::iterator itLowB;
    map<int, vector<pair<uint32_t, uint32_t> > >::iterator itUpB;

    itLowB = m_mpE.lower_bound(uiTS);
    itUpB = m_mpE.upper_bound(uiTE);
    while (itLowB != itUpB)
    {
        vDesEdges.insert(vDesEdges.end(), itLowB->second.begin(), itLowB->second.end());
        ++itLowB;
    }*/
    vector<ArrayEntry>::iterator itLowB;
    vector<ArrayEntry>::iterator itUpB;

    itLowB = upper_bound(m_vE.begin(), m_vE.end(), (ArrayEntry){uiTS, 0, 0}, cmp);
    itUpB = upper_bound(m_vE.begin(), m_vE.end(), (ArrayEntry){uiTE + 1, 0, 0}, cmp);

    while (itLowB != itUpB)
    {
        vDesEdges.push_back({itLowB->x, itLowB->y});
        ++itLowB;
    }
}

/**
    read data
**/
void Stream::readFromFile(char* pcFile)
{
    FILE* fp = fopen (pcFile, "rt");
    char acBuffer[100];
    ASSERT_MSG(NULL != fp, "invalid graph file");
    uint32_t x = 0;
    uint32_t y = 0;
    uint32_t t = 0;
    m_uiN = 0;
    m_uiM = 0;
    bool bFirst = true;
    while (!feof(fp))
    {
        char *pPos = fgets(acBuffer, 100, fp);
        if (NULL == pPos)
        {
            break;
        }

        int res = sscanf(acBuffer, "%d%d%d", &x, &y, &t)/*fscanf(acBuffer, "%d %d %d", &x, &y, &t) */;

        ASSERT_MSG(3 == res, "wrong file");


        //m_mpE[t].push_back({x, y});
        if(x != y)
        {
            m_vE.push_back({t, x, y});

            ++m_uiM;
            m_uiN = (m_uiN > x) ? m_uiN : x;
            m_uiN = (m_uiN > y) ? m_uiN : y;
            if (bFirst)
            {
                bFirst = false;
                m_uiTS = t;
            }
            else
            {
                m_uiTS = (m_uiTS < t) ? m_uiTS : t;
            }

            m_uiTE = (m_uiTE > t) ? m_uiTE : t;

            ASSERT_MSG(x != y, "loops exist in the graph" << x << y);
        }

    }

    fclose(fp);

    std::sort(m_vE.begin(), m_vE.end(), cmp);
}

/**
    save data
**/
void Stream::writeToFile(char* pcFile, vector<pair<uint32_t, uint32_t> > &vEdges) const
{
    FILE* fp = fopen (pcFile, "w+");
    ASSERT_MSG(NULL != fp, "invalid output file");
    for (const auto atE : vEdges)
    {
        fprintf(fp, "%d>%d\n", atE.first, atE.second);
    }

   fclose(fp);
}
