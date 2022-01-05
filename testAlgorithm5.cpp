#include <iostream>

#include "Lower_Envelope.h"
#include "Algorithms.h"

using namespace std;

int main()
{
    int n = 5;
    int start_node = 1;
    int end_node = n;
    vector<vector<vector<long double>>> arcs_length(n + 1, vector<vector<long double>>(n + 1));
    vector<vector<vector<int>>> arcs_speed(n + 1, vector<vector<int>>(n + 1));

    arcs_length[1][2].push_back(7);
    arcs_speed[1][2].push_back(0);

    arcs_length[1][2].push_back(5);
    arcs_speed[1][2].push_back(0);

    arcs_length[1][2].push_back(6.2);
    arcs_speed[1][2].push_back(0);

    arcs_length[1][2].push_back(7);
    arcs_speed[1][2].push_back(1);

    arcs_length[1][3].push_back(2.5);
    arcs_speed[1][3].push_back(0);

    arcs_length[1][3].push_back(4);
    arcs_speed[1][3].push_back(1);

    arcs_length[1][3].push_back(5);
    arcs_speed[1][3].push_back(1);

    arcs_length[2][4].push_back(15);
    arcs_speed[2][4].push_back(1);

    arcs_length[2][4].push_back(12.8);
    arcs_speed[2][4].push_back(0);

    arcs_length[3][4].push_back(16);
    arcs_speed[3][4].push_back(0);

    arcs_length[3][4].push_back(24);
    arcs_speed[3][4].push_back(1);

    arcs_length[4][5].push_back(5);
    arcs_speed[4][5].push_back(0);

    arcs_length[4][5].push_back(3);
    arcs_speed[4][5].push_back(1);

    vector<vector<long double>> time_intervals;

    time_intervals.push_back({0, 20});
    time_intervals.push_back({20, 100});
    time_intervals.push_back({100, 150});
    time_intervals.push_back({150, 1440});

    vector<vector<long double>> speed_profiles(2);

    speed_profiles[0].push_back(10 / 60.0);
    speed_profiles[0].push_back(20 / 60.0);
    speed_profiles[0].push_back(10 / 60.0);
    speed_profiles[0].push_back(20 / 60.0);

    speed_profiles[1].push_back(17 / 60.0);
    speed_profiles[1].push_back(25 / 60.0);
    speed_profiles[1].push_back(8 / 60.0);
    speed_profiles[1].push_back(22 / 60.0);

    vector<segment> res = Bellman_Ford(n, start_node, end_node, arcs_length, time_intervals, speed_profiles, arcs_speed);
    for (auto seg : res)
    {
        cout << "Segment " << seg.p.x << " " << seg.p.y << " " << seg.q.x << " " << seg.q.y << endl;
    }
    return 0;
}
