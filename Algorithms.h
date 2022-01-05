#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

long double T = 1440;

/*
* Function that compute the arrival time at a successor node j on arc(i,j)
* for any departure time from the predecessor node i.
* This is a c++ implementation of Algorithm 1 discussed in page 2.
* @param time_intervals is 2d dynamic array with shape (τ,2) where τ is the number
* of time intervals, every row i contains two long double values (starting time, ending time)
* of time interval number i.
* @param speed is a 1d dynamic array with size τ, containing the speed of every time interval.
* @param distance_ij is the total distance of arc(i,j) between node i and node j.
* @param departure_time is the departure time from node i.
* @return arrival_time is the arrival time at node j.
*/

long double compute_arr_time_at_suc_node(vector<vector<long double>> time_intervals,
                                         vector<long double> speed, long double distance_ij,
                                         long double departure_time)
{
    // STEP 1: t ← t0
    long double current_time = departure_time;
    // STEP 2: d ← d_ij
    long double total_distance = distance_ij;
    // STEP 3: τ ← 0
    int time_interval_index = 0;
    // STEP 4: find the time interval
    // time_intervals[time_interval_index][1] is the ending time for current interval.
    while (current_time >= time_intervals[time_interval_index][1])
    {
        // STEP 5: τ ← τ + 1
        time_interval_index++;
        if (time_interval_index >= time_intervals.size())
            return INF;
    }
    // STEP 7: ∆d ← 0 (total distance traveled up to period τ)
    long double total_distance_upto_now = 0;
    // STEP 8: δdτ ← ( ̄tτ − t ) * vτ_ij (cumulative distance traveled over & before period τ)
    long double cumulative_distance = (time_intervals[time_interval_index][1] - current_time) * speed[time_interval_index];
    // STEP 9: while δdτ < d_ij do
    while (cumulative_distance < total_distance)
    {
        // STEP 10: ∆d ← δdτ
        total_distance_upto_now = cumulative_distance;
        // STEP 11: t ← ̄tτ
        current_time = time_intervals[time_interval_index][1];
        // STEP 12: τ ← τ + 1
        time_interval_index++;
        if (time_interval_index >= time_intervals.size())
            return INF;

        // STEP 13: δdτ ← ∆d + ( ̄tτ − t ) * vτ_ij
        cumulative_distance = total_distance_upto_now + (time_intervals[time_interval_index][1] - current_time) * speed[time_interval_index];
    }
    // STEP 15: tarr_j ← t + ( dij − ∆d ) / vτ_ij
    long double arrival_time = current_time + (total_distance - total_distance_upto_now) / speed[time_interval_index];

    return arrival_time;
}

/*
* Function that computes the departure time at the predecessor node i on arc(i,j)
* for any arrival time from the successor node j.
* This is a c++ implementation of Algorithm 2 discussed in page 3.
* @param time_intervals is 2d dynamic array with shape (τ,2) where τ is the number
* of time intervals, every row i contains two long double values (starting time, ending time)
* of time interval number i.
* @param speed is a 1d dynamic array with size τ, containing the speed profiles of a given path
* within those intervals spτ.
* @param distance_ij is the total distance of arc(i,j) between node i and node j.
* @param arrival_time is the arrival time from the successor node j.
* @return departure_time is the departure time at the the predecessor node i.
*/

long double compute_dep_time_at_pre_node(vector<vector<long double>> time_intervals,
                                         vector<long double> speed_profile, long double distance_ij,
                                         long double arrival_time)
{
    // STEP 1: t ← tdep_j
    long double current_time = arrival_time;
    // STEP 2: d ← d_ij
    long double total_distance = distance_ij;
    // STEP 3: τ ← 0
    int time_interval_index = 0;
    // STEP 4: find the time interval
    // time_intervals[time_interval_index][1] is the ending time for current interval.
    while (arrival_time >= time_intervals[time_interval_index][1])
    {
        // STEP 5: τ ← τ + 1
        time_interval_index++;
        if (time_interval_index >= time_intervals.size())
            return INF;
    }
    // STEP 7: tdep_i ← t − ( d / v{τ,spτ}_ij )
    long double departure_time = current_time - total_distance / speed_profile[time_interval_index];
    // step 8: while tdep_i < tτ do
    // time_intervals[time_interval_index][0] is the starting time for current interval.
    while (departure_time < time_intervals[time_interval_index][0])
    {
        // STEP 9: d ← d − v{τ,spτ}_ij * ( t − tτ )
        total_distance -= (current_time - time_intervals[time_interval_index][0]) * speed_profile[time_interval_index];
        // STEP 10: t ← tτ
        current_time = time_intervals[time_interval_index][0];
        // STEP 11: τ ← τ − 1
        time_interval_index--;
        if (time_interval_index < 0)
            return -INF;
        // STEP 12: tdep_i ← t − ( d / v{τ,spτ}_ij )
        departure_time = current_time - total_distance / speed_profile[time_interval_index];
    }
    // STEP 14: return tdep_i
    return departure_time;
}

/*
* Function that find the breaking points of the piecewise-constant (step-wise) speed function.
* this function iterate over the range [0,T[ as a discrete array of time points
* with a small step (0.01 minute) and for every time point calculate the slope of the previous segment
* and the slope of the next segment and if the slopes values are different then this time point is a breaking point.
* @param time_intervals is 2d dynamic array with shape (τ,2) where τ is the number
* of time intervals, every row i contains two long double values (starting time, ending time)
* of time interval number i.
* @param speed is a 1d dynamic array with size τ, containing the speed profiles of a given path
* within those intervals spτ.
* @param distance_ij is the total distance of arc(i,j) between node i and node j.
* @return break_points is a dynamic array that contains the desired breaking points.
*/

vector<long double> compute_break_points(vector<vector<long double>> time_intervals,
                                         vector<long double> speed, long double distance_ij)
{
    // initializing the dynamic array of the breaking points.
    vector<long double> break_points;
    // add the first time 0 as a breaking point.
    break_points.push_back(0);
    // the number of time intervals
    int num_of_time_intervals = time_intervals.size();
    // iterate over the time intervals except the first and the last one.
    for (int time_interval_index = 1; time_interval_index < num_of_time_intervals - 1; time_interval_index++)
    {
        long double time_point = time_intervals[time_interval_index][1];
        // add the end time of the current time interval to the breaking points
        break_points.push_back(time_point);
        // add the departure time to the breaking points
        // when the arrival time is the end time of the current time interval
        break_points.push_back(compute_dep_time_at_pre_node(time_intervals, speed, distance_ij, time_point));
    }
    break_points.push_back(T - 1);

    // sort the breaking points
    sort(break_points.begin(), break_points.end());
    // return the dynamic array of the breaking points.
    return break_points;
}

/*
* Function that find Closed-Form Representation of the Arrival Times tarr_j(tdep_i).
* This is a c++ implementation of Algorithm 3 discussed in page 9.
* @param time_intervals is 2d dynamic array with shape (τ,2) where τ is the number
* of time intervals, every row i contains two long double values (starting time, ending time)
* of time interval number i.
* @param speed is a 1d dynamic array with size τ, containing the speed profiles of a given path
* within those intervals spτ.
* @param distance_ij is the total distance of arc(i,j) between node i and node j.
* @return A_Pieces is a simple dynamic array that store the Line Segments
* which discribe the linear function pieces.
*/

vector<segment> closed_form_construction_of_arr_time(vector<vector<long double>> time_intervals,
                                                     vector<long double> speed, long double distance_ij)
{
    // find the t1,t2,...,thij  a dynamic array that store the breaking points
    // of the piecewise-constant (step-wise) speed function νij : [0,T] → R+
    vector<long double> break_points = compute_break_points(time_intervals, speed, distance_ij);
    // h_ij is the number of breaking points.
    int h_ij = break_points.size();
    // the initial value of k = max { x | tx ≤ tdep_i(T) }
    int k = 0;
    // calculate the departure time at node i when the arrival time at node j is T-1.
    long double dep_time_at_arr_time_T = compute_dep_time_at_pre_node(time_intervals, speed, distance_ij, T - 1);
    // STEP 1: k ← k = max { x | tx ≤ tdep_i(T) }
    for (int i = 0; i < h_ij; i++)
    {
        if (break_points[i] <= dep_time_at_arr_time_T)
        {
            k = max(i, k);
        }
    }
    // the initial value of l ← min { x | tx ≥ tarr_j(0) }
    int l = INF;
    // calculate the arrival time at node j when the departure time at node i is 0.
    long double arr_time_at_dep_time_0 = compute_arr_time_at_suc_node(time_intervals, speed, distance_ij, 0);
    // STEP 2: l ← min { x | tx ≥ tarr_j(0) }
    for (int i = 0; i < h_ij; i++)
    {
        if (break_points[i] >= arr_time_at_dep_time_0)
        {
            l = min(i, l);
        }
    }
    // STEP 3: ABP = ∅
    vector<point> ABP;
    // STEP 4: ABP ← ( 0 , tarr_j(0) )
    ABP.push_back(point(0, compute_arr_time_at_suc_node(time_intervals, speed, distance_ij, 0)));
    // STEP 5: for x ∈ {1,...,k} do
    for (int x = 1; x <= k; x++)
    {
        // STEP 6: ABP ← ( tx , tarr_j(tx) )
        ABP.push_back(point(break_points[x], compute_arr_time_at_suc_node(time_intervals, speed, distance_ij, break_points[x])));
    }
    // STEP 8: for x ∈ {l,...,h_ij−1} do
    for (int x = l; x < h_ij; x++)
    {
        // STEP 9: ABP ← ( tdep_i(tx) , tx )
        ABP.push_back(point(compute_dep_time_at_pre_node(time_intervals, speed, distance_ij, break_points[x]), break_points[x]));
    }
    // STEP 11: A_Pieces = ∅
    vector<segment> A_Pieces;
    // STEP 12: Sort−and−Remove−Duplicates ( ABP )
    sort(ABP.begin(), ABP.end());
    ABP.erase(unique(ABP.begin(), ABP.end(), [](point BP_l, point BP_r)
                     { return (fabs(BP_l.x - BP_r.x) < EPS and fabs(BP_l.y - BP_r.y) < EPS); }),
              ABP.end());
    // STEP 13: for x ∈ {1,...,SIZE(ABP)−1} do
    for (int x = 0; x < ABP.size() - 1; x++)
    {
        // STEP 14: A_Pieces ← ( ABP[x] , ABP[x + 1] )
        A_Pieces.push_back(segment(ABP[x], ABP[x + 1]));
    }
    // STEP 16: return A_Pieces
    return A_Pieces;
}

/*
* Function that take a lower envelope and array of y coordinates.
* and for each y_coordinate it calculate it's x coordinate on the lower envelope.
*/
vector<point> get_x_on_LE(const vector<segment> &LE, vector<ld> y_coor)
{
    vector<point> ret;
    sort(y_coor.begin(), y_coor.end());
    int LE_index = 0;
    int sz_y = y_coor.size();
    int sz_LE = LE.size();
    for (int i = 0; i < sz_y; i++)
    {
        while (LE_index < sz_LE)
        {
            if (y_coor[i] <= LE[LE_index].q.y)
                break;
            LE_index++;
        }
        if (LE_index >= sz_LE)
            break;
        ld x_coor = get_x_on_seg(LE[LE_index], y_coor[i]);
        ret.push_back(point(x_coor, y_coor[i]));
    }
    return ret;
}

/*
* Function that compute f o g (x) = f(g(x))
* where f and g are piecwise linear function
* f is the compute_arr_time_at_suc_node funvtion from Algorithm1
* and g is the lower envelope (Epsilone in Bellman_Ford algorithm for some node)
*/
vector<segment> composite_f_o_g(vector<vector<long double>> time_intervals,
                                vector<long double> speed, long double distance_ij,
                                vector<segment> A_Pieces, vector<segment> LE)
{
    vector<segment> ret;
    vector<point> br_points;
    vector<point> x_coor;
    for (auto seg : LE)
    {
        x_coor.push_back(seg.p);
    }
    x_coor.push_back(LE.back().q);
    vector<ld> x_f;
    for (auto seg : A_Pieces)
    {
        x_f.push_back(seg.p.y);
    }
    x_f.push_back(A_Pieces.back().q.y);
    /// x for f are y for g
    vector<point> new_x = get_x_on_LE(LE, x_f);
    for (auto new_pnt : new_x)
    {
        x_coor.push_back(new_pnt);
    }
    sort(x_coor.begin(), x_coor.end());
    for (auto pnt : x_coor)
    {
        if (pnt.y < 0 or pnt.x < 0)
            continue;
        if (pnt.y > T or pnt.x > T)
            continue;
        ld y_val = compute_arr_time_at_suc_node(time_intervals, speed, distance_ij, pnt.y);
        if (y_val > T or y_val < 0)
            continue;
        br_points.push_back(point(pnt.x, y_val));
    }
    br_points.push_back(point(compute_dep_time_at_pre_node(time_intervals, speed, distance_ij, T - 1), T - 1));
    int sz_br = br_points.size();
    for (int i = 1; i < sz_br; i++)
    {
        ret.push_back(segment(br_points[i - 1], br_points[i]));
    }
    ret = remove_dup(ret);
    return ret;
}

/*
* Function that checks it two Lower Envelopes are equal or not.
*/
bool equal_LEs(vector<segment> a, vector<segment> b)
{
    a = remove_dup(a);
    b = remove_dup(b);
    if (a.size() != b.size())
        return false;
    for (int i = 0; i < a.size(); i++)
    {
        if (!equal_seg(a[i], b[i]))
            return false;
    }
    return true;
}

/*
* Implementation of Bellman Ford Algorithm on page 12.
* @param n the number of nodes in the network
* @param start_node the index of the start node to run bellman ford from.
* @param end_node the index of the end node to calculate shortest path to.
* @param arcs dynamic array of size [n][n][number of arcs between node i and j].
* contains for every arc(i,j) the distance of that arc.
* @param time intervals contain [tao][2] double value ( the start time and end time of each interval).
* @param speed_profile contains the speed value for each interval for each speed profile.
* @ arc_speed dynamic array of size [n][n][number of arcs between node i and j].
* contains for every arc(i,j) the index of speed profile for that arc.
* see test code for more understanding for the paeameters.
*/

vector<segment> Bellman_Ford(int n, int start_node, int end_node, vector<vector<vector<long double>>> arcs,
                             vector<vector<long double>> time_intervals, vector<vector<long double>> speed_profile,
                             vector<vector<vector<int>>> arcs_speed)
{
    queue<int> FIFO;                                                                                 // L
    vector<vector<segment>> weight(n + 1);                                                           // Ψarr_j
    vector<vector<segment>> temp_weight(n + 1);                                                      // Ψarr′_j
    vector<vector<vector<vector<segment>>>> A_Pieces(n + 1, vector<vector<vector<segment>>>(n + 1)); // t_arr(tdep_i)
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            for (int arc_index = 0; arc_index < arcs[i][j].size(); arc_index++)
            {
                int speed_profile_index = arcs_speed[i][j][arc_index];
                long double arc_distance = arcs[i][j][arc_index];
                // compute the path structure for every arc(i,j)
                A_Pieces[i][j].push_back(closed_form_construction_of_arr_time(time_intervals, speed_profile[speed_profile_index],
                                                                              arc_distance));
            }
        }
    }

    // for j∈V do
    // Ψarr′_j ← Ψarr_j ← { fID if i = j , f∞ if otherwise
    for (int i = 1; i <= n; i++)
    {
        if (i == start_node)
        {
            weight[i].push_back(segment(point(0, 0), point(INF, INF)));
            temp_weight[i].push_back(segment(point(0, 0), point(INF, INF)));
        }
        else
        {
            weight[i].push_back(segment(point(0, INF), point(INF, INF)));
            temp_weight[i].push_back(segment(point(0, INF), point(INF, INF)));
        }
    }
    // L ← { i }
    FIFO.push(start_node);
    // while L != ∅ do
    while (!FIFO.empty())
    {
        // for (x,y) ∈ A : x∈L do
        while (!FIFO.empty())
        {
            int x = FIFO.front();

            FIFO.pop();
            // Ψarr′_y(tdep_i) ← LowerEnvelope(Ψarr′_y(tdep_i),(tarr_y(tdep_i) ◦ Ψarr_x(tdep_i)))
            for (int y = 1; y <= n; y++)
            {
                for (int arc_index = 0; arc_index < arcs[x][y].size(); arc_index++)
                {
                    // comosite_LE is tarr_y(tdep_i)◦Ψarr_x(tdep_i)
                    vector<segment> comosite_LE = composite_f_o_g(
                        time_intervals, speed_profile[arcs_speed[x][y][arc_index]],
                        arcs[x][y][arc_index], A_Pieces[x][y][arc_index], weight[x]);
                    for (auto seg : temp_weight[y])
                        comosite_LE.push_back(seg);
                    // new_wight is LowerEnvelope(Ψarr′_y(tdep_i),(tarr_y(tdep_i) ◦ Ψarr_x(tdep_i)))
                    vector<segment> new_wight = LowerEnvelope(comosite_LE);
                    // Ψarr′_y(tdep_i) ← LowerEnvelope(Ψarr′_y(tdep_i),(tarr_y(tdep_i) ◦ Ψarr_x(tdep_i)))
                    temp_weight[y] = new_wight;
                }
            }
        }
        // for(y∈V) do
        for (int i = 1; i <= n; i++)
        {
            if (i == start_node)
                continue;
            // if Psiarr_y(tdep_i) != Psiarr′_y(tdep_i) then
            if (!equal_LEs(weight[i], temp_weight[i]))
            {
                // L ← L ∪ y
                FIFO.push(i);
            }
            // Psiarr_y(tdep_i) ← Psiarr′_y(tdep_i)
            weight[i] = temp_weight[i];
        }
    }
    return weight[end_node];
}

#endif // ALGORITHMS_H
