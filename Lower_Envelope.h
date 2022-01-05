#ifndef LOWER_ENVELOPE_H
#define LOWER_ENVELOPE_H

#include "Geometry.h"

/*
* Function that take a vector of segment (Lower Envelope) and merging 
* every two consecutive segments that can be merged into one segment. 
*/
vector<segment> remove_dup(const vector<segment> &segs)
{
    vector<segment> ret;
    for (auto seg : segs)
    {
        if (ret.size() == 0)
        {
            ret.push_back(seg);
        }
        else
        {
            segment lst = ret.back();
            ret.pop_back();
            if (can_merge_seg(lst, seg))
            {
                ret.push_back(merge_segment(lst, seg));
            }
            else
            {
                ret.push_back(lst);
                ret.push_back(seg);
            }
        }
    }
    return ret;
}
/*
* Implementation of Merging Two Lower Envelope Using Line Sweep Algorithm mentioned in Atallah paper.
* @param left_LE dynamic array of segment contains the returned Lower Envelope from the Left child of Divide and Conquer Algorithm.
* @param right_LE dynamic array of segment contains the returned Lower Envelope from the Right child of Divide and Conquer Algorithm.
* @return merged_LE dynamic array contains the Merged Lower Envelope.
*/
vector<segment> merge_LE(const vector<segment> &left_LE, const vector<segment> &right_LE)
{
    int left_index = 0;
    int right_index = 0;
    int sz_left_LE = left_LE.size();
    int sz_right_LE = right_LE.size();
    if (sz_left_LE == 0)
    {
        return right_LE;
    }
    if (sz_right_LE == 0)
    {
        return left_LE;
    }
    vector<segment> merged_LE;
    ld sweep_line_x = -INF;
    ld sweep_line_y = -INF;
    point temp_point = min(left_LE[left_index].p, right_LE[right_index].p);
    sweep_line_x = temp_point.x;
    sweep_line_y = temp_point.y;
    while (left_index < sz_left_LE and right_index < sz_right_LE)
    {
        ld x_next_event = get_x_next(left_LE[left_index], right_LE[right_index], sweep_line_x);
        ld y_next_event = get_y_event(left_LE[left_index], right_LE[right_index], x_next_event - 2 * EPS);
        merged_LE.push_back(segment(point(sweep_line_x, sweep_line_y), point(x_next_event, y_next_event)));
        sweep_line_x = x_next_event;

        bool change_segment = false;
        if (abs(sweep_line_x - left_LE[left_index].q.x) < EPS)
        {
            change_segment = true;
            left_index++;
        }
        if (abs(sweep_line_x - right_LE[right_index].q.x) < EPS)
        {
            change_segment = true;
            right_index++;
        }
        if (right_index >= sz_right_LE)
            break;
        if (left_index >= sz_left_LE)
            break;
        if (change_segment)
        {
            point temp_point = min(left_LE[left_index].p, right_LE[right_index].p);
            sweep_line_x = max(temp_point.x, x_next_event);
            sweep_line_y = get_y_event(left_LE[left_index], right_LE[right_index], sweep_line_x + 2 * EPS);
        }
        else
        {
            sweep_line_y = get_y_event(left_LE[left_index], right_LE[right_index], sweep_line_x + 2 * EPS);
        }
    }
    if (left_index < sz_left_LE and sweep_line_x >= left_LE[left_index].p.x)
    {
        sweep_line_y = get_y_on_seg(left_LE[left_index], sweep_line_x);
        merged_LE.push_back(segment(point(sweep_line_x, sweep_line_y), left_LE[left_index].q));
        left_index++;
    }
    while (left_index < sz_left_LE)
    {

        merged_LE.push_back(left_LE[left_index]);
        left_index++;
    }
    if (right_index < sz_right_LE and sweep_line_x >= right_LE[right_index].p.x)
    {
        sweep_line_y = get_y_on_seg(right_LE[right_index], sweep_line_x);
        merged_LE.push_back(segment(point(sweep_line_x, sweep_line_y), right_LE[right_index].q));
        right_index++;
    }

    while (right_index < sz_right_LE)
    {
        merged_LE.push_back(right_LE[right_index]);
        right_index++;
    }
    merged_LE = remove_dup(merged_LE);
    return merged_LE;
}

/*
* Implementation of Divide and Conquer Algorithm mentioned in Atallah paper.
* this function is used on the second phase of the Algorithm.
* @param segments dynamic array contains the segments to get the Lower Envelope for.
* @return LE dynamic array contains the the segments of the Lower Envelope.
*/
vector<segment> get_LE(const vector<segment> &segments)
{
    int sz_segs = segments.size();
    if (sz_segs < 2)
    {
        return segments;
    }
    int mid = sz_segs / 2;
    vector<segment> left_segments;
    vector<segment> right_segments;
    for (int i = 0; i < mid; i++)
    {
        left_segments.push_back(segments[i]);
    }
    for (int i = mid; i < sz_segs; i++)
    {
        right_segments.push_back(segments[i]);
    }
    // Divide
    vector<segment> left_LE = get_LE(left_segments);
    vector<segment> right_LE = get_LE(right_segments);
    // Conquer
    vector<segment> LE = merge_LE(left_LE, right_LE);

    return LE;
}
/*
* The first phase concludes by producing o(log n) subsets of S
* or each level of the binary tree; it combines all the sets S[i] for nodes i at that level
* @param A_Pieces the segments to patition
* @return dynamic array contains log n dynamic array each one represents a disjoint set.
*/
vector<vector<segment>> partition_phase(const vector<segment> &A_Pieces)
{
    // sorting 2n endpoint of the segments based on x coordinates.
    vector<pair<ld, int>> x_sorted_point;
    int n = A_Pieces.size();
    int segment_index = 1;
    for (auto A_Piece : A_Pieces)
    {
        point start_point = A_Piece.p;
        point end_point = A_Piece.q;
        ld x_coor = start_point.x;
        x_sorted_point.push_back({x_coor, segment_index});
        x_coor = end_point.x;
        x_sorted_point.push_back({x_coor, -segment_index});
        segment_index++;
    }
    sort(x_sorted_point.begin(), x_sorted_point.end());

    // generating the domain of each segment by an integer interva1 (i, j)
    //  such that (x_i, x_j) = (left(s), right(s)).
    vector<pair<int, int>> domain_segements;
    for (int i = 0; i < n; i++)
        domain_segements.push_back({-1, -1});
    int index = 1;
    for (pair<ld, int> temp : x_sorted_point)
    {
        ld x_coor = temp.first;
        segment_index = temp.second;
        if (segment_index < 0)
        {
            segment_index = -segment_index;
            domain_segements[segment_index - 1].second = index;
        }
        else
        {
            domain_segements[segment_index - 1].first = index;
        }
        index += 2;
    }

    // building a complete binary tree with 2n leaves
    // and label the nodes in symmetric order from 1 to 4n - 1
    queue<pair<int, int>> q;
    int mx = 0;
    for (int i = 1;; i += 2)
    {
        mx = i;
        if ((i & (i + 1)) == 0 and i >= 4 * n - 1)
        {
            break;
        }
    }
    vector<int> lvl(mx + 1);
    for (int i = 1;; i += 2)
    {
        q.push({i, 1});
        lvl[i] = 1;
        if ((i & (i + 1)) == 0 and i >= 4 * n - 1)
        {
            break;
        }
    }
    vector<int> tree_parents(mx + 1);
    int cur_lvl = 0;
    while (true)
    {
        if (q.size() < 2)
            break;
        int lft = q.front().first;
        cur_lvl = q.front().second + 1;
        q.pop();
        int rgt = q.front().first;
        q.pop();
        int mid = (lft + rgt) / 2;
        lvl[mid] = cur_lvl;
        tree_parents[lft] = mid;
        tree_parents[rgt] = mid;
        q.push({mid, cur_lvl});
    }
    // dynamic array contains O(log n) disjoint sets.
    vector<vector<segment>> paritioned_sets;
    paritioned_sets.resize(cur_lvl + 1);
    segment_index = 0;
    for (pair<int, int> domain_segment : domain_segements)
    {
        int a = domain_segment.first;
        int b = domain_segment.second;
        // getting the lca of the endpoint of the current domain segment.
        if (lvl[a] >= lvl[b])
            swap(a, b);
        while (lvl[a] < lvl[b])
            a = tree_parents[a];
        while (a != b)
        {
            a = tree_parents[a];
            b = tree_parents[b];
        }
        int lca = a;
        // adding the domain segment to the disjoint set at the index
        // of the level of the lca of it's integer endpoints.
        paritioned_sets[lvl[lca]].push_back(A_Pieces[segment_index]);
        segment_index++;
    }
    return paritioned_sets;
}

/*
* The third phase of the Algorithm.
* find the lower envelope of log n lower envelope using divide and conquer.
* @param LEs dynamic array contains log n lower envelopes.
* @return LE dynamic array contains the lower envelope of the given log n lower envelopes.
*/
vector<segment> third_phase(const vector<vector<segment>> &LEs)
{
    int sz = LEs.size();
    if (sz == 0)
    {
        vector<segment> ret;
        return ret;
    }
    if (sz == 1)
    {
        return LEs[0];
    }
    int mid = sz / 2;
    vector<vector<segment>> left_LEs;
    vector<vector<segment>> right_LEs;
    for (int i = 0; i < mid; i++)
    {
        left_LEs.push_back(LEs[i]);
    }
    for (int i = mid; i < sz; i++)
    {
        right_LEs.push_back(LEs[i]);
    }
    // Divide
    vector<segment> left_LE = third_phase(left_LEs);
    vector<segment> right_LE = third_phase(right_LEs);
    // Conquer
    vector<segment> LE = merge_LE(left_LE, right_LE);
    return LE;
}
/*
* Implemetaion of the O(n log n) Algorithm for finding the lower envelope
* of a given segments.
* @param segments dynamic array contains the n segments.
* @return LE dynamic array contains the lower envelope.
*/
vector<segment> LowerEnvelope(const vector<segment> &segments)
{

    vector<vector<segment>> disjoint_sets;

    vector<vector<segment>> LEs;

    vector<segment> LE;

    // the first phase of the Algorithm ( the parition phase).
    disjoint_sets = partition_phase(segments);

    // the second phase of the algorithm (finding the lower envelope of the log n disjoint sets).
    for (auto disjoint_set : disjoint_sets)
    {
        vector<segment> cur_LE = get_LE(disjoint_set);
        LEs.push_back(cur_LE);
    }

    // the third phase (finding the lower envelope of the log n lower envelopes).
    LE = third_phase(LEs);

    return LE;
}

#endif // LOWER_ENVELOPE_H
