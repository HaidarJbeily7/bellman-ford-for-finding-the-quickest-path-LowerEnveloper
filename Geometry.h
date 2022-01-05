#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>

using namespace std;

typedef long double ld;
const ld EPS = 1e-9;
const ld INF = 1e18;

/*
* Struct of point.
* has two members x,y the coordinates of the point.
* less than operator has been overloaded, first compare based on x values
* then if eqauls compare based on y values.
*/
struct point
{
    ld x, y;
    point() {}
    point(ld _x, ld _y)
    {
        x = _x;
        y = _y;
    }
    bool operator<(const point &p) const
    {
        return ((x < p.x - EPS) or (fabs(x - p.x) < EPS and y < p.y - EPS));
    }
};

/*
* Struct of Line.
* has three members (a,b,c) such that the equation of the line is a*x+b*y+c = 0
* has a construct that take two points and caluclate a,b,c.
* norm() is a method that normailze a,b,c.
* distance(point) takes a point and caluclate the distance between he Line and the given point.
*/

struct line
{
    ld a, b, c;
    line() {}
    line(point p, point q)
    {
        a = p.y - q.y;
        b = q.x - p.x;
        c = -a * p.x - b * p.y;
        norm();
    }
    void norm()
    {
        ld z = sqrt(a * a + b * b);
        if (abs(z) > EPS)
        {
            a /= z, b /= z, c /= z;
        }
    }
    ld distance(point p)
    {
        return a * p.x + b * p.y + c;
    }
};

/*
* det(a,b,c,d) caluclate the determinant of the matrix [[a,b],[c,d]]
*/
ld det(ld a, ld b, ld c, ld d)
{
    return a * d - b * c;
}

/*
* between(l,r,x) return true if a given x is between [l,r] , or false otherwise.
*/
inline bool between(ld l, ld r, ld x)
{
    return min(l, r) <= x + EPS and x <= max(l, r) + EPS;
}

/*
* intersect_1d(a,b,c,d) return true if two segments [a,b] ,[c,d] intersects in 1d
*/
inline bool intersect_1d(ld a, ld b, ld c, ld d)
{
    if (a > b)
        swap(a, b);
    if (c > d)
        swap(c, d);
    return max(a, c) <= min(b, d) + EPS;
}

/*
* Struct of segment.
* has two members p,q the two endpoints of the segment. 
*/
struct segment
{
    point p, q;
    segment(point a, point b)
    {
        p = point(a.x, a.y);
        q = point(b.x, b.y);
        if (q < p)
            swap(p, q);
    }
};

/*
* intersect_segment return true if two given segment intersects, 
* and store the left intersection point in the parameter left_int.
*/
bool intersect_segment(segment first_seg, segment second_seg, point &left_int)
{
    point a = first_seg.p;
    point b = first_seg.q;
    point c = second_seg.p;
    point d = second_seg.q;
    if (!intersect_1d(a.x, b.x, c.x, d.x) or
        !intersect_1d(a.y, b.y, c.y, d.y))
        return false;

    line m(a, b);
    line n(c, d);
    ld zn = det(m.a, m.b, n.a, n.b);
    if (abs(zn) < EPS)
    {
        if (abs(m.distance(c)) > EPS or abs(n.distance(a)) > EPS)
            return false;
        if (b < a)
            swap(a, b);
        if (d < c)
            swap(c, d);
        left_int = max(a, c);
        return true;
    }
    else
    {
        left_int.x = -det(m.c, m.b, n.c, n.b) / zn;
        left_int.y = -det(m.a, m.c, n.a, n.c) / zn;
        return between(a.x, b.x, left_int.x) and between(a.y, b.y, left_int.y) and between(c.x, d.x, left_int.x) and between(c.y, d.y, left_int.y);
    }
}

/*
* get_x_next return the x coordinate of the intersection of two segments
* to the right of the x coordinate of the sweep_line or return the leftmost
* x coordinate of the endpoints of the two segments. 
*/
ld get_x_next(segment left, segment right, ld sweep_line)
{
    point int_point = point(INF, INF);
    bool temp = intersect_segment(left, right, int_point);
    if (!temp)
        int_point = point(INF, INF);
    vector<ld> x_coor;
    x_coor.push_back(left.p.x);
    x_coor.push_back(left.q.x);
    x_coor.push_back(right.p.x);
    x_coor.push_back(right.q.x);
    x_coor.push_back(int_point.x);
    sort(x_coor.begin(), x_coor.end());
    ld ret = INF;
    for (ld x : x_coor)
    {
        if (x > sweep_line + EPS)
        {
            ret = x;
            break;
        }
    }
    return ret;
}

/*
* get_y_on_seg return the y coordinate of a given x coordinate 
* of a point on the segment 
* (the point of the intersection of the line x = x_coor with the given segment)
*/
ld get_y_on_seg(segment seg, ld x_coor)
{
    if (abs(seg.p.y - seg.q.y) < EPS)
        return seg.p.y;
    else
    {
        ld y_coor = seg.p.y + ((x_coor - seg.p.x) * (seg.p.y - seg.q.y)) / (seg.p.x - seg.q.x);
        return y_coor;
    }
}

/*
* get_x_on_seg return the x coordinate of a given y coordinate 
* of a point on the segment 
* (the point of the intersection of the line y = y_coor with the given segment)
*/
ld get_x_on_seg(segment seg, ld y_coor)
{
    if (abs(seg.p.x - seg.q.x) < EPS)
        return seg.p.x;
    else
    {
        ld x_coor = seg.p.x + ((y_coor - seg.p.y) * (seg.p.x - seg.q.x)) / (seg.p.y - seg.q.y);
        return x_coor;
    }
}

/*
* get_y_event return the y coordinate of the next event of the sweep line.
*/
ld get_y_event(segment left, segment right, ld x_coor)
{
    vector<ld> y_candidate;
    if (between(left.p.x, left.q.x, x_coor))
    {
        y_candidate.push_back(get_y_on_seg(left, x_coor));
    }
    if (between(right.p.x, right.q.x, x_coor))
    {
        y_candidate.push_back(get_y_on_seg(right, x_coor));
    }
    sort(y_candidate.begin(), y_candidate.end());
    if (y_candidate.size() == 0)
        return -INF;
    return y_candidate[0];
}

/*
* point_equal return true of two given points are eqauls, or false otherwise.
*/
bool point_equal(point a, point b)
{
    return abs(a.x - b.x) <= EPS and abs(a.y - b.y) <= EPS;
}

/*
* equal_seg return true of two given segments are eqauls, or false otherwise.
*/
bool equal_seg(segment a, segment b)
{
    return point_equal(a.p, b.p) and point_equal(a.q, b.q);
}

/*
* seg_len return the length of a given segment.
*/
ld seg_len(segment a)
{
    return sqrt((a.p.x - a.q.x) * (a.p.x - a.q.x) + (a.p.y - a.q.y) * (a.p.y - a.q.y));
}

/*
* can_merge_seg rturn true if two given segment could be merged in one segment, false oterwise.
*/
bool can_merge_seg(segment a, segment b)
{
    if (!point_equal(a.q, b.p))
        return false;
    segment c = segment(a.p, b.q);
    return abs(seg_len(a) + seg_len(b) - seg_len(c)) <= EPS;
}

/*
* merge_segment return the segment of merging two segments.
*/
segment merge_segment(segment a, segment b)
{
    return segment(a.p, b.q);
}

#endif // GEOMETRY_H
