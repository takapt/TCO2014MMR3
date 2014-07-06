#define NDEBUG

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <map>
#include <utility>
#include <set>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>
#include <complex>
#include <stack>
#include <queue>
#include <numeric>
#include <list>
#include <iomanip>
#include <fstream>
#include <bitset>

using namespace std;

#define foreach(it, c) for (__typeof__((c).begin()) it=(c).begin(); it != (c).end(); ++it)
template <typename T> void print_container(ostream& os, const T& c) { const char* _s = " "; if (!c.empty()) { __typeof__(c.begin()) last = --c.end(); foreach (it, c) { os << *it; if (it != last) os << _s; } } }
template <typename T> ostream& operator<<(ostream& os, const vector<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const set<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const multiset<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const deque<T>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const map<T, U>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const pair<T, U>& p) { os << "(" << p.first << ", " << p.second << ")"; return os; }

template <typename T> void print(T a, int n, const string& split = " ") { for (int i = 0; i < n; i++) { cerr << a[i]; if (i + 1 != n) cerr << split; } cerr << endl; }
template <typename T> void print2d(T a, int w, int h, int width = -1, int br = 0) { for (int i = 0; i < h; ++i) { for (int j = 0; j < w; ++j) { if (width != -1) cerr.width(width); cerr << a[i][j] << ' '; } cerr << endl; } while (br--) cerr << endl; }
template <typename T> void input(T& a, int n) { for (int i = 0; i < n; ++i) cin >> a[i]; }
#define dump(v) (cerr << #v << ": " << v << endl)

#define rep(i, n) for (int i = 0; i < (int)(n); ++i)
#define all(a) (a).begin(), (a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define clr(a, x) memset(a, x, sizeof(a))
#define sz(a) ((int)(a).size())
#define mp(a, b) make_pair(a, b)
#define ten(n) ((long long)(1e##n))

template <typename T, typename U> void upmin(T& a, const U& b) { a = min<T>(a, b); }
template <typename T, typename U> void upmax(T& a, const U& b) { a = max<T>(a, b); }
template <typename T> void uniq(T& a) { sort(a.begin(), a.end()); a.erase(unique(a.begin(), a.end()), a.end()); }
template <class T> string to_s(const T& a) { ostringstream os; os << a; return os.str(); }
template <class T> T to_T(const string& s) { istringstream is(s); T res; is >> res; return res; }
void fast_io() { cin.tie(0); ios::sync_with_stdio(false); }
bool in_rect(int x, int y, int w, int h) { return 0 <= x && x < w && 0 <= y && y < h; }


typedef pair<int, int> pint;
typedef long long ll;

const int DX[] = { 0, 1, 0, -1 };
const int DY[] = { 1, 0, -1, 0 };


#ifdef _MSC_VER
#include <Windows.h>
#else
#include <sys/time.h>
#endif
class Timer
{
    typedef double time_type;
    typedef unsigned int skip_type;

private:
    time_type start_time;
    time_type elapsed;

#ifdef _MSC_VER
    time_type get_ms() { return (time_type)GetTickCount64() / 1000; }
#else
    time_type get_ms() { struct timeval t; gettimeofday(&t, NULL); return (time_type)t.tv_sec * 1000 + (time_type)t.tv_usec / 1000; }
#endif

public:
    Timer() {}

    void start() { start_time = get_ms(); }
    time_type get_elapsed() { return elapsed = get_ms() - start_time; }
};

class Random
{
private:
    unsigned int  x, y, z, w;
public:
    Random(unsigned int x
             , unsigned int y
             , unsigned int z
             , unsigned int w)
        : x(x), y(y), z(z), w(w) { }
    Random() 
        : x(123456789), y(362436069), z(521288629), w(88675123) { }
    Random(unsigned int seed)
        : x(123456789), y(362436069), z(521288629), w(seed) { }

    unsigned int next()
    {
        unsigned int t = x ^ (x << 11);
        x = y;
        y = z;
        z = w;
        return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }

    int next_int() { return next(); }

    // [0, upper)
    int next_int(int upper) { return next() % upper; }

    // [low, high]
    int next_int(int low, int high) { return next_int(high - low + 1) + low; }

    double next_double(double upper) { return upper * next() / UINT_MAX; }
    double next_double(double low, double high) { return next_double(high - low) + low; }

    template <typename T>
    int select(const vector<T>& ratio)
    {
        T sum = accumulate(ratio.begin(), ratio.end(), (T)0);
        T v = next_double(sum) + (T)1e-6;
        for (int i = 0; i < (int)ratio.size(); ++i)
        {
            v -= ratio[i];
            if (v <= 0)
                return i;
        }
        return 0;
    }
};

template <typename FLOW, typename COST>
class PrimalDual
{
private:
    struct PrimalDualEdge
    {
        int to;
        FLOW cap;
        COST cost;
        int rev;
        PrimalDualEdge(int to, FLOW cap, COST cost, int rev)
            : to(to), cap(cap), cost(cost), rev(rev) { }
    };
public:
    int V;
    vector<vector<PrimalDualEdge> > g;

    PrimalDual(int V) : V(V), g(V) { }

    void add_edge(int from, int to, FLOW cap, COST cost)
    {
        g[from].push_back(PrimalDualEdge(to, cap, cost, g[to].size()));
        g[to].push_back(PrimalDualEdge(from, 0, -cost, g[from].size() - 1));
    }
    void add_undirected(int a, int b, FLOW cap, COST cost)
    {
        add_edge(a, b, cap, cost);
        add_edge(b, a, cap, cost);
    }

    COST min_cost_flow(int s, int t, FLOW f)
    {
        vector<COST> h(V);
        COST res = 0;
        int _f = f;
        while (_f > 0)
        {
            typedef pair<COST, int> _p;
            const COST _INF = (COST)((1LL << 60) | (1 << 29));
            priority_queue<_p, vector<_p>, greater<_p> > q;
            vector<COST> dis(V, _INF);
            vector<int> prevv(V), preve(V);
            dis[s] = 0;
            q.push(_p(0, s));
            while (!q.empty())
            {
                _p p = q.top(); q.pop();
                int v = p.second;
                COST cost = p.first;
                if (cost > dis[v])
                    continue;

                for (int i = 0; i < g[v].size(); ++i)
                {
                    PrimalDualEdge& e = g[v][i];
                    const COST _eps = 1e-10;
                    COST c = cost + e.cost + h[v] - h[e.to];
                    if (e.cap > 0 && c + _eps < dis[e.to])
                    {
                        dis[e.to] = c;
                        prevv[e.to] = v;
                        preve[e.to] = i;
                        q.push(_p(c, e.to));
                    }
                }
            }

            if (dis[t] == _INF)
            {
                // cant flow _f
                return -_INF;
            }

            for (int i = 0; i < V; ++i)
                h[i] += dis[i];

            FLOW d = _f;
            for (int i = t; i != s; i = prevv[i])
                d = min(d, g[prevv[i]][preve[i]].cap);
            _f -= d;
            res += d * h[t];
            for (int i = t; i != s; i = prevv[i])
            {
                PrimalDualEdge& e = g[prevv[i]][preve[i]];
                e.cap -= d;
                g[e.to][e.rev].cap += d;
            }
        }

        return res;
    }
};

enum Dir
{
    UP,
    RIGHT,
    DOWN,
    LEFT,

    NA,
};
const string dir_s[] = { "UP", "RIGHT", "DOWN", "LEFT", "NA" };
Dir to_dir(const string& s)
{
    int i = find(dir_s, dir_s + 4, s) - dir_s;
    assert(0 <= i && i < 4);
    return Dir(i);
}
Dir rev_dir(Dir dir)
{
    return Dir((dir + 2) % 4);
}
struct Pos
{
    int x, y;
    Pos(int x, int y)
        : x(x), y(y)
    {
    }
    Pos()
        : x(0), y(0)
    {
    }

    bool operator==(const Pos& other) const
    {
        return x == other.x && y == other.y;
    }
    bool operator !=(const Pos& other) const
    {
        return x != other.x || y != other.y;
    }

    void operator+=(const Pos& other)
    {
        x += other.x;
        y += other.y;
    }
    void operator-=(const Pos& other)
    {
        x -= other.x;
        y -= other.y;
    }

    Pos operator+(const Pos& other) const
    {
        Pos res = *this;
        res += other;
        return res;
    }
    Pos operator-(const Pos& other) const
    {
        Pos res = *this;
        res -= other;
        return res;
    }
    Pos operator*(int a) const
    {
        return Pos(x * a, y * a);
    }

    bool operator<(const Pos& other) const
    {
        if (x != other.x)
            return x < other.x;
        else
            return y < other.y;
    }
};
Pos operator*(int a, const Pos& pos)
{
    return pos * a;
}
Pos to_pos(Dir dir)
{
    assert(0 <= dir && dir < 4);
    return Pos(DX[dir], DY[dir]);
}



const int SOURCE_IMAGES = 200;


class Rect
{
public:
    Rect(Pos pos, int w, int h)
        : pos_(pos), w_(w), h_(h)
    {
    }
    Rect(int x, int y, int w, int h)
        : pos_(Pos(x, y)), w_(w), h_(h)
    {
    }
    Rect()
        : pos_(Pos(1919810, 114514)), w_(364364)
    {
    }

    Pos& pos() { return pos_; }
    const Pos& pos() const { return pos_; }
    int width() const { return w_; }
    int height() const { return h_; }

    int left() const { return pos_.x; }
    int right() const { return pos_.x + w_; }
    int top() const { return pos_.y; }
    int bottom() const { return pos_.y + height(); }

    void expand_right(int a)
    {
        w_ += a;
    }
    void expand_left(int a)
    {
        pos_.x -= a;
        w_ += a;
    }
    void expand_bottom(int a)
    {
        h_ += a;
    }
    void expand_top(int a)
    {
        pos_.y -= a;
        h_ += a;
    }
    void expand(Dir dir, int a)
    {
        switch (dir)
        {
            case LEFT:
                expand_left(a);
                break;
            case RIGHT:
                expand_right(a);
                break;
            case UP:
                expand_top(a);
                break;
            case DOWN:
                expand_bottom(a);
                break;
            default:
                abort();
        }
    }

    bool intersect(const Rect& other) const
    {
        return left() <= other.right() - 1 && right() - 1 >= other.left()
            && top() <= other.bottom() - 1 && bottom() - 1 >= other.top();
    }

    bool is_neighbor(const Rect& other) const
    {
        return (left() <= other.right() - 1 && right() - 1 >= other.left() &&
                (top() == other.bottom() || bottom() == other.top()))
            || (top() <= other.bottom() - 1 && bottom() - 1 >= other.top() &&
                (left() == other.right() || right() == other.left()));
    }

    bool valid(int w, int h) const
    {
        return width() > 0 && height() > 0
            && left() >= 0 && right() <= w
            && top() >= 0 && bottom() <= h;
    }
    bool valid() const
    {
        return width() > 0 && height() > 0;
    }

    bool operator==(const Rect& other) const
    {
        return pos_ == other.pos_ && w_ == other.w_ && h_ == other.h_;
    }
    bool operator !=(const Rect& other) const
    {
        return !(*this == other);
    }

    Pos pos_;
    int w_, h_;
};
const Rect invalid_rect(-114, -514, -1919, -810);

namespace std
{
    ostream& operator<<(ostream& os, Dir dir)
    {
        os << dir_s[dir];
        return os;
    }

    ostream& operator<<(ostream& os, const Pos& c)
    {
        char buf[256];
        sprintf(buf, "(%d, %d)", c.x, c.y);
        os << buf;
        return os;
    }

    ostream& operator<<(ostream& os, const Rect& r)
    {
        os << r.pos() << ", ";
        char buf[256];
        sprintf(buf, "%d * %d", r.width(), r.height());
        os << buf;
        return os;
    }
}

typedef int Pixel;
class Image
{
public:
    Image(int w, int h, const vector<int>& data, int& stream_i)
        : w_(w), h_(h)
    {
        a = vector<vector<Pixel>>(h, vector<Pixel>(w));
        rep(i, w * h)
            at(i % w, i / w) = data[stream_i++];
    }
    Image(int w, int h)
        : w_(w), h_(h)
    {
        a = vector<vector<Pixel>>(h, vector<Pixel>(w));
    }

    Image()
        : w_(-114514), h_(-1919810)
    {
    }

    int width() const { return w_; }
    int height() const { return h_; }

    Pixel& at(int x, int y)
    {
        assert(in_rect(x, y, width(), height()));
        return a[y][x];
    }
    Pixel& at(const Pos& pos)
    {
        return at(pos.x, pos.y);
    }

    Pixel at(int x, int y) const
    {
        assert(in_rect(x, y, width(), height()));
        return a[y][x];
    }
    Pixel at(const Pos& pos) const
    {
        return at(pos.x, pos.y);
    }

    Image trim(const Rect& rect) const
    {
        assert(in_rect(rect.left(), rect.top(), width(), height()));
        assert(in_rect(rect.right() - 1, rect.bottom() - 1, width(), height()));

        Image trimed(rect.width(), rect.height());
        rep(y, trimed.height()) rep(x, trimed.width())
            trimed.at(x, y) = at(rect.pos() + Pos(x, y));
        return trimed;
    }

    void replace(int lx, int ly, const Image& image)
    {
        rep(y, image.height()) rep(x, image.width())
            at(lx + x, ly + y) = image.at(x, y);
    }

    Image scale(int new_w, int new_h) const
    {
        assert(new_w > 0);
        assert(new_h > 0);

        static vector<int> ori_ys, new_ys, inter_ys;
        ori_ys.clear();
        new_ys.clear();
        inter_ys.clear();
        rep(i, height())
        {
            int y1 = i * new_h, y2 = y1 + new_h;
            rep(j, new_h)
            {
                int y3 = j * height(), y4 = y3 + height();
                int intr = min(y2, y4) - max(y1, y3);
                if (intr > 0)
                {
                    ori_ys.push_back(i);
                    new_ys.push_back(j);
                    inter_ys.push_back(intr);
                }
            }
        }
        static vector<int> ori_xs, new_xs, inter_xs;
        ori_xs.clear();
        new_xs.clear();
        inter_xs.clear();
        rep(i, width())
        {
            int x1 = i * new_w, x2 = x1 + new_w;
            rep(j, new_w)
            {
                int x3 = j * width(), x4 = x3 + width();
                int intr = min(x2, x4) - max(x1, x3);
                if (intr > 0)
                {
                    ori_xs.push_back(i);
                    new_xs.push_back(j);
                    inter_xs.push_back(intr);
                }
            }
        }

        Image scaled(new_w, new_h);
        rep(i, ori_ys.size())
        {
            int ori_y = ori_ys[i];
            int new_y = new_ys[i];
            int intr_y = inter_ys[i];
            rep(j, ori_xs.size())
            {
                int ori_x = ori_xs[j];
                int new_x = new_xs[j];
                int intr_x = inter_xs[j];
                scaled.at(new_x, new_y) += intr_y * intr_x * at(ori_x, ori_y);
            }
        }

        rep(y, scaled.height()) rep(x, scaled.width())
            scaled.at(x, y) = (2 * scaled.at(x, y) + width() * height()) / (2 * width() * height());
        return scaled;
    }

private:
    int w_, h_;
    vector<vector<Pixel>> a;
};

template <typename T>
class Array2D
{
public:
    Array2D(int w, int h)
        : w_(w), h_(h)
    {
    }

    Array2D(int w, int h, const T& init_val)
        : w_(w), h_(h)
    {
        clear(init_val);
    }

    Array2D()
        : w_(-114514), h_(-1919810)
    {
    }

    int width() const { return w_; }
    int height() const { return h_; }

    T& at(int x, int y)
    {
        assert(in_rect(x, y, width(), height()));
        return a[y][x];
    }
    T& at(const Pos& pos)
    {
        return at(pos.x, pos.y);
    }

    void clear(const T& val)
    {
        rep(y, height()) rep(x, width())
            at(x, y) = val;
    }

private:
    int w_, h_;
    T a[512][512];
};


ll sum_sq_diff(const Image& target, const Image& collage)
{
    assert(target.width() == collage.width());
    assert(target.height() == collage.height());

    ll sum = 0;
    rep(y, target.height()) rep(x, target.width())
    {
        int diff = target.at(x, y) - collage.at(x, y);
        sum += diff * diff;
    }
    return sum;
}

ll sum_sq_diff(const Image& target, const Rect& rect, const Image& collage)
{
    assert(rect.valid(target.width(), target.height()));
    assert(rect.width() == collage.width() && rect.height() == collage.height());

    ll sum = 0;
    rep(y, rect.height()) rep(x, rect.width())
    {
        int diff = target.at(rect.pos().x + x, rect.pos().y + y) - collage.at(x, y);
        sum += diff * diff;
    }
    return sum;
}

double score_collage(const Image& target, const Image& collage)
{
    return sqrt(double(sum_sq_diff(target, collage)) / (target.height() * target.width()));
}

Image make_collage(int w, int h, const vector<Image>& source, const vector<Rect>& target_rects)
{
    Image collage(w, h);
    rep(i, SOURCE_IMAGES)
    {
        if (target_rects[i].pos().x >= 0)
        {
            auto& r = target_rects[i];
            Image scaled = source[i].scale(r.width(), r.height());
            rep(y, scaled.height()) rep(x, scaled.width())
                collage.at(r.pos() + Pos(x, y)) = scaled.at(x, y);
        }
    }
    return collage;
}

vector<Rect> list_spaces(Array2D<bool>& used, int max_width, int max_height)
{
    static vector<Rect> spaces;
    spaces.clear();
    rep(ly, used.height()) rep(lx, used.width())
    {
        if (!used.at(lx, ly))
        {
            int hy = ly + 1;
            while (hy < used.height() && hy - ly < max_height && !used.at(lx, hy))
                ++hy;

            int hx = lx + 1;
            while (hx < used.width() && hx - lx < max_width)
            {
                bool ok = true;
                for (int y = ly; y < hy; ++y)
                    ok &= !used.at(hx, y);
                if (!ok)
                    break;

                ++hx;
            }

            assert(lx < hx);
            assert(ly < hy);
            spaces.push_back(Rect(Pos(lx, ly), hx - lx, hy - ly));
            assert(spaces.back().width() <= max_width);
            assert(spaces.back().height() <= max_height);

            for (int y = ly; y < hy; ++y)
            {
                for (int x = lx; x < hx; ++x)
                {
                    assert(!used.at(x, y));
                    used.at(x, y) = true;
                }
            }
        }
    }
#ifndef NDEBUG
    rep(y, used.height()) rep(x, used.width())
        assert(used.at(x, y));
#endif
    return spaces;
};

class Solution
{
public:
    Solution(Image& target, vector<Image>& source)
        : target(&target), source(&source), rects(vector<Rect>(SOURCE_IMAGES, invalid_rect))
    {
    }
    Solution()
        : target(nullptr), source(nullptr)
    {
    }

    Rect& rect(int i)
    {
        assert(0 <= i && i < (int)rects.size());
        return rects[i];
    }
    const Rect& rect(int i) const
    {
        assert(0 <= i && i < (int)rects.size());
        return rects[i];
    }

    vector<int> used_indices() const
    {
        vector<bool> used = used_source_table();
        vector<int> res;
        rep(i, SOURCE_IMAGES)
            if (used[i])
                res.push_back(i);
        return res;
    }

    Image make_collage() const
    {
        Image collage(target->width(), target->height());
        for (int i : used_indices())
        {
            const Rect& r = rect(i);
            Image scaled = (*source)[i].scale(r.width(), r.height());
            collage.replace(r.pos().x, r.pos().y, scaled);
        }
        return collage;
    }

    vector<Rect> list_spaces(int max_width, int max_height) const
    {
        Array2D<bool> used(target->width(), target->height(), false);
        for (auto& r : rects)
        {
            rep(y, r.height()) rep(x, r.width())
                used.at(r.pos().x + x, r.pos().y + y) = true;
        }
        return ::list_spaces(used, max_width, max_height);
    }

    vector<bool> used_source_table() const
    {
        vector<bool> used(SOURCE_IMAGES);
        rep(i, SOURCE_IMAGES)
            used[i] = rects[i].width() > 0;
        return used;
    }

    bool valid() const
    {
        int coverd[512][512];
        clr(coverd, 0);
        rep(i, SOURCE_IMAGES)
        {
            auto& r = rect(i);
            if (r.valid())
            {
                if (!r.valid(target->width(), target->height()))
                    return false;

                rep(y, r.height()) rep(x, r.width())
                    ++coverd[r.pos().y + y][r.pos().x + x];
            }
        }

        rep(y, target->height()) rep(x, target->width())
            if (coverd[y][x] != 1)
                return false;
        return true;
    }

    vector<Rect> used_rects() const
    {
        vector<Rect> res;
        for (int i : used_indices())
            res.push_back(rects[i]);
        return res;
    }

    vector<int> make_result() const
    {
        vector<int> res;
        rep(i, SOURCE_IMAGES)
        {
            if (rects[i].pos().x >= 0)
            {
                auto& r = rects[i];
                res.push_back(r.pos().y);
                res.push_back(r.pos().x);
                res.push_back(r.pos().y + r.height() - 1);
                res.push_back(r.pos().x + r.width() - 1);
            }
            else
            {
                rep(i, 4)
                    res.push_back(-1);
            }
        }
        return res;
    }

private:
    Image* target;
    vector<Image>* source;

    vector<Rect> rects;
};



// #define SOLUTION_LOG
#ifdef SOLUTION_LOG

string output_name;
#define SET_OUTPUT_NAME(name) output_name = name

vector<Solution> solution_logs;
#define ADD_SOLUTION_LOG(solution) solution_logs.push_back(solution)

#include <sys/stat.h>
void output_solution_logs()
{
    if (output_name.empty())
        return;

    int interval = max<int>(1, solution_logs.size() / 100);
    interval = 1;
    vector<int> use_i;
    for (int i = 0; i < (int)solution_logs.size(); i += interval)
        use_i.push_back(i);
    for (int i = max(0, (int)solution_logs.size() - 5); i < solution_logs.size(); ++i)
        use_i.push_back(i);
    uniq(use_i);

    cerr << "output solution logs" << endl;

    mkdir("solution_logs", 0755);
    mkdir(("solution_logs/" + output_name).c_str(), 0755);
    for (int i : use_i)
    {
        char filename[256];
        sprintf(filename, "solution_logs/%s/%06d", output_name.c_str(), i);
        ofstream fs(filename);
        for (int a : solution_logs[i].make_result())
            fs << a << endl;
    }
}
#define OUTPUT_SOLUTION_LOGS() output_solution_logs()

#else
#define SET_OUTPUT_NAME(name)
#define ADD_SOLUTION_LOG(solution)
#define OUTPUT_SOLUTION_LOGS()
#endif

#ifdef LOCAL
const double G_TLE = 10 * 1000;
#else
const double G_TLE = 9.8 * 1000;
#endif

Timer g_timer;

const int FIXED_SIZE = 16;

class Solver
{
public:
    Solver(Image& target, vector<Image>& source)
        : target(target), source(source)
    {
    }

    Solution match_images(const vector<Rect>& target_rects, Solution ori_solution)
    {
        PrimalDual<int, ll> pd(SOURCE_IMAGES + target_rects.size() + 2);
        const int target_begin = SOURCE_IMAGES;
        rep(j, target_rects.size())
        {
#ifndef NO_TIMER
            if (g_timer.get_elapsed() > G_TLE)
            {
                cerr << "match_images fail safe" << endl;
                return ori_solution;
            }
#endif

            Image tage = target.trim(target_rects[j]);
            rep(i, SOURCE_IMAGES)
            {
                if (source[i].width() < tage.width() || source[i].height() < tage.height())
                {
                    pd.add_edge(i, target_begin + j, 0, ten(16));
                }
                else
                {
                    Image scaled = source[i].scale(tage.width(), tage.height());
                    ll d = sum_sq_diff(tage, scaled);
                    pd.add_edge(i, target_begin + j, 1, d);
                }
            }
        }
        const int src = SOURCE_IMAGES + target_rects.size();
        const int sink = src + 1;
        rep(i, SOURCE_IMAGES)
            pd.add_edge(src, i, 1, 0);
        rep(j, target_rects.size())
            pd.add_edge(target_begin + j, sink, 1, 0);

        ll min_cost = pd.min_cost_flow(src, sink, target_rects.size());
        if (min_cost < 0 || min_cost >= ten(15))
            return empty_solution();
        assert(min_cost >= 0);
        assert(min_cost < ten(16));

        Solution solution = empty_solution();
        rep(i, SOURCE_IMAGES) rep(j, target_rects.size())
        {
            if (pd.g[pd.g[i][j].to][pd.g[i][j].rev].cap == 1)
                solution.rect(i) = target_rects[j];
        }
        assert(solution.used_indices().size() == target_rects.size());
        return solution;
    }

    Solution match_images_fast(const vector<Rect>& target_rects)
    {
        PrimalDual<int, ll> pd(SOURCE_IMAGES + target_rects.size() + 2);
        const int target_begin = SOURCE_IMAGES;
        rep(j, target_rects.size())
        {
            Image tage = target.trim(target_rects[j]).scale(FIXED_SIZE, FIXED_SIZE);
            const ll area = target_rects[j].width() * target_rects[j].height();
            assert(area > 0);
            rep(i, SOURCE_IMAGES)
            {
                if (source[i].width() < target_rects[j].width() || source[i].height() < target_rects[j].height())
                {
                    pd.add_edge(i, target_begin + j, 0, ten(15));
                }
                else
                {
                    const static ll FIXED_AREA = FIXED_SIZE * FIXED_SIZE;
                    ll cost = sum_sq_diff(tage, fixed_size_source[i]) * area / FIXED_AREA;
                    pd.add_edge(i, target_begin + j, 1, cost);
                }
            }
        }
        const int src = SOURCE_IMAGES + target_rects.size();
        const int sink = src + 1;
        rep(i, SOURCE_IMAGES)
            pd.add_edge(src, i, 1, 0);
        rep(j, target_rects.size())
            pd.add_edge(target_begin + j, sink, 1, 0);
        ll min_cost = pd.min_cost_flow(src, sink, target_rects.size());
        if (min_cost < 0 || min_cost >= ten(15))
            return empty_solution();
        assert(min_cost >= 0);
        assert(min_cost < ten(16));

        Solution solution = empty_solution();
        rep(i, SOURCE_IMAGES) rep(j, target_rects.size())
        {
            if (pd.g[pd.g[i][j].to][pd.g[i][j].rev].cap == 1)
                solution.rect(i) = target_rects[j];
        }
        assert(solution.used_indices().size() == target_rects.size());
        return solution;
    }

    vector<Rect> grid_rects(int width, int height, int rows, int cols)
    {
        vector<int> ys;
        rep(yi, rows)
            ys.push_back(height / rows * yi);
        ys.push_back(height);
        vector<int> xs;
        rep(xi, cols)
            xs.push_back(width / cols * xi);
        xs.push_back(width);

        vector<Rect> target_rects;
        rep(yi, rows) rep(xi, cols)
            target_rects.push_back(Rect(Pos(xs[xi], ys[yi]), xs[xi + 1] - xs[xi], ys[yi + 1] - ys[yi]));
        return target_rects;
    }

    Solution expand(Solution solution)
    {
        vector<int> used_i = solution.used_indices();
        int expand_i = used_i[rand() % used_i.size()];
        Rect& r = solution.rect(expand_i);
        Image& image = source[expand_i];

        const int rem_w = image.width() - r.width();
        const int rem_h = image.height() - r.height();
        assert(rem_w >= 0);
        assert(rem_h >= 0);

        Dir dir = Dir(rand() % 4);
        int div = (dir == LEFT || dir == RIGHT ? rem_w : rem_h);
        if (div == 0)
            return empty_solution();
        int len = 1 + rand() % div;
        r.expand(dir, len);
        if (r.width() > image.width() || r.height() > image.height() || !r.valid(target.width(), target.height()))
            return empty_solution();

        Solution res = empty_solution();
        for (int i : solution.used_indices())
        {
            if (i != expand_i && r.intersect(solution.rect(i)))
                solution.rect(i).expand(rev_dir(dir), -len);

            if (solution.rect(i).valid(target.width(), target.height()))
                res.rect(i) = solution.rect(i);
        }
        return res;
    }

    Solution shrink(Solution solution)
    {
        vector<int> used_i = solution.used_indices();
        const int shrink_i = used_i[rand() % used_i.size()];
        Rect& r = solution.rect(shrink_i);

        Dir dir = Dir(rand() % 4);
        int div = (dir == LEFT || dir == RIGHT ? r.width() : r.height());
        if (div == 0)
            return empty_solution();
        int len = 1 + rand() % div;
        r.expand(dir, -len);

        Solution res = empty_solution();
        for (int i : solution.used_indices())
        {
            if (solution.rect(i).valid(target.width(), target.height()))
                res.rect(i) = solution.rect(i);
        }
        return res;
    }

    Solution empty_solution() { return Solution(target, source); }

    Solution fill_space(Solution solution)
    {
        vector<bool> used = solution.used_source_table();
        vector<int> use_i;
        rep(i, SOURCE_IMAGES)
            if (!used[i])
                use_i.push_back(i);
        vector<Rect> spaces = solution.list_spaces(60, 60);
        PrimalDual<int, ll> pd(use_i.size() + spaces.size() + 2);
        if (use_i.size() < spaces.size())
            return empty_solution();

        const int target_begin = use_i.size();
        rep(j, spaces.size())
        {
            const ll area = spaces[j].width() * spaces[j].height();
            assert(area > 0);
            Image tage = target.trim(spaces[j]).scale(FIXED_SIZE, FIXED_SIZE);
            rep(i, use_i.size())
            {
                Image& src = source[use_i[i]];
                if (src.width() < spaces[j].width() || src.height() < spaces[j].height())
                {
                    pd.add_edge(i, target_begin + j, 0, ten(15));
                }
                else
                {
                    const static ll FIXED_AREA = FIXED_SIZE * FIXED_SIZE;
                    ll cost = sum_sq_diff(tage, fixed_size_source[use_i[i]]) * area / FIXED_AREA;
                    pd.add_edge(i, target_begin + j, 1, cost);
                }
            }
        }

        const int src = use_i.size() + spaces.size();
        const int sink = src + 1;
        rep(i, use_i.size())
            pd.add_edge(src, i, 1, 0);
        rep(j, spaces.size())
            pd.add_edge(target_begin + j, sink, 1, 0);
        ll min_cost = pd.min_cost_flow(src, sink, spaces.size());
        if (min_cost < 0 || min_cost >= ten(15))
            return empty_solution();

        rep(i, use_i.size()) rep(j, spaces.size())
        {
            if (pd.g[pd.g[i][j].to][pd.g[i][j].rev].cap == 1)
                solution.rect(use_i[i]) = spaces[j];
        }
        assert(solution.valid());
        return solution;
    }


    Solution improve(Solution solution, double tle)
    {
        Timer l_timer;
        l_timer.start();

        assert(solution.valid());
        ll cur_score = sum_sq_diff(target, solution.make_collage());
        ll best_score = cur_score;
        for (int loop = 0; ; ++loop)
        {
#ifdef NO_TIMER
            if (loop > 20000)
                break;
#else
            if (l_timer.get_elapsed() > tle)
                break;
#endif

            Solution nsol;
            int ra = rand() % 10000;
            if (ra < 15)
            {
//                 double prev = score_collage(target, solution.make_collage());
//                 nsol = match_images(solution.used_rects());
                nsol = match_images_fast(solution.used_rects());
//                 double cur = score_collage(target, nsol.make_collage());
//                 fprintf(stderr, "%.4f -> %.4f\n", prev, cur);
            }
            else
            if (ra < 9000)
                nsol = expand(solution);
            else
                nsol = shrink(solution);

            if (nsol.used_indices().empty())
                continue;

            nsol = fill_space(nsol);
            if (!nsol.used_indices().empty())
            {
                assert(nsol.valid());
                ll score = cur_score + diff_sum_sq_diff(nsol, solution);
//                 assert(score == sum_sq_diff(target, nsol.make_collage()));
                if (score < best_score)
                {
//                     fprintf(stderr, "%6d: %3d, %.5f\n", loop, (int)nsol.used_indices().size(), sqrt(double(score) / (target.width() * target.height())));
                    best_score = score;
                    cur_score = score;
                    solution = nsol;

                    ADD_SOLUTION_LOG(solution);
                }
            }
        }

        return solution;
    }

    ll diff_sum_sq_diff(Solution& cur, Solution& prev)
    {
        assert(cur.valid());
        assert(prev.valid());

        ll diff_ssd = 0;
        rep(i, SOURCE_IMAGES)
        {
            Rect& ar = cur.rect(i);
            Rect& br = prev.rect(i);
            if (ar != br)
            {
                ll a_ssd = ar.valid() ? sum_sq_diff(target, ar, source[i].scale(ar.width(), ar.height())) : 0;
                ll b_ssd = br.valid() ? sum_sq_diff(target, br, source[i].scale(br.width(), br.height())) : 0;
                diff_ssd += a_ssd - b_ssd;
            }
        }
        return diff_ssd;
    }

    Solution solve()
    {
        for (Image& image : source)
            fixed_size_source.push_back(image.scale(FIXED_SIZE, FIXED_SIZE));

        vector<Rect> target_rects = grid_rects(target.width(), target.height(), 7, 7);
        double match_time_cost = g_timer.get_elapsed();
//         Solution solution = match_images(target_rects);
        Solution solution = match_images_fast(target_rects);
        match_time_cost = g_timer.get_elapsed() - match_time_cost;
        assert(solution.valid());

        ADD_SOLUTION_LOG(solution);

        dump(match_time_cost / 1000);

        solution = improve(solution, G_TLE * 0.73 - match_time_cost);
        assert(solution.valid());

        ADD_SOLUTION_LOG(solution);

        dump(g_timer.get_elapsed());
        solution = match_images(solution.used_rects(), solution);
        dump(g_timer.get_elapsed());


        ADD_SOLUTION_LOG(solution);

        OUTPUT_SOLUTION_LOGS();

        return solution;
    }

private:
    Image target;
    vector<Image> source;

    vector<Image> fixed_size_source;
};

class CollageMaker
{
public:
    vector<int> compose(vector<int>& data)
    {
        g_timer.start();

        int stream_i = 0;
        Image target = input_image(data, stream_i);
        vector<Image> source;
        rep(i, SOURCE_IMAGES)
            source.push_back(input_image(data, stream_i));

        Solution solution = Solver(target, source).solve();
        dump(g_timer.get_elapsed() / 1000);

        return solution.make_result();
    }

private:

    Image input_image(const vector<int>& data, int& stream_i)
    {
        int h = data[stream_i++];
        int w = data[stream_i++];
        Image image(w, h, data, stream_i);
        assert(stream_i <= (int)data.size());
        return image;
    }
};


#ifdef LOCAL
int main(int argc, char** argv)
{
    if (argc > 1)
        SET_OUTPUT_NAME(argv[1]);

    int n;
    cin >> n;
    vector<int> data(n);
    input(data, n);

    vector<int> res = CollageMaker().compose(data);
    assert(res.size() == 800);
    for (int i : res)
        cout << i << endl;
    cout.flush();
}
#endif
