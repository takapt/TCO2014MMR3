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
    Rect()
        : pos_(Pos(1919810, 114514)), w_(364364)
    {
    }

    Pos& pos() { return pos_; }
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
        }
    }

    bool intersect(const Rect& other) const
    {
        return left() <= other.right() - 1 && right() - 1 >= other.left()
            && top() <= other.bottom() - 1 && bottom() - 1 >= other.top();
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

    ostream& operator<<(ostream& os, Rect& r)
    {
        os << r.pos() << ", ";
        char buf[256];
        sprintf(buf, "%d * %d", r.width(), r.height());
        os << buf;
        return os;
    }
}

class Image
{
public:
    Image(int w, int h, const vector<int>& data, int& stream_i)
        : w_(w), h_(h)
    {
        a = vector<vector<int>>(h, vector<int>(w));
        rep(i, w * h)
            at(i % w, i / w) = data[stream_i++];
    }
    Image(int w, int h)
        : w_(w), h_(h)
    {
        a = vector<vector<int>>(h, vector<int>(w));
    }

    Image()
        : w_(-114514), h_(-1919810)
    {
    }

    int width() const { return w_; }
    int height() const { return h_; }

    int& at(int x, int y)
    {
        assert(in_rect(x, y, width(), height()));
        return a[y][x];
    }
    int& at(const Pos& pos)
    {
        return at(pos.x, pos.y);
    }

    Image trim(Rect& rect)
    {
        assert(in_rect(rect.left(), rect.top(), width(), height()));
        assert(in_rect(rect.right() - 1, rect.bottom() - 1, width(), height()));

        Image trimed(rect.width(), rect.height());
        rep(y, trimed.height()) rep(x, trimed.width())
            trimed.at(x, y) = at(rect.pos() + Pos(x, y));
        return trimed;
    }

    void replace(int lx, int ly, Image& image)
    {
        rep(y, image.height()) rep(x, image.width())
            at(lx + x, ly + y) = image.at(x, y);
    }

    Image scale(int new_w, int new_h)
    {
//         assert(0 < new_w && new_w <= width());
//         assert(0 < new_h && new_h <= height());

        vector<int> ori_ys, new_ys, inter_ys;
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
        vector<int> ori_xs, new_xs, inter_xs;
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
//     int a[512][512];
    vector<vector<int>> a;
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


ll sum_sq_diff(Image& target, Image collage)
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

ll sum_sq_diff(Image& target, Rect& rect, Image collage)
{
    if (!(rect.valid(target.width(), target.height())))
        fprintf(stderr, "%d %d %d %d\n", rect.pos().x, rect.pos().y, rect.width(), rect.height());
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

double score_collage(Image& target, Image collage)
{
    return sqrt(double(sum_sq_diff(target, collage)) / (target.height() * target.width()));
}

Image make_collage(int w, int h, vector<Image>& source, vector<Rect>& target_rects)
{
    Image collage(w, h);
    rep(i, SOURCE_IMAGES)
    {
        if (target_rects[i].pos().x >= 0)
        {
            Rect& r = target_rects[i];
            Image scaled = source[i].scale(r.width(), r.height());
            rep(y, scaled.height()) rep(x, scaled.width())
                collage.at(r.pos() + Pos(x, y)) = scaled.at(x, y);
        }
    }
    return collage;
}

vector<Rect> list_spaces(Array2D<bool> used, int max_width = 50, int max_height = 50)
{
    vector<Rect> spaces;
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
        : target(&target), source(&source), rects(vector<Rect>(SOURCE_IMAGES, Rect(Pos(-1919, 810), -1, -1)))
    {
    }
    Solution()
        : target(nullptr), source(nullptr)
    {
    }

    Rect& rect(int i)
    {
        assert(0 <= i && i < rects.size());
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

    Image make_collage()
    {
        Image collage(target->width(), target->height());
        for (int i : used_indices())
        {
            Rect& r = rect(i);
            Image scaled = (*source)[i].scale(r.width(), r.height());
            collage.replace(r.pos().x, r.pos().y, scaled);
        }
        return collage;
    }

    vector<Rect> list_spaces(int max_width, int max_height)
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

    bool valid()
    {
        vector<int> ui = used_indices();
        rep(ii, ui.size()) rep(jj, ii)
        {
            int i = ui[ii], j = ui[jj];
            if (rects[i].intersect(rects[j]))
            {
                fprintf(stderr, "%3d %3d %3d %3d\n", rects[i].left(), rects[i].right(), rects[i].top(), rects[i].bottom());
                fprintf(stderr, "%3d %3d %3d %3d\n", rects[j].left(), rects[j].right(), rects[j].top(), rects[j].bottom());
                abort();
                return false;
            }
        }

        return true;
    }

    vector<Rect> used_rects() const
    {
        vector<Rect> res;
        for (int i : used_indices())
            res.push_back(rects[i]);
        return res;
    }

private:
    vector<Rect> rects;

    Image* target;
    vector<Image>* source;
};


#ifdef LOCAL
const double G_TLE = 10 * 1000;
#else
const double G_TLE = 9.8 * 1000;
#endif

Timer g_timer;

const int FIXED_SIZE = 32;

class Solver
{
public:
    Solver(Image& target, vector<Image>& source)
        : target(target), source(source)
    {
    }

    Solution match_images(vector<Rect> target_rects)
    {
        PrimalDual<int, ll> pd(source.size() + target_rects.size() + 2);
        const int target_begin = source.size();
        rep(j, target_rects.size())
        {
            Image tage = target.trim(target_rects[j]);
            rep(i, source.size())
            {
                if (source[i].width() < tage.width() || source[i].height() < tage.height())
                {
                    pd.add_edge(i, target_begin + j, 0, ten(15));
                }
                else
                {
                    Image scaled = source[i].scale(tage.width(), tage.height());
                    ll d = sum_sq_diff(tage, scaled);
                    pd.add_edge(i, target_begin + j, 1, d);
                }
            }
        }
        const int src = source.size() + target_rects.size();
        const int sink = src + 1;
        rep(i, source.size())
            pd.add_edge(src, i, 1, 0);
        rep(j, target_rects.size())
            pd.add_edge(target_begin + j, sink, 1, 0);
        ll min_cost = pd.min_cost_flow(src, sink, target_rects.size());
        assert(min_cost >= 0);
        assert(min_cost < ten(12));

        Solution solution = empty_solution();
        rep(i, source.size()) rep(j, target_rects.size())
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

    Solution expand(Solution solution, int expand_i)
    {
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

    Solution shrink(Solution solution, int shrink_i)
    {
        Rect& r = solution.rect(shrink_i);
        Image& image = source[shrink_i];

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
        vector<Rect> spaces = solution.list_spaces(60, 60);
        for (auto& space : spaces)
        {
            Image trimed_space = target.trim(space).scale(FIXED_SIZE, FIXED_SIZE);

            int best_i = -1;
            ll best_score = ten(18);
            rep(i, source.size())
            {
                if (!used[i] && source[i].width() >= space.width() && source[i].height() >= space.height())
                {
                    // naive
//                     Image scaled = source[i].scale(space.width(), space.height());
//                     ll score = sum_sq_diff(target, space, scaled);

                    // fast
                    ll score = sum_sq_diff(trimed_space, fixed_size_source[i]);

                    if (score < best_score)
                    {
                        best_score = score;
                        best_i = i;
                    }
                }
            }
            if (best_i == -1)
                return empty_solution();

            used[best_i] = true;
            solution.rect(best_i) = space;
        }
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
//             if (l_timer.get_elapsed() > tle)
//                 break;
            if (loop > 10000)
                break;

            const vector<int> ui = solution.used_indices();
            Solution nsol;
            int ra = rand() % 1000;
            if (ra == 0)
            {
                double prev = score_collage(target, solution.make_collage());
                nsol = match_images(solution.used_rects());
                double cur = score_collage(target, nsol.make_collage());
                fprintf(stderr, "%.4f -> %.4f\n", prev, cur);
            }
            else
            if (ra < 900)
            {
                int expand_i = rand() % ui.size();
                nsol = expand(solution, ui[expand_i]);
            }
            else
            {
                int shrink_i = rand() % ui.size();
                nsol = shrink(solution, ui[shrink_i]);
            }

            if (nsol.used_indices().empty())
                continue;

            nsol = fill_space(nsol);
            if (!nsol.used_indices().empty())
            {
                assert(nsol.valid());
                ll score = cur_score + diff_sum_sq_diff(nsol, solution);
                assert(score == sum_sq_diff(target, nsol.make_collage()));
                if (score < best_score)
//                 if (score < best_score * 0.999)
                {
                    fprintf(stderr, "%4d: %3d, %.5f\n", loop, (int)nsol.used_indices().size(), sqrt(double(score) / (target.width() * target.height())));
                    best_score = score;
                    cur_score = score;
                    solution = nsol;
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
        rep(i, source.size())
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
#ifndef NDEBUG
        ll a_ssd = sum_sq_diff(target, cur.make_collage());
        ll b_ssd = sum_sq_diff(target, prev.make_collage());
        assert(diff_ssd == a_ssd - b_ssd);
#endif
        return diff_ssd;
    }

    Solution solve()
    {
        for (Image& image : source)
            fixed_size_source.push_back(image.scale(FIXED_SIZE, FIXED_SIZE));

        vector<Rect> target_rects = grid_rects(target.width(), target.height(), 7, 7);
        double match_time_cost = g_timer.get_elapsed();
        Solution solution = match_images(target_rects);
        match_time_cost = g_timer.get_elapsed() - match_time_cost;
        assert(solution.valid());

        double rem_time = G_TLE - g_timer.get_elapsed();
        solution = improve(solution, rem_time - match_time_cost * 3);
        assert(solution.valid());

//         dump(score_collage(target, solution.make_collage()));
        vector<Rect> rects;
        for (int i : solution.used_indices())
            rects.push_back(solution.rect(i));
        solution = match_images(rects);
//         dump(score_collage(target, solution.make_collage()));

        return solution;
    }

private:
    Image target;
    vector<Image> source;

    vector<Image> fixed_size_source;
};

void analyze(vector<Image> images)
{
    int freq[256] = {};
    int ave_freq[256] = {};
    for (auto& image : images)
    {
        int ave = 0;
        rep(y, image.height()) rep(x, image.width())
        {
            assert(0 <= image.at(x, y) && image.at(x, y) < 256);
            ++freq[image.at(x, y)];

            ave += image.at(x, y);
        }
        ave /= image.width() * image.height();
        ++ave_freq[ave];
    }

    int ma = *max_element(ave_freq, ave_freq + 256);
    rep(i, 256)
        fprintf(stderr, "%3d| %s\n", i, string(300 * ave_freq[i] / ma, '*').c_str());
}

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
        vector<Rect> result_rects(SOURCE_IMAGES, Rect(Pos(-114514, 1919810), -1, -1));
        for (int i : solution.used_indices())
            result_rects[i] = solution.rect(i);

        dump(g_timer.get_elapsed() / 1000);

        return make_result(result_rects);
    }

private:
    vector<int> make_result(vector<Rect>& result_rects)
    {
        vector<int> res;
        rep(i, SOURCE_IMAGES)
        {
            if (result_rects[i].pos().x >= 0)
            {
                Rect& r = result_rects[i];
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

    Image input_image(const vector<int>& data, int& stream_i)
    {
        int h = data[stream_i++];
        int w = data[stream_i++];
        Image image(w, h, data, stream_i);
        assert(stream_i <= data.size());
        return image;
    }
};


#ifdef LOCAL
int main(int argc, char** argv)
{
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
