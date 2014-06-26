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

    Pos pos_;
    int w_, h_;
};

class Image
{
public:
    Image(int w, int h, const vector<int>& data, int& stream_i)
        : w_(w), h_(h)
    {
        rep(i, w * h)
            at(i % w, i / w) = data[stream_i++];
    }
    Image(int w, int h)
        : w_(w), h_(h)
    {
        rep(y, h) rep(x, w)
            at(x, y) = 0;
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

    Image scale(int new_w, int new_h)
    {
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
    int a[512][512];
};


double score_collage(Image& target, Image& collage)
{
    assert(target.width() == collage.width());
    assert(target.height() == collage.height());

    double score = 0;
    rep(y, target.height()) rep(x, target.width())
    {
        assert(0 <= collage.at(x, y) && collage.at(x, y) < 256);

        int diff = target.at(x, y) - collage.at(x, y);
        score += diff * diff;
    }

    return sqrt(score / (target.height() * target.width()));
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

class Solver
{
public:
    vector<Rect> solve(Image& target, vector<Image>& source)
    {
//         const int rows = 10;
//         const int cols = SOURCE_IMAGES / rows;
        int rows = 7, cols = 7;
        vector<int> ys;
        rep(yi, rows)
            ys.push_back(target.height() / rows * yi);
        ys.push_back(target.height());
        vector<int> xs;
        rep(xi, cols)
            xs.push_back(target.width() / cols * xi);
        xs.push_back(target.width());

        vector<Rect> target_rects;
        rep(i, rows * cols)
        {
            const int xi = i % cols;
            const int yi = i / cols;
            target_rects.push_back(Rect(Pos(xs[xi], ys[yi]), xs[xi + 1] - xs[xi], ys[yi + 1] - ys[yi]));
        }


        PrimalDual<int, double> pd(source.size() + target_rects.size() + 2);
        const int target_begin = source.size();
        rep(j, target_rects.size())
        {
            Image tage = target.trim(target_rects[j]);
            rep(i, source.size())
            {
                if (source[i].width() < tage.width() || source[i].height() < tage.height())
                {
                    pd.add_edge(i, target_begin + j, 1, 1e60);
                }
                else
                {
                    Image scaled = source[i].scale(tage.width(), tage.height());
                    double score = score_collage(tage, scaled);
                    pd.add_edge(i, target_begin + j, 1, score);
                }
            }
        }
        const int src = source.size() + target_rects.size();
        const int sink = src + 1;
        rep(i, source.size())
            pd.add_edge(src, i, 1, 0);
        rep(j, target_rects.size())
            pd.add_edge(target_begin + j, sink, 1, 0);
        double min_cost = pd.min_cost_flow(src, sink, target_rects.size());
        assert(min_cost >= 0);
        assert(min_cost < 1e50);

        vector<Rect> result_rects(SOURCE_IMAGES, Rect(Pos(-1919, 810), -1, -1));
        rep(i, source.size()) rep(j, target_rects.size())
        {
            if (pd.g[i][j].cap == 0)
                result_rects[i] = target_rects[j];
        }
        return result_rects;
    }
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
        int stream_i = 0;
        Image target = input_image(data, stream_i);
        vector<Image> source;
        rep(i, SOURCE_IMAGES)
            source.push_back(input_image(data, stream_i));

//         analyze(source);
//         exit(1);


        vector<Rect> result_rects = Solver().solve(target, source);

//         Image collage = make_collage(target.width(), target.height(), source, result_rects);
//         dump(score_collage(target, collage));


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
