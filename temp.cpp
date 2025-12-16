/* CHT_MAX
struct Line {
    ll a, b;
    Line (ll a, ll b) : a(a), b(b) {}

    ll compute(ll x) {
        return a * x + b;
    }

    ld intersect(Line other) {
        return (other.b - b) / (a - other.a);
    }

    friend bool bad(Line d1, Line d2, Line d3) {
        return d1.intersect(d2) <= d2.intersect(d3);
    }
};

struct CHT {
    vector <Line> hull;

    void add(Line d) {
        while(SZ(hull) >= 2 && bad(d, hull.back(), hull[SZ(hull) - 2])) hull.pop_back();
        hull.push_back(d);
    }

    ll query(ll x) {
        int l = 0, r = SZ(hull) - 1;

        while(l < r) {
            int mid = (l + r) >> 1;
            if(hull[mid].compute(x) <= hull[mid + 1].compute(x)) l = mid + 1;
            else r = mid;
        }

        return hull[l].compute(x);
    }
};
*/

/*  chtmin
struct Line {
    ll a, b;

    Line(ll a, ll b) : a(a), b(b) { }

    ll compute(ll x) {
        return a * x + b;
    }

    ld intersect(Line other) {
        return (other.b - b) / (a - other.a);
    }

    friend bool bad(Line d1, Line d2, Line d3) {
        return d1.intersect(d2) >= d2.intersect(d3);
    }
};

struct CHT {
    vector <Line> hull;

    void add(Line d) {
        while(SZ(hull) >= 2 && bad(d, hull.back(), hull[SZ(hull) - 2])) {
            hull.pop_back();
        }
        hull.push_back(d);
    }

    ll query(ll x) {
        int l = 0, r = SZ(hull) - 1;

        while(l < r) {
            int mid = (l + r) >> 1;
            if(hull[mid].compute(x) >= hull[mid + 1].compute(x)) l = mid + 1;
            else r = mid;
        }

        return hull[l].compute(x);
    }
};
*/
/* kmp
void build_kmp(void) {
    int k = kmp[1] = 0;
    for (int i = 2; i <= m; i++) {
        while (k > 0 && t[i] != t[k + 1]) k = kmp[k];
        kmp[i] = (t[i] == t[k + 1] ? ++ k : 0);
    }

    k = 0;
    for (int i = 1; i <= n; i++) {
        while (k > 0 && s[i] != t[k + 1]) k = kmp[k];
        match[i] = (s[i] == t[k + 1] ? ++ k : 0);
    }
}
*/
/*
int sa[N], tmp[N], pos[N];
int gap = 0;

bool cmp_sa(int i, int j){
    if (pos[i] != pos[j]) return pos[i] < pos[j];
    i+= gap; j+= gap;
    return i <= n && j <= n ? pos[i] < pos[j] : i > j;
}

void build_sa(){
    for (int i = 1; i <= n; i++){
        sa[i] = i;
        pos[i] = s[i];
    }
    for (gap = 1;; gap <<= 1){
        sort (sa + 1, sa + 1 + n, cmp_sa);

        for (int i = 2; i <= n; i++) tmp[i] = tmp[i - 1] + cmp_sa(sa[i - 1], sa[i]);
        for (int i = 1; i <= n; i++) pos[sa[i]] = tmp[i];
        if (tmp[n] == n - 1) break;
    }
}

inline pair <int, int> get(int l, int r){
    int t1 = (hsh[r].fi - 1LL * hsh[l - 1].fi * pw[r - l + 1].fi % sm1 + 1LL * sm1 * sm1) % sm1;
    int t2 = (hsh[r].se - 1LL * hsh[l - 1].se * pw[r - l + 1].se % sm2 + 1LL * sm2 * sm2) % sm2;
    return make_pair(t1, t2);
}

void build_hash(){
    pw[0] = {1, 1};
    for (int i = 1; i <= n; i++) {
        pw[i].fi = 1LL * pw[i - 1].fi * base % sm1;
        pw[i].se = 1LL * pw[i - 1].se * base % sm2;

        hsh[i].fi = 1LL * hsh[i - 1].fi * base % sm1 + (s[i] + 'a') % sm1;
        hsh[i].se = 1LL * hsh[i - 1].se * base % sm2 + (s[i] + 'a') % sm2;
    }
}
*/
/*
struct Points {
    ll x, y;

        A = y2 - y1;
        B = x1 - x2
        C = x1 * y2 - x2 * y1

        double det = A1 * B2 - A2 * B1;
        if (det == 0) {
            // Lines are parallel or coincident
            if (A1 * C2 == A2 * C1) {
                // Lines are coincident
            }
            else {
                // Lines are parallel
            }
        }
        else {
            double x = (B2 * C1 - B1 * C2) / det;
            double y = (A1 * C2 - A2 * C1) / det;
        }

} P[maxn];
/*

/*
ll ccw(Points a, Points b, Points c) {
    return (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x);
}

ll dttri(Points a, Points b, Points c) {
    /// ab x ac
    return abs((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

ll dist(Points a, Points b) {
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

const double eps = 1e-9;

int sign(double x) {
    if(x > eps) return 1;
    if(x < -eps) return -1;
    return 0;
}

double cross(Points AB, Points AC) {
    return AB.x * AC.y - AC.x * AB.y;
}

double dot(Points AB, Points AC) {
    return AB.x * AC.x + AB.y * AC.y;
}

Points operator -(Points a, Points b) {
    return {a.x - b.x, a.y - b.y};
}

bool intersect(Points A, Points B, Points C, Points D) {
    int ABxAC = sign(cross(B - A, C - A));
    int ABxAD = sign(cross(B - A, D - A));
    int CDxCA = sign(cross(D - C, A - C));
    int CDxCB = sign(cross(D - C, B - C));
    if (ABxAC == 0 | ABxAD == 0 || CDxCA == 0 || CDxCB == 0) {
        // C on segment AB if ABxAC = 0 and CA.CB <= 0
        if (ABxAC == 0 && sign(dot(A - C, B - C)) <= 0) return true;
        if (ABxAD == 0 && sign(dot(A - D, B - D)) <= 0) return true;
        if (CDxCA == 0 && sign(dot(C - A, D - A)) <= 0) return true;
        if (CDxCB == 0 && sign(dot(C - B, D - B)) <= 0) return true;
        return false;
    }
    return (ABxAC * ABxAD < 0 && CDxCA * CDxCB < 0);
}
*/

/*
struct Node{
	Node *link, *go[26], *child[26];
	vector <Node*> g;
	Node(){
		link = NULL;
		g.clear();
		FOR(i, 0, 25) go[i] = child[i] = NULL;
	}
};
Node *root = new Node();
inline void addString(const string &s){
	Node *p = root;
	for(const char &c : s){
		int j = c - 'a';
		if(p -> child[j] == NULL) p -> child[j] = new Node();
		p = p -> child[j];
	}
}
void AC_BFS(void){
    root -> link = root;
    FOR(i, 0, 25) root -> go[i] = root -> child[i] != NULL ? root -> child[i] : root;
    queue <Node*> que;
    que.push(root);
    while(!que.empty()){
        Node *p = que.front(); que.pop();
        FOR(i, 0, 25) if(p -> child[i] != NULL){
            Node *q = p -> child[i];
            q -> link = p == root ? p : p -> link -> go[i];
            q -> link -> g.push_back(q);
            FOR(j, 0, 25) q -> go[j] = q -> child[j] != NULL ? q -> child[j] : q -> link -> go[j];
            que.push(q);
        }
    }
}
*/
/*
#include <bits/stdc++.h>
#define ALL(A) (A).begin(), (A).end()
#define TIME  (1.0 * clock() / CLOCKS_PER_SEC)
#define file(Roxy) if(fopen(Roxy".inp", "r")){freopen(Roxy".inp", "r", stdin);freopen(Roxy".out", "w", stdout);}
#define fi first
#define se second
#define FOR(i, a, b) for (int i = (a); i <= (b); i++)
#define FOD(i, a, b) for (int i = (a); i >= (b); i--)
#define REP(i, n) for (int i = 0; i < (n); i++)
using namespace std;

const int MAXN = 2e5 + 5;
const string NAME = "";
const int TEST = 100;

mt19937 rd(time(0));

int Rand(int l, int r) {
    return l + rd()%(r - l + 1);
}

void create() {
    ofstream inp((NAME + ".inp").c_str());

    inp.close();
}

signed main() {
    srand(time(NULL));
    for(int i = 1; i <= TEST; i++) {
        create();
        system((NAME + "_trau.exe").c_str());
        system((NAME + ".exe").c_str());
        if(system(("fc " + NAME + ".out " + NAME + ".ans").c_str()) != 0) {
            cout << "WRONG !!";
            return 0;
        }
        cout << "CORRECT\n";
    }
    return 0;
}
*/
/*
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

template <class A, class B>
    bool maximize(A &a, const B b) {
        if (a < b) {
            a = b;
            return true;
        } return false;
    }

template <class A, class B>
    bool minimize(A &a, const B b) {
        if(a > b) {
            a = b;
            return true;
        } return false;
    }

using pII = pair <int, int>;
using vI = vector <int>;
using vL = vector <long long>;
using pLI = pair <long long, int>;

#define BIT(mask, i) ((mask >> (i)) & 1)
#define MASK(a) (1LL<<(a))
#define FOR(i, a, b) for(int i = a; i <= (int)b; i++)
#define FORD(i, a, b) for(int i = a; i >= (int)b; i--)
#define fi first
#define se second
#define pb push_back
#define all(a) a.begin(), a.end()
#define sz(a) (int)a.size()

const int mod = 1e9 + 7;
const int maxn = 2e5 + 5;

void process() {

}

int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);

    #define taskname "kieuoanh"
    if(fopen(taskname".inp", "r")) {
        freopen(taskname".inp", "r", stdin);
        freopen(taskname".out", "w", stdout);
    }

    int tc = 1; /// cin >> tc;
    while (tc--) process();

    return 0;
}
*/

