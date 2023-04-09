// Code By Arindam Dey
#include <bits/stdc++.h>
using namespace std;
#define mpp make_pair
#define pb push_back
#define em emplace_back
#define PI 3.1415926535897932384626433832795
#define mod 1000000007
#define ff first
#define ss second
#define endl "\n"
#define MOD 998244353
#define ain(a,n) for(ll i=0;i<n;i++)cin>>(a)[i]
#define all(str) (str).begin(), (str).end()
#define IOS                  \
    ios::sync_with_stdio(0); \
    cin.tie(0);              \
    cout.tie(0)
#define ll long long int
#define ld long double
typedef vector<ll> vi;
typedef vector<pair<ll,ll>> vii;
typedef vector<string> vs;
typedef vector<vector<ll>> vvi;
const ll INF=1e18;


// Code By Arindam Dey
bool isPrime(ll n){for(ll i=2;i*i<=n;i++){if(n%i==0)return false;}return true;}
void primeSieve(vector<bool> &sieve,ll maxi) {sieve[0] = sieve[1] = false; for (ll i = 2; i*i<= maxi; i++) { if (sieve[i]) {for (ll j = i*i; j <= maxi; j += i) sieve[j] = false; }}}
ll gcd(ll a,ll b){if(b==0)return a;return gcd(b,a%b);}
ll lcm(ll a,ll b){return (a / gcd(a, b)) * b;}
ll noofdigits(ll n) { return ((log10(n))+1); }
bool isPowerOf2(ll n) {return (n>0 && (n&(n-1))==0);}
map<ll,ll> primeFactorise(ll n){map<ll,ll> mpp;for(ll i=2;i*i<=n;i++){if(n%i==0){ll cnt=0;while(n%i==0){++cnt;n=n/i;}mpp[i]=cnt;}}if(n>1)mpp[n]=1;return mpp;}


ll GN=1e5+10;
vector<ll> fact(1e6+1,1);
// vector<bool> sieve(GN,1);
vector<ll> dsieve(GN,0);
 
/*---------------------------------------------------------------------------------------------------------------------------*/
 
// bool isPrime(ll n){for(ll i=2;i*i<=n;i++){if(n%i==0)return false;}return true;}
// void primeSieve(){sieve[0] = sieve[1] = false;for (ll i = 2; i*i<GN; i++) {if (sieve[i]) {for (ll j = i*i; j<GN; j += i)sieve[j] = false;}}}
void divSieve(){for(ll i=1;i<GN;i++){for(ll j=i;j<GN;j=j+i){++dsieve[j];}}}
void factSieve(ll m){for(ll i=1;i<GN;i++)fact[i]=(fact[i-1]%m*1LL*i%m)%m;}
//map<ll,> primeFactorise(ll n){mi mpp;for(ll i=2;i*i<=n;i++){if(n%i==0){ll cnt=0;while(n%i==0){++cnt;n=n/i;}mpp[i]=cnt;}}if(n>1)mpp[n]=1;return mpp;}
// ll gcd(ll a,ll b){if(b==0)return a;return gcd(b,a%b);}
// ll binexp(ll x,ll n){ll a=x;ll prod=1;while(n){if(n%2==1)prod=prod*a;a=a*a;n=n/2;}return prod;}
ll modexp(ll a,ll b,ll m) {a %= m;ll res = 1LL;while (b > 0) {if (b & 1)res = (res%m *1LL* a%m) % m;a = (a%m *1LL* a%m) % m;b >>= 1;}return res%m;}
ll modadd(ll a,ll b,ll m){return (a%m + b%m)%m;}
ll modsub(ll a,ll b,ll m){return (a%m - b%m+m)%m;}
ll modmul(ll a,ll b,ll m){return (a%m*1LL * b%m)%m;}
ll modinv(ll a,ll m){return modexp(a,m-2,m);}
ll moddiv(ll a,ll b,ll m){return (a%m*1LL*modinv(b,m)%m)%m;}
ll nCr(ll n,ll r,ll P){if(n<r||r<0)return 0;if(r==1)return n;ll nmr=fact[n]%P;ll dmr=(fact[r]%P*1LL*fact[n-r]%P)%P;return (nmr%P *1LL* modinv(dmr,P)%P)%P;}
  
//Segment Tree
// Start with v=1 , vi t(4*n), tl=0, tr=n-1
void buildtree(vi& a, ll v, ll tl, ll tr,vi& t)
{
    if (tl == tr)t[v] = a[tl];
    else {ll tm = (tl + tr) / 2;
    buildtree(a, v*2, tl, tm,t);
    buildtree(a, v*2+1, tm+1, tr,t);
    t[v] = t[v*2] + t[v*2+1];
    }
}
ll rangedsum(ll v, ll tl, ll tr, ll l, ll r,vi& t)
{
    if (l > r) return 0;
    if (l == tl && r == tr) {return t[v];}
    ll tm = (tl + tr) / 2;
    return rangedsum(v*2, tl, tm, l, min(r, tm),t)+ rangedsum(v*2+1, tm+1, tr, max(l, tm+1), r,t);
}
void update(ll v, ll tl, ll tr, ll pos, ll new_val,vi& t) 
{
    if (tl == tr) {t[v] = new_val;} 
    else {ll tm = (tl + tr) / 2;
    if (pos <= tm) update(v*2, tl, tm, pos, new_val,t);
    else update(v*2+1, tm+1, tr, pos, new_val,t);
    t[v] = t[v*2] + t[v*2+1];
    }
}
// Code By Arindam Dey
// Counting number of nodes in a subtree
// map<ll,ll> mp;
// ll dfs(ll node,ll par,vector<vector<ll>> &adj)
// {
//     ll c=0;
//     for(auto &x:adj[node])
//     {
//         if(x!=par)
//         {
//             c+=dfs(x,node,adj);
//         }
//     }
//     mp[node]=c+1;
//     return c+1;
// }


// Level Order Traversal On a Tree
// while(!q.empty())
// {
//   ll sz=q.size();
//   // ++level;
//   for(ll i=0;i<sz;i++)
//   {
//       ll top=q.front().first;
//       ll par=q.front().second;
//       q.pop();      
//       for(auto &x:adj[top])
//       {
//          if(x==par)continue;
//           q.push({x,top});
//       }
//   }
//   level++;
// }

ll findParent(ll node,vector<ll>& par)
{
    if(par[node]==node)return node;
    return par[node]=findParent(par[node],par);
}
void Union(ll u,ll v,vector<ll> &size,vector<ll> &par)
{
    ll pu=findParent(u,par);
    ll pv=findParent(v,par);
    if(pu==pv)return;
    if(size[pu]<=size[pv])
    {
        par[pu]=pv;
        size[pv]++;
    }
    else
    {
        par[pv]=pu;
        size[pu]++;
    }
}
// Code By Arindam Dey
//DP
ll KnapSack(vi val,vi wt, ll cap)
{
    ll n=val.size();
    vi dp(cap+1,0);
    for(ll i=0; i < n; i++)for(ll j=cap; j>=wt[i]; j--)
    dp[j] = max(dp[j] , val[i] + dp[j-wt[i]]);
    return dp[cap];
}

bool comp(pair<ll,ll> a,pair<ll,ll> b)
{
    if(a.ff==b.ff)
        return a.ss<b.ss;
    return a.ff<b.ss;
}
ll dfs(ll f,ll u,vvi &adj,vi &vis,multiset<pair<ll,ll>>& mst,map<ll,ll> &fr,ll parent)
{
    //cout<<u<<endl;
    ll flag=0;
    if(fr[u]==f && parent!=-1)
    {
        flag=1;
    }
    for(ll i=0;i<adj[u].size();i++)
    {
        if(vis[adj[u][i]]==0)
        {
            vis[adj[u][i]]=1;
            if(dfs(f,adj[u][i],adj,vis,mst,fr,u)==1)
            {
                flag=1;
                if(mst.find({u,adj[u][i]})!=mst.end())
                    mst.erase(mst.find({u,adj[u][i]}));
                else if(mst.find({adj[u][i],u})!=mst.end())
                    mst.erase(mst.find({adj[u][i],u}));         }
        }
    }
    return flag;
}
// Code By Arindam Dey
bool isValid(ll &x,ll &y,ll &n,ll &m)
{
    return (x>=0 && x<n && y>=0 && y<m);
}
long long binpow(long long a, long long b) {
    long long res = 1;
    while (b > 0) {
        if (b & 1)
            res = res * a;
        a = a * a;
        b >>= 1;
    }
    return res;
}

ll maxi=1e6;
vector<bool> prime(maxi,true);
void sieve(ll n)
{
    prime[0]=prime[1]=false;
    for (ll p = 2; p * p <= n; p++) {
        if (prime[p] == true) {
            for (ll i = p * p; i <= n; i += p)
                prime[i] = false;
        }
    }
}
// Code By Arindam Dey
int main() {
  
  #ifndef ONLINE_JUDGE
  freopen("input.txt", "r", stdin);
  freopen("output.txt", "w", stdout);
  #endif
  
  IOS;
 
  ll t=1;
  // cin>>t;
 
  while(t--)
  {
      ll n;
      cin>>n;
      vector<ll> v(n);
      ll sum=0;
      for(ll i=0;i<n;i++)
      {
        cin>>v[i];
        sum+=v[i];
      }
      if(sum%3!=0||n<=2)
      {
        cout<<0<<endl;
      }
      else
      {
        ll z=sum/3;
        vector<ll> pre(n),suf(n);
        ll curr=0;
        for(ll i=n-1;i>=0;i--)
        {
            curr+=v[i];
            if(curr==z)
            {
                suf[i]=1;
            }
        }
        for(ll i=n-2;i>=0;i--)
        {
            suf[i]+=suf[i+1];
        }
        curr=0;
        ll ans=0;
        for(ll i=0;i<=n-3;i++)
        {
            curr+=v[i];
            if(curr==z)
            {
                if(z==0)
                ans+=suf[i+2];
                else
                ans+=suf[i+1];
            }
        }
        cout<<ans<<endl;
      }

  }
  return 0;

}

 