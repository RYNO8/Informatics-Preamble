//BEGIN PACKED FROM/mnt/d/programming/C++/Preamble/tests/cht.test.cpp
//solution to Commando
//https://dmoj.ca/problem/apio10p1
//BEGIN PACKED FROM/mnt/d/programming/C++/Preamble/tests/../src/Constants.h
/*************************************************IMPORTS*************************************************/
#include<bits/stdc++.h>
//#include<unordered_map>
//#include<unordered_set>
//#include<functional>
//#include<algorithm>
//#include<iostream>
//#include<assert.h>
//#include<iterator>
//#include<utility>
//#include<numeric>
//#include<cstddef>
//#include<fstream>
//#include<iomanip>
//#include<sstream>
//#include<climits>
//#include<bitset>
//#include<random>
//#include<chrono>
//#include<math.h>
//#include<time.h>
//#include<string>
//#include<vector>
//#include<limits>
//#include<queue>
//#include<stack>
//#include<set>
//#include<map>
#ifdef _MSC_VER
#include<intrin.h>
#define __builtin_popcount __popcnt
#define __builtin_popcountll __popcntll
#endif
namespace DS{
/*************************************************CONSTANTS*************************************************/
typedef unsigned int uint;typedef long double ld;constexpr ld PI=3.1415926535897932384626433832795;constexpr ld E=2.7182818284590452353602874713527;constexpr ld PHI=1.6180339887498948482045868343656;constexpr ld GAMMA=0.5772156649015328606065120900824;typedef long long ll;constexpr ll SIXTYNINE=69;constexpr ll MAXR=1000;constexpr ll MAXC=1000;constexpr ll MAXN=1<<19;constexpr ll LOG_MAXN=19;constexpr ll SQRT_MAXN=724; 
//or 300
constexpr ll MOD=1000000007;constexpr ll MOD_POW=1; 
//highest value such that MOD-1 is divisible by 2^MOD_POW
typedef unsigned long long ull;std::vector<std::pair<int,int>>DIRS_RECTILINEAR={{0,1},{1,0},{0,-1},{-1,0}};std::vector<std::pair<int,int>>DIRS_DIAG={{1,1},{1,-1},{-1,1},{-1,-1}};std::vector<std::pair<int,int>>DIRS_ALL={{0,1},{1,0},{0,-1},{-1,0},{1,1},{1,-1},{-1,1},{-1,-1}};std::vector<std::pair<int,int>>DIRS=DIRS_RECTILINEAR;std::random_device rd;std::mt19937_64 rng(rd());std::uniform_int_distribution<int>int_dis(0,std::numeric_limits<int>::max());std::uniform_int_distribution<uint>uint_dis(0,std::numeric_limits<uint>::max());std::uniform_int_distribution<ll>ll_dis(0,std::numeric_limits<ll>::max());std::uniform_int_distribution<ull>ull_dis(0,std::numeric_limits<ull>::max());std::uniform_real_distribution<ld>prob_dist(0,1);
/*************************************************DISPLAY*************************************************/
//Displays a pair:
//(A,B)
template<typename A,typename B>std::ostream&operator<<(std::ostream&out,const std::pair<A,B>&c){out<< 
"("
 <<c.first<< 
","
 <<c.second<< 
")"
;return out;}
//Displays a vector:
//[T,T,...,T]
template<typename T>std::ostream&operator<<(std::ostream&out,const std::vector<T>&v){out<< 
"["
;for(int i=0;i<(int)v.size();++i){out<<v[i];if(i!=v.size()-1)out<< 
","
;}out<< 
"]"
;return out;}
//Displays an array:
//T T T...T
#define printArr0(a,N)for(int i=0;i<N;++i)std::cout<<a[i]<<" ";
#define printArr1(a,N)for(int i=1;i<=N;++i)std::cout<<a[i]<<" ";
#define printArr(a,N)printArr0(a,N);
//Displays a queue:
//<T,T,...,T>
template<typename T>std::ostream&operator<<(std::ostream&out,const std::queue<T>&q){out<< 
"<"
;for(int i=0;i<(int)q.size();++i){out<<q.front();q.push(q.front());q.pop();if(i!=q.size()-1)out<< 
","
;}out<< 
">"
;return out;}
//Displays a stack:
//<T,T,...,T>
template<typename T>std::ostream&operator<<(std::ostream&out,const std::stack<T>&s){out<< 
"<"
;for(int i=0;i<(int)s.size();++i){out<<s.top();s.push(s.top());s.pop();if(i!=s.size()-1)out<< 
","
;}out<< 
">"
;return out;}
//Displays a deque:
//<T,T,...,T>
template<typename T>std::ostream&operator<<(std::ostream&out,const std::deque<T>&q){out<< 
"<"
;for(int i=0;i<(int)q.size();++i){out<<q.front();q.push_back(q.front());q.pop_front();if(i!=q.size()-1)out<< 
","
;}out<< 
">"
;return out;}
//Displays a map:
//{A:B,A:B,...,A:B}
template<typename A,typename B>std::ostream&operator<<(std::ostream&out,const std::map<A,B>&m){out<< 
"{"
;for(auto it=m.begin();it!=m.end();){out<<it->first<< 
':'
 <<it->second;if(++it!=m.end())out<< 
","
;}out<< 
"}"
;return out;}
//Displays an unordered map:
//{A:B,A:B,...,A:B}
template<typename A,typename B>std::ostream&operator<<(std::ostream&out,const std::unordered_map<A,B>&m){out<< 
"{"
;for(auto it=m.begin();it!=m.end();){out<<it->first<< 
':'
 <<it->second;if(++it!=m.end())out<< 
","
;}out<< 
"}"
;return out;}
//Displays a set:
//{T,T,...,T}
template<typename T>std::ostream&operator<<(std::ostream&out,const std::set<T>&m){out<< 
"{"
;for(auto it=m.begin();it!=m.end();){out<<*it;if(++it!=m.end())out<< 
","
;}out<< 
"}"
;return out;}
//Displays an unordered set:
//{T,T,...,T}
template<typename T>std::ostream&operator<<(std::ostream&out,const std::unordered_set<T>&m){out<< 
"{"
;for(auto it=m.begin();it!=m.end();){out<<*it;if(++it!=m.end())out<< 
","
;}out<< 
"}"
;return out;}
//Displays a multiset:
//{T,T,...,T}
template<typename T>std::ostream&operator<<(std::ostream&out,const std::multiset<T>&m){out<< 
"{"
;for(auto it=m.begin();it!=m.end();){out<<*it;if(++it!=m.end())out<< 
","
;}out<< 
"}"
;return out;}
/*************************************************UTILITIES*************************************************/
//O(log a+log b)
//@note`a`and`b`should be non-negative
template<typename T>T gcd(T a,T b){if(b==T(0))return a;else return gcd(b,a%b);}
//O(log a+log b)
template<typename T>T gcdSlow(T x,T y){while(y!=0){T temp=y;y=x%temp;x=temp;}return x;}
//O(log a+log b)
//@note May overflow integer
template<typename T>T lcm(T a,T b){return(a/gcd(a,b))*b;}
//O(12)
//@returns the number of set bits in the binary represetnation of`x`
template<typename T>T popcount(T x){x-=(x>>1)&0x5555555555555555;x=(x&0x3333333333333333)+((x>>2)&0x3333333333333333);x=(x+(x>>4))&0x0f0f0f0f0f0f0f0f;return(x*0x0101010101010101)>>56;}
//O(1)
//@returns Whether the number is negative(-1),zero(0)or positive(+1)
template<typename T>T sgn(T x){return T(x>0)-T(x<0);}
//O(1)
//Sets`a`to the minimum of`a`and`b`
//@note`a`is provided by reference
template<typename T>void pMin(T&a,T b){if(b<a)a=b;}
//O(1)
//Sets`a`to the maximum of`a`and`b`
//@note`a`is provided by reference
template<typename T>void pMax(T&a,T b){if(b>a)a=b;}
//O(N)
//@returns The number of characters of the printed representation of`obj`
template<typename T>int reprLen(T obj){std::stringstream s;s<<obj;return s.str().size();}};
//END PACKED FROM/mnt/d/programming/C++/Preamble/tests/../src/Constants.h
//BEGIN PACKED FROM/mnt/d/programming/C++/Preamble/tests/../src/CHT.h
namespace DS{ 
//Representation of an infinite line using intersection gradient form
 struct CHTLine{ll m,b; 
//O(1)
//Displays the line
//@param`out`The string representation of the graph is piped to this output stream
//@param`newLine`Indicates whether to end with a trailing`\\n`
 void print(std::ostream&out=std::cout,bool newLine=true)const{out<< 
"y="
 <<m<< 
" x+"
 <<b;if(newLine)out<< 
'\n'
;} 
//O(1)
//@returns the x value of the intersection
 double intersect(CHTLine a)const{return(double)(b-a.b)/(a.m-m);} 
//O(1)
//@returns Evalute the y value when x=x
//@note May overflow
 ll eval(ll x)const{return m*x+b;}};std::ostream&operator<<(std::ostream&out,CHTLine line){line.print(out);return out;} 
//Convex Hull Trick
//"infinte" range and domain,but added lines need to have non-decreasing gradient
 class CHT:public std::vector<CHTLine>{public: 
//amortised O(1)
//Include a unique line into the convex hull,where`m`is non-decreasing
 void addLine(CHTLine l){while(size()>=2&&back().intersect(*++rbegin())>=back().intersect(l)){pop_back();}while(!empty()&&back().m==l.m&&back().b>=l.b){pop_back();}push_back(l);} 
//O(log n)
//@returns Minimum y value when x=`x`
 ll getMinima(ll x)const{if(empty())return std::numeric_limits<ll>::max();int l=-1,r=size()-1;while(l+1<r){int m=(l+r)/2;if((*this)[m].intersect((*this)[m+1])>=x)r=m;else l=m;}return(*this)[r].eval(x);}}; 
//Li Chao tree
//Lines can be inserted and queried in any order,but range is limited to`MAXN`
//Initially contains the line 0x+0
 class LiChaoTree{private:CHTLine tree[2*MAXN];public: 
//O(1)
//Initialisation
 LiChaoTree(){} 
//O(log MAXN)
//Inserts a line
 void addLine(CHTLine line){int node=1,l=0,r=MAXN-1;do{int m=(l+r)/2;bool lGood=line.eval(l)>=tree[node].eval(l);bool mGood=line.eval(m)>=tree[node].eval(m);if(mGood)tree[node]=line;if(lGood!=mGood){node=node*2;r=m;}else{node=node*2+1;l=m+1;}}while(l!=r);} 
//O(log MAXN)
//@returns Maximum y value when x=`x`
 ll getLine(int x){ 
//implemented iteratively
 int node=1,l=0,r=MAXN-1;ll ans=LLONG_MIN;assert(l<=x&&x<=r&& 
"x out of range"
);while(true){ans=std::max(ans,tree[node].eval(x));int m=(l+r)/2;if(l==r)return ans;else if(x<m){node=node*2;r=m;}else{node=node*2+1;l=m+1;}}}};};
//END PACKED FROM/mnt/d/programming/C++/Preamble/tests/../src/CHT.h
using namespace std;using namespace DS;ll N,A,B,C;ll arr[1000005],prefix[1000005],dp[1000005];int main(){cin.tie(0);ios::sync_with_stdio(0); 
//ifstream cin{"in.txt"};
 cin>>N>>A>>B>>C;for(ll i=1;i<=N;++i)cin>>arr[i];for(ll i=1;i<=N;++i)prefix[i]=prefix[i-1]+arr[i];CHT cht;dp[0]=0;for(ll i=0;i<=N;++i){if(i==0)dp[i]=0;else dp[i]=cht.getMinima(prefix[i])+A*prefix[i]*prefix[i]+B*prefix[i]+C; 
//Ax^2+Bx;
 CHTLine l={.m=-2*A*prefix[i],.b=A*prefix[i]*prefix[i]-B*prefix[i]+dp[i]};cht.addLine(l);}cout<<dp[N]<< 
"\n"
;}
//END PACKED FROM/mnt/d/programming/C++/Preamble/tests/cht.test.cpp
