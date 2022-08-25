int suspect(long long b, int t, long long u, long long n)
{
long long prod=1;
while(u)
{
if(u&1) prod=((prod*b)%n);
b=((b*b)%n);
u/=2;
}
if(prod == 1) return 1;
for(int i = 1; i <= t; i++)
{
if(prod == n-1) return 1;
prod = (prod * prod) % n;
}
return 0;
}
int isprime(unsigned int n)
{
long long k = n - 1;
int t = 0;
while(!(k%2)) { t++; k/=2; }
if(n>2 && n%2==0) return 0;
if(n>3 && n%3==0) return 0;
if(n>5 && n%5==0) return 0;
if(n>7 && n%7==0) return 0;
if(suspect(61, t, k, n) && suspect(7, t, k, n) && suspect(2, t, k, n))
return 1;
return 0;
}