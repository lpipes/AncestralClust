#include "math.h"

double LnGamma (double alpha)
{
  /* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
     Stirling's formula is used for the central polynomial part of the procedure.
     Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
     Communications of the Association for Computing Machinery, 9:684
     */
  double x=alpha, f=0, z;

  if (x<7) {
    f=1;  z=x-1;
    while (++z<7)  f*=z;
    x=z;   f=-log(f);
  }
  z = 1/(x*x);
  return  f + (x-0.5)*log(x) - x + .918938533204673
    + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
        +.083333333333333)/x;
}

double IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
  /* returns the incomplete gamma ratio I(x,alpha) where x is the upper 
     limit of the integration and alpha is the shape parameter.
     returns (-1) if in error
     ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
     (1) series expansion     if (alpha>x || x<=1)
     (2) continued fraction   otherwise
     RATNEST FORTRAN by
     Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
19: 285-287 (AS32)
*/
  int i;
  double p=alpha, g=ln_gamma_alpha;
  double accurate=1e-20, overflow=1e30;
  double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

  if (x==0) return (0);
  if (x<0 || p<=0) return (-1);

  factor=exp(p*log(x)-x-g);   
  if (x>1 && x>=p) goto l30;
  /* (1) series expansion */
  gin=1;  term=1;  rn=p;
l20:
  rn++;
  term*=x/rn;   gin+=term;

  if (term > accurate) goto l20;
  gin*=factor/p;
  goto l50;
l30:
  /* (2) continued fraction */
  a=1-p;   b=a+x+1;  term=0;
  pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
  gin=pn[2]/pn[3];
l32:
  a++;  b+=2;  term++;   an=a*term;
  for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
  if (pn[5] == 0) goto l35;
  rn=pn[4]/pn[5];   dif=fabs(gin-rn);
  if (dif>accurate) goto l34;
  if (dif<=accurate*rn) goto l42;
l34:
  gin=rn;
l35:
  for (i=0; i<4; i++) pn[i]=pn[i+2];
  if (fabs(pn[4]) < overflow) goto l32;
  for (i=0; i<4; i++) pn[i]/=overflow;
  goto l32;
l42:
  gin=1-factor*gin;

l50:
  /*printf("Incompletegamma got %f %f %f and returned %f\n",x,  alpha,  ln_gamma_alpha,gin);*/
  printf("");
  return (gin);
}

double PointNormal (double prob)
{
  /* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
     returns (-9999) if in error
     Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
     Applied Statistics 22: 96-97 (AS70)

     Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
     normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
     points of the normal distribution.  26: 118-121.

*/
  double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
  double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
  double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
  double y, z=0, p=prob, p1;

  p1 = (p<0.5 ? p : 1-p);
  if (p1<1e-20) return (-9999);

  y = sqrt (log(1/(p1*p1)));
  z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
  return (p<0.5 ? -z : z);
}



double PointChi2 (double prob, double v)
{
  /* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
     returns -1 if in error.   0.000002<prob<0.999998
     RATNEST FORTRAN by
     Best DJ & Roberts DE (1975) The percentage points of the
     Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
     Converted into C by Ziheng Yang, Oct. 1993.
     */
  double e=.5e-6, aa=.6931471805, p=prob, g;
  double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

  if (p<.000002 || p>.999998 || v<=0) return (-1);

  g = LnGamma (v/2);
  xx=v/2;   c=xx-1;
  if (v >= -1.24*log(p)) goto l1;

  ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
  if (ch-e<0) return (ch);
  goto l4;
l1:
  if (v>.32) goto l3;
  ch=0.4;   a=log(1-p);
l2:
  q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
  t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
  ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
  if (fabs(q/ch-1)-.01 <= 0) goto l4;
  else                       goto l2;

l3:
  x=PointNormal (p);
  p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
  if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
  q=ch;   p1=.5*ch;
  if ((t=IncompleteGamma (p1, xx, g))<0) {
    printf ("\nerr IncompleteGamma");
    return (-1);
  }
  p2=p-t;
  t=p2*exp(xx*aa+g+p1-c*log(ch));
  b=t/ch;  a=0.5*t-b*c;

  s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
  s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
  s3=(210+a*(462+a*(707+932*a)))/2520;
  s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
  s5=(84+264*a+c*(175+606*a))/2520;
  s6=(120+c*(346+127*c))/5040;
  ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
  if (fabs(q/ch-1) > e) goto l4;

  return (ch);
}



#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))


double CDFfunGamma(double x, double par[2])

{
  return IncompleteGamma(par[1]*x,par[0],LnGamma(par[0]));
}

void definegammaquantiles(int k, double par[2])
{
  int i;
  double mean=0.0;

  // printf("Theta values (alpha = %lf, beta = %lf:",par[0],par[1]);
  for (i=0; i<k; i++){
    statevector[i] = PointGamma(((double)i+0.5)/(double)k, par[0],par[1]);
    //mean+=statevector[i];
  }
}

void initlogfactorial()
{
  int i;

  Logfactorial[0]=0.0; 
  Logfactorial[1]=0.0;
  for (i=2; i<MAXNUMBEROFINDINSPECIES+1; i++)
    Logfactorial[i] = Logfactorial[i-1] + log((double)i);
}




/***********************************************************
 *  This eigen() works for eigenvalue/vector analysis
 *         for real general square matrix A
 *         A will be destroyed
 *         rr,ri are vectors containing eigenvalues
 *         vr,vi are matrices containing (right) eigenvectors
 *
 *              A*[vr+vi*i] = [vr+vi*i] * diag{rr+ri*i}
 *
 *  Algorithm: Handbook for Automatic Computation, vol 2
 *             by Wilkinson and Reinsch, 1971
 *             most of source codes were taken from a public domain
 *             solftware called MATCALC.
 *  Credits:   to the authors of MATCALC
 *
 *  return     -1 not converged
 *              0 no complex eigenvalues/vectors
 *              1 complex eigenvalues/vectors
 *  Tianlin Wang at University of Illinois
 *  Thu May  6 15:22:31 CDT 1993
 ***************************************************************/

#define BASE        2    /* base of floating point arithmetic */
#define DIGITS     40    /* no. of digits to the base BASE in the fraction */
#define MAXITER    30    /* max. no. of iterations to converge */

#define pos(i,j,n)      ((i)*(n)+(j))

/*A is the matrix, rr = root real(nx1), ri = root imaginary(nx1), vr = real part of eigenvector(nxn), w = working space (size of 2n), set job to 1 (decides if both eigenvectors and eigenvalues should be calculated*/

int eigen(int job, double A[], int n, double rr[], double ri[],
    double vr[], double vi[], double w[]);
void balance(double mat[], int n, int *low, int *hi, double scale[]);
void unbalance(int n, double vr[], double vi[], int low, int hi,
    double scale[]);
int realeig(int job, double mat[], int n,int low, int hi, double valr[],
    double vali[], double vr[], double vi[]);
void elemhess(int job, double mat[], int n, int low, int hi,
    double vr[], double vi[], int work[]);

int eigen(int job, double A[], int n, double rr[], double ri[],
    double vr[], double vi[], double work[])
{
  /*  double work[n*2]: working space
  */
  int low,hi,i,j,k, it, istate=0;
  double tiny=sqrt(pow((double)BASE,(double)(1-DIGITS))), t;

  balance(A,n,&low,&hi,work);
  elemhess(job,A,n,low,hi,vr,vi, (int*)(work+n));
  if (-1 == realeig(job,A,n,low,hi,rr,ri,vr,vi)) return (-1);
  if (job) unbalance(n,vr,vi,low,hi,work);
  /* sort, added by Z. Yang */
  for (i=0; i<n; i++) {
    for (j=i+1,it=i,t=rr[i]; j<n; j++)
      if (t<rr[j]) { t=rr[j]; it=j; }
    rr[it]=rr[i];   rr[i]=t;
    t=ri[it];       ri[it]=ri[i];  ri[i]=t;
    for (k=0; k<n; k++) {
      t=vr[k*n+it];
	  vr[k*n+it]=vr[k*n+i];
	  vr[k*n+i]=t;

      t=vi[k*n+it];
	  vi[k*n+it]=vi[k*n+i];
	  vi[k*n+i]=t;
    }
    if (fabs(ri[i])>tiny) istate=1;
  }

  return (istate) ;
}

/* complex funcctions
*/

complex compl (double re,double im)
{
  complex r;

  r.re = re;
  r.im = im;
  return(r);
}

complex conj (complex a)
{
  a.im = -a.im;
  return(a);
}

#define csize(a) (fabs(a.re)+fabs(a.im))

complex cplus (complex a, complex b)
{
  complex c;
  c.re = a.re+b.re;
  c.im = a.im+b.im;
  return (c);
}

complex cminus (complex a, complex b)
{
  complex c;
  c.re = a.re-b.re;
  c.im = a.im-b.im;
  return (c);
}

complex cby (complex a, complex b)
{
  complex c;
  c.re = a.re*b.re-a.im*b.im ;
  c.im = a.re*b.im+a.im*b.re ;
  return (c);
}

complex cdiv (complex a,complex b)
{
  double ratio, den;
  complex c;

  if (fabs(b.re) <= fabs(b.im)) {
    ratio = b.re / b.im;
    den = b.im * (1 + ratio * ratio);
    c.re = (a.re * ratio + a.im) / den;
    c.im = (a.im * ratio - a.re) / den;
  }
  else {
    ratio = b.im / b.re;
    den = b.re * (1 + ratio * ratio);
    c.re = (a.re + a.im * ratio) / den;
    c.im = (a.im - a.re * ratio) / den;
  }
  return(c);
}

complex cexp (complex a)
{
  complex c;
  c.re = exp(a.re);
  if (fabs(a.im)==0) c.im = 0;
  else  { c.im = c.re*sin(a.im); c.re*=cos(a.im); }
  return (c);
}

complex cfactor (complex x, double a)
{
  complex c;
  c.re = a*x.re;
  c.im = a*x.im;
  return (c);
}

int cxtoy (complex x[], complex y[], int n)
{
  int i;
  FOR (i,n) y[i]=x[i];
  return (0);
}

int cmatby (complex a[], complex b[], complex c[], int n,int m,int k)
  /* a[n*m], b[m*k], c[n*k]  ......  c = a*b
  */
{
  int i,j,i1;
  complex t;

  FOR (i,n)  FOR(j,k) {
    for (i1=0,t=compl(0,0); i1<m; i1++)
      t = cplus (t, cby(a[i*m+i1],b[i1*k+j]));
    c[i*k+j] = t;
  }
  return (0);
}

int cmatout (FILE * fout, complex x[], int n, int m)
{
  int i,j;
  for (i=0,FPN(fout); i<n; i++,FPN(fout))
    FOR(j,m) fprintf(fout, "%7.3f%7.3f  ", x[i*m+j].re, x[i*m+j].im);
  return (0);
}

int cmatinv( complex x[], int n, int m, double space[])
{
  /* x[n*m]  ... m>=n
  */
  int i,j,k, *irow=(int*) space;
  double xmaxsize, ee=1e-20;
  complex xmax, t,t1;

  FOR(i,n)  {
    xmaxsize = 0.;
    for (j=i; j<n; j++) {
      if ( xmaxsize < csize (x[j*m+i]))  {
        xmaxsize = csize (x[j*m+i]);
        xmax = x[j*m+i];
        irow[i] = j;
      }
    }
    if (xmaxsize < ee)   {
      printf("\nDet goes to zero at %8d!\t\n", i+1);
      return(-1);
    }
    if (irow[i] != i) {
      FOR(j,m) {
        t = x[i*m+j];
        x[i*m+j] = x[irow[i]*m+j];
        x[ irow[i]*m+j] = t;
      }
    }
    t = cdiv (compl(1,0), x[i*m+i]);
    FOR(j,n) {
      if (j == i) continue;
      t1 = cby (t,x[j*m+i]);
      FOR(k,m)  x[j*m+k] = cminus (x[j*m+k], cby(t1,x[i*m+k]));
      x[j*m+i] = cfactor (t1, -1);
    }
    FOR(j,m)   x[i*m+j] = cby (x[i*m+j], t);
    x[i*m+i] = t;
  }
  for (i=n-1; i>=0; i--) {
    if (irow[i] == i) continue;
    FOR(j,n)  {
      t = x[j*m+i];
      x[j*m+i] = x[j*m+irow[i]];
      x[ j*m+irow[i]] = t;
    }
  }
  return (0);
}


void balance(double mat[], int n,int *low, int *hi, double scale[])
{
  /* Balance a matrix for calculation of eigenvalues and eigenvectors
  */
  double c,f, g, r,s;
  int i,j,k,l,done;
  /* search for rows isolating an eigenvalue and push them down */
  for (k = n - 1; k >= 0; k--) {
    for (j = k; j >= 0; j--) {
      for (i = 0; i <= k; i++) {
        if (i != j && fabs(mat[pos(j,i,n)]) != 0) break;
      }

      if (i > k) {
        scale[k] = j;

        if (j != k) {
          for (i = 0; i <= k; i++) {
            c = mat[pos(i,j,n)];
            mat[pos(i,j,n)] = mat[pos(i,k,n)];
            mat[pos(i,k,n)] = c;
          }

          for (i = 0; i < n; i++) {
            c = mat[pos(j,i,n)];
            mat[pos(j,i,n)] = mat[pos(k,i,n)];
            mat[pos(k,i,n)] = c;
          }
        }
        break;
      }
    }
    if (j < 0) break;
  }

  /* search for columns isolating an eigenvalue and push them left */

  for (l = 0; l <= k; l++) {
    for (j = l; j <= k; j++) {
      for (i = l; i <= k; i++) {
        if (i != j && fabs(mat[pos(i,j,n)]) != 0) break;
      }
      if (i > k) {
        scale[l] = j;
        if (j != l) {
          for (i = 0; i <= k; i++) {
            c = mat[pos(i,j,n)];
            mat[pos(i,j,n)] = mat[pos(i,l,n)];
            mat[pos(i,l,n)] = c;
          }

          for (i = l; i < n; i++) {
            c = mat[pos(j,i,n)];
            mat[pos(j,i,n)] = mat[pos(l,i,n)];
            mat[pos(l,i,n)] = c;
          }
        }

        break;
      }
    }

    if (j > k) break;
  }

  *hi = k;
  *low = l;

  /* balance the submatrix in rows l through k */

  for (i = l; i <= k; i++) {
    scale[i] = 1;
  }

  do {
    for (done = 1,i = l; i <= k; i++) {
      for (c = 0,r = 0,j = l; j <= k; j++) {
        if (j != i) {
          c += fabs(mat[pos(j,i,n)]);
          r += fabs(mat[pos(i,j,n)]);
        }
      }

      if (c != 0 && r != 0) {
        g = r / BASE;
        f = 1;
        s = c + r;

        while (c < g) {
          f *= BASE;
          c *= BASE * BASE;
        }

        g = r * BASE;

        while (c >= g) {
          f /= BASE;
          c /= BASE * BASE;
        }

        if ((c + r) / f < 0.95 * s) {
          done = 0;
          g = 1 / f;
          scale[i] *= f;

          for (j = l; j < n; j++) {
            mat[pos(i,j,n)] *= g;
          }

          for (j = 0; j <= k; j++) {
            mat[pos(j,i,n)] *= f;
          }
        }
      }
    }
  } while (!done);
}


/*
 * Transform back eigenvectors of a balanced matrix
 * into the eigenvectors of the original matrix
 */
void unbalance(int n,double vr[],double vi[], int low, int hi, double scale[])
{
  int i,j,k;
  double tmp;

  for (i = low; i <= hi; i++) {
    for (j = 0; j < n; j++) {
      vr[pos(i,j,n)] *= scale[i];
      vi[pos(i,j,n)] *= scale[i];
    }
  }

  for (i = low - 1; i >= 0; i--) {
    if ((k = (int)scale[i]) != i) {
      for (j = 0; j < n; j++) {
        tmp = vr[pos(i,j,n)];
        vr[pos(i,j,n)] = vr[pos(k,j,n)];
        vr[pos(k,j,n)] = tmp;

        tmp = vi[pos(i,j,n)];
        vi[pos(i,j,n)] = vi[pos(k,j,n)];
        vi[pos(k,j,n)] = tmp;
      }
    }
  }

  for (i = hi + 1; i < n; i++) {
    if ((k = (int)scale[i]) != i) {
      for (j = 0; j < n; j++) {
        tmp = vr[pos(i,j,n)];
        vr[pos(i,j,n)] = vr[pos(k,j,n)];
        vr[pos(k,j,n)] = tmp;

        tmp = vi[pos(i,j,n)];
        vi[pos(i,j,n)] = vi[pos(k,j,n)];
        vi[pos(k,j,n)] = tmp;
      }
    }
  }
}

/*
 * Reduce the submatrix in rows and columns low through hi of real matrix mat to
 * Hessenberg form by elementary similarity transformations
 */
void elemhess(int job,double mat[],int n,int low,int hi, double vr[],
    double vi[], int work[])
{
  /* work[n] */
  int i,j,m;
  double x,y;

  for (m = low + 1; m < hi; m++) {
    for (x = 0,i = m,j = m; j <= hi; j++) {
      if (fabs(mat[pos(j,m-1,n)]) > fabs(x)) {
        x = mat[pos(j,m-1,n)];
        i = j;
      }
    }

    if ((work[m] = i) != m) {
      for (j = m - 1; j < n; j++) {
        y = mat[pos(i,j,n)];
        mat[pos(i,j,n)] = mat[pos(m,j,n)];
        mat[pos(m,j,n)] = y;
      }

      for (j = 0; j <= hi; j++) {
        y = mat[pos(j,i,n)];
        mat[pos(j,i,n)] = mat[pos(j,m,n)];
        mat[pos(j,m,n)] = y;
      }
    }

    if (x != 0) {
      for (i = m + 1; i <= hi; i++) {
        if ((y = mat[pos(i,m-1,n)]) != 0) {
          y = mat[pos(i,m-1,n)] = y / x;

          for (j = m; j < n; j++) {
            mat[pos(i,j,n)] -= y * mat[pos(m,j,n)];
          }

          for (j = 0; j <= hi; j++) {
            mat[pos(j,m,n)] += y * mat[pos(j,i,n)];
          }
        }
      }
    }
  }
  if (job) {
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
        vr[pos(i,j,n)] = 0.0; vi[pos(i,j,n)] = 0.0;
      }
      vr[pos(i,i,n)] = 1.0;
    }

    for (m = hi - 1; m > low; m--) {
      for (i = m + 1; i <= hi; i++) {
        vr[pos(i,m,n)] = mat[pos(i,m-1,n)];
      }

      if ((i = work[m]) != m) {
        for (j = m; j <= hi; j++) {
          vr[pos(m,j,n)] = vr[pos(i,j,n)];
          vr[pos(i,j,n)] = 0.0;
        }
        vr[pos(i,m,n)] = 1.0;
      }
    }
  }
}

/*
 * Calculate eigenvalues and eigenvectors of a real upper Hessenberg matrix
 * Return 1 if converges successfully and 0 otherwise
 */

int realeig(int job,double mat[],int n,int low, int hi, double valr[],
    double vali[], double vr[],double vi[])
{
  complex v;
  double p=0,q=0,r=0,s=0,t,w,x,y,z=0,ra,sa,norm,eps;
  int niter,en,i,j,k,l,m;
  double precision  = pow((double)BASE,(double)(1-DIGITS));

  eps = precision;
  for (i=0; i<n; i++) {
    valr[i]=0.0;
    vali[i]=0.0;
  }
  /* store isolated roots and calculate norm */
  for (norm = 0,i = 0; i < n; i++) {
    for (j = max(0,i-1); j < n; j++) {
      norm += fabs(mat[pos(i,j,n)]);
    }
    if (i < low || i > hi) valr[i] = mat[pos(i,i,n)];
  }
  t = 0;
  en = hi;

  while (en >= low) {
    niter = 0;
    for (;;) {

      /* look for single small subdiagonal element */

      for (l = en; l > low; l--) {
        s = fabs(mat[pos(l-1,l-1,n)]) + fabs(mat[pos(l,l,n)]);
        if (s == 0) s = norm;
        if (fabs(mat[pos(l,l-1,n)]) <= eps * s) break;
      }

      /* form shift */

      x = mat[pos(en,en,n)];

      if (l == en) {             /* one root found */
        valr[en] = x + t;
        if (job) mat[pos(en,en,n)] = x + t;
        en--;
        break;
      }

      y = mat[pos(en-1,en-1,n)];
      w = mat[pos(en,en-1,n)] * mat[pos(en-1,en,n)];

      if (l == en - 1) {                /* two roots found */
        p = (y - x) / 2;
        q = p * p + w;
        z = sqrt(fabs(q));
        x += t;
        if (job) {
          mat[pos(en,en,n)] = x;
          mat[pos(en-1,en-1,n)] = y + t;
        }
        if (q < 0) {                /* complex pair */
          valr[en-1] = x+p;
          vali[en-1] = z;
          valr[en] = x+p;
          vali[en] = -z;
        }
        else {                      /* real pair */
          z = (p < 0) ? p - z : p + z;
          valr[en-1] = x + z;
          valr[en] = (z == 0) ? x + z : x - w / z;
          if (job) {
            x = mat[pos(en,en-1,n)];
            s = fabs(x) + fabs(z);
            p = x / s;
            q = z / s;
            r = sqrt(p*p+q*q);
            p /= r;
            q /= r;
            for (j = en - 1; j < n; j++) {
              z = mat[pos(en-1,j,n)];
              mat[pos(en-1,j,n)] = q * z + p *
                mat[pos(en,j,n)];
              mat[pos(en,j,n)] = q * mat[pos(en,j,n)] - p*z;
            }
            for (i = 0; i <= en; i++) {
              z = mat[pos(i,en-1,n)];
              mat[pos(i,en-1,n)] = q * z + p * mat[pos(i,en,n)];
              mat[pos(i,en,n)] = q * mat[pos(i,en,n)] - p*z;
            }
            for (i = low; i <= hi; i++) {
              z = vr[pos(i,en-1,n)];
              vr[pos(i,en-1,n)] = q*z + p*vr[pos(i,en,n)];
              vr[pos(i,en,n)] = q*vr[pos(i,en,n)] - p*z;
            }
          }
        }
        en -= 2;
        break;
      }
      if (niter == MAXITER) return(-1);
      if (niter != 0 && niter % 10 == 0) {
        t += x;
        for (i = low; i <= en; i++) mat[pos(i,i,n)] -= x;
        s = fabs(mat[pos(en,en-1,n)]) + fabs(mat[pos(en-1,en-2,n)]);
        x = y = 0.75 * s;
        w = -0.4375 * s * s;
      }
      niter++;
      /* look for two consecutive small subdiagonal elements */
      for (m = en - 2; m >= l; m--) {
        z = mat[pos(m,m,n)];
        r = x - z;
        s = y - z;
        p = (r * s - w) / mat[pos(m+1,m,n)] + mat[pos(m,m+1,n)];
        q = mat[pos(m+1,m+1,n)] - z - r - s;
        r = mat[pos(m+2,m+1,n)];
        s = fabs(p) + fabs(q) + fabs(r);
        p /= s;
        q /= s;
        r /= s;
        if (m == l || fabs(mat[pos(m,m-1,n)]) * (fabs(q)+fabs(r)) <=
            eps * (fabs(mat[pos(m-1,m-1,n)]) + fabs(z) +
              fabs(mat[pos(m+1,m+1,n)])) * fabs(p)) break;
      }
      for (i = m + 2; i <= en; i++) mat[pos(i,i-2,n)] = 0;
      for (i = m + 3; i <= en; i++) mat[pos(i,i-3,n)] = 0;
      /* double QR step involving rows l to en and columns m to en */
      for (k = m; k < en; k++) {
        if (k != m) {
          p = mat[pos(k,k-1,n)];
          q = mat[pos(k+1,k-1,n)];
          r = (k == en - 1) ? 0 : mat[pos(k+2,k-1,n)];
          if ((x = fabs(p) + fabs(q) + fabs(r)) == 0) continue;
          p /= x;
          q /= x;
          r /= x;
        }
        s = sqrt(p*p+q*q+r*r);
        if (p < 0) s = -s;
        if (k != m) {
          mat[pos(k,k-1,n)] = -s * x;
        }
        else if (l != m) {
          mat[pos(k,k-1,n)] = -mat[pos(k,k-1,n)];
        }
        p += s;
        x = p / s;
        y = q / s;
        z = r / s;
        q /= p;
        r /= p;
        /* row modification */
        for (j = k; j <= (!job ? en : n-1); j++){
          p = mat[pos(k,j,n)] + q * mat[pos(k+1,j,n)];
          if (k != en - 1) {
            p += r * mat[pos(k+2,j,n)];
            mat[pos(k+2,j,n)] -= p * z;
          }
          mat[pos(k+1,j,n)] -= p * y;
          mat[pos(k,j,n)] -= p * x;
        }
        j = min(en,k+3);
        /* column modification */
        for (i = (!job ? l : 0); i <= j; i++) {
          p = x * mat[pos(i,k,n)] + y * mat[pos(i,k+1,n)];
          if (k != en - 1) {
            p += z * mat[pos(i,k+2,n)];
            mat[pos(i,k+2,n)] -= p*r;
          }
          mat[pos(i,k+1,n)] -= p*q;
          mat[pos(i,k,n)] -= p;
        }
        if (job) {             /* accumulate transformations */
          for (i = low; i <= hi; i++) {
            p = x * vr[pos(i,k,n)] + y * vr[pos(i,k+1,n)];
            if (k != en - 1) {
              p += z * vr[pos(i,k+2,n)];
              vr[pos(i,k+2,n)] -= p*r;
            }
            vr[pos(i,k+1,n)] -= p*q;
            vr[pos(i,k,n)] -= p;
          }
        }
      }
    }
  }

  if (!job) return(0);
  if (norm != 0) {
    /* back substitute to find vectors of upper triangular form */
    for (en = n-1; en >= 0; en--) {
      p = valr[en];
      if ((q = vali[en]) < 0) {            /* complex vector */
        m = en - 1;
        if (fabs(mat[pos(en,en-1,n)]) > fabs(mat[pos(en-1,en,n)])) {
          mat[pos(en-1,en-1,n)] = q / mat[pos(en,en-1,n)];
          mat[pos(en-1,en,n)] = (p - mat[pos(en,en,n)]) /
            mat[pos(en,en-1,n)];
        }
        else {
          v = cdiv(compl(0.0,-mat[pos(en-1,en,n)]),
              compl(mat[pos(en-1,en-1,n)]-p,q));
          mat[pos(en-1,en-1,n)] = v.re;
          mat[pos(en-1,en,n)] = v.im;
        }
        mat[pos(en,en-1,n)] = 0;
        mat[pos(en,en,n)] = 1;
        for (i = en - 2; i >= 0; i--) {
          w = mat[pos(i,i,n)] - p;
          ra = 0;
          sa = mat[pos(i,en,n)];
          for (j = m; j < en; j++) {
            ra += mat[pos(i,j,n)] * mat[pos(j,en-1,n)];
            sa += mat[pos(i,j,n)] * mat[pos(j,en,n)];
          }
          if (vali[i] < 0) {
            z = w;
            r = ra;
            s = sa;
          }
          else {
            m = i;
            if (vali[i] == 0) {
              v = cdiv(compl(-ra,-sa),compl(w,q));
              mat[pos(i,en-1,n)] = v.re;
              mat[pos(i,en,n)] = v.im;
            }
            else {                      /* solve complex equations */
              x = mat[pos(i,i+1,n)];
              y = mat[pos(i+1,i,n)];
              v.re = (valr[i]- p)*(valr[i]-p) + vali[i]*vali[i] - q*q;
              v.im = (valr[i] - p)*2*q;
              if ((fabs(v.re) + fabs(v.im)) == 0) {
                v.re = eps * norm * (fabs(w) +
                    fabs(q) + fabs(x) + fabs(y) + fabs(z));
              }
              v = cdiv(compl(x*r-z*ra+q*sa,x*s-z*sa-q*ra),v);
              mat[pos(i,en-1,n)] = v.re;
              mat[pos(i,en,n)] = v.im;
              if (fabs(x) > fabs(z) + fabs(q)) {
                mat[pos(i+1,en-1,n)] =
                  (-ra - w * mat[pos(i,en-1,n)] +
                   q * mat[pos(i,en,n)]) / x;
                mat[pos(i+1,en,n)] = (-sa - w * mat[pos(i,en,n)] -
                    q * mat[pos(i,en-1,n)]) / x;
              }
              else {
                v = cdiv(compl(-r-y*mat[pos(i,en-1,n)],
                      -s-y*mat[pos(i,en,n)]),compl(z,q));
                mat[pos(i+1,en-1,n)] = v.re;
                mat[pos(i+1,en,n)] = v.im;
              }
            }
          }
        }
      }
      else if (q == 0) {                             /* real vector */
        m = en;
        mat[pos(en,en,n)] = 1;
        for (i = en - 1; i >= 0; i--) {
          w = mat[pos(i,i,n)] - p;
          r = mat[pos(i,en,n)];
          for (j = m; j < en; j++) {
            r += mat[pos(i,j,n)] * mat[pos(j,en,n)];
          }
          if (vali[i] < 0) {
            z = w;
            s = r;
          }
          else {
            m = i;
            if (vali[i] == 0) {
              if ((t = w) == 0) t = eps * norm;
              mat[pos(i,en,n)] = -r / t;
            }
            else {            /* solve real equations */
              x = mat[pos(i,i+1,n)];
              y = mat[pos(i+1,i,n)];
              q = (valr[i] - p) * (valr[i] - p) + vali[i]*vali[i];
              t = (x * s - z * r) / q;
              mat[pos(i,en,n)] = t;
              if (fabs(x) <= fabs(z)) {
                mat[pos(i+1,en,n)] = (-s - y * t) / z;
              }
              else {
                mat[pos(i+1,en,n)] = (-r - w * t) / x;
              }
            }
          }
        }
      }
    }
    /* vectors of isolated roots */
    for (i = 0; i < n; i++) {
      if (i < low || i > hi) {
        for (j = i; j < n; j++) {
          vr[pos(i,j,n)] = mat[pos(i,j,n)];
        }
      }
    }
    /* multiply by transformation matrix */

    for (j = n-1; j >= low; j--) {
      m = min(j,hi);
      for (i = low; i <= hi; i++) {
        for (z = 0,k = low; k <= m; k++) {
          z += vr[pos(i,k,n)] * mat[pos(k,j,n)];
        }
        vr[pos(i,j,n)] = z;
      }
    }
  }
  /* rearrange complex eigenvectors */
  for (j = 0; j < n; j++) {
    if (vali[j] != 0) {
      for (i = 0; i < n; i++) {
        vi[pos(i,j,n)] = vr[pos(i,j+1,n)];
        vr[pos(i,j+1,n)] = vr[pos(i,j,n)];
        vi[pos(i,j+1,n)] = -vi[pos(i,j,n)];
      }
      j++;
    }
  }
  return(0);
}

int matinv( double x[], int n, int m, double space[])
{
  /* x[n*m]  ... m>=n
  */
  register int i,j,k;
  int *irow=(int*) space;
  double ee=1.0e-20, t,t1,xmax;
  double det=1.0;

  FOR (i,n)  {
    xmax = 0.;
    for (j=i; j<n; j++) {
      if (xmax < fabs(x[j*m+i]))  {
        xmax = fabs( x[j*m+i] );
        irow[i] = j;
      }
    }
    det *= xmax;
    if (xmax < ee)   {
      printf("\nDet becomes zero at %3d!\t\n", i+1);
      return(-1);
    }
    if (irow[i] != i) {
      FOR (j,m) {
        t = x[i*m+j];
        x[i*m+j] = x[irow[i] * m + j];
        x[ irow[i] * m + j] = t;
      }
    }
    t = 1./x[i*m+i];
    FOR (j,n) {
      if (j == i) continue;
      t1 = t*x[j*m+i];
      FOR(k,m)  x[j*m+k] -= t1*x[i*m+k];
      x[j*m+i] = -t1;
    }
    FOR(j,m)   x[i*m+j] *= t;
    x[i*m+i] = t;
  }                            /* i  */
  for (i=n-1; i>=0; i--) {
    if (irow[i] == i) continue;
    FOR(j,n)  {
      t = x[j*m+i];
      x[j*m+i] = x[ j*m + irow[i] ];
      x[ j*m + irow[i] ] = t;
    }
  }
  return (0);
}



void inittransitionmatrix(double lambda)

{
  int i, j;
  double RIVAL[STATESPACE], RIVEC[STATESPACE][STATESPACE],  A[STATESPACE][STATESPACE], workspace[2*STATESPACE];

  for (i=0; i<STATESPACE; i++)
    for (j=0; j<STATESPACE; j++)
      A[i][j] = 0.0;
  A[0][0]=-lambda;
  A[0][1]=lambda;
  A[STATESPACE-1][STATESPACE-1]=-lambda;
  A[STATESPACE-1][STATESPACE-2]=lambda;
  for (i=1; i<STATESPACE-1; i++){
    A[i][i] = -2.0*lambda;
    A[i][i-1] = A[i][i+1] = lambda;
  }
  /*  for (i=0; i<STATESPACE; i++){
      for (j=0; j<STATESPACE; j++)
      printf("%lf ",A[i][j]);
      printf("\n");
      }*/
  if (eigen(1, A[0], STATESPACE, RRVAL, RIVAL, RRVEC[0], RIVEC[0], workspace) != 0)
  {
    printf("Transitions matrix did not converge or contained non-real values!\n");
    exit(-1);
  }
  for (i=0; i<STATESPACE; i++)
    for (j=0; j<STATESPACE; j++)
      LRVEC[i][j] = RRVEC[i][j];
  if (matinv(RRVEC[0],STATESPACE, STATESPACE, workspace) != 0)
    printf("Could not invert matrix!\nResults may not be reliable!\n");
}

void maketransitionmatrix(int matnum, double t)
{
  int i, j, k;
  double EXPOS[STATESPACE];

  for (k=0; k<STATESPACE; k++)
    EXPOS[k] = exp(t*RRVAL[k]);
  for (i=0; i<STATESPACE; i++)
  {
    for (j=0; j<STATESPACE; j++)
    {
      if (matnum==0) PMAT1[i][j] = 0.0;
      else PMAT2[i][j] = 0.0;
      for (k=0; k<STATESPACE; k++)
      {
        if (matnum==0) PMAT1[i][j] =  PMAT1[i][j] + RRVEC[k][j]*LRVEC[i][k]*EXPOS[k];
        else PMAT2[i][j] =  PMAT2[i][j] + RRVEC[k][j]*LRVEC[i][k]*EXPOS[k];
      }
    }
  }
  /*    if (matnum==1){
        printf("PMAT2 (t=%lf)\n",t);
        for (i=0; i<STATESPACE; i++){
        for (j=0; j<STATESPACE; j++)
        printf("%lf ",PMAT2[i][j]);
        printf("\n");} exit(-1);
        }*/
}

/*this function needs to be checked*/
double pdata(double theta, int s, int nt, int nb)
{
  int k, n;
  double p, totp=0.0;

  n = nt/nb; /*might want to do something more clever here*/
  for (k=1; k<n; k++)
  {
    p = exp(Logfactorial[n-2] - Logfactorial[k-1] - Logfactorial[n-2-k+1] + (double)(s+1)*log(theta/((double)k+theta)));
    if ((k/2)*k==k) totp = totp - p;
    else totp = totp + p;
  }
  //    printf("Pdata (%f): %lf\n",theta,totp*(double)(n-1)/theta);
  return totp*(double)(n-1)/theta;
}
