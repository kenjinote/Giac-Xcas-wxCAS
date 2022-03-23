// -*- mode:C++ ; compile-command: "g++-3.4 -I.. -g -c ifactor.cc -DHAVE_CONFIG_H -DIN_GIAC" -*-
#include "giacPCH.h"

#include "path.h"
/*
 *  Copyright (C) 2003,7 R. De Graeve & B. Parisse, 
 *  Institut Fourier, 38402 St Martin d'Heres
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

using namespace std;
#include <fstream>
//#include <unistd.h> // For reading arguments from file
#include "ifactor.h"
#include "pari.h"
#include "usual.h"
#include "sym2poly.h"
#include "rpn.h"
#include "prog.h"
#include "misc.h"
#include "giacintl.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  //#undef HAVE_LIBPARI 

#if 0
  static inline gen f(gen x,gen k,gen n){
    return (x*x+k)%n;
  }
  
  static gen pollar(gen n, gen k){
    gen g,m,m1,a,a1,c,x,x1,y,y1,p,j;
    g=1;
    m=2;
    a=1;
    x=2;
    y=2;
    x1=2;
    p=1;
    c=0;
    
    // cout<<"k=(1,-1,3)"<<endl;
    //cin>>k;
    while (g==1) {
      a=2*a+1;//a=2^(e+1)-1=2*l(m)-1 
      while (is_strictly_greater(a,m,context0) && (g==1)) { // ok
	x=f(x,k,n);
	m=m+1;
	p=(p*(x1-x))%n;
	//if (is_strictly_positive(-p)) {p=p+n;}
	c=c+1;
	if (c==32) {
	  g=gcd(abs(p,context0),n); // ok
	  if (g==1) {
	    y=x;
	    y1=x1;
	    a1=a;
	    m1=m;
	    c=0;
	    p=1;
	  }
	}
      }//m=a=2^e-1=l(m)
      if (g==1) {
	x1=x;//x1=x_m=x_l(m)-1
	j=3*iquo(a+1,2);
	for (gen i=m+1;is_positive(j-i,0);i=i+1){
	  x=f(x,k,n);//x_i de i=2^e=l(m) jusque 3*l(m)/2
	}
	m=j;
      }
    }
    //g<>1 ds le paquet de 32
    
    g=1;
    a=iquo(a1-1,2);
    m=m1;
    x=y;
    x1=y1;
    while (g==1) {
      a=2*a+1;
      while (is_strictly_greater(a,m,context0) && (g==1)) { // ok
	x=f(x,k,n);
	m=m+1;
	p=(x1-x)%n;
	//if (is_strictly_positive(-p)) {p=p+n;}
	//c=c+1;
	g=gcd(abs(p,context0),n);  // ok
      }
      if (g==1) {
	x1=x;
	j=3*iquo(a+1,2);
	for (gen i=m+1;is_positive(j-i,0);i=i+1){
	  x=f(x,k,n);//x_i
	}
	m=j;
	//x=f(x1,k,n);
	//m=m+1;
      }
    }
    
    
    if (g==n) {
      if (k==1) 
	return(pollar(n,-1)); 
      else {
	if (k*k==1)
	  return(pollar(n,3));
	else
	  return(g); 
      }
    } 
    else {
      //n=iquo(n,g);
      //cout<<g<<endl;
      return (g);
    }  
  }

#else
  static inline gen f(const gen &x,const gen &k,const gen &n,gen &q){
    return irem(x*x+k,n,q);
  }
  
  static gen pollar_old(gen n, gen k){
    gen g,x,x1,y,y1,p,q;
    long m,m1,a,a1,j;
    g=1;
    m=2;
    a=1;
    x=2;
    y=2;
    x1=2;
    p=1;
    int c=0;
    while (g==1) {
      a=2*a+1;//a=2^(e+1)-1=2*l(m)-1 
      while (g==1 && a>m) { // ok
	x=f(x,k,n,q);
	m += 1;
	if (m > long(1)<<30)
	  return undeferr(gettext("Integer too large"));
	p=irem(p*(x1-x),n,q);
	c += 1;
	if (c==32) {
	  g=gcd(abs(p,context0),n); // ok
	  if (g==1) {
	    y=x;
	    y1=x1;
	    p=1;
	    a1=a;
	    m1=m;
	    c=0;
	  }
	}
      }//m=a=2^e-1=l(m)
      if (g==1) {
	x1=x;//x1=x_m=x_l(m)-1
	j=3*(a+1)/2; // j=3*iquo(a+1,2);
	for (long i=m+1;i<=j;i++)
	  x=f(x,k,n,q);
	m=j;
      }
    }
    //g<>1 ds le paquet de 32
    
    g=1;
    a=(a1-1)/2; // a=iquo(a1-1,2);
    m=m1;
    x=y;
    x1=y1;
    while (g==1) {
      a=2*a+1;
      while (g==1 && a>m) { // ok
	x=f(x,k,n,q);
	m += 1;
	if (m > long(1)<<30)
	  return undeferr(gettext("Integer too large"));
	p=irem(x1-x,n,q);
	g=gcd(abs(p,context0),n);  // ok
      }
      if (g==1) {
	x1=x;
	j=3*(a+1)/2; // j=3*iquo(a+1,2);
	for (long i=m+1;j>=i;i++){
	  x=f(x,k,n,q);//x_i
	}
	m=j;
      }
    }
    
    if (g==n) {
      if (k==1) 
	return(pollar_old(n,-1)); 
      else {
	if (k*k==1)
	  return(pollar_old(n,3));
	else
	  return(g); 
      }
    } 
    else 
      return (g);
  }

  static gen pollar(gen n, gen k){
    k.uncoerce();
    n.uncoerce();
    int m,m1,a,a1,j;
    m1=m=2;
    a1=a=1;
    int c=0;
    mpz_t g,x,x1,x2,x2k,y,y1,p,q;
    mpz_init_set_si(g,1); // ? mp_init_size to specify size
    mpz_init_set_si(x,2);
    mpz_init_set_si(x1,2);
    mpz_init_set_si(y,2);
    mpz_init(y1);
    mpz_init(x2);
    mpz_init(x2k);
    mpz_init_set_si(p,1);
    mpz_init(q);
    while (!ctrl_c && mpz_cmp_si(g,1)==0) {
      a=2*a+1;//a=2^(e+1)-1=2*l(m)-1 
      while (!ctrl_c && mpz_cmp_si(g,1)==0 && a>m) { // ok
	// x=f(x,k,n,q);
#ifdef USE_GMP_REPLACEMENTS
	mp_sqr(&x,&x2);
#else 
	mpz_mul(x2,x,x);
#endif
	mpz_add(x2k,x2,*k._ZINTptr);
	mpz_tdiv_r(x,x2k,*n._ZINTptr);
	m += 1;
	if (m > (1<<20) )
	  return undeferr(gettext("Pollard-rho too many tries"));
	// p=irem(p*(x1-x),n,q);
	mpz_sub(q,x1,x);
	mpz_mul(x2,p,q);
	mpz_tdiv_r(p,x2,*n._ZINTptr);
	c += 1;
	if (c==32) {
	  // g=gcd(abs(p,context0),n); 
	  mpz_abs(q,p);
	  mpz_gcd(g,q,*n._ZINTptr);
	  if (mpz_cmp_si(g,1)==0) {
	    mpz_set(y,x); // y=x;
	    mpz_set(y1,x1); // y1=x1;
	    mpz_set_si(p,1); // p=1;
	    a1=a;
	    m1=m;
	    c=0;
	  }
	}
      }//m=a=2^e-1=l(m)
      if (mpz_cmp_si(g,1)==0) {
	mpz_set(x1,x); // x1=x;//x1=x_m=x_l(m)-1
	j=3*(a+1)/2; // j=3*iquo(a+1,2);
	for (long i=m+1;i<=j;i++){
	  // x=f(x,k,n,q);
	  mpz_mul(x2,x,x);
	  mpz_add(x2k,x2,*k._ZINTptr);
	  mpz_tdiv_r(x,x2k,*n._ZINTptr);
	}
	m=j;
      }
    }
    //g<>1 ds le paquet de 32
    
    mpz_set(x,y); // x=y;
    mpz_set(x1,y1); // x1=y1;
    mpz_set_si(g,1); // g=1;
    a=(a1-1)/2; // a=iquo(a1-1,2);
    m=m1;
    while (!ctrl_c && mpz_cmp_si(g,1)==0) {
      a=2*a+1;
      while (!ctrl_c && mpz_cmp_si(g,1)==0 && a>m) { // ok
	// x=f(x,k,n,q);
	mpz_mul(x2,x,x);
	mpz_add(x2k,x2,*k._ZINTptr);
	mpz_tdiv_r(x,x2k,*n._ZINTptr);
	m += 1;
	if (m > (1<<20) )
	  return undeferr(gettext("Pollard-rho too many tries"));
	// p=irem(x1-x,n,q);
	mpz_sub(q,x1,x);
	mpz_tdiv_r(p,q,*n._ZINTptr);
	// g=gcd(abs(p,context0),n);  // ok
	mpz_abs(q,p);
	mpz_gcd(g,q,*n._ZINTptr);
      }
      if (mpz_cmp_si(g,1)==0) {
	mpz_set(x1,x); // x1=x;
	j=3*(a+1)/2; // j=3*iquo(a+1,2);
	for (long i=m+1;j>=i;i++){
	  // x=f(x,k,n,q);
	  mpz_mul(x2,x,x);
	  mpz_add(x2k,x2,*k._ZINTptr);
	  mpz_tdiv_r(x,x2k,*n._ZINTptr);
	}
	m=j;
      }
    }
    
    mpz_clear(x);
    mpz_clear(x1);
    mpz_clear(x2);
    mpz_clear(x2k);
    mpz_clear(y);
    mpz_clear(y1);
    mpz_clear(p);
    mpz_clear(q);
    if (ctrl_c){
      mpz_clear(g);
      return 0;
    }
    if (mpz_cmp(g,*n._ZINTptr)==0) {
      if (k==1) {
	mpz_clear(g);
	return(pollar(n,-1)); 
      }
      else {
	if (k*k==1){
	  mpz_clear(g);
	  return(pollar(n,3));
	}
	else {
	  ref_mpz_t * ptr=new ref_mpz_t;
	  mpz_init_set(ptr->z,g);
	  mpz_clear(g);
	  return ptr;
	} 
      }
    } 
    ref_mpz_t * ptr=new ref_mpz_t;
    mpz_init_set(ptr->z,g);
    mpz_clear(g);
    return ptr;
  }

#endif
  // const short int giac_primes[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997};
  const short int giac_primes[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,1601,1607,1609,1613,1619,1621,1627,1637,1657,1663,1667,1669,1693,1697,1699,1709,1721,1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,1823,1831,1847,1861,1867,1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,1949,1951,1973,1979,1987,1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063,2069,2081,2083,2087,2089,2099,2111,2113,2129,2131,2137,2141,2143,2153,2161,2179,2203,2207,2213,2221,2237,2239,2243,2251,2267,2269,2273,2281,2287,2293,2297,2309,2311,2333,2339,2341,2347,2351,2357,2371,2377,2381,2383,2389,2393,2399,2411,2417,2423,2437,2441,2447,2459,2467,2473,2477,2503,2521,2531,2539,2543,2549,2551,2557,2579,2591,2593,2609,2617,2621,2633,2647,2657,2659,2663,2671,2677,2683,2687,2689,2693,2699,2707,2711,2713,2719,2729,2731,2741,2749,2753,2767,2777,2789,2791,2797,2801,2803,2819,2833,2837,2843,2851,2857,2861,2879,2887,2897,2903,2909,2917,2927,2939,2953,2957,2963,2969,2971,2999,3001,3011,3019,3023,3037,3041,3049,3061,3067,3079,3083,3089,3109,3119,3121,3137,3163,3167,3169,3181,3187,3191,3203,3209,3217,3221,3229,3251,3253,3257,3259,3271,3299,3301,3307,3313,3319,3323,3329,3331,3343,3347,3359,3361,3371,3373,3389,3391,3407,3413,3433,3449,3457,3461,3463,3467,3469,3491,3499,3511,3517,3527,3529,3533,3539,3541,3547,3557,3559,3571,3581,3583,3593,3607,3613,3617,3623,3631,3637,3643,3659,3671,3673,3677,3691,3697,3701,3709,3719,3727,3733,3739,3761,3767,3769,3779,3793,3797,3803,3821,3823,3833,3847,3851,3853,3863,3877,3881,3889,3907,3911,3917,3919,3923,3929,3931,3943,3947,3967,3989,4001,4003,4007,4013,4019,4021,4027,4049,4051,4057,4073,4079,4091,4093,4099,4111,4127,4129,4133,4139,4153,4157,4159,4177,4201,4211,4217,4219,4229,4231,4241,4243,4253,4259,4261,4271,4273,4283,4289,4297,4327,4337,4339,4349,4357,4363,4373,4391,4397,4409,4421,4423,4441,4447,4451,4457,4463,4481,4483,4493,4507,4513,4517,4519,4523,4547,4549,4561,4567,4583,4591,4597,4603,4621,4637,4639,4643,4649,4651,4657,4663,4673,4679,4691,4703,4721,4723,4729,4733,4751,4759,4783,4787,4789,4793,4799,4801,4813,4817,4831,4861,4871,4877,4889,4903,4909,4919,4931,4933,4937,4943,4951,4957,4967,4969,4973,4987,4993,4999,5003,5009,5011,5021,5023,5039,5051,5059,5077,5081,5087,5099,5101,5107,5113,5119,5147,5153,5167,5171,5179,5189,5197,5209,5227,5231,5233,5237,5261,5273,5279,5281,5297,5303,5309,5323,5333,5347,5351,5381,5387,5393,5399,5407,5413,5417,5419,5431,5437,5441,5443,5449,5471,5477,5479,5483,5501,5503,5507,5519,5521,5527,5531,5557,5563,5569,5573,5581,5591,5623,5639,5641,5647,5651,5653,5657,5659,5669,5683,5689,5693,5701,5711,5717,5737,5741,5743,5749,5779,5783,5791,5801,5807,5813,5821,5827,5839,5843,5849,5851,5857,5861,5867,5869,5879,5881,5897,5903,5923,5927,5939,5953,5981,5987,6007,6011,6029,6037,6043,6047,6053,6067,6073,6079,6089,6091,6101,6113,6121,6131,6133,6143,6151,6163,6173,6197,6199,6203,6211,6217,6221,6229,6247,6257,6263,6269,6271,6277,6287,6299,6301,6311,6317,6323,6329,6337,6343,6353,6359,6361,6367,6373,6379,6389,6397,6421,6427,6449,6451,6469,6473,6481,6491,6521,6529,6547,6551,6553,6563,6569,6571,6577,6581,6599,6607,6619,6637,6653,6659,6661,6673,6679,6689,6691,6701,6703,6709,6719,6733,6737,6761,6763,6779,6781,6791,6793,6803,6823,6827,6829,6833,6841,6857,6863,6869,6871,6883,6899,6907,6911,6917,6947,6949,6959,6961,6967,6971,6977,6983,6991,6997,7001,7013,7019,7027,7039,7043,7057,7069,7079,7103,7109,7121,7127,7129,7151,7159,7177,7187,7193,7207,7211,7213,7219,7229,7237,7243,7247,7253,7283,7297,7307,7309,7321,7331,7333,7349,7351,7369,7393,7411,7417,7433,7451,7457,7459,7477,7481,7487,7489,7499,7507,7517,7523,7529,7537,7541,7547,7549,7559,7561,7573,7577,7583,7589,7591,7603,7607,7621,7639,7643,7649,7669,7673,7681,7687,7691,7699,7703,7717,7723,7727,7741,7753,7757,7759,7789,7793,7817,7823,7829,7841,7853,7867,7873,7877,7879,7883,7901,7907,7919,7927,7933,7937,7949,7951,7963,7993,8009,8011,8017,8039,8053,8059,8069,8081,8087,8089,8093,8101,8111,8117,8123,8147,8161,8167,8171,8179,8191,8209,8219,8221,8231,8233,8237,8243,8263,8269,8273,8287,8291,8293,8297,8311,8317,8329,8353,8363,8369,8377,8387,8389,8419,8423,8429,8431,8443,8447,8461,8467,8501,8513,8521,8527,8537,8539,8543,8563,8573,8581,8597,8599,8609,8623,8627,8629,8641,8647,8663,8669,8677,8681,8689,8693,8699,8707,8713,8719,8731,8737,8741,8747,8753,8761,8779,8783,8803,8807,8819,8821,8831,8837,8839,8849,8861,8863,8867,8887,8893,8923,8929,8933,8941,8951,8963,8969,8971,8999,9001,9007,9011,9013,9029,9041,9043,9049,9059,9067,9091,9103,9109,9127,9133,9137,9151,9157,9161,9173,9181,9187,9199,9203,9209,9221,9227,9239,9241,9257,9277,9281,9283,9293,9311,9319,9323,9337,9341,9343,9349,9371,9377,9391,9397,9403,9413,9419,9421,9431,9433,9437,9439,9461,9463,9467,9473,9479,9491,9497,9511,9521,9533,9539,9547,9551,9587,9601,9613,9619,9623,9629,9631,9643,9649,9661,9677,9679,9689,9697,9719,9721,9733,9739,9743,9749,9767,9769,9781,9787,9791,9803,9811,9817,9829,9833,9839,9851,9857,9859,9871,9883,9887,9901,9907,9923,9929,9931,9941,9949,9967,9973};
  static const int giac_last_prime=giac_primes[sizeof(giac_primes)/sizeof(short)-1];

  static const char _ithprime_s []="ithprime";
  static symbolic symb_ithprime(const gen & args){
    return symbolic(at_ithprime,args);
  }
  static gen ithprime(const gen & g,GIAC_CONTEXT){
    if (g.type!=_INT_)
      return symb_ithprime(g);
    int i=g.val;
    if (i<0)
      return gendimerr(contextptr);
    if (i==0)
      return 1;
    if (i<=sizeof(giac_primes)/sizeof(short int))
      return giac_primes[i-1];
    return symb_ithprime(g);
  }
  gen _ithprime(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type==_VECT)
      return apply(args,_ithprime,contextptr);
    return ithprime(args,contextptr);
  }
  static define_unary_function_eval (__ithprime,&giac::_ithprime,_ithprime_s);
  define_unary_function_ptr5( at_ithprime ,alias_at_ithprime,&__ithprime,0,true);

  bool is_divisible_by(const gen & n,unsigned long a){
    if (n.type==_ZINT){
#ifdef USE_GMP_REPLACEMENTS
      mp_digit c;
      mp_mod_d(n._ZINTptr, a, &c);
      return c==0;
#else
      return mpz_divisible_ui_p(*n._ZINTptr,a);
#endif
    }
    return n.val%a==0;
  }

  // find trivial factors of n, 
  // if add_last is true the remainder is put in the vecteur,
  // otherwise n contains the remainder
  vecteur pfacprem(gen & n,bool add_last,GIAC_CONTEXT){
    gen a;
    gen q;
    int p,i,prime;
    vecteur v(2);
    vecteur u;
    if (is_zero(n))
      return u;
    if (n.type==_ZINT){
      ref_mpz_t * cur = new ref_mpz_t;
      mpz_t div,q,r;
      mpz_set(cur->z,*n._ZINTptr);
      mpz_init_set(q,*n._ZINTptr);
      mpz_init(r);
      mpz_init(div);
      for (i=0;i<sizeof(giac_primes)/sizeof(short int);++i){
	if (mpz_cmp_si(cur->z,1)==0) 
	  break;
	prime=giac_primes[i];
#ifdef USE_GMP_REPLACEMENTS
	mp_digit c;
	mp_mod_d(&cur->z, prime, &c);
	if (c==0){
	  mpz_set_ui(div,prime);
	  for (p=0;;p++){
	    mpz_tdiv_qr(q,r,cur->z,div);
	    if (mpz_cmp_si(r,0))
	      break;
	    mp_exch(&cur->z,&q);
	  }
	  // *logptr(contextptr) << "Factor " << prime << " " << p << endl;
	  u.push_back(prime);
	  u.push_back(p);
	}
#else
	if (mpz_divisible_ui_p(cur->z,prime)){
	  mpz_set_ui(div,prime);
	  for (p=0;;p++){
	    mpz_tdiv_qr(q,r,cur->z,div);
	    if (mpz_cmp_si(r,0))
	      break;
	    mpz_swap(cur->z,q);
	  }
	  // *logptr(contextptr) << "Factor " << prime << " " << p << endl;
	  u.push_back(prime);
	  u.push_back(p);
	}
#endif
      } // end for on smal primes
      mpz_clear(div); mpz_clear(r); mpz_clear(q);
      n=cur;
    }
    else {
      for (i=0;i<sizeof(giac_primes)/sizeof(short int);++i){
	if (n==1) 
	  break;
	a.val=giac_primes[i];
	p=0;
	while (is_divisible_by(n,a.val)){ // while (irem(n,a,q)==0){
	  n=iquo(n,a); 
	  p=p+1;
	}
	if (p!=0){
	  // *logptr(contextptr) << "Factor " << a << " " << p << endl;
	  u.push_back(a);
	  u.push_back(p);
	}
      }
    }
    if (add_last && i==1229 && !is_one(n)){
      u.push_back(n);
      u.push_back(1);
      n=1;
    }
    //v[0]=n;
    //v[1]=u;
    
    return(u);
  }

  static vecteur facprem(gen & n){
    gen a,b,q;
    int p;
    vecteur v(2);    
    if (n==1) {return vecteur(0);}
    if ( (n.type==_INT_ && n.val<giac_last_prime*giac_last_prime) || is_probab_prime_p(n)) {
      v[0]=n;v[1]=1;
      n=1;
      return v;
    }
    b=n;
    if (debug_infolevel>5)
      cerr << "Pollard begin " << clock() << endl;
    a=pollar(b,1);
    if (is_zero(a))
      return makevecteur(gensizeerr(gettext("Stopped by user interruption")));
    p=0;
    while ( a!=b && is_greater(a,giac_last_prime*giac_last_prime,context0) && is_probab_prime_p(a)==0 ) {
      b=a;
      a=pollar(b,1);
      if (is_zero(a))
	return makevecteur(gensizeerr(gettext("Stopped by user interruption")));
    }
    if (debug_infolevel>5)
      cerr << "Pollard end " << clock() << endl;
    //n=iquo(n,a);
    
    while (irem(n,a,q)==0){
      n=iquo(n,a); // n=q does not work because irem assumes n and q does not point to the same _ZINT
      p=p+1;
    }
    v[0]=a;
    if (a==b){v[1]=-p;} else {v[1]=p;}
    return v;
  }

  vecteur ifactors(const gen & n0){
    if (is_zero(n0))
      return vecteur(1,gensizeerr(gettext("ifactors")));
    if (is_one(n0))
      return vecteur(0);
#ifdef HAVE_LIBPARI
    gen g(pari_ifactor(n0),context0); // FIXME GIAC_CONTEXT
    if (g.type==_VECT){
      matrice m(mtran(*g._VECTptr));
      vecteur res;
      const_iterateur it=m.begin(),itend=m.end();
      for (;it!=itend;++it){
	if (it->type!=_VECT) return vecteur(1,gensizeerr(gettext("ifactor.cc/ifactors")));
	res.push_back(it->_VECTptr->front());
	res.push_back(it->_VECTptr->back());
      }
      return res;
    }
#else // LIBPARI
    gen n(n0);
    vecteur f;
    vecteur g;
    vecteur u;
    f=pfacprem(n,false,context0);
    //cout<<n<<" "<<f<<endl;
    while (n!=1) {
      g=facprem(n);
      if (is_undef(g))
	return g;
      //cout<<n<<" "<<g<<endl;
      u=mergevecteur(u,g);
    }
    
    g=mergevecteur(f,u);
    return g;
#endif // LIBPARI
    return 0;
  }

  vecteur ifactors(const gen & r,const gen & i,const gen & ri){
    gen norm=r*r+i*i;
    gen reste(ri);
    const vecteur & facto = ifactors(norm);
    if (is_undef(facto))
      return facto;
    int l=facto.size()/2;
    vecteur res;
    for (int i=0;i<l;++i){
      gen prime=facto[2*i];
      int mult=facto[2*i+1].val,multp=0;
      int n=smod(prime,4).val;
      if (n==2){
	res.push_back(1+cst_i);
	res.push_back(mult);
	reste=reste/pow(1+cst_i,mult,context0);
	continue;
      }
      if (n==-1){
	res.push_back(prime);
	res.push_back(mult/2);
	reste=reste/pow(prime,mult/2,context0);
	continue;
      }
      prime=pa2b2(prime);
      prime=gen(prime[0],prime[1]);
      for (;mult>0;--mult,++multp){
	if (!is_zero(reste % prime))
	  break;
	reste=reste/prime;
      }
      if (multp){
	res.push_back(prime);
	res.push_back(multp);
      }
      if (mult){
	prime=conj(prime,context0);
	res.push_back(prime);
	res.push_back(mult);
	reste=reste/pow(prime,mult,context0);
      }
    }
    if (!is_one(reste)){
      res.insert(res.begin(),1);
      res.insert(res.begin(),reste);
    }
    return res;
  }

  gen ifactors(const gen & args,int maplemode){
    if ( (args.type==_INT_) || (args.type==_ZINT)){
      if (is_zero(args)){
	if (maplemode==1)
	  return makevecteur(args,vecteur(0));
	else
	  return makevecteur(args);
      }
      vecteur v(ifactors(abs(args,context0))); // ok
      if (!v.empty() && is_undef(v.front()))
	return v.front();
      if (maplemode!=1){
	if (is_positive(args,context0))
	  return v;
	return mergevecteur(makevecteur(minus_one,plus_one),v);
      }
      vecteur res;
      const_iterateur it=v.begin(),itend=v.end();
      for (;it!=itend;it+=2){
	res.push_back(makevecteur(*it,*(it+1)));
      }
      if (is_positive(args,context0))
	return makevecteur(plus_one,res);
      else
	return makevecteur(minus_one,res);	
    }
    if (args.type==_CPLX && is_integer(*args._CPLXptr) && is_integer(*(args._CPLXptr+1)))
      return ifactors(*args._CPLXptr,*(args._CPLXptr+1),args);
    return gentypeerr(gettext("ifactors"));
  }
  static symbolic symb_ifactors(const gen & args){
    return symbolic(at_ifactors,args);
  }
  gen _ifactors(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type==_VECT)
      return apply(args,_ifactors,contextptr);
    gen g(args);
    if (!is_integral(g))
      return gensizeerr(contextptr);
    return ifactors(g,0);
  }
  static const char _ifactors_s []="ifactors";
  static define_unary_function_eval (__ifactors,&giac::_ifactors,_ifactors_s);
  define_unary_function_ptr5( at_ifactors ,alias_at_ifactors,&__ifactors,0,true);

  static const char _facteurs_premiers_s []="facteurs_premiers";
  static define_unary_function_eval (__facteurs_premiers,&giac::_ifactors,_facteurs_premiers_s);
  define_unary_function_ptr5( at_facteurs_premiers ,alias_at_facteurs_premiers,&__facteurs_premiers,0,true);

  static const char _maple_ifactors_s []="maple_ifactors";
  static symbolic symb_maple_ifactors(const gen & args){
    return symbolic(at_maple_ifactors,args);
  }
  gen _maple_ifactors(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type==_VECT)
      return apply(args,_maple_ifactors,contextptr);
    return ifactors(args,1);
  }
  static define_unary_function_eval (__maple_ifactors,&giac::_maple_ifactors,_maple_ifactors_s);
  define_unary_function_ptr5( at_maple_ifactors ,alias_at_maple_ifactors,&__maple_ifactors,0,true);

  static vecteur in_factors(const gen & gf){
    if (gf.type!=_SYMB)
      return makevecteur(gf,plus_one);
    unary_function_ptr & u=gf._SYMBptr->sommet;
    if (u==at_inv){
      vecteur v=in_factors(gf._SYMBptr->feuille);
      iterateur it=v.begin(),itend=v.end();
      for (;it!=itend;it+=2)
	*(it+1)=-*(it+1);
      return v;
    }
    if (u==at_neg){
      vecteur v=in_factors(gf._SYMBptr->feuille);
      v.push_back(minus_one);
      v.push_back(plus_one);
      return v;
    }
    if ( (u==at_pow) && (gf._SYMBptr->feuille._VECTptr->back().type==_INT_) ){
      vecteur v=in_factors(gf._SYMBptr->feuille._VECTptr->front());
      gen k=gf._SYMBptr->feuille._VECTptr->back();
      iterateur it=v.begin(),itend=v.end();
      for (;it!=itend;it+=2)
	*(it+1)=k* *(it+1);
      return v;
    }
    if (u!=at_prod)
      return makevecteur(gf,plus_one);
    vecteur res;
    const_iterateur it=gf._SYMBptr->feuille._VECTptr->begin(),itend=gf._SYMBptr->feuille._VECTptr->end();
    for (;it!=itend;++it){
      res=mergevecteur(res,in_factors(*it));
    }
    return res;
  }
  vecteur factors(const gen & g,const gen & x,GIAC_CONTEXT){
    gen gf=factor(g,x,false,contextptr);
    vecteur res=in_factors(gf);
    if (xcas_mode(contextptr)!=1)
      return res;
    gen coeff(1);
    vecteur v;
    const_iterateur it=res.begin(),itend=res.end();
    for (;it!=itend;it+=2){
      if (lidnt(*it).empty())
	coeff=coeff*(pow(*it,*(it+1),contextptr));
      else
	v.push_back(makevecteur(*it,*(it+1)));
    }
    return makevecteur(coeff,v);
  }
  static const char _factors_s []="factors";
  gen _factors(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type==_VECT)
      return apply(args,_factors,contextptr);
    return factors(args,vx_var,contextptr);
  }
  static define_unary_function_eval (__factors,&giac::_factors,_factors_s);
  define_unary_function_ptr5( at_factors ,alias_at_factors,&__factors,0,true);

  gen ifactors2ifactor(const vecteur & l){
    int s;
    s=l.size();
    gen r;
    vecteur v(s/2);
    for (int j=0;j<s;j=j+2){
      if (!is_one(l[j+1]))
	v[j/2]=symbolic(at_pow,gen(makevecteur(l[j],l[j+1]),_SEQ__VECT));
      else
	v[j/2]=l[j];
    }
    if (v.size()==1){
#ifdef GIAC_HAS_STO_38
      return symb_quote(v.front());
#else
      return v.front();
#endif
    }
    r=symbolic(at_prod,v);
#ifdef GIAC_HAS_STO_38
    r=symb_quote(r);
#endif
    return(r);
  }
  gen ifactor(const gen & n){
    vecteur l;
    l=ifactors(n);
    if (!l.empty() && is_undef(l.front())) return l.front();
    return ifactors2ifactor(l);
  }
  gen _ifactor(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type==_CPLX && is_integer(*args._CPLXptr) && is_integer(*(args._CPLXptr+1))){
      const vecteur & v=ifactors(*args._CPLXptr,*(args._CPLXptr+1),args);
      return ifactors2ifactor(v);
    }
    gen n=args;
    if (!is_integral(n))
      return gensizeerr(contextptr);
    if (is_strictly_positive(-n,0))
      return -_ifactor(-n,contextptr);
    if (n.type==_INT_ && n.val<=3)
      return n;
    return ifactor(n);
  }
  static const char _ifactor_s []="ifactor";
  static define_unary_function_eval (__ifactor,&_ifactor,_ifactor_s);
  define_unary_function_ptr5( at_ifactor ,alias_at_ifactor,&__ifactor,0,true);

  static const char _factoriser_entier_s []="factoriser_entier";
  static define_unary_function_eval (__factoriser_entier,&_ifactor,_factoriser_entier_s);
  define_unary_function_ptr5( at_factoriser_entier ,alias_at_factoriser_entier,&__factoriser_entier,0,true);

  static vecteur divis(const vecteur & l3){
    vecteur l1(1);
    gen d,e;
    int s=l3.size();
    l1[0]=1;//l3.push_back(..);
    for (int k=0;k<s;k=k+2) {
      vecteur l2;
      int s1;
      s1=l1.size();
      vecteur l4(s1);
      d=l3[k];
      e=l3[k+1];
      int ei;
      if (e.type==_INT_){
	ei=e.val;
      }
      else
	return vecteur(1,gensizeerr(gettext("Integer too large")));
      for (int j=1;j<=ei;j++){
	for (int l=0;l<s1;l++){ 
	  l4[l]=l1[l]*pow(d,j);
	}
	l2=mergevecteur(l2,l4);
      }
      l1=mergevecteur(l1,l2);
    }
    return(l1); 
  }
  gen idivis(const gen & n){
    vecteur l3(ifactors(n));
    if (!l3.empty() && is_undef(l3.front())) return l3.front();
    return divis(l3);
  }
  gen _idivis(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type==_VECT)
      return apply(args,_idivis,contextptr);
    gen n=args;
    if (!is_integral(n) && !is_integer(n)) 
      return gensizeerr(contextptr);
    return idivis(abs(n,contextptr));
  }
  static const char _idivis_s []="idivis";
  static define_unary_function_eval (__idivis,&_idivis,_idivis_s);
  define_unary_function_ptr5( at_idivis ,alias_at_idivis,&__idivis,0,true);

  gen _divis(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type==_VECT)
      return apply(args,_divis,contextptr);
    return divis(factors(args,vx_var,contextptr));
  }
  static const char _divis_s []="divis";
  static define_unary_function_eval (__divis,&_divis,_divis_s);
  define_unary_function_ptr5( at_divis ,alias_at_divis,&__divis,0,true);

  /*
  gen ichinreme(const vecteur & a,const vecteur & b){
    vecteur r(2);
    gen p=a[1],q=b[1],u,v,d;
    egcd(p,q,u,v,d);
    if (d!=1)  return gensizeerr(contextptr);
    r[0]=(u*p*b[0]+v*q*a[0]%p*q);
    r[1]=p*q;
    return(r);
  }
  gen _ichinreme(const gen & args){
  if ( args){
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=4) )
      return gensizeerr(contextptr);
    vecteur a(2).type==_STRNG && args.subtype==-1{
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=4) )
      return gensizeerr(contextptr);
    vecteur a(2))) return  args){
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=4) )
      return gensizeerr(contextptr);
    vecteur a(2);
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=4) )
      return gensizeerr(contextptr);
    vecteur a(2),b(2);
    a[0]=args[0];
    a[1]=args[1];
    b[0]=args[2];
    b[1]=args[3];
    //gen a=args[0],p=args[1], b=args[2],q=args[3];
    return ichinreme(a,b);
  }
  static const char _ichinreme_s []="ichinreme";
  static define_unary_function_eval (__ichinreme,&_ichinreme,_ichinreme_s);
  define_unary_function_ptr5( at_ichinreme ,alias_at_ichinreme,&__ichinreme,0,true); 
  */

  gen euler(const gen & e){
    if (e==0)
      return e;
    vecteur v(ifactors(e));
    if (!v.empty() && is_undef(v.front())) return v.front();
    const_iterateur it=v.begin(),itend=v.end();
    for (gen res(plus_one);;){
      if (it==itend)
	return res;
      gen p=*it;
      ++it;
      int n=it->val;
      res = res * (p-plus_one)*pow(p,n-1);
      ++it;
    }
  }
  static symbolic symb_euler(const gen & args){
    return symbolic(at_euler,args);
  }
  static const char _euler_s []="euler";
  gen _euler(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type==_VECT)
      return apply(args,_euler,contextptr);
    if ( (args.type==_INT_) || (args.type==_ZINT))
      return euler(args);
    return gentypeerr(contextptr);
  }
  static define_unary_function_eval (__euler,&giac::_euler,_euler_s);
  define_unary_function_ptr5( at_euler ,alias_at_euler,&__euler,0,true);

  gen pa2b2(const gen & p){
    if ((p%4)!=1) return gensizeerr(gettext("pa2b2"));// car p!=1 mod 4
    gen q=(p-1)/4;
    gen a=2;
    gen ra;
    ra=powmod(a,q,p);
    //on cherche ra^2=-1 mod p avec ra!=1 et ra !=p-1
    while ((a!=p-1) && ((ra==1)|| (ra==p-1))){
      a=a+1;
      ra=powmod(a,q,p);
    }
    if ((ra==1)||(ra==p-1))  return gensizeerr(gettext("pa2b2"));//car p n'est pas premier
    gen ux=1,uy=ra,vx=0,vy=p,wx,wy; 
    gen m=1;
    while(m!=0){
      if (is_positive(vx*vx+vy*vy-ux*ux+uy*uy,0)){
	//on echange u et v
	wx=vx;
	wy=vy;
	vx=ux;
	vy=uy;
	ux=wx;
	uy=wy;
      }
      gen alpha=inv(2,context0)-(ux*vx+uy*vy)*inv(vx*vx+vy*vy,context0);
      //m=partie entiere de alpha (-v.v/2<(u+mv).v<=v.v/2)
      m=_floor(alpha,0);
      ux=ux+m*vx;
      uy=uy+m*vy;
    }
    vecteur v(2);
    //v repond a la question
    v[0]=abs(vx,context0); // ok
    v[1]=abs(vy,context0); // ok
    if (vx*vx+vy*vy!=p)
      return gensizeerr(gettext("pa2b2"));
    return v;
  }
  gen _pa2b2(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (!is_integer(args)) 
      return gensizeerr(contextptr);
    gen n=args;
    return pa2b2(n);
  }
  static const char _pa2b2_s []="pa2b2";
  static define_unary_function_eval (__pa2b2,&_pa2b2,_pa2b2_s);
  define_unary_function_ptr5( at_pa2b2 ,alias_at_pa2b2,&__pa2b2,0,true);

  static gen ipropfrac(const gen & a,const gen & b){
    gen r=a%b;
    gen q=(a-r)/b;
    gen d=gcd(r,b);
    r=r/d;
    gen b1=b/d;
    gen v;
    v=symbolic(at_division,makevecteur(r,b1));
    gen w;
    w=symbolic(at_plus,makevecteur(q,v));    
    return w;
  }
  gen _propfrac(const gen & arg,GIAC_CONTEXT){
    if ( arg.type==_STRNG && arg.subtype==-1) return  arg;
    gen args(arg);
    vecteur v;
    if (arg.type==_VECT && arg._VECTptr->size()==2){
      v=vecteur(1,arg._VECTptr->back());
      args=arg._VECTptr->front();
      lvar(args,v);
    }
    else
      v=lvar(arg);
    gen g=e2r(args,v,contextptr);
    gen a,b;
    fxnd(g,a,b);
    if (v.empty())
      return ipropfrac(a,b);
    else {
      gen d=r2e(b,v,contextptr);
      g=_quorem(makevecteur(r2e(a,v,contextptr),d,v.front()),contextptr);
      if (is_undef(g)) return g;
      vecteur &v=*g._VECTptr;
      return v[0]+rdiv(v[1],d);
    }
  }
  static const char _propfrac_s []="propfrac";
  static define_unary_function_eval (__propfrac,&_propfrac,_propfrac_s);
  define_unary_function_ptr5( at_propfrac ,alias_at_propfrac,&__propfrac,0,true);

  gen iabcuv(const gen & a,const gen & b,const gen & c){
    gen d=gcd(a,b);
    if (c%d!=0)  return gensizeerr(gettext("iabcuv"));
    gen a1=a/d,b1=b/d,c1=c/d;
    gen u,v,w;
    egcd(a1,b1,u,v,w);
    vecteur r(2);
    r[0]=smod(u*c1,b);
    r[1]=iquo(c-r[0]*a,b);
    return r;
  }
  gen _iabcuv(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=3) )
      return gensizeerr(contextptr);
    gen a=args[0],b=args[1],c=args[2];
    return iabcuv(a,b,c);
  }
  static const char _iabcuv_s []="iabcuv";
  static define_unary_function_eval (__iabcuv,&_iabcuv,_iabcuv_s);
  define_unary_function_ptr5( at_iabcuv ,alias_at_iabcuv,&__iabcuv,0,true);

  gen abcuv(const gen & a,const gen & b,const gen & c,const gen & x,GIAC_CONTEXT){
    gen g=_egcd(makevecteur(a,b,x),contextptr);
    if (is_undef(g)) return g;
    vecteur & v=*g._VECTptr;
    gen h=_quorem(makevecteur(c,v[2],x),contextptr);
    if (is_undef(h)) return h;
    vecteur & w=*h._VECTptr;
    if (!is_zero(w[1]))
      return gensizeerr(gettext("No solution in ring"));
    gen U=v[0]*w[0],V=v[1]*w[0];
    if (_degree(makevecteur(c,x),contextptr).val<_degree(makevecteur(a,x),contextptr).val+_degree(makevecteur(b,x),contextptr).val ){
      U=_rem(makevecteur(U,b,x),contextptr);
      V=_rem(makevecteur(V,a,x),contextptr);
    }
    return makevecteur(U,V);
  }
  gen _abcuv(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if ( (args.type!=_VECT) || (args._VECTptr->size()<3) )
      return gensizeerr(contextptr);
    vecteur & v =*args._VECTptr;
    if (v.size()>3)
      return abcuv(v[0],v[1],v[2],v[3],contextptr);
    return abcuv(v[0],v[1],v[2],vx_var,contextptr);
  }
  static const char _abcuv_s []="abcuv";
  static define_unary_function_eval (__abcuv,&_abcuv,_abcuv_s);
  define_unary_function_ptr5( at_abcuv ,alias_at_abcuv,&__abcuv,0,true);

  gen simp2(const gen & a,const gen & b,GIAC_CONTEXT){
    vecteur r(2);
    gen d=gcd(a,b);
    r[0]=normal(a/d,contextptr);
    r[1]=normal(b/d,contextptr);
    return r;
  }
  gen _simp2(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2) )
      return gensizeerr(contextptr);
    gen a=args[0],b=args[1];
    if ( (a.type==_VECT) || (b.type==_VECT) )
      return gensizeerr(contextptr);
    return simp2(a,b,contextptr);
  }
  static const char _simp2_s []="simp2";
  static define_unary_function_eval (__simp2,&_simp2,_simp2_s);
  define_unary_function_ptr5( at_simp2 ,alias_at_simp2,&__simp2,0,true);
 
  gen fxnd(const gen & a){
    vecteur v(lvar(a));
    gen g=e2r(a,v,context0); // ok
    gen n,d;
    fxnd(g,n,d);
    return makevecteur(r2e(n,v,context0),r2e(d,v,context0)); // ok
  }
  gen _fxnd(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type==_VECT)
      return apply(args,fxnd);
    return fxnd(args);
  }
  static const char _fxnd_s []="fxnd";
  static define_unary_function_eval (__fxnd,&_fxnd,_fxnd_s);
  define_unary_function_ptr5( at_fxnd ,alias_at_fxnd,&__fxnd,0,true);
 
#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
