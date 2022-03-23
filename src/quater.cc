// -*- mode:C++ ; compile-command: "g++-3.4 -I.. -g -c quater.cc -DHAVE_CONFIG_H -DIN_GIAC" -*-
#include "giacPCH.h"
/*
 *  Copyright (C) 2001,2007 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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

#include "quater.h"
#include "unary.h"
#include "sym2poly.h"
#include "usual.h"
#include "intg.h"
#include "subst.h"
#include "derive.h"
#include "lin.h"
#include "vecteur.h"
#include "gausspol.h"
#include "plot.h"
#include "prog.h"
#include "modpoly.h"
#include "series.h"
#include "tex.h"
#include "ifactor.h"
#include "risch.h"
#include "solve.h"
#include "giacintl.h"
using namespace std;

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  static const char _quaternion_s []="quaternion";
  gen _quaternion(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type!=_VECT)
      return quaternion(args);
    vecteur v(*args._VECTptr);
    if (v.size()==1)
      return quaternion(v.front());
    if (v.size()!=4)
      return gensizeerr(gettext("Quaternion has 1 or 4 arguments"));
    return quaternion(v[0],v[1],v[2],v[3]);
  }
  static define_unary_function_eval (__quaternion,&giac::_quaternion,_quaternion_s);
  define_unary_function_ptr5( at_quaternion ,alias_at_quaternion,&__quaternion,0,true); // auto-register
  
  string quaternion::print(GIAC_CONTEXT) const {
    return string(_quaternion_s)+"("+r.print()+","+i.print()+","+j.print()+","+k.print()+")";
    this->dbgprint();
  }

  quaternion::quaternion(const gen & g){
    if (g.type==_USER){
      const quaternion  * q =dynamic_cast<const quaternion *>(g._USERptr);
      if (q)
	*this=*q;
    }
    else {
      r=g;
      i=zero;
      j=zero;
      k=zero;
    }
  };

  static const char _galois_field_s []="GF";

  void lrdm(modpoly & p,int n); // form intg.cc

  // Is the polynomial v irreducible and primitive modulo p?
  // If it is only irreducible, returns 2 and sets vmin
  int is_irreducible_primitive(const modpoly & v,const gen & p,modpoly & vmin){
    vmin=v;
    int m=v.size()-1;
    if (m<2)
      return 0; // setsizeerr(gettext("irreducibility: degree too short")+gen(v).print());
    gen gpm=pow(p,m);
    if (gpm.type!=_INT_)
      return 0; // setsizeerr(gettext("Characteristic or degree too large p^m=")+gpm.print());
    int pm=gpm.val;
    environment * env=new environment;
    env->modulo=p;
    env->pn=env->modulo;
    env->moduloon=true;
    vecteur polyx(2),g;
    polyx[0]=1;
    vecteur test(p.val+1);
    test[0]=1;
    // Irreducible: v must be prime with x^(p^k)-x, for k<=m/2
    for (int k=1;k<=m/2;k++){
      gcdmodpoly(operator_minus(test,polyx,env),v,env,g);
      if (!is_one(g)){
	delete env;
	return 0;
      }
      test=powmod(test,p,v,env);
    }
    // Primi: must not divide x^[(p^m-1)/d]-1 for any prime divisor d of p^m-1
    gen tmp=pm-1;
    vecteur vp(pfacprem(tmp,true,context0));
    int ntest=vp.size();
    for (int i=0;i<ntest;i+=2){
      // Compute x^[(p^m-1)/d] mod v, mod p, is it 1?
      int pm_d=(pm-1)/vp[i].val;
      test=powmod(polyx,pm_d,v,env);
      if (is_one(test)){
	vecteur cyclic;
	// Find a cyclic element in GF, [1,0] does not work
	for (int k=p.val+1;k<pm;k++){
	  cyclic=vecteur(m);
	  // decompose k in base p
	  for (int j=0,k1=k;j<m;++j){
	    cyclic[j]=k1%p.val;
	    k1/=p.val;
	  }
	  cyclic=trim(cyclic,0);
	  // ?cyclic
	  for (int i=0;i<ntest;i+=2){
	    int pm_d=(pm-1)/vp[i].val;
	    test=powmod(cyclic,pm_d,v,env);
	    if (is_one(test))
	      break; // not cyclic
	  }
	  if (!is_one(test)) // cyclic! 
	    break;
	}
	// cyclic is cyclic, find it's minimal polynomial
	// Compute 1,cyclic, ..., cyclic^m and find kernel 
	matrice minmat(m+1);
	minmat[0]=vecteur(1,1);
	for (int i=1;i<=m;++i)
	  minmat[i]=operator_mod(operator_times(cyclic,*minmat[i-1]._VECTptr,env),v,env);
	for (int i=0;i<=m;++i)
	  lrdm(*minmat[i]._VECTptr,m-1);
	minmat=mtran(minmat);
	matrice minred,pivots; gen det;
	if (!modrref(minmat,minred,pivots,det,0,m,0,m+1,true,0,p,0))
	  return 0;
	// Extract kernel from last column
	vmin=vecteur(m+1,1);
	for (int i=1;i<=m;++i)
	  vmin[i]=-minred[m-i][m];
	// vecteur tmpv;
	// cout << is_irreducible_primitive(vmin,p,tmpv) << endl;
	delete env;
	return 2;
      }
      /* vecteur test(pm_d+1);
	 test[0]=1;
	 test[pm_d]=-1; 
      if (is_zero(operator_mod(test,v,env))){
	delete env;
	return false;
      }
      */
    }
    delete env;
    return 1;
  }

  vecteur find_irreducible_primitive(int p,int m){
    // Now test all possible coeffs for test[k] until it's irreducible
    int pm=int(std::pow(double(p),double(m)));
    for (int k=0;k<pm;k++){
      vecteur test(m+1),test2;
      test[0]=1;
      // decompose k in base p
      for (int j=1,k1=k;j<=m;++j){
	test[j]=k1%p;
	k1/=p;
      }
      if (is_irreducible_primitive(test,p,test2))
	return test2;
    }
    return vecteur(1,gensizeerr(gettext("No irreducible primitive polynomial found")));
  }
  gen _galois_field(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (is_integer(args)){ // must be a power of a prime
      gen pm=abs(args,context0); // ok
      vecteur u(pfacprem(pm,true,contextptr));
      if (u.size()!=2)
	return gensizeerr(gettext("Not a power of a prime"));
      return galois_field(u);
    }
    if (args.type!=_VECT)
      return galois_field(args);
    vecteur v(*args._VECTptr);
    int s=v.size();
    if (s==3 && v[1].type!=_INT_){
      v.push_back(undef);
      ++s;
    }
    if (s==2 || s==3)
      return galois_field(args);
    if (s!=4)
      return gensizeerr(gettext("galois_field has 1 or 4 arguments (charac p, irred poly P, var name x, value as a poly of x or as a vector)"));
    vecteur a,P,vmin;
    gen & x= v[2];
    gen xid(x);
    if (x.type==_VECT && !x._VECTptr->empty())
      xid=x._VECTptr->front();
    if (!is_undef(v[3]) && v[3].type!=_VECT)
      v[3]=_e2r(makevecteur(v[3],xid),context0); // ok
    if (v[1].type!=_VECT)
      v[1]=_e2r(makevecteur(v[1],xid),context0); // ok
    if (v[1].type!=_VECT)
      return gensizeerr();
    if (!is_irreducible_primitive(*v[1]._VECTptr,v[0],vmin))
      return gensizeerr(gettext("Not irreducible or not primitive polynomial")+args.print());
    return galois_field(v[0],gen(vmin,_POLY1__VECT),v[2],v[3]);
  }
  static define_unary_function_eval (__galois_field,&giac::_galois_field,_galois_field_s);
  define_unary_function_ptr5( at_galois_field ,alias_at_galois_field,&__galois_field,0,true); // auto-register
  
  string galois_field::print(GIAC_CONTEXT) const {
    gen xid(x);
    if (x.type==_VECT && x._VECTptr->size()>=2){
      xid=x._VECTptr->front();
      if (!is_undef(a))
	return x._VECTptr->back().print()+"("+r2e(a,xid,contextptr).print()+")";      
    }
    return string(_galois_field_s)+"("+p.print()+","+r2e(P,xid,contextptr).print()+","+x.print()+","+r2e(a,xid,contextptr).print()+")";
    this->dbgprint(); // not reached, it's for the debugger
    return "";
  }

  galois_field::galois_field(const gen & g){
    if (g.type==_USER){
      const galois_field  * q =dynamic_cast<const galois_field *>(g._USERptr);
      if (q)
	*this=*q;
      else {
	P=gensizeerr(gettext("Unable to convert user type to galois field"));
      }
    }
    else {
      if (g.type!=_VECT || g._VECTptr->size()<2 || g._VECTptr->front().type!=_INT_ || (*g._VECTptr)[1].type!=_INT_)
	P=gensizeerr(gettext("Expecting characteristic p, integer m"));
      else {
	int p0=g._VECTptr->front().val; // max(absint(),2);
	if (p0<2)
	  P=gensizeerr(gettext("Bad characteristic: ")+print_INT_(p0));
	else {
	  int m0=(*g._VECTptr)[1].val; // max(absint(),2);
	  if (m0<2)
	    P=gensizeerr(gettext("Exponent must be >=2: ")+print_INT_(m0));
	  else {
	    p=p0;
	    P=find_irreducible_primitive(p0,m0);
	    x=g._VECTptr->size()>2?(*g._VECTptr)[2]:vx_var;
	    a=undef;
	  }
	}
      }
    }
  };
  
  void galois_field::reduce(){
    if (!is_undef(a)){
      a = smod(a,p);
      if (a.type!=_VECT)
	a=gen(vecteur(1,a),_POLY1__VECT);
    }
  }

  galois_field::galois_field(const gen p_,const gen & P_,const gen & x_,const gen & a_):p(p_),P(P_),x(x_),a(a_) {
    reduce();
  }

  galois_field::galois_field(const galois_field & q):p(q.p),P(q.P),x(q.x),a(q.a) { 
    reduce();
  }

  gen galois_field::operator + (const gen & g) const { 
    if (is_integer(g))
      return galois_field(p,P,x,a+g);
    if (g.type!=_USER)
      return sym_add(*this,g,context0); // ok symbolic(at_plus,makevecteur(g,*this));
    if (galois_field * gptr=dynamic_cast<galois_field *>(g._USERptr)){
      if (gptr->p!=p || gptr->P!=P || is_undef(P) || is_undef(gptr->P))
	return gensizeerr();
      if (a.type==_VECT && gptr->a.type==_VECT){
	vecteur res;
	environment * env=new environment;
	env->modulo=p;
	env->pn=env->modulo;
	env->moduloon=true;
	addmodpoly(*a._VECTptr,*gptr->a._VECTptr,env,res);
	delete env;
	return galois_field(p,P,x,res);
      }
      return galois_field(p,P,x,a+gptr->a);
    }
    else
      return gensizeerr();
  }

  gen galois_field::operator - (const gen & g) const { 
    if (is_integer(g))
      return galois_field(p,P,x,a-g);
    if (g.type!=_USER)
      return sym_add(*this,-g,context0); // ok symbolic(at_plus,makevecteur(-g,*this));
    if (galois_field * gptr=dynamic_cast<galois_field *>(g._USERptr)){
      if (gptr->p!=p || gptr->P!=P || is_undef(P) || is_undef(gptr->P))
	return gensizeerr();
      if (a.type==_VECT && gptr->a.type==_VECT){
	vecteur res;
	environment * env=new environment;
	env->modulo=p;
	env->pn=env->modulo;
	env->moduloon=true;
	submodpoly(*a._VECTptr,*gptr->a._VECTptr,env,res);
	delete env;
	return galois_field(p,P,x,res);
      }
      return galois_field(p,P,x,a-gptr->a);
    }
    else
      return gensizeerr();
  }

  gen galois_field::operator - () const { 
    return galois_field(p,P,x,-a);
  }

  gen galois_field::operator * (const gen & g) const { 
    if (is_integer(g)){
      gen tmp=smod(g,p);
      if (giac::is_zero(tmp))
	return zero;
      return galois_field(p,P,x,g*a);
    }
    if (g.type!=_USER)
      return sym_mult(*this,g,context0); // ok symbolic(at_prod,makevecteur(g,*this));
    if (galois_field * gptr=dynamic_cast<galois_field *>(g._USERptr)){
      if (gptr->p!=p || gptr->P!=P || P.type!=_VECT || is_undef(P) || is_undef(gptr->P))
	return gensizeerr();
      if (a.type==_VECT && gptr->a.type==_VECT){
	vecteur res;
	environment * env=new environment;
	env->modulo=p;
	env->pn=env->modulo;
	env->moduloon=true;
	mulmodpoly(*a._VECTptr,*gptr->a._VECTptr,env,res);
	res=operator_mod(res,*P._VECTptr,env),
	delete env;
	return galois_field(p,P,x,res);
      }
      return galois_field(p,P,x,a*gptr->a);
    }
    else
      return gensizeerr();
  }

  gen galois_field::inv () const {
    if (a.type!=_VECT || P.type!=_VECT || is_undef(P))
      return gensizeerr(gettext("galois field inv"));
    vecteur & A = *a._VECTptr;
    if (A.empty())
      return galois_field(p,P,x,undef);
    modpoly u,v,d;
    environment * env=new environment;
    env->modulo=p;
    env->pn=env->modulo;
    env->moduloon=true;
    egcd(A,*P._VECTptr,env,u,v,d);
    delete env;
    // d should be [1]
    if (d!=vecteur(1,1))
      return gensizeerr(gettext("GF inv internal bug"));
    return galois_field(p,P,x,u);
  }

  bool galois_field::operator == (const gen & g) const {
    if (is_zero())
      return giac::is_zero(g);
    if (g.type!=_USER)
      return a==vecteur(1,g);
    if (galois_field * gptr=dynamic_cast<galois_field *>(g._USERptr)){
      if (gptr->p!=p || gptr->P!=P)
	return false;
      return gptr->a==a;
    }
    return false;
  }

  bool galois_field::is_zero () const {
    return a.type==_VECT && ( a._VECTptr->empty() || (a._VECTptr->size()==1 && a._VECTptr->front()==0) );
  }

  bool galois_field::is_one () const {
    return a.type==_VECT && a._VECTptr->size()==1 && a._VECTptr->front()==1;
  }

  bool galois_field::is_minus_one () const {
    return a.type==_VECT && a._VECTptr->size()==1 && smod(a._VECTptr->front(),p)==-1;
  }

  gen galois_field::operator () (const gen & g,GIAC_CONTEXT) const {
    if (is_undef(a)){
      gen res;
      if (g.type==_VECT)
	res=g;
      else {
	gen xid(x);
	if (x.type==_VECT && !x._VECTptr->empty())
	  xid=x._VECTptr->front();
	res=_e2r(makevecteur(g,xid),contextptr);
      }
      if (res.type==_VECT){
	environment env;
	env.modulo=p;
	env.pn=env.modulo;
	env.moduloon=true;
	res=smod(res,p);
	res=operator_mod(*res._VECTptr,*P._VECTptr,&env);
      }
      return galois_field(p,P,x,res);
    }
    return *this;
  }

  gen galois_field::operator [] (const gen & g) {
    if (g.type==_INT_){
      int i= g.val;
      if (xcas_mode(context0)) --i;
      switch (i){
      case 0:
	return p;
      case 1:
	return P;
      case 2:
	return x;
      case 3:
	return a;
      }
    }
    return undef;
  }

  gen galois_field::operator >(const gen & g) const {
    if (g.type!=_USER)
      return undef;
    galois_field * gf=dynamic_cast<galois_field *>(g._USERptr);
    if (!gf)
      return undef;
    return is_strictly_positive(p-gf->p,context0); // ok
  }

  gen galois_field::operator <(const gen & g) const {
    if (g.type!=_USER)
      return undef;
    galois_field * gf=dynamic_cast<galois_field *>(g._USERptr);
    if (!gf)
      return undef;
    return is_strictly_positive(gf->p-p,context0); // ok
  }

  gen galois_field::operator <=(const gen & g) const {
    if (g.type!=_USER)
      return undef;
    galois_field * gf=dynamic_cast<galois_field *>(g._USERptr);
    if (!gf)
      return undef;
    return is_positive(gf->p-p,context0); // ok
  }

  gen galois_field::operator >=(const gen & g) const {
    if (g.type!=_USER)
      return undef;
    galois_field * gf=dynamic_cast<galois_field *>(g._USERptr);
    if (!gf)
      return undef;
    return is_positive(p-gf->p,0);
  }

  polynome galois_field::poly_reduce(const polynome & q) const {
    polynome res(q.dim);
    vector< monomial<gen> >::const_iterator it=q.coord.begin(),itend=q.coord.end();
    for (;it!=itend;++it){
      gen g=it->value;
      if (is_integer(g)){
	gen tmp=smod(g,p);
	if (giac::is_zero(tmp))
	  continue;
	g=galois_field(p,P,x,g*a);
      }
      res.coord.push_back(monomial<gen>(g,it->index));
    }
    return res;
  }

  gen galois_field::polygcd(const polynome & p,const polynome & q,polynome & res) const {
    res=Tgcdpsr(poly_reduce(p),poly_reduce(q));
    if (!res.coord.empty())
      res=res/res.coord.front().value;
    return 0;
  }

  gen galois_field::makegen(int i) const {
    if (P.type!=_VECT || p.type!=_INT_)
      return gendimerr();
    unsigned n=P._VECTptr->size()-1;
    //    i += pow(p,int(n)).val;
    vecteur res;
    for (unsigned j=0;j<n;++j){
      if (!i)
	break;
      res.push_back(gen(i%p.val));
      i=i/p.val;
    }
    reverse(res.begin(),res.end());
    return galois_field(p,P,x,res);
  }

  gen galois_field::polyfactor (const polynome & p0,factorization & f) const {
    f.clear();
    if (p0.coord.empty())
      return 0;
    polynome p=p0.coord.front().value.inverse(context0)*p0;
    if (p.dim!=1)
      return gendimerr(gettext("Multivariate GF factorization not yet implemented"));
    if (P.type!=_VECT)
      return gensizeerr(gettext("GF polyfactor"));
    environment env;
    env.moduloon=false;
    env.coeff=*this;
    env.modulo=this->p.to_int();
    int exposant=int(this->P._VECTptr->size())-1;
    env.pn=giac::pow(this->p,exposant);
    factorization sqff_f(squarefree_fp(p,env.modulo.val,exposant));
    if (!sqff_ffield_factor(sqff_f,env.modulo.val,&env,f))
      return gensizeerr(gettext("GF polyfactor"));
    f.push_back(facteur<polynome>(
				  polynome(
					   monomial<gen>(p0.coord.front().value,0,p.dim)
					   ),
				  1));
    return 0;
  }

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
