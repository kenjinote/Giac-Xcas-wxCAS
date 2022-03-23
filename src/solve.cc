// -*- mode:C++ ; compile-command: "g++-3.4 -I.. -g -c -DHAVE_CONFIG_H -DIN_GIAC solve.cc" -*-
#include "giacPCH.h"

/*
 *  Copyright (C) 2001,7 B. Parisse, R. De Graeve
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
#include <stdexcept>
#include <cmath>
#include "gen.h"
#include "solve.h"
#include "modpoly.h"
#include "unary.h"
#include "symbolic.h"
#include "usual.h"
#include "sym2poly.h"
#include "subst.h"
#include "derive.h"
#include "plot.h"
#include "prog.h"
#include "series.h"
#include "alg_ext.h"
#include "intg.h"
#include "rpn.h"
#include "lin.h"
#include "misc.h"
#include "cocoa.h"
#include "ti89.h"
#include "maple.h"
#include "giacintl.h"
#ifdef HAVE_LIBGSL
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  // FIXME intvar_counter should be contextized
  static int intvar_counter=0;
  static int realvar_counter=0;
  string print_intvar_counter(GIAC_CONTEXT){
    if (intvar_counter<0)
      return print_INT_(-intvar_counter);
    string res=print_INT_(intvar_counter);
    ++intvar_counter;
    return res;
  }

  string print_realvar_counter(GIAC_CONTEXT){
    if (realvar_counter<0)
      return print_INT_(int(-realvar_counter));
    string res=print_INT_(int(realvar_counter));
    ++realvar_counter;
    return res;
  }

  gen _reset_solve_counter(const gen & args,const context * contextptr){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (is_zero(args)){
      intvar_counter=0;
      return 1;
    }
    if (is_one(args)){
      realvar_counter=0;
      return 1;
    }
    if (args.type==_VECT && args._VECTptr->size()==2){
      intvar_counter=int(evalf_double(args._VECTptr->front(),1,contextptr)._DOUBLE_val);
      realvar_counter=int(evalf_double(args._VECTptr->back(),1,contextptr)._DOUBLE_val);      
    }
    else {
      intvar_counter=0;
      realvar_counter=0;
    }
    return 1;
  }
  static const char _reset_solve_counter_s []="reset_solve_counter";
  static define_unary_function_eval (___reset_solve_counter,&_reset_solve_counter,_reset_solve_counter_s);
  define_unary_function_ptr5( at_reset_solve_counter ,alias_at_reset_solve_counter,&___reset_solve_counter,0,true);

  void set_merge(vecteur & v,const vecteur & w){
    if (is_undef(w)){
      v=w;
      return;
    }
    const_iterateur it=w.begin(),itend=w.end();
    for (;it!=itend;++it)
      if (!equalposcomp(v,*it))
	v.push_back(*it);
  }

  static gen one_tour(GIAC_CONTEXT){
    if (angle_radian(contextptr)) 
      return cst_two_pi;
    else
      return 360;
  }
  static gen one_half_tour(GIAC_CONTEXT){
    if (angle_radian(contextptr)) 
      return cst_pi;
    else
      return 180;
  }
  static gen isolate_exp(const gen & e,int isolate_mode,GIAC_CONTEXT){
    if (isolate_mode &1)
      return ln(e,contextptr);
    if (e.type!=_VECT){
      if (is_strictly_positive(-e,contextptr))
	return vecteur(0);
      else
	return ln(e,contextptr);
    }
    // check in real mode for negative ln
    const_iterateur it=e._VECTptr->begin(),itend=e._VECTptr->end();
    vecteur res;
    for (;it!=itend;++it){
      if (!is_strictly_positive(-*it,contextptr))
	res.push_back(ln(*it,contextptr));
    }
    return res;
  }
  static gen isolate_ln(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return simplify(exp(e,contextptr),contextptr);
  }
  static gen isolate_sin(const gen & e,int isolate_mode,GIAC_CONTEXT){
    gen asine=asin(e,contextptr);
    if (!(isolate_mode & 2))
      return makevecteur(asine,one_half_tour(contextptr)-asine);
    identificateur * x=new identificateur(string("n_")+print_intvar_counter(contextptr));
    return makevecteur(asine+(*x)*one_tour(contextptr),one_half_tour(contextptr)-asine+(*x)*one_tour(contextptr));
  }
  static gen isolate_cos(const gen & e,int isolate_mode,GIAC_CONTEXT){
    gen acose=acos(e,contextptr);
    if (!(isolate_mode & 2))
      return makevecteur(acose,-acose);
    identificateur * x=new identificateur(string("n_")+print_intvar_counter(contextptr));
    return makevecteur(acose+(*x)*one_tour(contextptr),-acose+(*x)*one_tour(contextptr));
  }
  static gen isolate_tan(const gen & e,int isolate_mode,GIAC_CONTEXT){
    if (!(isolate_mode & 2))
      return atan(e,contextptr);
    identificateur * x=new identificateur(string("n_")+print_intvar_counter(contextptr));
    return atan(e,contextptr)+(*x)*one_half_tour(contextptr);
  }
  static gen isolate_asin(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return sin(e,contextptr);
  }
  static gen isolate_acos(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return cos(e,contextptr);
  }
  static gen isolate_atan(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return tan(e,contextptr);
  }

  static gen isolate_asinh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return sinh(e,contextptr);
  }
  static gen isolate_acosh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return cosh(e,contextptr);
  }
  static gen isolate_atanh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return tanh(e,contextptr);
  }

  static gen isolate_sinh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    gen asine= asinh(e,contextptr);
    if (!(isolate_mode & 2))
      return asine;
    identificateur * x=new identificateur(string("n_")+print_intvar_counter(contextptr));
    return makevecteur(asine+(*x)*one_tour(contextptr)*cst_i,(one_half_tour(contextptr)+(*x)*one_tour(contextptr))*cst_i-asine);
  }
  static gen isolate_cosh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    gen acose=acosh(e,contextptr);
    if (!(isolate_mode & 2))
      return makevecteur(acose,-acose);
    identificateur * x=new identificateur(string("n_")+print_intvar_counter(contextptr));
    return makevecteur(acose+(*x)*one_tour(contextptr)*cst_i,-acose+(*x)*one_tour(contextptr)*cst_i);
  }
  static gen isolate_tanh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    if (!(isolate_mode & 2))
      return atanh(e,contextptr);
    identificateur * x=new identificateur(string("n_")+print_intvar_counter(contextptr));
    return atanh(e,contextptr)+(*x)*one_half_tour(contextptr)*cst_i;
  }

  static gen (* const isolate_fcns[] ) (const gen &,int,GIAC_CONTEXT) = { isolate_exp,isolate_ln,isolate_sin,isolate_cos,isolate_tan,isolate_asin,isolate_acos,isolate_atan,isolate_sinh,isolate_cosh,isolate_tanh,isolate_asinh,isolate_acosh,isolate_atanh};

  static vecteur find_excluded(const gen & g,GIAC_CONTEXT){
    if (g.type!=_IDNT)
      return vecteur(0);
    gen g2=g._IDNTptr->eval(eval_level(contextptr),g,contextptr);
    if ((g2.type==_VECT) && (g2.subtype==_ASSUME__VECT)){
      vecteur v=*g2._VECTptr;
      if ( (v.size()==3) && (v.front()==_DOUBLE_) && (v[2].type==_VECT)){
	return *v[2]._VECTptr;
      }
    }
    return vecteur(0);
  }

  vecteur protect_sort(const vecteur & res,GIAC_CONTEXT){
#ifndef NO_STDEXCEPT
    try {
#endif
      gen tmp=_sort(res,contextptr);
      if (tmp.type==_VECT){
	vecteur w=*tmp._VECTptr,res;
	iterateur it=w.begin(),itend=w.end();
	for (;it!=itend;++it){
	  if (res.empty() || *it!=res.back())
	    res.push_back(*it);
	}
	return res;
      }
#ifndef NO_STDEXCEPT
    }
    catch (std::runtime_error & e){
      cerr << e.what() << endl;
    }
#endif
    return res;
  }

  vecteur find_singularities(const gen & e,const identificateur & x,int cplxmode,GIAC_CONTEXT){
    vecteur lv(lvarxpow(e,x));
    vecteur res;
    vecteur l(lvar(e));
    gen p=e2r(e,l,contextptr),n,d;
    vecteur pv=gen2vecteur(p);
    const_iterateur jt=pv.begin(),jtend=pv.end();
    for (;jt!=jtend;++jt){
      fxnd(*jt,n,d);
      if (d.type==_POLY)
	res=solve(r2e(d,l,contextptr),x,cplxmode,contextptr);
    }
    const_iterateur it=lv.begin(),itend=lv.end();
    for (;it!=itend;++it){
      if (it->type!=_SYMB)
	continue;
      const unary_function_ptr & u=it->_SYMBptr->sommet;
      gen & f=it->_SYMBptr->feuille;
      res=mergevecteur(res,find_singularities(f,x,cplxmode,contextptr));
      if (u==at_ln || u==at_sign)
	res=mergevecteur(res,solve(f,x,cplxmode,contextptr));
      if (u==at_pow && f.type==_VECT && f._VECTptr->size()==2)
	res=mergevecteur(res,solve(f._VECTptr->front(),x,cplxmode,contextptr));	
      if (u==at_tan)
	res=mergevecteur(res,solve(cos(f,contextptr),x,cplxmode,contextptr));
      if (u==at_piecewise && f.type==_VECT){
	vecteur & v = *f._VECTptr;
	int s=v.size();
	for (int i=0;i<s-1;i+=2){
	  gen & e =v[i];
	  if (e.is_symb_of_sommet(at_superieur_strict) ||
	      e.is_symb_of_sommet(at_inferieur_strict) ||
	      e.is_symb_of_sommet(at_superieur_egal) ||
	      e.is_symb_of_sommet(at_inferieur_egal)){
	    vecteur tmp=solve(e._SYMBptr->feuille._VECTptr->front()-e._SYMBptr->feuille._VECTptr->back(),x,cplxmode,contextptr);
	    // is *it continuous at tmp
	    const_iterateur jt=tmp.begin(),jtend=tmp.end();
	    for (;jt!=jtend;++jt){
	      if (!is_zero(limit(*it,x,*jt,1,contextptr)-limit(*it,x,*jt,-1,contextptr)))
		res.push_back(*jt);
	    }
	  }
	}
      }
    }
    if (cplxmode)
      return res;
    return protect_sort(res,contextptr);
  }

  static void solve_ckrange(const identificateur & x,vecteur & v,int isolate_mode,GIAC_CONTEXT){
    vecteur w,excluded(find_excluded(x,contextptr));
    // assumption on x, either range or integer
    if (find_range(x,w,contextptr)==2){
      int s=v.size();
      for (int i=0;i<s;++i){
	if (is_integer(v[i]))
	  w.push_back(v[i]);
      }
      v=w;
      return;
    }
    if (w.size()!=1 || w.front().type!=_VECT)
      return;
    w=*w.front()._VECTptr;
    if (w.size()!=2)
      return;
    gen l(w.front()),m(w.back());
    vecteur newv;
    iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      *it=simplifier(*it,contextptr);
      if (equalposcomp(excluded,*it))
	continue;
      gen sol=*it;
      if (l!=minus_inf && sign(l-sol,contextptr)==1)
	continue;
      if (m!=plus_inf && sign(sol-m,contextptr)==1)
	continue;
      sol=evalf(sol,eval_level(contextptr),contextptr);
      if (sol.type==_CPLX && !(isolate_mode &1) && !is_zero(im(sol,contextptr)))
	continue;
      if (sol.type!=_DOUBLE_){ // check for trig solutions
	newv.push_back(*it);
	vecteur lv(lvar(sol));
	if (lv.size()!=1 || l==minus_inf || m==plus_inf)
	  continue;
	gen n(lv.front()),a,b;
	// check linearity
	if (n.type!=_IDNT || n.print(contextptr).substr(0,2)!="n_" || !is_linear_wrt(*it,n,a,b,contextptr))
	  continue;
	newv.pop_back();
	a=normal(a,contextptr);
	b=normal(b,contextptr);
	if (!is_positive(a,contextptr))
	  swapgen(l,m);
	int n0(_ceil(evalf((l-b)/a,eval_level(contextptr),contextptr),contextptr).val),n1(_floor(evalf((m-b)/a,eval_level(contextptr),contextptr),contextptr).val);
	for (;n0<=n1;++n0)
	  newv.push_back(n0*a+b);
      }
      else {
	if (is_strictly_greater(l,sol,contextptr))
	  continue;
	if (is_strictly_greater(sol,m,contextptr))
	  continue;
	newv.push_back(*it);
      }
    }
    v=newv;
  }

  // Helper for the solver, make a translation using x^(n-1) coeff
  // and find gcd of deg, return true if non-trivial gcd found
  static bool translate_gcddeg(const vecteur & v,vecteur & v_translated, gen & x_translation,int & gcddeg){
    int s=v.size();
    if (s<4)
      return false;
    x_translation=-v[1]/((s-1)*v[0]);
    v_translated=taylor(v,x_translation,0);
    gcddeg=0;
    for (int i=1;i<s;++i){
      if (!is_zero(v_translated[i]))
	gcddeg=gcd(gcddeg,i);
    }
    if (gcddeg<=1)
      return false;
    int newdeg=(s-1)/gcddeg+1;
    // compress v_translated, keep only terms with index multiple of gcddeg
    for (int i=1;i<newdeg;++i){
      v_translated[i]=v_translated[i*gcddeg];
    }
    v_translated=vecteur(v_translated.begin(),v_translated.begin()+newdeg);
    return true;
  }

  static vecteur solve_inequation(const gen & e0,const identificateur & x,int direction,GIAC_CONTEXT);

  static vecteur solve_piecewise(const gen & args,const gen & value,const identificateur & x,int isolate_mode,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return vecteur(1,gensizeerr(contextptr));
    vecteur & piece_args=*args._VECTptr;
    vecteur failtest; // all these tests must fail to keep solution
    gen successtest,equation; // this test must succeed
    int s=piece_args.size();
    vecteur res;
    for (int i=0;i<s;i+=2){
      if (i)
	failtest.push_back(successtest);
      if (i+1==s){
	successtest=1;
	equation=piece_args[i];
      }
      else {
	successtest=piece_args[i];
	equation=piece_args[i+1];
      }
      int fails=failtest.size();
      vecteur sol=solve(equation-value,x,isolate_mode,contextptr);
      // now test whether solutions in sol are acceptable
      const_iterateur it=sol.begin(),itend=sol.end();
      for (;it!=itend;++it){
	const gen & g=*it;
	if (g==x){
	  if (fails){
	    gen tmp=symb_not(symbolic(at_ou,gen(failtest,_SEQ__VECT)));
	    res.push_back(is_one(successtest)?tmp:symb_and(tmp,successtest));
	  }
	  else
	    res.push_back(is_one(successtest)?g:successtest);
	  continue;
	}
	if (!is_zero(derive(g,x,contextptr))){
	  if (fails){
	    gen tmp=symb_not(symbolic(at_ou,gen(failtest,_SEQ__VECT)));
	    tmp=is_one(successtest)?tmp:symb_and(tmp,successtest);
	    res.push_back(symb_and(tmp,g));
	  }
	  else
	    res.push_back(is_one(successtest)?g:symb_and(successtest,g));
	  continue;
	}
	int j;
	for (j=0;j<fails;++j){
	  if (is_one(subst(failtest[j],x,g,false,contextptr)))
	    break;
	}
	if (j==fails && is_one(subst(successtest,x,g,false,contextptr)))
	  res.push_back(g);
      }
    }
    return res;
  }

  // inner solver
  void in_solve(const gen & e,const identificateur & x,vecteur &v,int isolate_mode,GIAC_CONTEXT){
    if (has_op(e,*at_equal)){
      v=vecteur(1,gensizeerr(gettext("Bad equal in")+e.print(contextptr)));
      return;
    }
    bool complexmode=isolate_mode & 1;
    vecteur lv(lvarx(e,x));
    int s=lv.size();
    if (!s)
      return;
    if (s>1){
      for (int i=0;i<s;++i){
	gen xvar=lv[i];
	if (xvar._SYMBptr->sommet==at_sign){
	  gen new_e=subst(e,xvar,1,false,contextptr);
	  vecteur vplus;
	  in_solve(new_e,x,vplus,isolate_mode,contextptr);
	  const_iterateur it=vplus.begin(),itend=vplus.end();
	  for (;it!=itend;++it){
	    if (is_one(subst(xvar,x,*it,false,contextptr)))
	      v.push_back(*it);
	  }
	  new_e=subst(e,xvar,-1,false,contextptr);
	  vecteur vminus;
	  in_solve(new_e,x,vminus,isolate_mode,contextptr);
	  it=vminus.begin();
	  itend=vminus.end();
	  for (;it!=itend;++it){
	    if (is_one(-subst(xvar,x,*it,false,contextptr)))
	      v.push_back(*it);
	  }
	  return;
	}
      }
#ifndef NO_STDEXCEPT
      throw(std::runtime_error("Unable to isolate "+string(x.print(contextptr))+" in "+e.print(contextptr)));
#endif
      v=vecteur(1,undeferr(gettext("Unable to isolate ")+string(x.print(contextptr))+" in "+e.print(contextptr)));
      return;
    }
    gen xvar(lv.front());
    if (xvar!=x){ // xvar must be a unary function of x
      if (xvar.type!=_SYMB){
	v=vecteur(1,gentypeerr(contextptr));
	return;
      }
      if (xvar._SYMBptr->sommet!=at_piecewise && xvar._SYMBptr->feuille.type==_VECT){
#ifndef NO_STDEXCEPT
	throw(std::runtime_error("Unable to isolate "+string(x.print(contextptr))+" in "+xvar.print(contextptr)));
#endif
	v=vecteur(1,undeferr(gettext("Unable to isolate ")+string(x.print(contextptr))+" in "+xvar.print(contextptr)));
	return;
      }
      if (xvar._SYMBptr->sommet==at_sign){
	gen new_e=subst(e,xvar,1,false,contextptr);
	if (is_zero(new_e)){
	  v=solve_inequation(symbolic(at_superieur_strict,makevecteur(xvar._SYMBptr->feuille,0)),x,1,contextptr);
	}
	else {
	  new_e=subst(e,xvar,-1,false,contextptr);
	  if (is_zero(new_e)){
	    v=solve_inequation(symbolic(at_inferieur_strict,makevecteur(xvar._SYMBptr->feuille,0)),x,-1,contextptr);
	  }
	}
	return;
      }
      int pos=equalposcomp(solve_fcns_tab,xvar._SYMBptr->sommet);
      if (xvar._SYMBptr->sommet==at_piecewise)
	pos=-1;
      if (!pos){
#ifndef NO_STDEXCEPT
	throw(std::runtime_error(string("Unable to isolate function ")+xvar._SYMBptr->sommet.ptr()->print(contextptr)));
#endif
	v=vecteur(1,undeferr(string("Unable to isolate function ")+xvar._SYMBptr->sommet.ptr()->print(contextptr)));
	return;
      }
      // solve with respect to xvar
      identificateur localt(" t");
      // ck_parameter_t();
      gen new_e=subst(e,xvar,localt,false,contextptr);
      vecteur new_v=solve(new_e,localt,isolate_mode,contextptr);
      const_iterateur it=new_v.begin(),itend=new_v.end();
      for (;it!=itend;++it){
	if (pos==-1){
	  // solve piecewise()==*it
	  set_merge(v,solve_piecewise(xvar._SYMBptr->feuille,*it,x,isolate_mode,contextptr));
	  if (is_undef(v)) return;
	  continue;
	}
	gen res=isolate_fcns[pos-1](*it,isolate_mode,contextptr);
	if (res.type!=_VECT)
	  set_merge(v,solve(xvar._SYMBptr->feuille-res,x,isolate_mode,contextptr));
	else {
	  const_iterateur it=res._VECTptr->begin(),itend=res._VECTptr->end();
	  for (;it!=itend;++it)
	    set_merge(v,solve(xvar._SYMBptr->feuille-*it,x,isolate_mode,contextptr));
	}
      }
      solve_ckrange(x,v,isolate_mode,contextptr);
      return;
    } // end xvar!=x
    // rewrite e as a univariate polynomial, first add other vars to x
    vecteur newv;
    lv=vecteur(1,vecteur(1,x));
    alg_lvar(e,lv);
    vecteur lvrat(1,x);
    lvar(e,lvrat);
    if (lvrat==lv.front())
      lv=lvrat;
    // int lv_size=lv.size();
    gen num,den,f;
    f=e2r(e,lv,contextptr);
    fxnd(f,num,den);
    if (num.type!=_POLY)
      return;
    vecteur w=polynome2poly1(*num._POLYptr,1);
    if (lv.front().type==_VECT){
      lv.front()=vecteur(lv.front()._VECTptr->begin()+1,lv.front()._VECTptr->end());
      if (lv.front()._VECTptr->empty())
	lv.erase(lv.begin()); // remove x from lv (CDR_VECT)	
    }
    else
      lv.erase(lv.begin()); // remove x from lv (CDR_VECT)
    int deg;
    vecteur w_translated;
    gen delta_x;
    if (translate_gcddeg(w,w_translated,delta_x,deg)){
      // composite polynomials
      gen invdeg=inv(deg,contextptr);
      gen newe=symb_horner(*r2sym(w_translated,lv,contextptr)._VECTptr,x);
      delta_x=r2sym(delta_x,lv,contextptr);
      vecteur vtmp;
      in_solve(newe,x,vtmp,isolate_mode,contextptr);
      vecteur unitroot(1,plus_one),munitroot;
      if (complexmode){
	for (int k=1;k<deg;++k)
	  unitroot.push_back(exp(2*k*cst_pi/deg*cst_i,contextptr));
	for (int k=0;k<deg;++k)
	  munitroot.push_back(exp((1+2*k)*cst_pi/deg*cst_i,contextptr));
      }
      const_iterateur it=vtmp.begin(),itend=vtmp.end();
      for (;it!=itend;++it){
	bool negatif=is_strictly_positive(-*it,contextptr);
	gen tmp=pow((negatif?-*it:*it),invdeg,contextptr);
	if (complexmode){
	  const_iterateur jt,jtend;
	  if (!negatif){
	    jt=unitroot.begin();
	    jtend=unitroot.end();
	  }
	  else {
	    jt=munitroot.begin();
	    jtend=munitroot.end();
	  }
	  for (;jt!=jtend;++jt)
	    newv.push_back(delta_x + (*jt) * tmp);
	}
	else {
	  if (deg%2)
	    newv.push_back(delta_x + (negatif?-tmp:tmp));
	  else {
	    if (!negatif){
	      newv.push_back(delta_x + tmp);
	      newv.push_back(delta_x - tmp);
	    }
	  }
	}
      }
      solve_ckrange(x,newv,isolate_mode,contextptr);
      v=mergevecteur(v,newv);
      return;
    }
    // if degree(w)=0, 1 or 2 solve it, otherwise error (should return ext)
    int d=w.size()-1;
    if (!d)
      return;
    if (d==1){
      gen tmp=rdiv(-r2sym(w.back(),lv,contextptr),r2sym(w.front(),lv,contextptr));
      if (!complexmode && has_i(tmp))
	return;
      newv.push_back(tmp);
      solve_ckrange(x,newv,isolate_mode,contextptr);
      v=mergevecteur(v,newv);
      return;
    }
    if (d>2){
      if (has_num_coeff(w)){
	if (complexmode)
	  newv=proot(w,epsilon(contextptr));
	else
	  newv=real_proot(w,epsilon(contextptr),contextptr);
	solve_ckrange(x,newv,isolate_mode,contextptr);
	v=mergevecteur(v,newv);
	return;
      }
      int n=is_cyclotomic(w,epsilon(contextptr));
      if (!n){
	*logptr(contextptr) << "Warning! Algebraic extension not implemented yet for poly " << r2sym(w,lv,contextptr) << endl;
	w=*evalf(w,1,contextptr)._VECTptr;
	if (has_num_coeff(w)){
#ifndef NO_STDEXCEPT
	  try {
#endif
	    if (complexmode)
	      newv=proot(w,epsilon(contextptr));
	    else
	      newv=real_proot(w,epsilon(contextptr),contextptr);
	    solve_ckrange(x,newv,isolate_mode,contextptr);
	    v=mergevecteur(v,newv);
#ifndef NO_STDEXCEPT
	  }
	  catch (std::runtime_error & ){
	  }
#endif
	  return;
	}
	return;
      }
      if (complexmode){
	for (int j=1;j<=n/2;++j){
	  if (gcd(j,n)==1){
	    if (n%2){
	      newv.push_back(exp(rdiv(gen(2*j)*cst_pi*cst_i,n),contextptr));
	      newv.push_back(exp(rdiv(gen(-2*j)*cst_pi*cst_i,n),contextptr));
	    }
	    else {
	      newv.push_back(exp(rdiv(gen(j)*cst_pi*cst_i,n/2),contextptr));
	      newv.push_back(exp(rdiv(gen(-j)*cst_pi*cst_i,n/2),contextptr));
	    }
	  }
	}
      }
      solve_ckrange(x,newv,isolate_mode,contextptr);
      v=mergevecteur(v,newv);
      return ;
    }
    gen b_over_2=rdiv(w[1],plus_two);
    if (b_over_2.type!=_FRAC){
      gen a=r2sym(w.front(),lv,contextptr);
      gen minus_b_over_2=r2sym(-b_over_2,lv,contextptr);
      gen delta_prime=r2sym(pow(b_over_2,2,contextptr)-w.front()*w.back(),lv,contextptr);
      if (!complexmode && is_strictly_positive(-delta_prime,contextptr))
	return;
      newv.push_back(rdiv(minus_b_over_2+sqrt(delta_prime,contextptr),a));
      if (!is_zero(delta_prime))
	newv.push_back(rdiv(minus_b_over_2-sqrt(delta_prime,contextptr),a));
    }
    else {
      gen two_a=r2sym(plus_two*w.front(),lv,contextptr);
      gen minus_b=r2sym(-w[1],lv,contextptr);
      gen delta=r2sym(w[1]*w[1]-gen(4)*w.front()*w.back(),lv,contextptr);
      if (!complexmode && is_positive(-delta,contextptr))
	return;
      newv.push_back(rdiv(minus_b+sqrt(delta,contextptr),two_a));
      newv.push_back(rdiv(minus_b-sqrt(delta,contextptr),two_a));
    }
    solve_ckrange(x,newv,isolate_mode,contextptr);
    v=mergevecteur(v,newv);
  }

  // v assumed to represent an irreducible dense 1-d poly
  vecteur solve(const vecteur & v,bool complexmode,GIAC_CONTEXT){
    vecteur res;
    int d=v.size()-1;
    if (d<1)
      return res;
    if (d==1){
      res.push_back(rdiv(-v.back(),v.front()));
      return res;
    }
    if (!is_one(v.front())){
      // if v is not monic, set Y=a*X
      gen a(v.front()),puissance(plus_one);
      vecteur w;
      w.reserve(d+1);
      for (int i=0;i<=d;++i,puissance=puissance*a)
	w.push_back(v[i]*puissance);
      return divvecteur(solve(divvecteur(w,a),complex_mode(contextptr),contextptr),a);
    }
    // should call sym2rxroot for extensions of extensions
    vecteur tmp(2,zero);
    tmp.front()=plus_one;
    if (d==2){
      gen b(v[1]),c(v[2]);
      gen bprime(rdiv(b,plus_two));
      if (!has_denominator(bprime)){
	gen delta(bprime*bprime-c);
	if (!complexmode && is_positive(-delta,contextptr))
	  return res;
	vecteur w(3,zero);
	w.front()=plus_one;
	w.back()=-delta;
	tmp.back()=-bprime;
	res.push_back(algebraic_EXTension(tmp,w));
	tmp.front()=minus_one;
	tmp.back()=-bprime;
	res.push_back(algebraic_EXTension(tmp,w));	
      }
      else {
	if (!complexmode && is_positive(4*c-b*b,contextptr))
	  return res;
	tmp.back()=zero;
	res.push_back(algebraic_EXTension(tmp,v));
	tmp.front()=minus_one;
	tmp.back()=-b;
	res.push_back(algebraic_EXTension(tmp,v));
      }
      return res;
    }
    // should return a list of d algebraic extension with order number
    res.push_back(algebraic_EXTension(tmp,v));
    return res;
  }

  static vecteur in_solve_inequation(const gen & e0,const gen &e,const identificateur & x,int direction,const gen & rangeg,const vecteur & veq_excluded,const vecteur & veq_not_singu,const vecteur & excluded_not_singu,const vecteur & singu,GIAC_CONTEXT){
    if (rangeg.type!=_VECT)
      return vecteur(0);
    vecteur rangev = *rangeg._VECTptr,range=rangev;
    if (rangev.size()==2){
      gen &a=rangev.front();
      gen & b=rangev.back();
      // keep only values inside a,b
      range=vecteur(1,a);
      const_iterateur it=veq_excluded.begin(),itend=veq_excluded.end();
      for (;it!=itend;++it){
	if (is_strictly_greater(*it,a,contextptr))
	  break;
      }
      for (;it!=itend;++it){
	if (is_greater(*it,b,contextptr))
	  break;
	range.push_back(*it);
      }
      range.push_back(b);
    }
    else {
      range=mergevecteur(rangev,veq_excluded);
      range=protect_sort(range,contextptr);
    }
    vecteur res;
    int s=range.size();
    if (s<2)
      return vecteur(1,gensizeerr(contextptr));
    if (s==2 && range[0]==minus_inf && range[1]==plus_inf){
      gen test=sign(subst(e,x,0,false,contextptr),contextptr);
      if (direction<0)
	test=-test;
      if (is_one(test))
	return vecteur(1,x);
      if (is_one(-test))
	return vecteur(0);
      return vecteur(1,gensizeerr(gettext("Unable to check sign ")+test.print()));
    }
    vecteur add_eq,already_added;
    for (int i=0;i<s-1;++i){
      gen l=range[i],m=range[i+1];
      if (l==m)
	continue;
      gen test;
      if (l==minus_inf)
	test=m-1;
      else {
	if (m==plus_inf)
	  test=l+1;
	else
	  test=(l+m)/2;
      }
      test=eval(subst(e0,x,test,false,contextptr),eval_level(contextptr),contextptr);
      if (test!=1){
	if (!equalposcomp(already_added,l) && equalposcomp(veq_not_singu,l))
	  add_eq.push_back(l);
	continue;
      }
      already_added.push_back(m);
      gen symb_sup,symb_inf;
      if (equalposcomp(excluded_not_singu,l) || equalposcomp(singu,l) ||
	  ( !(direction %2) && equalposcomp(veq_not_singu,l)) )
	symb_inf=symb_superieur_strict(x,l);
      else
	symb_inf=symb_superieur_egal(x,l);
      if (equalposcomp(excluded_not_singu,m) || equalposcomp(singu,m) ||
	  ( !(direction %2) && equalposcomp(veq_not_singu,m)) )
	symb_sup=symb_inferieur_strict(x,m);
      else
	symb_sup=symb_inferieur_egal(x,m);
      if (l==minus_inf)
	res.push_back(symb_sup);
      else {
	if (m==plus_inf)
	  res.push_back(symb_inf);
	else
	  res.push_back(symbolic(at_and,makevecteur(symb_inf,symb_sup))); 
      }
    }
    if (direction % 2)
      res=mergevecteur(add_eq,res);
    return res;
  }

  static bool ck_sorted(const vecteur & v,GIAC_CONTEXT){
    int vs=v.size();
    for (int i=1;i<vs;++i){
      if (!ck_is_greater(v[i],v[i-1],contextptr))
	return false;
    }
    return true;
  }

  // works for continuous functions only
  static vecteur solve_inequation(const gen & e0,const identificateur & x,int direction,GIAC_CONTEXT){
    if (has_num_coeff(e0))
      return vecteur(1,gensizeerr(gettext("Unable to solve inequations with approx coeffs ")+e0.print(contextptr)));
    gen e(e0._SYMBptr->feuille._VECTptr->front()-e0._SYMBptr->feuille._VECTptr->back());
    if (e.is_symb_of_sommet(at_superieur_strict) ||
	e.is_symb_of_sommet(at_inferieur_strict) ||
	e.is_symb_of_sommet(at_superieur_egal) ||
	e.is_symb_of_sommet(at_inferieur_egal))
      return vecteur(1,gensizeerr(gettext("Inequation inside inequation not implemented ")+e.print()));
    if (is_zero(ratnormal(derive(e,x,contextptr))))
      return vecteur(1,gensizeerr(gettext("Inequation is constant with respect to ")+string(x.print(contextptr))));
    vecteur veq_not_singu,veq,singu;
    singu=find_singularities(e,x,0,contextptr);
    veq_not_singu=solve(e,x,0,contextptr);
    veq=mergevecteur(veq_not_singu,singu);
    vecteur range,excluded_not_singu(find_excluded(x,contextptr));
    vecteur excluded=mergevecteur(excluded_not_singu,singu);
    vecteur veq_excluded=mergevecteur(excluded,veq);
    veq_excluded=protect_sort(veq_excluded,contextptr);
    if (!ck_sorted(veq_excluded,contextptr))
      return vecteur(1,gensizeerr(gettext("Unable to sort ")+gen(veq_excluded).print(contextptr)));
    // From the lower bound of range to the higher bound
    // find the sign 
    vecteur res;
    if (!find_range(x,range,contextptr))
      return res;
    for (unsigned i=0;i<range.size();i++){
      res=mergevecteur(res,in_solve_inequation(e0,e,x,direction,range[i],veq_excluded,veq_not_singu,excluded_not_singu,singu,contextptr));
    }
    return res;
  }

  static bool modsolve(const gen & e,const identificateur & x,const gen & modulo,vecteur &v,GIAC_CONTEXT){
    if (modulo.type!=_INT_)
      return false; // setdimerr(gettext("Modular equation with modulo too large"));
    int m=modulo.val;
    for (int i=0;i<m;++i){
      gen tmp=subst(e,x,i,false,contextptr);
      if (is_zero(tmp.eval(eval_level(contextptr),contextptr)))
	v.push_back(i);
    }
    return true;
  }
  
  static void clean(gen & e,const identificateur & x,GIAC_CONTEXT){
    if (e.type!=_SYMB)
      return;
    if (e._SYMBptr->sommet==at_inv || (e._SYMBptr->sommet==at_pow && ck_is_positive(-e._SYMBptr->feuille._VECTptr->back(),contextptr))){
      gen ef=e._SYMBptr->feuille;
      if (e._SYMBptr->sommet==at_pow)
	ef=ef._VECTptr->front();
      // search for a tan in the variables
      vecteur lv(lvarx(e,x));
      if (lv.size()!=1)
	e=1;
      else {
	gen xvar(lv.front());
	if (!xvar.is_symb_of_sommet(at_tan))
	  e=1;
      }
      return;
    }
    if (e._SYMBptr->sommet==at_prod){
      gen ef=e._SYMBptr->feuille;
      if (ef.type!=_VECT)
	return;
      vecteur v=*ef._VECTptr;
      int vs=v.size();
      for (int i=0;i<vs;++i)
	clean(v[i],x,contextptr);
      ef=gen(v,ef.subtype);
      e=symbolic(at_prod,ef);
    }
  }
  
  // detect product and powers
  void solve(const gen & e,const identificateur & x,vecteur &v,int isolate_mode,GIAC_CONTEXT){
    if (is_zero(e)){
      v.push_back(x);
      return;
    }
    switch (e.type){
    case _IDNT:
      if (*e._IDNTptr==x && !equalposcomp(find_excluded(x,contextptr),zero))
	addtolvar(zero,v);
      return;
    case _SYMB:
      if ( e._SYMBptr->sommet==at_pow && ck_is_strictly_positive(e._SYMBptr->feuille._VECTptr->back(),contextptr) ){
	vecteur tmpv;
	solve(e._SYMBptr->feuille._VECTptr->front(),x,tmpv,isolate_mode,contextptr);
	int ncopy=1;
	// make copies of the answer (xcas_mode(contextptr)==1 compatibility)
	if (xcas_mode(contextptr)==1 && e._SYMBptr->feuille._VECTptr->back().type==_INT_)
	  ncopy=e._SYMBptr->feuille._VECTptr->back().val;
	const_iterateur it=tmpv.begin(),itend=tmpv.end();
	for (;it!=itend;++it){
	  for (int i=0;i<ncopy;++i)
	    v.push_back(*it);
	}
	return;
      }
      if (e._SYMBptr->sommet==at_prod){
	const_iterateur it=e._SYMBptr->feuille._VECTptr->begin(),itend=e._SYMBptr->feuille._VECTptr->end();
	for (;it!=itend;++it){
	  solve(*it,x,v,isolate_mode,contextptr);
	  if (is_undef(v)) break;
	}
	return;
      }
      if (e._SYMBptr->sommet==at_neg){
	solve(e._SYMBptr->feuille,x,v,isolate_mode,contextptr);
	return;
      }
      if (!(isolate_mode & 2) && 
	  (e._SYMBptr->sommet==at_inv || (e._SYMBptr->sommet==at_pow && ck_is_positive(-e._SYMBptr->feuille._VECTptr->back(),contextptr)))
	  ){
	gen ef=e._SYMBptr->feuille;
	if (e._SYMBptr->sommet==at_pow)
	  ef=ef._VECTptr->front();
	// search for a tan in the variables
	vecteur lv(lvarx(e,x));
	if (lv.size()!=1)
	  return;
	gen xvar(lv.front());
	if (!xvar.is_symb_of_sommet(at_tan))
	  return;
	gen arg=xvar._SYMBptr->feuille;
	// solve arg=pi/2[pi]
	in_solve(arg-isolate_tan(plus_inf,isolate_mode,contextptr),x,v,isolate_mode,contextptr);
	return;
      }
      in_solve(e,x,v,isolate_mode,contextptr);
      break;
    default:
      return;
    }
  }

  // find the arguments of sqrt inside expression e
  static vecteur lvarfracpow(const gen & e){
    vecteur l0=lvar(e),l;
    const_iterateur it=l0.begin(),itend=l0.end();
    for (;it!=itend;++it){
      if (it->_SYMBptr->sommet!=at_pow)
	continue;
      vecteur & arg=*it->_SYMBptr->feuille._VECTptr;
      gen g=arg[1],expnum,expden;
      if (g.type==_FRAC){
	expnum=g._FRACptr->num;
	expden=g._FRACptr->den;
      }
      else {
	if ( (g.type!=_SYMB) || (g._SYMBptr->sommet!=at_prod) )
	  continue;
	gen & arg1=g._SYMBptr->feuille;
	if (arg1.type!=_VECT)
	  continue;
	vecteur & v=*arg1._VECTptr;
	if ( (v.size()!=2) || (v[1].type!=_SYMB) || (v[1]._SYMBptr->sommet==at_inv) )
	  continue;
	expnum=v[0];
	expden=v[1]._SYMBptr->feuille;
      }
      if (expden.type!=_INT_)
	continue;
      l.push_back(arg[0]);
      l.push_back(expden.val);
      l.push_back(*it);
    }
    return l;
  }

  static vecteur lvarfracpow(const gen & g,const identificateur & x,GIAC_CONTEXT){
    vecteur l0=lvarfracpow(g),l;
    const_iterateur it=l0.begin(),itend=l0.end();
    for (;it!=itend;++it){
      if (!is_zero(derive(*it,x,contextptr))){
	l.push_back(*it);
	++it;
	l.push_back(*it);
	++it;
	l.push_back(*it);
      }
      else
	it+=2;
    }
    return l;
  }

  static void solve_fracpow(const gen & e,const identificateur & x,const vecteur & eq,const vecteur & listvars,vecteur & fullres,int isolate_mode,GIAC_CONTEXT){
    vecteur equations(eq);
    if (e.type==_IDNT){
      if (*e._IDNTptr==x && !equalposcomp(find_excluded(x,contextptr),zero)){
	addtolvar(zero,fullres);
	return;
      }
    }
    if (e.type==_SYMB){
      if ( (e._SYMBptr->sommet==at_pow) && (ck_is_positive(e._SYMBptr->feuille._VECTptr->back(),contextptr)) ){
	solve_fracpow(e._SYMBptr->feuille._VECTptr->front(),x,equations,listvars,fullres,isolate_mode,contextptr);
	return;
      }
      if (e._SYMBptr->sommet==at_prod){
	const_iterateur it=e._SYMBptr->feuille._VECTptr->begin(),itend=e._SYMBptr->feuille._VECTptr->end();
	for (;it!=itend;++it)
	  solve_fracpow(*it,x,equations,listvars,fullres,isolate_mode,contextptr);
	return;
      }
      if (e._SYMBptr->sommet==at_neg){
	solve_fracpow(e._SYMBptr->feuille,x,equations,listvars,fullres,isolate_mode,contextptr);
	return;
      }
    } // end if (e.type==_SYMB)
    /*
      // new code with resultant in all var except the first one
      // disadvantage: does not check that listvars[i] are admissible
      // example assume(M<0); solve(sqrt(x)=M);
    gen expr(e);
    int s=listvars.size();
    for (int i=1;i<s;++i){
      // expr must be rationnal wrt listvars[i]
      vecteur vtmp(1,listvars[i]);
      if (listvars[i].type!=_IDNT)
	setsizeerr();
      rlvarx(expr,listvars[i],vtmp);
      // IMPROVE: maybe a function applied to expr is rationnal
      if (vtmp.size()!=1)
	setsizeerr(gettext("Solve with fractional power:")+expr.print(contextptr)+" is not rationnal w.r.t. "+listvars[i].print(contextptr));
      if (!is_zero(derive(expr,listvars[i],contextptr)))
	expr=_resultant(makevecteur(expr,equations[i-1],listvars[i]),contextptr);
    }
    expr=factor(expr,false,contextptr);
    if (is_zero(derive(expr,x,contextptr)))
      return;
    solve(expr,x,fullres,isolate_mode,contextptr);
    return;
    */
    // old code with Groebner basis
    equations.push_back(e);      
    vecteur res=gsolve(equations,listvars,complex_mode(contextptr),contextptr);
    iterateur it=res.begin(),itend=res.end();
    for (;it!=itend;++it)
      *it=(*it)[0];
    _purge(vecteur(listvars.begin()+1,listvars.end()),contextptr);
    if (listvars[0].type==_IDNT){
      fullres=mergevecteur(res,fullres);
      return;
    }
    // recursive call to solve composevar=*it with respect to x
    for (it=res.begin();it!=itend;++it){
      fullres=mergevecteur(fullres,solve(*it-listvars[0],x,isolate_mode,contextptr));
    }
  }

  static vecteur solve_cleaned(const gen & e,const identificateur & x,int isolate_mode,GIAC_CONTEXT){
    gen expr(e);
    if (expr.is_symb_of_sommet(at_prod)){
      vecteur v=gen2vecteur(expr._SYMBptr->feuille),res;
      for (unsigned i=0;i<v.size();++i){
	res=mergevecteur(res,solve_cleaned(v[i],x,isolate_mode,contextptr));
      }
      return res;
    }
    if (expr.is_symb_of_sommet(at_pow)){
      gen & f =expr._SYMBptr->feuille;
      if (f.type==_VECT && f._VECTptr->size()==2 && is_strictly_positive(f._VECTptr->back(),contextptr))
	return solve_cleaned(f._VECTptr->front(),x,isolate_mode,contextptr);
    }
    // Check for re/im/conj in complexmode
    bool complexmode=isolate_mode & 1;
    if (complexmode){
      vecteur lc=mergevecteur(lop(expr,at_conj),mergevecteur(lop(expr,at_re),lop(expr,at_im)));
      int s=lc.size();
      for (int i=0;i<s;++i){
	gen f=lc[i]._SYMBptr->feuille;
	if (!is_zero(derive(f,x,contextptr))){
	  identificateur xrei(" x"),ximi(" y");
	  gen xre(xrei),xim(ximi);
	  bool savec=complex_mode(contextptr);
	  bool savecv=complex_variables(contextptr);
	  complex_mode(false,contextptr);
	  complex_variables(false,contextptr);
	  gen tmp=subst(e,x,xre+cst_i*xim,false,contextptr);
	  vecteur res=gsolve(makevecteur(re(tmp,contextptr),im(tmp,contextptr)),makevecteur(xre,xim),false,contextptr);
	  complex_mode(savec,contextptr);
	  complex_variables(savecv,contextptr);
	  s=res.size();
	  for (int j=0;j<s;++j){
	    if (res[j].type==_VECT && res[j]._VECTptr->size()==2){
	      gen a=res[j]._VECTptr->front();
	      gen b=res[j]._VECTptr->back();
	      if (is_zero(a))
		res[j]=cst_i*b;
	      else {
		if (is_zero(b))
		  res[j]=a;
		else
		  res[j]=symbolic(at_plus,gen(makevecteur(a,cst_i*b),_SEQ__VECT));
	      }
	    }
	  }
	  return res;
	}
      }
    }
    if ( (approx_mode(contextptr) || has_num_coeff(e)) && lidnt(e)==makevecteur(x)){
      vecteur vtmp=makevecteur(e,x);
      vecteur res=gen2vecteur(in_fsolve(vtmp,contextptr));
      solve_ckrange(x,res,isolate_mode,contextptr);
      return res;
    }
    // should rewrite e in terms of a minimal number of variables
    // first factorization of e
    // Checking for abs
    vecteur la;
    if (!complexmode)
      la=lop(expr,at_abs);
    const_iterateur itla=la.begin(),itlaend=la.end();
    for (;itla!=itlaend;++itla){
      gen g=itla->_SYMBptr->feuille;
      if (is_zero(derive(g,x,contextptr)))
	continue;
      vecteur res;
      gen ee=subst(expr,*itla,g,false,contextptr);
      vecteur v1=solve(ee,x,isolate_mode,contextptr);
      const_iterateur it=v1.begin(),itend=v1.end();
      for (;it!=itend;++it){
	if (*it==x){
	  res=mergevecteur(res,solve_inequation(symbolic(at_superieur_egal,gen(makevecteur(g,0),_SEQ__VECT)),x,1,contextptr));
	  continue;
	}
	if (contains(*it,x)){
	  *logptr(contextptr) << "Warning, trying to solve " << g << ">=0 with " << *it << endl;
	  gen tmp=symbolic(at_solve,gen(makevecteur(symbolic(at_superieur_egal,gen(makevecteur(g,0),_SEQ__VECT)),x),_SEQ__VECT));
	  tmp=_tilocal(gen(makevecteur(tmp,*it),_SEQ__VECT),contextptr);
	  if (tmp.type==_VECT)
	    res=mergevecteur(res,*tmp._VECTptr);
	  continue;
	}
	gen g1=subst(g,x,*it,false,contextptr);
	if (normal(abs(g1,contextptr)-g1,contextptr)==0) // was ratnormal, but insufficient
	  res.push_back(*it);
      }
      ee=subst(expr,*itla,-g,false,contextptr);
      v1=solve(ee,x,isolate_mode,contextptr);
      it=v1.begin(); itend=v1.end();
      for (;it!=itend;++it){
	if (*it==x){
	  res=mergevecteur(res,solve_inequation(symbolic(at_inferieur_egal,gen(makevecteur(g,0),_SEQ__VECT)),x,1,contextptr));
	  continue;
	}
	if (contains(*it,x)){
	  *logptr(contextptr) << "Warning, trying to solve " << g << "<=0 with " << *it << endl;
	  gen tmp=symbolic(at_solve,gen(makevecteur(symbolic(at_inferieur_egal,gen(makevecteur(g,0),_SEQ__VECT)),x),_SEQ__VECT));
	  tmp=_tilocal(gen(tmp,_SEQ__VECT),contextptr);
	  if (tmp.type==_VECT)
	    res=mergevecteur(res,*tmp._VECTptr);
	  continue;
	}
	gen g1=subst(g,x,*it,false,contextptr);
	if (normal(abs(g1,contextptr)+g1,contextptr)==0) 
	  res.push_back(*it);
      }
      return res;
    }
    vecteur lv(lvarx(expr,x));
    int s=lv.size();
    if (s>1){
      expr=halftan_hyp2exp(expr,contextptr);
      lv=lvarx(expr,x);
      s=lv.size();
      if (s>1){
	gen tmp=_texpand(expr,contextptr);
	if (lvarx(tmp,x).size()==1){
	  expr=tmp;
	}
      }
    }
    // Checking for fractional power
    // Remark: algebraic extension could also be solved using resultant
    vecteur ls(lvarfracpow(expr,x,contextptr));
    if (!ls.empty()){ // Use auxiliary variables
      int s=ls.size()/3;
      vecteur substin,substout,equations,listvars(lvarx(expr,x,true));
      // remove ls from listvars, add aux var instead
      for (int i=0;i<s;++i){
	gen lsvar=ls[3*i+2];
	int j=equalposcomp(listvars,ls);
	if (j)
	  listvars.erase(listvars.begin()+j-1);
      }
      if (listvars.size()!=1)
	return vecteur(1,gensizeerr(gettext("unable to isolate ")+gen(listvars).print(contextptr)));
      for (int i=0;i<s;++i){
	gen lsvar=ls[3*i+2];
	substin.push_back(lsvar);
	gen tmp("c__"+print_INT_(i),contextptr);
	if (!(ls[3*i+1].val %2))
	  assumesymbolic(symb_superieur_egal(tmp,0),0,contextptr); 
	listvars.push_back(tmp);
	substout.push_back(tmp);
	equations.push_back(pow(tmp,ls[3*i+1],contextptr)-ls[3*i]);
      }
      gen expr1=subst(expr,substin,substout,false,contextptr);
      expr1=factor(expr1,false,contextptr);
      if (is_undef(expr1))
	return vecteur(1,expr1);
      vecteur fullres;
      solve_fracpow(expr1,x,equations,listvars,fullres,isolate_mode,contextptr);
      // Check that expr at x=fullres is 0
      // Only if expr1 does not depend on other variables than x
      vecteur othervar(1,x),res;
      lidnt(expr1,othervar);
      if (othervar.size()==listvars.size()){
	const_iterateur it=fullres.begin(),itend=fullres.end();
	for (;it!=itend;++it){
	  vecteur algv=alg_lvar(*it);
	  if ( (!algv.empty() && algv.front().type==_VECT && !algv.front()._VECTptr->empty()) || is_zero(limit(expr,x,*it,0,contextptr)))
	    res.push_back(*it);
	}
      }
      else
	res=fullres;
      return res;
    }
    lv=lvarx(expr,x);
    if (lv.size()>1){
      gen tmp=factor(simplify(expr,contextptr),false,contextptr);
      if (is_undef(tmp))
	return vecteur(1,tmp);
      int lvs=0;
      if (tmp.is_symb_of_sommet(at_prod) && tmp._SYMBptr->feuille.type==_VECT){
	vecteur & f=*tmp._SYMBptr->feuille._VECTptr;
	int fs=f.size();
	for (int i=0;i<fs;++i){
	  lvs=lvarx(f[i],x).size();
	  if (lvs>1)
	    break;
	}
      }
      else
	lvs=lvarx(tmp,x).size();
      if (lvs<2)
	expr=tmp;
    }
    // -> exp/ln
    expr=pow2expln(expr,x,contextptr);
    bool setcplx=complexmode && complex_mode(contextptr)==false;
    if (setcplx)
      complex_mode(true,contextptr);
    expr=factor(expr,false,contextptr); // factor in complex or real mode
    if (is_undef(expr))
      return vecteur(1,expr);
    if (setcplx)
      complex_mode(false,contextptr);
    lv=lvarx(expr,x);
    s=lv.size();
    vecteur v;
    if (!s){
      if (is_zero(expr))
	v.push_back(x);
      return v;
    }
    solve(expr,x,v,isolate_mode,contextptr);
    if (!(isolate_mode & 2) && !has_num_coeff(v)){
      // check solutions if there is a tan inside
      for (int i=0;i<s;++i){
	if (lv[i].is_symb_of_sommet(at_tan)){
	  vecteur res;
	  const_iterateur it=v.begin(),itend=v.end();
	  for (;it!=itend;++it){
	    if (is_zero(recursive_normal(limit(_tan2sincos2(expr,contextptr),x,*it,0,contextptr),contextptr)) || is_zero(recursive_normal(limit(expr,x,*it,0,contextptr),contextptr)))
	      res.push_back(*it);
	  }
	  return res;
	}
      }
    }
    return v;
  }

  vecteur solve(const gen & e,const identificateur & x,int isolate_mode,GIAC_CONTEXT){
    if (is_undef(e)) return vecteur(0);
    gen expr(e);
    gen modulo;
    if (has_mod_coeff(expr,modulo)){
      vecteur v;
      if (!modsolve(expr,x,modulo,v,contextptr))
	return vecteur(1,gensizeerr(gettext("Modulo too large")));
      return v;
    }
    // Inequation?
    if (e.type==_SYMB){ 
      if (e._SYMBptr->sommet==at_inferieur_strict)
	return solve_inequation(e,x,-2,contextptr);
      if (e._SYMBptr->sommet==at_inferieur_egal)
	return solve_inequation(e,x,-1,contextptr);
      if (e._SYMBptr->sommet==at_superieur_strict)
	return solve_inequation(e,x,2,contextptr);
      if (e._SYMBptr->sommet==at_superieur_egal)
	return solve_inequation(e,x,1,contextptr);
      if (e._SYMBptr->sommet==at_equal ||e._SYMBptr->sommet==at_same)
	expr = e._SYMBptr->feuille._VECTptr->front()-e._SYMBptr->feuille._VECTptr->back();
    }
    clean(expr,x,contextptr);
    return solve_cleaned(expr,x,isolate_mode,contextptr);
  }

  gen remove_and(const gen & g,const unary_function_ptr * u){
    if (g.type==_VECT){
      vecteur res;
      const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
      for (;it!=itend;++it){
	gen tmp=remove_and(*it,u);
	if (tmp.type!=_VECT){
	  tmp=remove_and(*it,at_and);
	  res.push_back(tmp);
	}
	else
	  res=mergevecteur(res,*tmp._VECTptr);
      }
      return res;
    }
    if (!g.is_symb_of_sommet(u))
      return g;
    return remove_and(g._SYMBptr->feuille,u);
  }

  vecteur solve(const gen & e,const gen & x,int isolate_mode,GIAC_CONTEXT){
    bool complexmode=isolate_mode & 1;
    vecteur res;
    if (x.type!=_IDNT){
      if (x.type==_VECT && x._VECTptr->size()==1 && e.type==_VECT && e._VECTptr->size()==1){
	vecteur res=solve(e._VECTptr->front(),x._VECTptr->front(),isolate_mode,contextptr);
	iterateur it=res.begin(),itend=res.end();
	for (;it!=itend;++it)
	  *it=vecteur(1,*it);	
	return res;
      }
      if ( (x.type==_VECT) && (e.type==_VECT) )
	return gsolve(*e._VECTptr,*x._VECTptr,complexmode,contextptr);
      identificateur xx("x");
      res=solve(subst(e,x,xx,false,contextptr),xx,isolate_mode,contextptr);
      return res;
    }
    if (e.type==_VECT){
      const_iterateur it=e._VECTptr->begin(),itend=e._VECTptr->end();
      gen curx=x._IDNTptr->eval(1,x,contextptr);
      res=vecteur(1,x); // everything is solution up to now
      for (;it!=itend;++it){
	if (res==vecteur(1,x))
	  res=solve(*it,*x._IDNTptr,isolate_mode,contextptr);
	else { // check every element of res
	  vecteur newres;
	  const_iterateur jt=res.begin(),jtend=res.end();
	  for (;jt!=jtend;++jt){
	    if (jt->is_symb_of_sommet(at_superieur_strict) ||
		jt->is_symb_of_sommet(at_inferieur_strict) ||
		jt->is_symb_of_sommet(at_superieur_egal) ||
		jt->is_symb_of_sommet(at_inferieur_egal) ||
		jt->is_symb_of_sommet(at_and)){
	      assumesymbolic(*jt,0,contextptr); // assume and solve next equation
	      newres=mergevecteur(newres,solve(*it,*x._IDNTptr,isolate_mode,contextptr));
	      _purge(x,contextptr);
	    }
	    else { 
	      if (is_zero(normal(subst(*it,x,*jt,true,contextptr),1,contextptr)))
		newres.push_back(*jt);
	    }
	  } // end for (;jt!=jtend;++jt) loop on previous solutions
	  res=newres;
	} // end else
      } // end for (;it!=itend;++it) loop on equations
      if (curx!=x)
	sto(curx,x,contextptr);
      return res;
    }
    else
      res=solve(e,*x._IDNTptr,isolate_mode,contextptr);
    return res;
  }

  static gen symb_solution(const gen & g,const gen & var,GIAC_CONTEXT){
    if (var.type!=_VECT){
      if (var.type==_IDNT && !lvarx(g,var).empty())
	return g;
      else
	return symbolic(at_equal,makevecteur(var,g));
    }
    vecteur v=*var._VECTptr;
    if (g.type!=_VECT || g._VECTptr->size()!=v.size())
      return gensizeerr(contextptr);
    iterateur it=v.begin(),itend=v.end(),jt=g._VECTptr->begin();
    for (;it!=itend;++it,++jt)
      *it=symbolic(at_equal,makevecteur(*it,*jt));
    if (xcas_mode(contextptr)==3)
      return symbolic(at_and,v);
    else
      return v;
  }

  static gen quote_inferieur_strict(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_inferieur_strict,eval(g,eval_level(contextptr),contextptr)));
  }

  static gen quote_superieur_strict(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_superieur_strict,eval(g,eval_level(contextptr),contextptr)));
  }

  static gen quote_inferieur_egal(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_inferieur_egal,eval(g,eval_level(contextptr),contextptr)));
  }

  static gen quote_superieur_egal(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_superieur_egal,eval(g,eval_level(contextptr),contextptr)));
  }

  static gen quote_conj(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_conj,eval(g,eval_level(contextptr),contextptr)));
  }

  static gen quote_re(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_re,eval(g,eval_level(contextptr),contextptr)));
  }

  static gen quote_im(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_im,eval(g,eval_level(contextptr),contextptr)));
  }

  vecteur solvepreprocess(const gen & args,bool complexmode,GIAC_CONTEXT){
    gen g(args);
    if (g.type==_VECT && !g._VECTptr->empty() && g._VECTptr->front().is_symb_of_sommet(at_and)){
      vecteur v(*g._VECTptr);
      v.front()=remove_and(v.front(),at_and);
      g=gen(v,g.subtype);
    }
    // quote < <= > and >=
    vector<const unary_function_ptr *> quote_inf;
    quote_inf.push_back(at_inferieur_strict);
    quote_inf.push_back(at_inferieur_egal);
    quote_inf.push_back(at_superieur_strict);
    quote_inf.push_back(at_superieur_egal);
    if (complexmode){
      quote_inf.push_back(at_conj);
      quote_inf.push_back(at_re);
      quote_inf.push_back(at_im);
    }
    vector< gen_op_context > quote_inf_v;
    quote_inf_v.push_back(quote_inferieur_strict);
    quote_inf_v.push_back(quote_inferieur_egal);
    quote_inf_v.push_back(quote_superieur_strict);
    quote_inf_v.push_back(quote_superieur_egal);
    if (complexmode){
      quote_inf_v.push_back(quote_conj);
      quote_inf_v.push_back(quote_re);
      quote_inf_v.push_back(quote_im);
    }
    g=subst(g,quote_inf,quote_inf_v,true,contextptr);
    return plotpreprocess(g,contextptr);
  }

  gen solvepostprocess(const gen & g,const gen & x,GIAC_CONTEXT){
    if (g.type!=_VECT)
      return g;
    vecteur res=*g._VECTptr;
    // convert solution to an expression
    iterateur it=res.begin(),itend=res.end();
    if (it==itend)
      return res;
    if (x.type==_VECT || xcas_mode(contextptr)==3){
      for (;it!=itend;++it)
	*it=symb_solution(*it,x,contextptr);
    }
    if (xcas_mode(contextptr)==3)
      return symbolic(at_ou,res);
    if (xcas_mode(contextptr)==2)
      return gen(res,_SET__VECT);
    return gen(res,_SEQ__VECT);
  }

  gen _solve(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    int isolate_mode=int(complex_mode(contextptr)) | int(int(all_trig_sol(contextptr)) << 1);
    vecteur v(solvepreprocess(args,complex_mode(contextptr),contextptr));
    int s=v.size();
    if (s==2 && v[1].is_symb_of_sommet(at_equal))
      return _fsolve(gen(makevecteur(v[0],v[1]._SYMBptr->feuille[0],v[1]._SYMBptr->feuille[1]),_SEQ__VECT),contextptr);
    if (s>2)
      return _fsolve(args,contextptr);
    gen arg1(v.front());
    if (arg1.type==_VECT){ // Flatten equations which are list of equations
      vecteur w;
      const_iterateur it=arg1._VECTptr->begin(),itend=arg1._VECTptr->end();
      for (;it!=itend;++it){
	gen tmp=equal2diff(*it);
	if (tmp.type==_VECT){
	  const_iterateur jt=tmp._VECTptr->begin(),jtend=tmp._VECTptr->end();
	  for (;jt!=jtend;++jt)
	    w.push_back(*jt);
	}
	else
	  w.push_back(tmp);
      }
      arg1=w;
    }
    vecteur res=solve(apply(arg1,equal2diff),v.back(),isolate_mode,contextptr);
    // if (is_fully_numeric(res))
    if (v.back().type!=_VECT && lidnt(res).empty() && is_zero(im(res,contextptr)))
      res=protect_sort(res,contextptr);
    if (!xcas_mode(contextptr))
      return res;
    return solvepostprocess(res,v[1],contextptr);
  }
  static const char _solve_s []="solve";
  static define_unary_function_eval_quoted (__solve,&_solve,_solve_s);
  define_unary_function_ptr5( at_solve ,alias_at_solve,&__solve,_QUOTE_ARGUMENTS,true);

  // also sets iszero to -2 if endpoints have same sign, -1 if err or undef
  // 1 if zero found, 2 if sign reversal (no undef),
  vecteur bisection_solver(const gen & equation,const gen & var,const gen & a0,const gen &b0,int & iszero,GIAC_CONTEXT){
    iszero=0;
    gen a(evalf_double(a0,1,contextptr)),b(evalf_double(b0,1,contextptr));
    gen fa,fb;
#ifndef NO_STDEXCEPT
    try {
#endif
      fa=subst(equation,var,a,false,contextptr);
      fb=subst(equation,var,b,false,contextptr);
#ifndef NO_STDEXCEPT
    } catch (std::runtime_error & ){
      iszero=-1;
      return vecteur(0);
    }
#endif
    int ntries=40;
    gen ab=(b-a)/ntries;
    if (fb.type!=_DOUBLE_){      
      for (int i=0;i<ntries;++i){
	b -= ab;
	fb=subst(equation,var,b,false,contextptr);
	if (fb.type==_DOUBLE_)
	  break;
      }
    }
    ab=(b-a)/ntries;
    if (fb.type==_DOUBLE_ && fa.type!=_DOUBLE_){
      for (int i=0;i<ntries;++i){
	a += ab;
	fa=subst(equation,var,a,false,contextptr);
	if (fa.type==_DOUBLE_)
	  break;
      }
    }
    if (fa.type!=_DOUBLE_ || fb.type!=_DOUBLE_){
      iszero=-1;
      return vecteur(0);
    }
    double faorig=fa._DOUBLE_val,fborig=fb._DOUBLE_val;
    if (fa._DOUBLE_val*fb._DOUBLE_val>0){
      bool test1=fa._DOUBLE_val>0;
      bool found=false;
      gen b0=b;
      // discretization of [a,b] searching a sign reversal
      for (int i=1;i<=6;i++){
	int ntest=1 << (i-1);
	gen decal=(b0-a)/gen(1 << i);
	b=a+decal;
	double bd=b._DOUBLE_val;
	decal=2*decal;
	for (int j=0;j<ntest;j++,b+=decal){
#ifndef NO_STDEXCEPT
	  try {
#endif
	    fb=subst(equation,var,b,false,contextptr);
#ifndef NO_STDEXCEPT
	  } catch (std::runtime_error & ){
	    iszero=-1;
	    return vecteur(0);
	  }
#endif
	  if (fb.type!=_DOUBLE_){
	    iszero=-1;
	    return vecteur(0);
	  }
	  double fbd=fb._DOUBLE_val;
	  bool test2=fbd>0;
	  if (test1 ^ test2){
	    found=true;
	    break;
	  }
	}
	if (found)
	  break;
      }
      if (!found){
	iszero=-2;
	return vecteur(0);
      }
    }
    if (is_exactly_zero(fa)){
      iszero=1;
      return vecteur(1,a);
    }
    if (is_exactly_zero(fb)){
      iszero=1;
      return vecteur(1,b);
    }
    // sign change in [a,b]
    // number of steps
    gen n=ln(abs(b-a,contextptr),contextptr)-ln(max(abs(b,contextptr),abs(a,contextptr),contextptr),contextptr)+53;
    n=_floor(n/0.69,contextptr);
    for (int i=0;i<n.val;i++){
      gen c=(a+b)/2,fc;
#ifndef NO_STDEXCEPT
      try {
#endif
	fc=subst(equation,var,c,false,contextptr);
#ifndef NO_STDEXCEPT
      } catch (std::runtime_error & ){
	iszero=-1;
	return vecteur(0);
      }
#endif
      if (fc.type!=_DOUBLE_){
	iszero=-1;
	if (fa.type==_DOUBLE_){
	  b=c; fb=fc;
	  continue;
	}
	if (fb.type==_DOUBLE_){
	  a=c; fa=fc;
	  continue;
	}
	return vecteur(1,c);
      }
      if (is_exactly_zero(fc)){
	iszero=1;
	return vecteur(1,c);
      }
      if (fa._DOUBLE_val*fc._DOUBLE_val>0){
	a=c;
	fa=fc;
      }
      else {
	b=c;
	fb=fc;
      }
    }
    iszero=2;
    if (fa.type==_DOUBLE_ && fb.type==_DOUBLE_ && fabs(fa._DOUBLE_val*fb._DOUBLE_val/faorig/fborig)<1e-10)
      iszero=1;
    return vecteur(1,(a+b)/2);
  }

  static void set_nearest_first(const gen & guess,vecteur & res,GIAC_CONTEXT){
    int s=res.size();
    if (s<2)
      return;
    int pos=0,i;
    gen minabs=evalf_double(abs(res[0]-guess,contextptr),1,contextptr);
    for (i=1;i<s;++i){
      gen curabs=evalf_double(abs(res[i]-guess,contextptr),1,contextptr);
      if (is_strictly_greater(minabs,curabs,contextptr)){
	minabs=curabs;
	pos=i;
      }
    }
    if (pos){
      minabs=res[0];
      res[0]=res[pos];
      res[pos]=minabs;
    }
  }

  bool is_idnt_function38(const gen & g){
    if (g.type!=_IDNT)
      return false;
    const char * ch = g._IDNTptr->id_name;
    if (strlen(ch)==2 && ch[1]>='0' && ch[1]<='9'){
      switch (ch[0]){
      case 'F': case 'X': case 'Y': case 'R':
	return true;
      }
    }
    return false;
  }

  void lidnt_solve(const gen &g,vecteur & res){
    vecteur v=lidnt(g);
    for (int i=0;i<v.size();++i){
      if (!is_idnt_function38(v[i]))
	res.push_back(v[i]);
    }
  }

  vecteur lidnt_solve(const gen & g){
    vecteur res;
    lidnt_solve(g,res);
    return res;
  }

  void lidnt_function38(const gen &g,vecteur & res){
    vecteur v=lidnt(g);
    for (int i=0;i<v.size();++i){
      if (is_idnt_function38(v[i]))
	res.push_back(v[i]);
    }
  }

  vecteur lidnt_function38(const gen & g){
    vecteur res;
    lidnt_function38(g,res);
    return res;
  }

  // Find zero or extrema of equation for variable near guess in real mode
  // For polynomial input, returns all zeros or extrema
  // On entry type=0 for zeros, =1 for extrema
  //  guess might be a single value or vecteur with 2 values (an interval)
  //  bisection is used if guess is an interval
  //  if guess is a single value, guess is checked to be in [xmin,xmax]
  // returns 0 if zero(s) were found, 1 if extrema found, 2 if sign reversal found
  vecteur solve_zero_extremum(const gen & equation0,const gen & variable0,const gen & guess,double xmin, double xmax,int & type,GIAC_CONTEXT){
    if (variable0.type!=_IDNT)
      return vecteur(1,gentypeerr(contextptr));
    vecteur l0(1,variable0);
    lidnt(equation0,l0);
    vecteur l1=gen2vecteur(eval(l0,1,contextptr));
    identificateur id_solve("aspen_x");
    gen variable(id_solve);
    l1.front()=variable;
    gen eq0=subst(equation0,l0,l1,false,contextptr),eq;
    // ofstream of("log"); of << equation0 << endl << eq0 << endl << l0 << endl << l1 << endl; of.close();
    eq0=remove_equal(eval(eq0,1,contextptr));
    vecteur res;
    if (is_undef(eq0) || is_inf(eq0)){
      type=-2;
      return res;
    }
    gen a,b;
    if (is_linear_wrt(eq0,variable,a,b,contextptr)){
      a=ratnormal(a);
      if (is_zero(a)){
	type=-1;
	return 0;
      }
      type=0;
      a=-b/a;
      b=im(a,contextptr);
      if (is_zero(b))
	res=vecteur(1,re(a,contextptr));
      else
	res=vecteur(1,undef);
      return res;
    }
    bool interval=false;
    a=xmin;b=xmax;
    if (guess.type==_VECT){
      if (guess._VECTptr->size()!=2)
	return vecteur(1,gendimerr(contextptr));
      // Find in [a,b]
      interval=true;
      a=guess._VECTptr->front();
      b=guess._VECTptr->back();
    }
    else {
      gen tmp=evalf_double(guess,1,contextptr);
      if (tmp.type==_DOUBLE_){
	if (tmp._DOUBLE_val>xmax)
	  b=tmp;
	if (tmp._DOUBLE_val<xmin)
	  a=tmp;
      }
    }
    // Check if equation is smooth, if not, find an interval for solving
#ifndef NO_STDEXCEPT
    try {
#endif
      eq=derive(eq0,variable,contextptr); 
      if (is_undef(eq))
	interval=true;
#ifndef NO_STDEXCEPT
    } catch (std::runtime_error &){
      eq=undef;
      interval=true;
    }
#endif
    // ofstream of("log"); of << eq << " " << diffeq << endl; of.close();
    if (is_zero(ratnormal(eq))){
      type=-1;
      return res;
    }
    if (type==0){ // Find zero
      if (interval){
	int iszero;
	res=bisection_solver(eq0,*variable._IDNTptr,a,b,iszero,contextptr);
	if (iszero<=0)
	  res.clear();
	if (iszero==2)
	  type=2;
      }
      else {
#ifndef NO_STDEXCEPT
	try {
#endif
	  if (lvar(eq0)==vecteur(1,variable)){
	    res=solve(eq0,*variable._IDNTptr,0,contextptr);
	    set_nearest_first(guess,res,contextptr);
	    if (res.empty())
	      type=1;
	  }
#ifndef NO_STDEXCEPT
	}
	catch (std::runtime_error & ){
	  res.clear();
	}
#endif
      }
      if (!res.empty() && is_undef(res.front()))
	res.clear();
      if (type==0 && res.empty()){
	gen sol=_fsolve(gen(makevecteur(evalf(eq0,1,contextptr),symbolic(at_equal,makevecteur(variable,guess))),_SEQ__VECT),contextptr);
	sol=evalf2bcd_nock(sol,1,contextptr);
	if (sol.type==_VECT){
	  res=*sol._VECTptr;
	  set_nearest_first(guess,res,contextptr);
	}
	else {
	  if (sol.type==_FLOAT_)
	    res=vecteur(1,sol);
	}
      }
      if (!res.empty() && !is_undef(res))
	return *eval(res,1,contextptr)._VECTptr;
    }
    if (type==0)
      type=1;
    if (type==1 && !is_undef(eq)){ // Find extremum
      if (interval){
	int iszero;
	res=bisection_solver(eq,variable,a,b,iszero,contextptr);
	if (iszero<=0)
	  res.clear();
      }
      else {
#ifndef NO_STDEXCEPT
	try {
#endif
	  if (lvar(eq)==vecteur(1,variable)){
	    res=solve(eq,*variable._IDNTptr,0,contextptr);
	    if (!res.empty() && is_undef(res.front()))
	      res.clear();
	    if (res.empty())
	      type=2;
	    else
	      set_nearest_first(guess,res,contextptr);
	  }
#ifndef NO_STDEXCEPT
	}
	catch (std::runtime_error & ){
	  res.clear();
	}
#endif
      }
      if (!res.empty() && is_undef(res.front()))
	res.clear();
      if (type==1 && res.empty()){
	gen sol=_fsolve(gen(makevecteur(evalf(eq,1,contextptr),symbolic(at_equal,makevecteur(variable,guess))),_SEQ__VECT),contextptr);
	sol=evalf2bcd_nock(sol,1,contextptr);
	if (sol.type==_VECT){
	  res=*sol._VECTptr;
	  set_nearest_first(guess,res,contextptr);
	}
	else {
	  if (sol.type==_FLOAT_)
	    res=vecteur(1,sol);
	}
      }
      if (!res.empty() && !is_undef(res))
	return *eval(res,1,contextptr)._VECTptr;
    }
    type=2; // Find singularities
    res=find_singularities(eq0,*variable._IDNTptr,0,contextptr);
    if (res.empty()) type=3;
    return *eval(res,1,contextptr)._VECTptr;
  }
  vecteur solve_zero_extremum(const gen & equation0,const gen & variable,const gen & guess,int & type,GIAC_CONTEXT){
#ifndef NO_STDEXCEPT
    try {
#endif
      return solve_zero_extremum(equation0,variable,guess,gnuplot_xmin,gnuplot_xmax,type,contextptr);
#ifndef NO_STDEXCEPT
    } catch(std::runtime_error & ){
      type=-2;
      return vecteur(1,undef);
    }
#endif
  }
  gen _solve_zero_extremum(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    vecteur v(solvepreprocess(args,complex_mode(contextptr),contextptr));
    int s=v.size();
    if (s<3 || v[1].type!=_IDNT)
      return gensizeerr(contextptr);
    int type=0;
    if (s==4 && v[3].type==_INT_)
      type=v[3].val;
    vecteur res=solve_zero_extremum(v[0],v[1],v[2],type,contextptr);
    return makevecteur(type,res);
  }
  static const char _solve_zero_extremum_s []="solve_zero_extremum";
  static define_unary_function_eval_quoted (__solve_zero_extremum,&_solve_zero_extremum,_solve_zero_extremum_s);
  define_unary_function_ptr5( at_solve_zero_extremum ,alias_at_solve_zero_extremum,&__solve_zero_extremum,_QUOTE_ARGUMENTS,true);

  double nan(){
    double x=0.0;
    return 0.0/x;
  }
#ifdef HAVE_LIBGSL
  // p should point a vector with elements the expression f(x) and x
  // OR with f(x), f'(x) and x
  static double my_f (double x0, void * p) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()<2))
      setsizeerr(gettext("solve.cc/my_f"));
#endif	// DEBUG_SUPPORT
    gen & f=params._VECTptr->front();
    gen & x=params._VECTptr->back();
    gen res=evalf(subst(f,x,x0,false,context0),1,context0);
    if (res.type==_REAL)
      res=evalf_double(res,1,context0);
#ifdef DEBUG_SUPPORT
    if (res.type>_IDNT)
      setsizeerr();
#endif
    if (res.type!=_DOUBLE_)
      return nan();
    else
      return res._DOUBLE_val;
  }

  static double my_df (double x0, void * p) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()!=3))
      setsizeerr(gettext("solve.cc/my_df"));
#endif	// DEBUG_SUPPORT
    gen & f=(*params._VECTptr)[1];
    gen & x=params._VECTptr->back();
    gen res=evalf_double(subst(f,x,x0,false,context0),1,context0);
#ifdef DEBUG_SUPPORT
    if (res.type>_IDNT)
      setsizeerr();
#endif
    if (res.type!=_DOUBLE_)
      return nan();
    else
      return res._DOUBLE_val;
  }

  static void my_fdf (double x0, void * p,double * fx,double * dfx) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()!=3))
      setsizeerr(gettext("solve.cc/my_fdf"));
#endif	// DEBUG_SUPPORT
    gen & f=params._VECTptr->front();
    gen & df=(*params._VECTptr)[1];
    gen & x=params._VECTptr->back();
    gen res=evalf_double(subst(f,x,x0,false,context0),1,context0);
    if (res.type!=_DOUBLE_)
      *fx=nan();
    else
      *fx=res._DOUBLE_val;
    res=evalf_double(subst(df,x,x0,false,context0),1,context0);
    if (res.type!=_DOUBLE_)
      *dfx=nan();
    else
      *dfx=res._DOUBLE_val;
  }

  static int my_F (const gsl_vector * x0, void * p,gsl_vector * F) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()<2))
      setsizeerr(gettext("solve.cc/my_F"));
#endif	// DEBUG_SUPPORT
    gen & f=params._VECTptr->front();
    gen & x=params._VECTptr->back();
    gen res=evalf_double(subst(f,x,gsl_vector2vecteur(x0),false,context0),1,context0);
    if (res.type!=_VECT)
      return !GSL_SUCCESS;
    return vecteur2gsl_vector(*res._VECTptr,F,context0);
  }

  static int my_dF (const gsl_vector *x0, void * p,gsl_matrix * J) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()!=3))
      setsizeerr(gettext("solve.cc/my_dF"));
#endif	// DEBUG_SUPPORT
    gen & f=(*params._VECTptr)[1];
    gen & x=params._VECTptr->back();
    gen res=evalf_double(subst(f,x,gsl_vector2vecteur(x0),false,context0),1,context0);
    if (res.type!=_VECT)
      return !GSL_SUCCESS;
    else
      return matrice2gsl_matrix(*res._VECTptr,J,context0);
  }

  static int my_FdF (const gsl_vector * x0, void * p,gsl_vector * fx,gsl_matrix * dfx) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()!=3))
      setsizeerr(gettext("solve.cc/my_FdF"));
#endif	// DEBUG_SUPPORT
    gen & f=params._VECTptr->front();
    gen & df=(*params._VECTptr)[1];
    gen & x=params._VECTptr->back();
    gen g0=gsl_vector2vecteur(x0);
    gen res=evalf_double(subst(f,x,g0,false,context0),1,context0);
    if (res.type!=_VECT)
      return !GSL_SUCCESS;
    int ires=vecteur2gsl_vector(*res._VECTptr,fx,context0);
    if (ires!=GSL_SUCCESS)
      return !GSL_SUCCESS;
    res=evalf_double(subst(df,x,g0,false,context0),1,context0);
    if (res.type!=_VECT)
      return !GSL_SUCCESS;
    return matrice2gsl_matrix(*res._VECTptr,dfx,context0);
  }

  gen msolve(const gen & f,const vecteur & vars,const vecteur & g,int method,double eps,GIAC_CONTEXT){
    vecteur guess(g);
    bool with_derivative=false;
    int dim=vars.size();
    switch (method){
    case _NEWTONJ_SOLVER: case _HYBRIDSJ_SOLVER: case _HYBRIDJ_SOLVER:
      with_derivative=true;
      break;
    case _DNEWTON_SOLVER: case _HYBRIDS_SOLVER: case _HYBRID_SOLVER:
      with_derivative=false;
      break;
    }
    if (with_derivative){
      gen difff=derive(f,vars,contextptr);
      if (is_undef(difff) || difff.type!=_VECT)
	return vecteur(vars.size(),undef);
      gen params(makevecteur(f,mtran(*difff._VECTptr),vars));
      gsl_multiroot_function_fdf FDF;
      FDF.f=&my_F;
      FDF.df=&my_dF;
      FDF.fdf=&my_FdF;
      FDF.n=dim;
      FDF.params=&params;
      const gsl_multiroot_fdfsolver_type * T=0;
      switch (method){
      case _NEWTONJ_SOLVER: 
	T=gsl_multiroot_fdfsolver_gnewton;
	break;
      case _HYBRIDSJ_SOLVER:
	T=gsl_multiroot_fdfsolver_hybridsj;
	break;
      case _HYBRIDJ_SOLVER:
	T=gsl_multiroot_fdfsolver_hybridj;
	break;
      }
      gsl_multiroot_fdfsolver * s= gsl_multiroot_fdfsolver_alloc (T, dim);
      gsl_vector * X=vecteur2gsl_vector(guess,contextptr);
      gsl_multiroot_fdfsolver_set (s, &FDF, X);
      int maxiter=SOLVER_MAX_ITERATE,res=0;
      vecteur oldguess;
      for (;maxiter;--maxiter){
	oldguess=guess;
	res=gsl_multiroot_fdfsolver_iterate(s);
	if ( (res==GSL_EBADFUNC) || (res==GSL_ENOPROG) )
	  break;
	guess=gsl_vector2vecteur(gsl_multiroot_fdfsolver_root(s));
	if (is_strictly_greater(eps,abs(guess-oldguess,contextptr),contextptr))
	  break;
      }
      gsl_multiroot_fdfsolver_free(s);
      if ( (res==GSL_EBADFUNC) || (res==GSL_ENOPROG) )
	return vecteur(dim,gensizeerr(contextptr));
      return guess;
    }
    else {
      gen params(makevecteur(f,vars));
      gsl_multiroot_function F;
      F.f=&my_F;
      F.n=dim;
      F.params=&params;
      const gsl_multiroot_fsolver_type * T=0;
      switch (method){
      case _DNEWTON_SOLVER: 
	T=gsl_multiroot_fsolver_dnewton;
	break;
      case _HYBRIDS_SOLVER:
	T=gsl_multiroot_fsolver_hybrids;
	break;
      case _HYBRID_SOLVER:
	T=gsl_multiroot_fsolver_hybrid;
	break;
      }
      gsl_multiroot_fsolver * s= gsl_multiroot_fsolver_alloc (T, dim);
      gsl_vector * X=vecteur2gsl_vector(guess,contextptr);
      gsl_multiroot_fsolver_set (s, &F, X);
      int maxiter=SOLVER_MAX_ITERATE,res=0;
      vecteur oldguess;
      for (;maxiter;--maxiter){
	oldguess=guess;
	res=gsl_multiroot_fsolver_iterate(s);
	if ( (res==GSL_EBADFUNC) || (res==GSL_ENOPROG) )
	  break;
	guess=gsl_vector2vecteur(gsl_multiroot_fsolver_root(s));
	if (is_strictly_greater(eps,abs(guess-oldguess,contextptr),contextptr))
	  break;
      }
      gsl_multiroot_fsolver_free(s);
      if ( (res==GSL_EBADFUNC) || (res==GSL_ENOPROG) )
	return vecteur(1,gensizeerr(contextptr));
      return guess;
    }
  }
#endif // HAVE_LIBGSL

  // fsolve(expr,var[,interval/guess,method])
  gen _fsolve(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    vecteur v(plotpreprocess(args,contextptr));
    return in_fsolve(v,contextptr);
  }

  gen in_fsolve(vecteur & v,GIAC_CONTEXT){
    if (is_undef(v))
      return v;    
    double gsl_eps=abs_calc_mode(contextptr)==38?1e-5:epsilon(contextptr);
    int s=v.size();
    if (s<2)
      return gentoofewargs("fsolve");
    gen v0=remove_equal(v[0]);
    if (s==2 && v[1].type==_IDNT){ 
      gen v00= evalf(v0,1,contextptr);
      // no initial guess, check for poly-like equation
      vecteur lv(lvar(v00));
      int lvs=lv.size();
      if (lv==vecteur(1,v[1])){
	gen tmp=_e2r(makevecteur(v00,v[1]),contextptr);
	if (tmp.type==_FRAC)
	  tmp=tmp._FRACptr->num;
	tmp=evalf(tmp,eval_level(contextptr),contextptr);
	if (tmp.type==_VECT)
	  return complex_mode(contextptr)?proot(*tmp._VECTptr,epsilon(contextptr)):real_proot(*tmp._VECTptr,epsilon(contextptr),contextptr);
      }
      if (lvs>1)
	v00=halftan_hyp2exp(v00,contextptr);
      lv=lvar(v00);
      lvs=lv.size();
      if (lvs==1 && lv[0].type==_SYMB && lv[0]._SYMBptr->feuille.type!=_VECT){
	int pos=equalposcomp(solve_fcns_tab,lv[0]._SYMBptr->sommet);
	if (pos){
	  gen tmp=_e2r(makevecteur(v00,lv[0]),contextptr);
	  if (tmp.type==_FRAC)
	    tmp=tmp._FRACptr->num;
	  tmp=evalf(tmp,eval_level(contextptr),contextptr);
	  if (tmp.type==_VECT){
	    vecteur res0=complex_mode(contextptr)?proot(*tmp._VECTptr,epsilon(contextptr)):real_proot(*tmp._VECTptr,epsilon(contextptr),contextptr);
	    vecteur res;
	    const_iterateur it=res0.begin(),itend=res0.end();
	    for (;it!=itend;++it){
	      vecteur res0val=gen2vecteur(isolate_fcns[pos-1](*it,complex_mode(contextptr),contextptr));
	      const_iterateur jt=res0val.begin(),jtend=res0val.end();
	      for (;jt!=jtend;++jt){
		gen fs=_fsolve(gen(makevecteur(lv[0]._SYMBptr->feuille-*jt,v[1]),_SEQ__VECT),contextptr);
		if (fs.type==_VECT)
		  res=mergevecteur(res,*fs._VECTptr);
		else
		  res.push_back(fs);
	      }
	    }
	    return res;
	  }
	}
      }
    }
    gen gguess;
    if (v[1].is_symb_of_sommet(at_equal)){
      gguess=v[1]._SYMBptr->feuille[1];
      v[1]=v[1]._SYMBptr->feuille[0];
      v.insert(v.begin()+2,gguess);
      ++s;
    }
    if (s>=3)
      gguess=v[2];
    if (gguess.is_symb_of_sommet(at_equal))
      return gensizeerr(contextptr);
    // check method
    int method=_NEWTON_SOLVER;
    //int method=0;
    if ( (s>=5) && (v[4].type==_DOUBLE_) )
      gsl_eps=v[4]._DOUBLE_val;
    if ( (s>=4) && (v[3].type==_INT_) )
      method=v[3].val;
    if (v[1].type==_VECT){
      int dim=v[1]._VECTptr->size();
      if (!dim)
	return gensizeerr(contextptr);
      if (s>=3){
	if (gguess.type!=_VECT)
	  return gensizeerr(contextptr);
	if (gguess._VECTptr->size()!=unsigned(dim))
	  return gensizeerr(contextptr);
      }
      else {
	gguess=vecteur(dim);
	gguess[0]=(gnuplot_xmin+gnuplot_xmax)/2;
	if (dim>1)
	  gguess[1]=(gnuplot_ymin+gnuplot_ymax)/2;
	if (dim>2)
	  gguess[2]=(gnuplot_zmin+gnuplot_zmax)/2;
	if (dim>3)
	  gguess[3]=(gnuplot_tmin+gnuplot_tmax)/2;
      }
#ifdef HAVE_LIBGSL
      if (method!=_NEWTON_SOLVER)
	return msolve(v0,*v[1]._VECTptr,*gguess._VECTptr,method,gsl_eps,contextptr);
#endif
    }
#ifndef HAVE_LIBGSL
    if (v[1].type!=_VECT && gguess.type==_VECT && gguess._VECTptr->size()==2){
      int iszero;
      vecteur res= bisection_solver(v0,v[1],gguess[0],gguess[1],iszero,contextptr);
      if (!res.empty() && iszero!=1)
	*logptr(contextptr) << (iszero==-1?"Warning: undefined":"Warning: sign reversal") << endl;
      return res;
    }
#endif
#ifdef HAVE_LIBGSL
    if (method!=_NEWTON_SOLVER){
      bool with_derivative=false;
      switch (method){
      case _BISECTION_SOLVER: case _FALSEPOS_SOLVER: case _BRENT_SOLVER:
	with_derivative=false;
	break;
      case _NEWTON_SOLVER: case _SECANT_SOLVER: case _STEFFENSON_SOLVER:
	with_derivative=true;
	break;
      }
      gen params;
      if (with_derivative){
	gen dv0=derive(v0,v[1],contextptr);
	if (is_undef(dv0))
	  return dv0;
	params= makevecteur(v0,dv0,v[1]);
	double guess((gnuplot_xmin+gnuplot_xmax)/2),oldguess;
	if (s>=3){
	  gen g=evalf(gguess,eval_level(contextptr),contextptr);
	  if (gguess.type==_DOUBLE_)
	    guess=gguess._DOUBLE_val;
	}
	gsl_function_fdf FDF ;     
	FDF.f = &my_f ;
	FDF.df = &my_df ;
	FDF.fdf = &my_fdf ;
	FDF.params = &params ;
	gsl_root_fdfsolver * slv=0;
	switch (method){
	case _NEWTON_SOLVER:
	  slv=gsl_root_fdfsolver_alloc (gsl_root_fdfsolver_newton);
	  break;
	case _SECANT_SOLVER: 
	  slv=gsl_root_fdfsolver_alloc (gsl_root_fdfsolver_secant);
	  break;
	case _STEFFENSON_SOLVER:
	  slv=gsl_root_fdfsolver_alloc (gsl_root_fdfsolver_steffenson);
	  break;
	}
	if (!slv)
	  return gensizeerr(contextptr);
	gsl_root_fdfsolver_set(slv,&FDF,guess);
	int maxiter=SOLVER_MAX_ITERATE,res=0;
	for (;maxiter;--maxiter){
	  oldguess=guess;
	  res=gsl_root_fdfsolver_iterate(slv);
	  guess=gsl_root_fdfsolver_root(slv);
	  if ( (res==GSL_EBADFUNC) || (res==GSL_EZERODIV) )
	    break;
	  if (fabs(guess-oldguess)<gsl_eps)
	    break;
	}
	gsl_root_fdfsolver_free(slv);
	if (!maxiter)
	  return gensizeerr(contextptr);
	if ( (res==GSL_EBADFUNC) || (res==GSL_EZERODIV) )
	  return undef;
	else
	  return guess;
      }
      else {
	params= makevecteur(v0,v[1]);
	double x_low,x_high;
	if (s>=3) {
	  vecteur w;
	  if (gguess.type==_VECT)
	    w=*gguess._VECTptr;
	  else {
	    if ( (gguess.type==_SYMB) && (gguess._SYMBptr->sommet==at_interval) )
	      w=*gguess._SYMBptr->feuille._VECTptr;
	  }
	  if (w.size()!=2)
	    return gentypeerr(contextptr);
	  gen low=w[0].evalf(eval_level(contextptr),contextptr);
	  gen high=w[1].evalf(eval_level(contextptr),contextptr);
	  if ( (low.type!=_DOUBLE_) || (high.type!=_DOUBLE_) )
	    return gensizeerr(contextptr);
	  x_low=low._DOUBLE_val;
	  x_high=high._DOUBLE_val;
	}
	else {
	  x_low=gnuplot_xmin;
	  x_high=gnuplot_xmax;
	}
	if (x_low>x_high){
	  double tmp=x_low;
	  x_low=x_high;
	  x_high=tmp;
	}
	gsl_function F ;
	F.function=&my_f;
	F.params = &params ;
	gsl_root_fsolver * slv =0 ;
	switch (method){
	case  _BISECTION_SOLVER:
	  slv=gsl_root_fsolver_alloc (gsl_root_fsolver_bisection);
	  break;
	case _FALSEPOS_SOLVER: 
	  slv=gsl_root_fsolver_alloc (gsl_root_fsolver_falsepos);
	  break;
	case _BRENT_SOLVER:
	  slv=gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
	  break;
	}
	if (!slv)
	  return gensizeerr(contextptr);
	gsl_root_fsolver_set (slv,&F,x_low,x_high);
	int res=0;
	int maxiter=SOLVER_MAX_ITERATE;
	for (;maxiter && (x_high-x_low>gsl_eps);--maxiter){
	  res=gsl_root_fsolver_iterate(slv);
	  if (res==GSL_EBADFUNC)
	    break;
	  x_low=gsl_root_fsolver_x_lower(slv);
	  x_high= gsl_root_fsolver_x_upper(slv);
	}
	gsl_root_fsolver_free (slv);
	if (res==GSL_EBADFUNC)
	  return undef;
	return makevecteur(x_low,x_high);
      } // end if derivative
    }
#else // HAVE_LIBGSL
    if (method!=_NEWTON_SOLVER)
      return gensizeerr(gettext("Not linked with GSL"));
#endif    // HAVE_LIBGSL
    else  {// newton method, call newton
      gguess=newton(v0,v[1],gguess,NEWTON_DEFAULT_ITERATION,gsl_eps,1e-12,contextptr);
      if (is_greater(1e-8,im(gguess,contextptr)/re(gguess,contextptr),contextptr))
	return re(gguess,contextptr);
      return gguess;
    }
    return undef;
  } // end f_solve
  static const char _fsolve_s []="fsolve";
  static define_unary_function_eval_quoted (__fsolve,&_fsolve,_fsolve_s);
  define_unary_function_ptr5( at_fsolve ,alias_at_fsolve,&__fsolve,_QUOTE_ARGUMENTS,true);

  vecteur sxa(const vecteur & sl_orig,const vecteur & x,GIAC_CONTEXT){
    vecteur sl(sl_orig);
    int d;
    d=x.size();
    int de;
    de=sl.size();
    for (int i=0;i<de;i++){
      //gen e:
      //e=sl[i];    
      if ( (sl[i].type==_SYMB) && ((*sl[i]._SYMBptr).sommet==at_equal || (*sl[i]._SYMBptr).sommet==at_same)){
	sl[i]=(*sl[i]._SYMBptr).feuille[0]-(*sl[i]._SYMBptr).feuille[1];
      }
    }
    vecteur A;
    for (int i=0;i<de;i++){
      vecteur li(d+1);
      gen lo=sl[i];
      for (int j=0;j<d;j++){
	lo=subst(lo,x[j],0,false,contextptr);
	li[j]=derive(sl[i],x[j],contextptr);
      }
      li[d]=lo;
      A.push_back(li);
    }
    return(A);
  }

  vecteur linsolve(const vecteur & sl,const vecteur & x,GIAC_CONTEXT){
    vecteur A; 
    if (ckmatrix(sl)){
      unsigned int n=sl.size();
      A=mtran(sl);
      if (ckmatrix(x)){
	if (x.size()==1){
	  if (x.front()._VECTptr->size()!=n)
	    return vecteur(1,gendimerr(contextptr));
	  A.push_back(-x.front());
	}
	else {
	  if (x.size()!=n)
	    return vecteur(1,gendimerr(contextptr));
	  matrice xm=mtran(x);
	  if (xm.size()!=1)
	    return vecteur(1,gensizeerr(contextptr));
	  A.push_back(-xm.front());
	}
      }
      else {
	if (x.size()!=n)
	  return vecteur(1,gendimerr(contextptr));
	A.push_back(-x);
      }
      A=mtran(A);
      vecteur B=-mker(A,contextptr);
      if (is_undef(B) || B.empty())
	return B;
      // The last element of B must have a non-zero last component
      vecteur Bend=*B.back()._VECTptr;
      gen last=Bend.back();
      if (is_zero(last))
	return vecteur(0);
      gen R=Bend/last;
      // The solution is sum(B[k]*Ck+Blast/last)
      int s=B.size();
      for (int k=0;k<s-1;k++)
	R=R+gen("C_"+print_INT_(k),contextptr)*B[k];
      vecteur res=*R._VECTptr;
      res.pop_back();
      return res;
    }
    A=sxa(sl,x,contextptr);
    vecteur B,R(x);
    gen rep;
    B=mrref(A,contextptr);
    //cout<<B<<endl;
    int d=x.size();
    int de=sl.size();
    for (int i=0; i<de;i++){
      vecteur li(d+1);
      for(int k=0;k<d+1;k++){
	li[k]=B[i][k];
      }
      int j;
      j=i;
      while (li[j]==0 && j<d){
	j=j+1;
      }
      if (j==d && !is_zero(li[d])){
	return vecteur(0);
      } 
      else {
	if (j<d){
	  rep=-li[d];
	  for (int k=j+1;k<d;k++){
	    rep=rep-li[k]*x[k];
	  }
	  rep=rdiv(rep,li[j]);
	  R[j]=rep;
	}
      }
    }
    return R;
  }

  gen equal2diff(const gen & g){
    if ( (g.type==_SYMB) && (g._SYMBptr->sommet==at_equal || g._SYMBptr->sommet==at_same) ){
      vecteur & v=*g._SYMBptr->feuille._VECTptr;
      return v[0]-v[1];
    }
    else
      return g;
  }

  gen symb_linsolve(const gen & syst,const gen & vars){
    return symbolic(at_linsolve,makevecteur(syst,vars));
  }
 
  gen linsolve(const gen & syst,const gen & vars,GIAC_CONTEXT){
    if ((syst.type!=_VECT)||(vars.type!=_VECT))
      return symb_linsolve(syst,vars);
    gen res=linsolve(*syst._VECTptr,*vars._VECTptr,contextptr);
    if (!has_i(syst) && has_i(res))
      res=_evalc(res,contextptr);
    else
      res=normal(res,contextptr);
    return res;
  }
  
  gen _linsolve(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    vecteur v(plotpreprocess(args,contextptr));
    if (is_undef(v))
      return v;
    int s=v.size();
    if (s!=2)
      return gentoomanyargs("linsolve");
    if (v[1].type==_IDNT)
      v[1]=eval(v[1],eval_level(contextptr),contextptr);
    gen syst=apply(v[0],equal2diff),vars=v[1];
    return linsolve(syst,v[1],contextptr);
  }
  static const char _linsolve_s []="linsolve";
  static define_unary_function_eval_quoted (__linsolve,&_linsolve,_linsolve_s);
  define_unary_function_ptr5( at_linsolve ,alias_at_linsolve,&__linsolve,_QUOTE_ARGUMENTS,true);

  static const char _resoudre_systeme_lineaire_s []="resoudre_systeme_lineaire";
  static define_unary_function_eval_quoted (__resoudre_systeme_lineaire,&_linsolve,_resoudre_systeme_lineaire_s);
  define_unary_function_ptr5( at_resoudre_systeme_lineaire ,alias_at_resoudre_systeme_lineaire,&__resoudre_systeme_lineaire,_QUOTE_ARGUMENTS,true);

  /*
  gen iter(const gen & f, const gen & x,const gen & arg,int maxiter,double eps, int & b){
    gen a=arg;
    complex<double> olda;
    complex<double> ad;
    b=0;
    ad=gen2complex_d(a);
    //cout<<"a"<<a<<endl;
    //cout<<"ad"<<ad<<endl;
    for (int j=1;j<=maxiter;j++){
      olda=ad;    
      // cout << f << " " << x << " " << a << endl;
      a=subst(f,x,a).evalf();
      // cout<<"a"<<a<<endl;
      //a=a.evalf();
      //ad=a._DOUBLE_val;
      ad=gen2complex_d(a);
      // cout<<"a"<<a<<endl;
      // cout<<"ad"<<ad<<endl;
      // cout<<"j"<<j<<"abs"<<abs(ad-olda)<<endl;
      if (eps>abs(ad-olda)) {
	b=1; return(a);
      }
    } 
    return(a); 
  }
  
  gen newtona(const gen & f, const gen & x, const gen & arg,int niter1, int niter2, double eps1,double eps2,double prefact1,double prefact2, int & b){
    if (x.type!=_IDNT)
      settypeerr(gettext("2nd arg must be an identifier"));
    //cout<<a<<endl;
    gen a=arg;
    gen g1;
    gen g;
    g1=x-gen(prefact1)*rdiv(f,derive(f,x));
    // sym_sub(x,sym_mult(rdiv(9,10),rdiv(f,derive(f,x))));
    try {
      a= iter(g1,x,a,niter1,eps1,b);
      g=x-gen(prefact2)*rdiv(f,derive(f,x));
      a= iter(g,x,a,niter2,eps2,b); 
    }
    catch (std::runtime_error & err){
      b=0;
    }
    return a;
  }
 
  gen newton(const gen & f, const gen & x,const gen & guess,int niter1,int niter2,double eps1,double eps2,double prefact1,double prefact2){
    bool guess_first=is_undef(guess);
    for (int j=1;j<5;j++,niter2 *=2, niter1 *=2){ 
      gen a;
      int b;
      //on prend un départ au hasard (a=x0=un _DOUBLE_)
      // a=gen(2.0);
      if (guess_first)
	a=j*4*(rand()/(RAND_MAX+1.0)-0.5);
      else {
	a=guess;
	guess_first=true;
      }
      // cout<<j<<"j"<<a<<endl; 
      gen e;
      e=newtona(f, x, a,niter1,niter2,eps1,eps2,prefact1,prefact2,b);
      if (b==1) return e;
      gen c;
      c=j*4*(rand()/(RAND_MAX+1.0)-0.5);
      // cout<<j<<"j"<<c<<endl;
      // g=x-gen(0.5)*rdiv(f,derive(f,x));
      gen ao(gen(0.0),c);
      // cout<<"ao"<<ao<<endl;
      gen e0= newtona(f, x, ao,niter1,niter2,eps1,eps2,prefact1,prefact2,b);
      if (b==1) 
	return(e0);
      gen a1(a,c);
      // cout<<j<<"j,a1"<<a1<<endl;
      e0= newtona(f, x, a1,niter1,niter2,eps1,eps2,prefact1,prefact2,b);
      if (b==1) 
	return(e0);
    }
    setsizeerr(gettext("nontrouve"));
    return(0);
  }
  */

  static gen newton_rand(int j,GIAC_CONTEXT){
    gen a=evalf(j*4*(gen(rand())/(gen(RAND_MAX)+1)-plus_one_half),1,contextptr);
    if (j>2 && complex_mode(contextptr))
      a=a+cst_i*evalf(j*4*(gen(rand())/(gen(RAND_MAX)+1)-plus_one_half),1,contextptr);
    return a;
  }

  gen newton(const gen & f0, const gen & x,const gen & guess_,int niter,double eps1, double eps2,GIAC_CONTEXT){
    bool out=niter!=NEWTON_DEFAULT_ITERATION;
    gen guess(guess_);
    // ofstream of("log"); of << f0 << endl << x << endl << guess << endl << niter ; 
    gen f(eval(f0,1,context0));  // eval of f wrt context0 is intentionnal, replace UTPN by erf
    if (guess.is_symb_of_sommet(at_interval))
      guess=(guess._SYMBptr->feuille[0]+guess._SYMBptr->feuille[1])/2;
    gen a,b,d,fa,fb,invdf=inv(derive(f,x,contextptr),contextptr),epsg1(eps1),epsg2(eps2);
    if (is_undef(invdf))
      return invdf;
    if (ckmatrix(invdf))
      invdf=mtran(*invdf._VECTptr);
    bool guess_first=is_undef(guess);
    // Main loop with random initialization
    int j=1;
    for (;j<=5 ;j++,niter += 5){ 
      if (guess_first){
	if (f.type==_VECT){
	  int s=f._VECTptr->size();
	  vecteur v(s);
	  for (int i=0;i<s;++i)
	    v[i]=(newton_rand(j,contextptr));
	  a=v;
	}
	else
	  a=newton_rand(j,contextptr);
      }
      else {
	a=guess;
	guess_first=true;
      }
#ifndef NO_STDEXCEPT
      try {
#endif
	fa=evalf(eval(subst(f,x,a,false,contextptr),eval_level(contextptr),contextptr),1,contextptr); 
	// First loop to localize the solution with prefactor
	gen lambda(1);
	int k;
	for (k=0;k<niter;++k){
	  if (ctrl_c) { 
	    interrupted = true; 
	    return gensizeerr(gettext("Stopped by user interruption.")); 
	  }
	  d=subst(invdf,x,a,false,contextptr);
	  // of << k << " " << d << " " << invdf << " " << x << " " << a << " ";
	  d=eval(d,eval_level(contextptr),contextptr);
	  // of << d << " " << fa << " ";
	  d=-evalf(d*fa,1,contextptr);
	  if (d.type!=_FLOAT_ && d.type!=_DOUBLE_ && d.type!=_CPLX && d.type!=_REAL && !is_undef(d) && !is_inf(d))
	    return gensizeerr(contextptr);
	  if (k==0 && is_zero(d) && !is_zero(fa)){
	    a=newton_rand(j,contextptr);
	    continue;
	  }
	  // of << d << " " << endl;
	  // of << k << " " << invdf << " " << " " << f << " " << x << " " << a << " " << fa << " " << d << " " << epsg1 << endl;
	  // cerr << k << " " << invdf << " " << " " << f << " " << x << " " << a << " " << fa << " " << d << " " << epsg1 << endl;
	  b=a+lambda*d;
	  gen babs=_l2norm(b,contextptr);
	  if (is_inf(babs) || is_undef(babs)){
	    guess_first=true;
	    k=niter;
	    break;
	  }
	  if (is_positive(epsg1-_l2norm(d,contextptr),contextptr)){
	    a=b;
	    break;
	  }
	  fb=evalf(eval(subst(f,x,b,false,contextptr),eval_level(contextptr),contextptr),1,contextptr);
	  if (is_positive(_l2norm(fb,contextptr)-_l2norm(fa,contextptr),contextptr)){
	    // Decrease prefactor and try again
	    lambda=evalf(plus_one_half,1,contextptr)*lambda;
	  }
	  else {
	    // Save new value of a and increase the prefactor slightly
	    if (is_positive(lambda-0.9,contextptr))
	      lambda=1;
	    else
	      lambda=evalf(gen(12)/gen(10),1,contextptr)*lambda;
	    a=b;
	    fa=fb;
	  }
	} // end for (k<niter)
	if (k==niter){
	  if (out)
	    return a;
	  continue;
	}
	// Second loop to improve precision (prefactor 1)
	for (k=0;k<niter;++k){
	  if (ctrl_c) { 
	    interrupted = true; 
	    return gensizeerr(gettext("Stopped by user interruption.")); 
	  }
	  d=-evalf(subst(invdf,x,a,false,contextptr)*subst(f,x,a,false,contextptr),1,contextptr);
	  a=a+d;
	  if (is_positive(epsg2-_l2norm(d,contextptr),contextptr))
	    break;
	}
	if (k!=niter || is_positive(epsg1-_l2norm(d,contextptr),contextptr))
	  break;
#ifndef NO_STDEXCEPT
      } catch (std::runtime_error & ){
	continue; // start with a new initial point
      }
#endif
    } // end for
    if (j>5)
      return undef;
    return a;
  }

  gen _newton(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    double gsl_eps=epsilon(contextptr);
    if (args.type!=_VECT)
      return newton(args,vx_var,undef,NEWTON_DEFAULT_ITERATION,gsl_eps,1e-12,contextptr);
    vecteur v=*args._VECTptr;
    int s=v.size();
    if (s<2)
      return gensizeerr(contextptr);
    if (s==2){
      if (v[1].is_symb_of_sommet(at_equal))
	return newton(v[0],v[1]._SYMBptr->feuille[0],v[1]._SYMBptr->feuille[1],NEWTON_DEFAULT_ITERATION,gsl_eps,1e-12,contextptr);
      return newton(v[0],v[1],undef,NEWTON_DEFAULT_ITERATION,gsl_eps,1e-12,contextptr);
    }
    int niter=NEWTON_DEFAULT_ITERATION;
    double eps=epsilon(contextptr);
    for (int j=3;j<s;++j){
      if (v[j].type==_INT_)
	niter=v[j].val;
      else {
	gen tmp=evalf_double(v[j],1,contextptr);
	if (tmp.type==_DOUBLE_)
	  eps=tmp._DOUBLE_val;
      }
    }
    gen res=newton(v[0],v[1],v[2],niter,1e-10,eps,contextptr);
    if (debug_infolevel)
      *logptr(contextptr) << res << endl;
    return res;
    return gentoomanyargs("newton");
  }
  static const char _newton_s []="newton";
  static define_unary_function_eval (__newton,&_newton,_newton_s);
  define_unary_function_ptr5( at_newton ,alias_at_newton,&__newton,0,true);
  
  bool has_num_coeff(const vecteur & v){
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (has_num_coeff(*it))
	return true;
    }
    return false;
  }
  
  bool has_num_coeff(const polynome & p){
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (;it!=itend;++it){
      if (has_num_coeff(it->value))
	return true;
    }
    return false;
  }

  bool has_num_coeff(const gen & e){
    switch (e.type){
    case _ZINT: case _INT_: case _IDNT: case _USER:
      return false;
    case _DOUBLE_: case _REAL: case _FLOAT_:
      return true;
    case _CPLX:
      return (e._CPLXptr->type==_DOUBLE_) || ((e._CPLXptr+1)->type==_DOUBLE_);
    case _SYMB:
      return has_num_coeff(e._SYMBptr->feuille);
    case _VECT:
      return has_num_coeff(*e._VECTptr);
    case _POLY:
      return has_num_coeff(*e._POLYptr);
    case _FRAC:
      return has_num_coeff(e._FRACptr->num) || has_num_coeff(e._FRACptr->den);
    default:
      return false;
    }
    return 0;
  }

  bool has_mod_coeff(const vecteur & v,gen & modulo){
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (has_mod_coeff(*it,modulo))
	return true;
    }
    return false;
  }

  bool has_mod_coeff(const polynome & p,gen & modulo){
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (;it!=itend;++it){
      if (has_mod_coeff(it->value,modulo))
	return true;
    }
    return false;
  }

  bool has_mod_coeff(const gen & e,gen & modulo){
    switch (e.type){
    case _MOD:
      modulo = *(e._MODptr+1);
      return true;
    case _SYMB:
      return has_mod_coeff(e._SYMBptr->feuille,modulo);
    case _VECT:
      return has_mod_coeff(*e._VECTptr,modulo);
    case _POLY:
      return has_mod_coeff(*e._POLYptr,modulo);
    default:
      return false;
    }
  }

  polynome spoly(const polynome & p,const polynome & q,environment * env){
    if (p.coord.empty())
      return q;
    if (q.coord.empty())
      return p;
    const index_t & pi = p.coord.front().index.iref();
    const index_t & qi = q.coord.front().index.iref();
    index_t lcm= index_lcm(pi,qi);
    polynome tmp=p.shift(lcm-pi,q.coord.front().value)-q.shift(lcm-qi,p.coord.front().value);
    return (env && env->moduloon)?smod(tmp,env->modulo):tmp;
  }


  static bool GroebnerDiv (const polynome & celuici,const polynome & other, polynome & quo, polynome & rem, bool allowrational = true ) {  
    if (celuici.coord.empty()){
      quo=celuici;
      rem=celuici; 
      return true;
    }
    assert(!other.coord.empty());
    quo.coord.clear();
    rem=celuici;
    const index_m & b_max = other.coord.front().index;
    const gen & b=other.coord.front().value;
    std::vector< monomial<gen> >::const_iterator it,itend;
    for (;;){
      itend=rem.coord.end();
      it=rem.coord.begin();
      // look in rem for a monomial >= b_max
      for (; it!=itend ;++it){
	if (it->index>=b_max)
	  break;
      }
      if (it==itend) // no monomial of rem are divisible by LT(b): finished
	break;
      gen q=rdiv(it->value,b);
      if (!allowrational){ // quotient not allowed
	if ( has_denominator(q) || (!is_zero(q*b - it->value)) )
	  return false;
      }
      quo.coord.push_back(monomial<gen>(q,it->index-b_max));
      polynome temp=other.shift(it->index-b_max,q);
      rem = rem-temp;
    } // end for
    return(true);    
  }


  template <class T>
  static void polynome2map(const polynome & p,map<index_t,gen,T> & pmap){
    std::vector< monomial<gen> >::const_iterator pt=p.coord.begin(),ptend=p.coord.end();
    pmap.clear();
    for (;pt!=ptend;++pt){
      pmap[pt->index.iref()]=pt->value;
    }
  }

  // find integer content and divide
  template <class T>
  static gen mapppz(map<index_t,gen,T> & pmap){
    typename std::map< index_t,gen,T>::iterator it=pmap.begin(),itend=pmap.end();
    gen d;
    for (;it!=itend;++it){
      d=gcd(d,it->second);
      if (is_one(d))
	return plus_one;
    }
    for (it=pmap.begin();it!=itend;++it)
      it->second=it->second/d;
    return d;
  }

  static polynome reduce2(const polynome & p,vectpoly::const_iterator it0,vectpoly::const_iterator itend){
    if (p.coord.empty())
      return p;
    typedef std::map< index_t,gen,const std::pointer_to_binary_function < const index_t &, const index_t &, bool> > application;
    application rem(ptr_fun(p.is_strictly_greater));
    polynome2map(p,rem);
    application::const_iterator pt,ptend;
    vectpoly::const_iterator it;
    for (;;){
      ptend=rem.end();
      // look in rem for a monomial >= to a monomial in it0, then it0+1 
      for (it=it0; it!=itend ;++it){
	for (pt=rem.begin();pt!=ptend;++pt){
	  if (pt->first>=it->coord.front().index)
	    break;
	}
	if (pt!=ptend)
	  break;
      }
      if (it==itend) // no monomial of rem are divisible by LT(b): finished
	break;
      gen a(pt->second),b(it->coord.front().value) ;
      simplify(a,b);
      if (is_minus_one(b)){
	a=-a;
	b=1;
      }
      polynome temp=it->shift(pt->first-it->coord.front().index,a);
      application::iterator pit;
      if (!is_one(b)){
	for (pit=rem.begin();pit!=ptend;++pit)
	  pit->second = b*pit->second;
      }
      // substract temp from rem
      vector< monomial<gen> >::const_iterator jt=temp.coord.begin(),jtend=temp.coord.end();
      for (;jt!=jtend;++jt){
	pit=rem.find(jt->index.iref());
	if (pit==ptend)
	  rem[jt->index.iref()]=-jt->value;
	else {
	  pit->second -= jt->value;
	  if (is_zero(pit->second)){
	    rem.erase(pit->first);
	    ptend=rem.end();
	  }
	}
      }
      if (!is_one(b))
	mapppz(rem);
    }
    // convert back
    polynome res(p.dim,p);
    pt=rem.begin();
    res.coord.reserve(rem.size());
    for (;pt!=ptend;++pt){
      res.coord.push_back(monomial<gen>(pt->second,pt->first));
    }
    gen d=ppz(res);
    if (!is_one(d)) cerr << d << endl;
    return res;
  }

  static polynome reducegb(const polynome & p,vectpoly::const_iterator it0,vectpoly::const_iterator itend,environment * env){
    if (p.coord.empty())
      return p;
    polynome rem(p),quo,tmp;
    for (;it0!=itend;++it0){
      if (!GroebnerDiv(rem,*it0,quo,tmp,true)){
#ifdef NO_STDEXCEPT
	setsizeerr();
#endif
	return rem;
      }
      rem=tmp;
      ppz(rem);
    }
    return rem;
  }

  polynome reduce(const polynome & p,vectpoly::const_iterator it0,vectpoly::const_iterator itend,environment * env){
    if (p.coord.empty())
      return p;
    polynome rem(p);
    std::vector< monomial<gen> >::const_iterator pt,ptend;
    vectpoly::const_iterator it;
    for (;;){
      ptend=rem.coord.end();
      // look in rem for a monomial >= to a monomial in it0, then it0+1 
      for (it=it0; it!=itend ;++it){
	for (pt=rem.coord.begin();pt!=ptend;++pt){
	  if (pt->index>=it->coord.front().index)
	    break;
	}
	if (pt!=ptend)
	  break;
      }
      if (it==itend) // no monomial of rem are divisible by LT(b): finished
	break;
      gen a(pt->value),b(it->coord.front().value) ;
      if (env && env->moduloon){
	polynome temp=it->shift(pt->index-it->coord.front().index,a*invmod(b,env->modulo));
	rem = smod(rem - temp,env->modulo) ; // FIXME: improve!
      }
      else {
	simplify(a,b);
	polynome temp=it->shift(pt->index-it->coord.front().index,a);
	if (is_one(b))
	  rem = rem-temp;
	else {
	  rem = b*rem - temp;
	  ppz(rem);
	}
      }
    }
    if (env && env->moduloon)
      ;
    else
      ppz(rem);
    return rem;
  }

  polynome reduce(const polynome & p,const vectpoly & v,environment * env){
    vectpoly::const_iterator it=v.begin(),itend=v.end();
    return reduce(p,it,itend,env);
  }

  void reduce(vectpoly & res,environment * env){
    if (res.empty())
      return;
    sort(res.begin(),res.end(),tensor_is_greater<gen>);
    // reduce res
    int s=res.size();
    for (int i=0;i<s-1;){
      polynome & p=res[i];
      polynome pred=reduce(p,res.begin()+i+1,res.end(),env);
      if (pred.coord.empty()){
	res.erase(res.begin()+i);
	--s;
	continue;
      }
      if (pred.coord.size()==p.coord.size() && pred*p.coord.front().value==p*pred.coord.front().value){
	++i;
	continue;
      }
      res[i]=pred;
      sort(res.begin()+i,res.end(),tensor_is_greater<gen>);
      i=0;
    }
  }

  vectpoly gbasis(const vectpoly & v,const gen & order,bool with_cocoa,bool with_f5,environment * env){
    if (v.size()==1){
      return v;
    }
    vectpoly res(v);
#ifndef NO_STDEXCEPT
    try {
#endif
      if (with_cocoa){
	bool ok=with_f5?f5(res,order):cocoa_gbasis(res,order);
	if (ok){
	  if (debug_infolevel)
	    cerr << res << endl;
	  return res;
	}
      }
#ifndef NO_STDEXCEPT
    } catch (...){
      cerr << "Unable to compute gbasis with CoCoA" << endl;
    }
#endif
    reduce(res,env);
    bool notfound=true;
    for (;notfound;){
      // cerr << res << endl;
      notfound=false;
      vectpoly::const_iterator it=res.begin(),itend=res.end(),jt;
      vectpoly newres(res);
      for (;it!=itend;++it){
	for (jt=it+1;jt!=itend;++jt){
	  polynome toadd(spoly(*it,*jt,env));
	  toadd=reduce(toadd,newres,env);
	  if (!toadd.coord.empty()){
	    newres.push_back(toadd); // should be at the right place
	    notfound=true;
	  }
	}
      }
      reduce(newres,env);
      res=newres;
    }
    return res;
  }

  static gen in_ideal(const vectpoly & r,const vectpoly & v,const gen & order,bool with_cocoa,bool with_f5,environment * env){
#ifndef NO_STDEXCEPT
    try {
      if (with_cocoa){
	return cocoa_in_ideal(r,v,order);
      }
    } catch (...){
     return -1;
    }
#endif
    return -1;
  }

  gen remove_equal(const gen & f){
    if ( (f.type==_SYMB) && (f._SYMBptr->sommet==at_equal || f._SYMBptr->sommet==at_same ) ){
      vecteur & v=*f._SYMBptr->feuille._VECTptr;
      return v.front()-v.back();
    }
    if (f.type==_VECT)
      return apply(f,remove_equal);
    return f;
  }

  vecteur remove_equal(const_iterateur it,const_iterateur itend){
    vecteur conditions;
    conditions.reserve(itend-it);
    for (;it!=itend;++it){
	conditions.push_back(remove_equal(*it));
    }
    return conditions;
  }

  bool vecteur2vector_polynome(const vecteur & eq_in,const vecteur & l,vectpoly & eqp){
    // remove all denominators
    const_iterateur it=eq_in.begin(),itend=eq_in.end();
    for (;it!=itend;++it){
      gen n,d;
      fxnd(*it,n,d);
      if (n.type==_POLY){
	// should reordre n with total degree+revlex order here
	eqp.push_back(*n._POLYptr);
	continue;
      }
      if (!is_zero(n))
	return false;
    }
    return true;
  }

  vecteur gsolve(const vecteur & eq_orig,const vecteur & var_orig,bool complexmode,GIAC_CONTEXT){
    // replace variables in var_orig by true identificators
    vecteur var(var_orig);
    iterateur it=var.begin(),itend=var.end();
    int s=itend-it; // # of unknowns
    bool need_subst=false;
    vector<identificateur> tab_idnt(s);
    for (int i=0;it!=itend;++it,++i){
      if (it->type!=_IDNT){
	*it=tab_idnt[i]; 
	need_subst=true;
      }
    }
    vecteur eq(remove_equal(eq_orig.begin(),eq_orig.end()));
    if (need_subst)
      eq=subst(eq,var_orig,var,false,contextptr);
    if (approx_mode(contextptr)){
#ifdef HAVE_LIBGSL
      return makevecteur(msolve(eq,var,multvecteur(zero,var),_HYBRID_SOLVER,epsilon(contextptr),contextptr));
#else
      return vecteur(1,undef);
#endif
    }
    bool convertapprox=has_num_coeff(eq);
    if (convertapprox)
      eq=*exact(evalf(eq,1,contextptr),contextptr)._VECTptr;
    // check rational
    for (it=var.begin();it!=itend;++it){
      if (it->type!=_IDNT) // should not occur!
	return vecteur(1,gensizeerr(gettext("Bad var ")+it->print(contextptr)));
      vecteur l(rlvarx(eq,*it));
      if (l.size()>1)
	return vecteur(1,string2gen(gen(l).print(contextptr)+" is not rational w.r.t. "+it->print(contextptr),false));
    }
    vecteur l(1,var);
    alg_lvar(eq,l);
    // convert eq to polynomial
    vecteur eq_in(*e2r(eq,l,contextptr)._VECTptr);
    vectpoly eqp;
    // remove all denominators
    it=eq_in.begin();
    itend=eq_in.end();
    for (;it!=itend;++it){
      gen n,d;
      fxnd(*it,n,d);
      if (n.type==_POLY){
	// should reordre n with total degree+revlex order here
	eqp.push_back(*n._POLYptr);
	continue;
      }
      if (!is_zero(n))
	return vecteur(0); // no solution since cst equation
    }
    vectpoly eqpr(gbasis(eqp,_PLEX_ORDER));
    // should reorder eqpr with lex order here
    // solve from right to left
    reverse(eqpr.begin(),eqpr.end());
    vecteur sols(1,vecteur(0)); // sols=[ [] ]
    vectpoly::const_iterator jt=eqpr.begin(),jtend=eqpr.end();
    for (;jt!=jtend;++jt){
      // the # of found vars is the size of sols.front()
      if (sols.empty())
	break;
      vecteur newsols;
      gen g(r2e(*jt,l,contextptr));
      const_iterateur st=sols.begin(),stend=sols.end();
      for (;st!=stend;++st){
	int foundvars=st->_VECTptr->size();
	vecteur current=*st->_VECTptr;
	gen curg=ratnormal(ratnormal(subst(g,vecteur(var.end()-foundvars,var.end()),*st,false,contextptr)));
	gen x;
	int xpos=0;
	// First search in current an identifier curg depends on
	for (;xpos<foundvars;++xpos){
	  x=current[xpos];
	  if (x==var[s-foundvars+xpos] && !is_zero(derive(curg,x,contextptr)) )
	    break;
	}
	if (xpos==foundvars){
	  xpos=0;
	  // find next var g depends on 
	  for (;foundvars<s;++foundvars){
	    x=var[s-foundvars-1];
	    current.insert(current.begin(),x);
	    if (!is_zero(derive(curg,x,contextptr)))
	      break;
	  }
	  if (s==foundvars){
	    if (is_zero(curg))
	      newsols.push_back(current);
	    continue;
	  }
	}
	// solve
	vecteur xsol(solve(curg,*x._IDNTptr,complexmode,contextptr));
	const_iterateur xt=xsol.begin(),xtend=xsol.end();
	for (;xt!=xtend;++xt){
	  // current[xpos]=*xt;
	  newsols.push_back(subst(current,*x._IDNTptr,*xt,false,contextptr));
	}
      } // end for (;st!=stend;)
      sols=newsols;
    }
    // Add var at the beginning of each solution of sols if needed
    it=sols.begin(); 
    itend=sols.end();
    for (;it!=itend;++it){
      int ss=it->_VECTptr->size();
      if (ss<s)
	*it=mergevecteur(vecteur(var.begin(),var.begin()+s-ss),*it->_VECTptr);
    }
    if (need_subst)
      sols=subst(sols,var,var_orig,false,contextptr);
    if (convertapprox)
      sols=*evalf_VECT(sols,0,1,contextptr)._VECTptr;
    return sols;
  }

  static void read_gbargs(const vecteur & v,int start,int s,gen & order,bool & with_cocoa,bool & with_f5){
    for (int i=start;i<s;++i){
      if (v[i].is_symb_of_sommet(at_equal)){
	gen & tmp=v[i]._SYMBptr->feuille;
	if (tmp.type==_VECT && tmp._VECTptr->front().type==_INT_ && tmp._VECTptr->back().type==_INT_){
	  switch (tmp._VECTptr->front().val){
	  case _WITH_COCOA:
	    with_cocoa=tmp._VECTptr->back().val!=0;
	    break;
	  case _WITH_F5:
	    with_f5=tmp._VECTptr->back().val!=0;
	    break;
	  }
	}
      }
      if (v[i].type==_INT_ && v[i].subtype==_INT_GROEBNER){
	switch (v[i].val){
	case _WITH_COCOA:
	  with_cocoa=true;
	  break;
	case _WITH_F5:
	  with_f5=true;
	  break;
	default:
	  order=v[i].val;
	}
      }
    }
  }


  void change_monomial_order(polynome & p,const gen & order){
    switch (order.val){
      // should be strict, but does not matter since monomials are !=
    case _REVLEX_ORDER: 
      p.is_strictly_greater=total_revlex_is_greater<deg_t>;
      p.m_is_greater=std::ptr_fun(m_total_revlex_is_greater<gen>);
      break;
    case _TDEG_ORDER:
      p.is_strictly_greater=total_lex_is_greater<deg_t>;
      p.m_is_greater=std::ptr_fun(m_total_lex_is_greater<gen>);
      break;
    }
    p.tsort();
  }

  static void change_monomial_order(vectpoly & eqp,const gen & order){
    // change polynomial order
    if (order.type==_INT_ && order.val){
      vectpoly::iterator it=eqp.begin(),itend=eqp.end();
      for (;it!=itend;++it){
	change_monomial_order(*it,order);
      }
    }
  }

  // gbasis([Pi],[vars]) -> [Pi']
  gen _gbasis(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type!=_VECT)
      return symbolic(at_gbasis,args);
    vecteur & v = *args._VECTptr;
    int s=v.size();
    if (s<2)
      return gentoofewargs("gbasis");
    if ( (v[0].type!=_VECT) || (v[1].type!=_VECT) )
      return gensizeerr(contextptr);
    vecteur l=vecteur(1,v[1]);
    alg_lvar(v[0],l);
    // v[2] will serve for ordering
    gen order=_PLEX_ORDER; // _REVLEX_ORDER;
    bool with_f5=false,with_cocoa=true;
    read_gbargs(v,2,s,order,with_cocoa,with_f5);
    // convert eq to polynomial
    vecteur eq_in(*e2r(v[0],l,contextptr)._VECTptr);
    vectpoly eqp;
    if (!vecteur2vector_polynome(eq_in,l,eqp))
      return vecteur(1,plus_one);
    gen coeff;
    environment env ;
    if (!eqp.empty() && coefftype(eqp.front(),coeff)==_MOD){
      with_cocoa = false;
      env.moduloon = true;
      env.modulo = *(coeff._MODptr+1);
      env.pn=env.modulo;
      vectpoly::iterator it=eqp.begin(),itend=eqp.end();
      for (;it!=itend;++it)
	*it=unmodularize(*it);
    }
    else
      env.moduloon = false;
    if (!with_cocoa)
      change_monomial_order(eqp,order);
    vectpoly eqpr(gbasis(eqp,order,with_cocoa,with_f5,&env));
    vecteur res;
    vectpoly::const_iterator it=eqpr.begin(),itend=eqpr.end();
    res.reserve(itend-it);
    for (;it!=itend;++it)
      res.push_back(r2e(*it,l,contextptr));
    return res;
  }
  static const char _gbasis_s []="gbasis";
  static define_unary_function_eval (__gbasis,&_gbasis,_gbasis_s);
  define_unary_function_ptr5( at_gbasis ,alias_at_gbasis,&__gbasis,0,true);
  
  // greduce(P,[gbasis],[vars])
  gen _greduce(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type!=_VECT)
      return symbolic(at_gbasis,args);
    vecteur v = *args._VECTptr;
    int s=v.size();
    if (s<2)
      return gentoofewargs("greduce");
    if (s<3)
      v.push_back(lidnt(v[1]));
    if (v[1].type!=_VECT) 
      v[1]=vecteur(1,v[1]);
    v[1]=remove_equal(v[1]);
    if (v[2].type!=_VECT)
      return gensizeerr(contextptr);
    vecteur l=vecteur(1,v[2]);
    alg_lvar(v[0],l);
    alg_lvar(v[1],l);
    // v[3] will serve for ordering
    gen order=_PLEX_ORDER; // _REVLEX_ORDER;
    bool with_f5=false,with_cocoa=true;
    read_gbargs(v,3,s,order,with_cocoa,with_f5);
    gen eq(e2r(v[0],l,contextptr));
    if (eq.type!=_POLY)
      return v[0];
    gen coeff;
    environment env ;
    if (coefftype(*eq._POLYptr,coeff)==_MOD){
      with_cocoa = false;
      env.moduloon = true;
      env.modulo = *(coeff._MODptr+1);
      env.pn=env.modulo;
    }
    else
      env.moduloon = false;
    vecteur eq_in(*e2r(v[1],l,contextptr)._VECTptr);
    vectpoly eqp;
    if (!vecteur2vector_polynome(eq_in,l,eqp))
      return zero;
    change_monomial_order(eqp,order);
    polynome p(*eq._POLYptr);
    change_monomial_order(p,order);
    vectpoly rescocoa;
    if (!env.moduloon && cocoa_greduce(vectpoly(1,p),eqp,order,rescocoa))
      return r2e(rescocoa.front(),l,contextptr);
    gen C(eq._POLYptr->constant_term());
    // FIXME: get constant term, substract one to get the correct constant
    eq=eq-C+plus_one;
    // polynome res(env.moduloon?reduce(p,eqp.begin(),eqp.end(),&env):reducegb(p,eqp.begin(),eqp.end(),&env));
    polynome res(reduce(p,eqp.begin(),eqp.end(),&env));
    gen C1(res.constant_term());
    if (env.moduloon){
      res=invmod(C1,env.modulo)*res;
      modularize(res,env.modulo);
    }
    else
      res=res/C1;
    return r2e(res-plus_one,l,contextptr)+C;
  }
  static const char _greduce_s []="greduce";
  static define_unary_function_eval (__greduce,&_greduce,_greduce_s);
  define_unary_function_ptr5( at_greduce ,alias_at_greduce,&__greduce,0,true);
  
  // in_ideal([Pi],[gb],[vars]) -> true/false
  gen _in_ideal(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if (args.type!=_VECT)
      return gensizeerr(contextptr);
    vecteur & v = *args._VECTptr;
    int s=v.size();
    if (s<3)
      return gentoofewargs("in_ideal");
    if ( v[1].type!=_VECT || v[2].type!=_VECT )
      return gensizeerr(contextptr);
    vecteur atester=gen2vecteur(v[0]);
    vecteur l=vecteur(1,v[2]);
    alg_lvar(v[1],l);
    alg_lvar(v[0],l);
    gen order=_PLEX_ORDER; // _REVLEX_ORDER;
    bool with_f5=false,with_cocoa=true;
    read_gbargs(v,3,s,order,with_cocoa,with_f5);
    // convert eq to polynomial
    vecteur eq_in(*e2r(v[1],l,contextptr)._VECTptr);
    vecteur r(*e2r(atester,l,contextptr)._VECTptr);
    vectpoly eqp,eqr;
    if (!vecteur2vector_polynome(eq_in,l,eqp) || !vecteur2vector_polynome(r,l,eqr))
      return gensizeerr(contextptr);
    gen coeff;
    environment env ;
    if (!eqp.empty() && coefftype(eqp.front(),coeff)==_MOD){
      with_cocoa = false;
      env.moduloon = true;
      env.modulo = *(coeff._MODptr+1);
      env.pn=env.modulo;
      vectpoly::iterator it=eqp.begin(),itend=eqp.end();
      for (;it!=itend;++it)
	*it=unmodularize(*it);
    }
    else
      env.moduloon = false;
    if (!with_cocoa){
      change_monomial_order(eqp,order);
      change_monomial_order(eqr,order);
    }
    // is r in ideal eqp?
    gen res=in_ideal(eqr,eqp,order,with_cocoa,with_f5,&env);
    if (res.type==_VECT && !res._VECTptr->size()==1 && v[0].type!=_VECT)
      return res._VECTptr->front();
    return res;
  }
  static const char _in_ideal_s []="in_ideal";
  static define_unary_function_eval (__in_ideal,&_in_ideal,_in_ideal_s);
  define_unary_function_ptr5( at_in_ideal ,alias_at_in_ideal,&__in_ideal,0,true);

  // returns 0 for 0 solution, 1 for 1 solution, 2 for infinity solution
  // -1 on error
  int aspen_linsolve(const matrice & m,GIAC_CONTEXT){
    gen k=_ker(exact(gen(m),contextptr),contextptr);
    if (is_undef(k) || k.type!=_VECT) return -1;
    if (k._VECTptr->empty()) return 0;
    if (is_zero(k._VECTptr->back()._VECTptr->back()))
      return 0;
    if (k._VECTptr->size()==1)
      return 1;
    return 2;
  }
  // returns 0 for 0 solution, 1 for 1 solution, 2 for infinity solution
  // -1 on error
  int aspen_linsolve_2x2(const gen & a,const gen &b,const gen &c,
			  const gen &d,const gen & e,const gen & f,GIAC_CONTEXT){
    matrice m(makevecteur(makevecteur(a,b,c),makevecteur(d,e,f)));
    return aspen_linsolve(m,contextptr);
  }
  // returns 0 for 0 solution, 1 for 1 solution, 2 for infinity solution
  // -1 on error
  int aspen_linsolve_3x3(const gen & a,const gen &b,const gen &c,const gen &d,
			  const gen & e,const gen &f,const gen & g,const gen &h,
			  const gen & i,const gen & j,const gen &k,const gen &l,GIAC_CONTEXT){
    matrice m(makevecteur(makevecteur(a,b,c,d),makevecteur(e,f,g,h),makevecteur(i,j,k,l)));
    return aspen_linsolve(m,contextptr);
  }

#if defined(GIAC_HAS_STO_38) || defined(ConnectivityKit)
  gen fmin_cobyla(const gen & f,const vecteur & constraints,const vecteur & variables,const vecteur & guess,const gen & eps0,const gen & maxiter0,GIAC_CONTEXT){
    return gensizeerr(contextptr);
  }
#else // GIAC_HAS_STO_38
  struct gen_context {
    gen g; //  should be a vector [function,conditions,variables]
    const context * contextptr;
  };
  // state is a pointer of type gen_context
  int cobyla_giac_function(int n, int m, double *x, double *f, double *con,void *state){
    gen_context * gptr=(gen_context *)state;
    if (gptr->g.type!=_VECT || gptr->g._VECTptr->size()!=3)
      return 1; //error
    gen F=(*gptr->g._VECTptr)[0];
    vecteur conditions=gen2vecteur((*gptr->g._VECTptr)[1]);
    vecteur variables=gen2vecteur((*gptr->g._VECTptr)[2]);
    if (conditions.size()!=m || variables.size()!=n)
      return 1;
    vecteur values(n);
    for (int i=0;i<n;++i)
      values[i]=x[i];
    gen Fx=subst(F,variables,values,false,gptr->contextptr);
    Fx=evalf_double(Fx,1,gptr->contextptr);
    if (Fx.type!=_DOUBLE_)
      return 1;
    *f=Fx._DOUBLE_val;
    gen conditionsx=subst(conditions,variables,values,false,gptr->contextptr);
    if (conditionsx.type!=_VECT || conditionsx._VECTptr->size()!=m)
      return 1;
    vecteur & conditionsv=*conditionsx._VECTptr;
    for (int i=0;i<m;++i){
      gen cx=evalf_double(conditionsv[i],1,gptr->contextptr);
      if (cx.type!=_DOUBLE_)
	return 1;
      con[i]=cx._DOUBLE_val;
    }
    return 0;
  }

  // COBYLA will try to make all the values of the constraints positive.
  // So if you want to input a constraint j such as x[i] <= MAX, set:
  // con[j] = MAX - x[i]
  gen fmin_cobyla(const gen & f,const vecteur & constraints,const vecteur & variables,const vecteur & guess,const gen & eps0,const gen & maxiter0,GIAC_CONTEXT){
    vecteur con(constraints);
    iterateur it=con.begin(),itend=con.end();
    for (;it!=itend;++it){
      if (it->type!=_SYMB || it->_SYMBptr->feuille.type!=_VECT || it->_SYMBptr->feuille._VECTptr->size()!=2)
	continue;
      if (it->_SYMBptr->sommet==at_superieur_strict || it->_SYMBptr->sommet==at_superieur_egal)
	*it=it->_SYMBptr->feuille._VECTptr->front()-it->_SYMBptr->feuille._VECTptr->back();
      if (it->_SYMBptr->sommet==at_inferieur_strict || it->_SYMBptr->sommet==at_inferieur_egal)
	*it=it->_SYMBptr->feuille._VECTptr->back()-it->_SYMBptr->feuille._VECTptr->front();
    }
    gen fcv=makevecteur(f,con,variables);
    gen_context gc={fcv,contextptr};
    int n=variables.size(),m=constraints.size(),message(3),maxfun(1000);
    gen maxiter(maxiter0);
    gen eps0d=evalf_double(eps0,1,contextptr);
    if (is_greater(1,maxiter,contextptr) || is_greater(eps0d,1,contextptr))
      swapgen(maxiter,eps0d);
    if (is_integral(maxiter))
      maxfun=maxiter.val;
    if (guess.size()!=n)
      return gendimerr(contextptr);
    double x[n];
    for (int i=0;i<n;++i){
      gen tmp=evalf_double(guess[i],1,contextptr);
      if (tmp.type!=_DOUBLE_)
	return gensizeerr(contextptr);
      x[i]=tmp._DOUBLE_val;
    }
    double eps;
    if (eps0d.type==_DOUBLE_)
      eps=eps0d._DOUBLE_val;
    else {
      eps=(x[0]>1e-10?std::abs(x[0]):1)*epsilon(contextptr);
      if (eps<1e-13)
	eps=1e-13;
    }
    int cres=cobyla(n,m,x,x[0]/100,eps,message,&maxfun,cobyla_giac_function,&gc);
    if (cres==0){
      vecteur res(n);
      for (int i=0;i<n;++i)
	res[i]=x[i];
      return res;
    }
    return undef;
  }

#endif // GIAC_HAS_STO_38

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#if !defined(GIAC_HAS_STO_38) && !defined(ConnectivityKit)
/* cobyla : contrained optimization by linear approximation */

/*
 * Copyright (c) 1992, Michael J. D. Powell (M.J.D.Powell@damtp.cam.ac.uk)
 * Copyright (c) 2004, Jean-Sebastien Roy (js@jeannot.org)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * This software is a C version of COBYLA2, a contrained optimization by linear
 * approximation package developed by Michael J. D. Powell in Fortran.
 * 
 * The original source code can be found at :
 * http://plato.la.asu.edu/topics/problems/nlores.html
 */

static char const rcsid[] =
  "@(#) $Jeannot: cobyla.c,v 1.11 2004/04/18 09:51:36 js Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define abs(x) ((x) >= 0 ? (x) : -(x))

/*
 * Return code strings
 */
const char *cobyla_rc_string[6] =
{
  "N<0 or M<0",
  "Memory allocation failed",
  "Normal return from cobyla",
  "Maximum number of function evaluations reached",
  "Rounding errors are becoming damaging",
  "User requested end of minimization"
};

static int cobylb(int *n, int *m, int *mpp, double *x, double *rhobeg,
  double *rhoend, int *iprint, int *maxfun, double *con, double *sim,
  double *simi, double *datmat, double *a, double *vsig, double *veta,
  double *sigbar, double *dx, double *w, int *iact, cobyla_function *calcfc,
  void *state);
static int trstlp(int *n, int *m, double *a, double *b, double *rho,
  double *dx, int *ifull, int *iact, double *z__, double *zdota, double *vmultc,
  double *sdirn, double *dxnew, double *vmultd);

/* ------------------------------------------------------------------------ */

int cobyla(int n, int m, double *x, double rhobeg, double rhoend, int iprint,
  int *maxfun, cobyla_function *calcfc, void *state)
{
  int icon, isim, isigb, idatm, iveta, isimi, ivsig, iwork, ia, idx, mpp, rc;
  int *iact;
  double *w;

/*
 * This subroutine minimizes an objective function F(X) subject to M
 * inequality constraints on X, where X is a vector of variables that has 
 * N components. The algorithm employs linear approximations to the 
 * objective and constraint functions, the approximations being formed by 
 * linear interpolation at N+1 points in the space of the variables. 
 * We regard these interpolation points as vertices of a simplex. The 
 * parameter RHO controls the size of the simplex and it is reduced 
 * automatically from RHOBEG to RHOEND. For each RHO the subroutine tries 
 * to achieve a good vector of variables for the current size, and then 
 * RHO is reduced until the value RHOEND is reached. Therefore RHOBEG and 
 * RHOEND should be set to reasonable initial changes to and the required 
 * accuracy in the variables respectively, but this accuracy should be 
 * viewed as a subject for experimentation because it is not guaranteed. 
 * The subroutine has an advantage over many of its competitors, however, 
 * which is that it treats each constraint individually when calculating 
 * a change to the variables, instead of lumping the constraints together 
 * into a single penalty function. The name of the subroutine is derived 
 * from the phrase Constrained Optimization BY Linear Approximations. 
 *
 * The user must set the values of N, M, RHOBEG and RHOEND, and must 
 * provide an initial vector of variables in X. Further, the value of 
 * IPRINT should be set to 0, 1, 2 or 3, which controls the amount of 
 * printing during the calculation. Specifically, there is no output if 
 * IPRINT=0 and there is output only at the end of the calculation if 
 * IPRINT=1. Otherwise each new value of RHO and SIGMA is printed. 
 * Further, the vector of variables and some function information are 
 * given either when RHO is reduced or when each new value of F(X) is 
 * computed in the cases IPRINT=2 or IPRINT=3 respectively. Here SIGMA 
 * is a penalty parameter, it being assumed that a change to X is an 
 * improvement if it reduces the merit function 
 *      F(X)+SIGMA*MAX(0.0,-C1(X),-C2(X),...,-CM(X)), 
 * where C1,C2,...,CM denote the constraint functions that should become 
 * nonnegative eventually, at least to the precision of RHOEND. In the 
 * printed output the displayed term that is multiplied by SIGMA is 
 * called MAXCV, which stands for 'MAXimum Constraint Violation'. The 
 * argument MAXFUN is an int variable that must be set by the user to a 
 * limit on the number of calls of CALCFC, the purpose of this routine being 
 * given below. The value of MAXFUN will be altered to the number of calls 
 * of CALCFC that are made. The arguments W and IACT provide real and 
 * int arrays that are used as working space. Their lengths must be at 
 * least N*(3*N+2*M+11)+4*M+6 and M+1 respectively. 
 *
 * In order to define the objective and constraint functions, we require 
 * a subroutine that has the name and arguments 
 *      SUBROUTINE CALCFC (N,M,X,F,CON) 
 *      DIMENSION X(*),CON(*)  . 
 * The values of N and M are fixed and have been defined already, while 
 * X is now the current vector of variables. The subroutine should return 
 * the objective and constraint functions at X in F and CON(1),CON(2), 
 * ...,CON(M). Note that we are trying to adjust X so that F(X) is as 
 * small as possible subject to the constraint functions being nonnegative. 
 *
 * Partition the working space array W to provide the storage that is needed 
 * for the main calculation.
 */

  if (n == 0)
  {
    if (iprint>=1) fprintf(stderr, "cobyla: N==0.\n");
    *maxfun = 0;
    return 0;
  }

  if (n < 0 || m < 0)
  {
    if (iprint>=1) fprintf(stderr, "cobyla: N<0 or M<0.\n");
    *maxfun = 0;
    return -2;
  }

  /* workspace allocation */
  w = (double *) malloc((n*(3*n+2*m+11)+4*m+6)*sizeof(*w));
  if (w == NULL)
  {
    if (iprint>=1) fprintf(stderr, "cobyla: memory allocation error.\n");
    *maxfun = 0;
    return -1;
  }
  iact = (int *) malloc((m+1)*sizeof(*iact));
  if (iact == NULL)
  {
    if (iprint>=1) fprintf(stderr, "cobyla: memory allocation error.\n");
    free(w);
    *maxfun = 0;
    return -1;
  }
  
  /* Parameter adjustments */
  --iact;
  --w;
  --x;

  /* Function Body */
  mpp = m + 2;
  icon = 1;
  isim = icon + mpp;
  isimi = isim + n * n + n;
  idatm = isimi + n * n;
  ia = idatm + n * mpp + mpp;
  ivsig = ia + m * n + n;
  iveta = ivsig + n;
  isigb = iveta + n;
  idx = isigb + n;
  iwork = idx + n;
  rc = cobylb(&n, &m, &mpp, &x[1], &rhobeg, &rhoend, &iprint, maxfun,
      &w[icon], &w[isim], &w[isimi], &w[idatm], &w[ia], &w[ivsig], &w[iveta],
      &w[isigb], &w[idx], &w[iwork], &iact[1], calcfc, state);

  /* Parameter adjustments (reverse) */
  ++iact;
  ++w;

  free(w);
  free(iact);
  
  return rc;
} /* cobyla */

/* ------------------------------------------------------------------------- */
int cobylb(int *n, int *m, int *mpp, double 
    *x, double *rhobeg, double *rhoend, int *iprint, int *
    maxfun, double *con, double *sim, double *simi, 
    double *datmat, double *a, double *vsig, double *veta,
     double *sigbar, double *dx, double *w, int *iact, cobyla_function *calcfc,
     void *state)
{
  /* System generated locals */
  int sim_dim1, sim_offset, simi_dim1, simi_offset, datmat_dim1, 
      datmat_offset, a_dim1, a_offset, i__1, i__2, i__3;
  double d__1, d__2;

  /* Local variables */
  double alpha, delta, denom, tempa, barmu;
  double beta, cmin = 0.0, cmax = 0.0;
  double cvmaxm, dxsign, prerem = 0.0;
  double edgmax, pareta, prerec = 0.0, phimin, parsig = 0.0;
  double gamma;
  double phi, rho, sum = 0.0;
  double ratio, vmold, parmu, error, vmnew;
  double resmax, cvmaxp;
  double resnew, trured;
  double temp, wsig, f;
  double weta;
  int i__, j, k, l;
  int idxnew;
  int iflag = 0;
  int iptemp;
  int isdirn, nfvals, izdota;
  int ivmc;
  int ivmd;
  int mp, np, iz, ibrnch;
  int nbest, ifull, iptem, jdrop;
  int rc = 0;

/* Set the initial values of some parameters. The last column of SIM holds */
/* the optimal vertex of the current simplex, and the preceding N columns */
/* hold the displacements from the optimal vertex to the other vertices. */
/* Further, SIMI holds the inverse of the matrix that is contained in the */
/* first N columns of SIM. */

  /* Parameter adjustments */
  a_dim1 = *n;
  a_offset = 1 + a_dim1 * 1;
  a -= a_offset;
  simi_dim1 = *n;
  simi_offset = 1 + simi_dim1 * 1;
  simi -= simi_offset;
  sim_dim1 = *n;
  sim_offset = 1 + sim_dim1 * 1;
  sim -= sim_offset;
  datmat_dim1 = *mpp;
  datmat_offset = 1 + datmat_dim1 * 1;
  datmat -= datmat_offset;
  --x;
  --con;
  --vsig;
  --veta;
  --sigbar;
  --dx;
  --w;
  --iact;

  /* Function Body */
  iptem = min(*n,4);
  iptemp = iptem + 1;
  np = *n + 1;
  mp = *m + 1;
  alpha = .25;
  beta = 2.1;
  gamma = .5;
  delta = 1.1;
  rho = *rhobeg;
  parmu = 0.;
  if (*iprint >= 2) {
    fprintf(stderr,
      "cobyla: the initial value of RHO is %12.6E and PARMU is set to zero.\n",
      rho);
  }
  nfvals = 0;
  temp = 1. / rho;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sim[i__ + np * sim_dim1] = x[i__];
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
      sim[i__ + j * sim_dim1] = 0.;
      simi[i__ + j * simi_dim1] = 0.;
    }
    sim[i__ + i__ * sim_dim1] = rho;
    simi[i__ + i__ * simi_dim1] = temp;
  }
  jdrop = np;
  ibrnch = 0;

/* Make the next call of the user-supplied subroutine CALCFC. These */
/* instructions are also used for calling CALCFC during the iterations of */
/* the algorithm. */

L40:
  if (nfvals >= *maxfun && nfvals > 0) {
    if (*iprint >= 1) {
      fprintf(stderr,
        "cobyla: maximum number of function evaluations reach.\n");
    }
    rc = 1;
    goto L600;
  }
  ++nfvals;
  if (calcfc(*n, *m, &x[1], &f, &con[1], state))
  {
    if (*iprint >= 1) {
      fprintf(stderr, "cobyla: user requested end of minimization.\n");
    }
    rc = 3;
    goto L600;
  }
  resmax = 0.;
  if (*m > 0) {
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
      d__1 = resmax, d__2 = -con[k];
      resmax = max(d__1,d__2);
    }
  }
  if (nfvals == *iprint - 1 || *iprint == 3) {
    fprintf(stderr, "cobyla: NFVALS = %4d, F =%13.6E, MAXCV =%13.6E\n",
      nfvals, f, resmax);
    i__1 = iptem;
    fprintf(stderr, "cobyla: X =");
    for (i__ = 1; i__ <= i__1; ++i__) {
      if (i__>1) fprintf(stderr, "  ");
      fprintf(stderr, "%13.6E", x[i__]);
    }
    if (iptem < *n) {
      i__1 = *n;
      for (i__ = iptemp; i__ <= i__1; ++i__) {
        if (!((i__-1) % 4)) fprintf(stderr, "\ncobyla:  ");
        fprintf(stderr, "%15.6E", x[i__]);
      }
    }
    fprintf(stderr, "\n");
  }
  con[mp] = f;
  con[*mpp] = resmax;
  if (ibrnch == 1) {
    goto L440;
  }

/* Set the recently calculated function values in a column of DATMAT. This */
/* array has a column for each vertex of the current simplex, the entries of */
/* each column being the values of the constraint functions (if any) */
/* followed by the objective function and the greatest constraint violation */
/* at the vertex. */

  i__1 = *mpp;
  for (k = 1; k <= i__1; ++k) {
    datmat[k + jdrop * datmat_dim1] = con[k];
  }
  if (nfvals > np) {
    goto L130;
  }

/* Exchange the new vertex of the initial simplex with the optimal vertex if */
/* necessary. Then, if the initial simplex is not complete, pick its next */
/* vertex and calculate the function values there. */

  if (jdrop <= *n) {
    if (datmat[mp + np * datmat_dim1] <= f) {
      x[jdrop] = sim[jdrop + np * sim_dim1];
    } else {
      sim[jdrop + np * sim_dim1] = x[jdrop];
      i__1 = *mpp;
      for (k = 1; k <= i__1; ++k) {
        datmat[k + jdrop * datmat_dim1] = datmat[k + np * datmat_dim1]
            ;
        datmat[k + np * datmat_dim1] = con[k];
      }
      i__1 = jdrop;
      for (k = 1; k <= i__1; ++k) {
        sim[jdrop + k * sim_dim1] = -rho;
        temp = 0.f;
        i__2 = jdrop;
        for (i__ = k; i__ <= i__2; ++i__) {
          temp -= simi[i__ + k * simi_dim1];
        }
        simi[jdrop + k * simi_dim1] = temp;
      }
    }
  }
  if (nfvals <= *n) {
    jdrop = nfvals;
    x[jdrop] += rho;
    goto L40;
  }
L130:
  ibrnch = 1;

/* Identify the optimal vertex of the current simplex. */

L140:
  phimin = datmat[mp + np * datmat_dim1] + parmu * datmat[*mpp + np * 
      datmat_dim1];
  nbest = np;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    temp = datmat[mp + j * datmat_dim1] + parmu * datmat[*mpp + j * 
        datmat_dim1];
    if (temp < phimin) {
      nbest = j;
      phimin = temp;
    } else if (temp == phimin && parmu == 0.) {
      if (datmat[*mpp + j * datmat_dim1] < datmat[*mpp + nbest * 
          datmat_dim1]) {
        nbest = j;
      }
    }
  }

/* Switch the best vertex into pole position if it is not there already, */
/* and also update SIM, SIMI and DATMAT. */

  if (nbest <= *n) {
    i__1 = *mpp;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = datmat[i__ + np * datmat_dim1];
      datmat[i__ + np * datmat_dim1] = datmat[i__ + nbest * datmat_dim1]
          ;
      datmat[i__ + nbest * datmat_dim1] = temp;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = sim[i__ + nbest * sim_dim1];
      sim[i__ + nbest * sim_dim1] = 0.;
      sim[i__ + np * sim_dim1] += temp;
      tempa = 0.;
      i__2 = *n;
      for (k = 1; k <= i__2; ++k) {
        sim[i__ + k * sim_dim1] -= temp;
        tempa -= simi[k + i__ * simi_dim1];
      }
      simi[nbest + i__ * simi_dim1] = tempa;
    }
  }

/* Make an error return if SIGI is a poor approximation to the inverse of */
/* the leading N by N submatrix of SIG. */

  error = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
      temp = 0.;
      if (i__ == j) {
        temp += -1.;
      }
      i__3 = *n;
      for (k = 1; k <= i__3; ++k) {
        temp += simi[i__ + k * simi_dim1] * sim[k + j * sim_dim1];
      }
      d__1 = error, d__2 = abs(temp);
      error = max(d__1,d__2);
    }
  }
  if (error > .1) {
    if (*iprint >= 1) {
      fprintf(stderr, "cobyla: rounding errors are becoming damaging.\n");
    }
    rc = 2;
    goto L600;
  }

/* Calculate the coefficients of the linear approximations to the objective */
/* and constraint functions, placing minus the objective function gradient */
/* after the constraint gradients in the array A. The vector W is used for */
/* working space. */

  i__2 = mp;
  for (k = 1; k <= i__2; ++k) {
    con[k] = -datmat[k + np * datmat_dim1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      w[j] = datmat[k + j * datmat_dim1] + con[k];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = 0.;
      i__3 = *n;
      for (j = 1; j <= i__3; ++j) {
        temp += w[j] * simi[j + i__ * simi_dim1];
      }
      if (k == mp) {
        temp = -temp;
      }
      a[i__ + k * a_dim1] = temp;
    }
  }

/* Calculate the values of sigma and eta, and set IFLAG=0 if the current */
/* simplex is not acceptable. */

  iflag = 1;
  parsig = alpha * rho;
  pareta = beta * rho;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    wsig = 0.;
    weta = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      d__1 = simi[j + i__ * simi_dim1];
      wsig += d__1 * d__1;
      d__1 = sim[i__ + j * sim_dim1];
      weta += d__1 * d__1;
    }
    vsig[j] = 1. / sqrt(wsig);
    veta[j] = sqrt(weta);
    if (vsig[j] < parsig || veta[j] > pareta) {
      iflag = 0;
    }
  }

/* If a new vertex is needed to improve acceptability, then decide which */
/* vertex to drop from the simplex. */

  if (ibrnch == 1 || iflag == 1) {
    goto L370;
  }
  jdrop = 0;
  temp = pareta;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    if (veta[j] > temp) {
      jdrop = j;
      temp = veta[j];
    }
  }
  if (jdrop == 0) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      if (vsig[j] < temp) {
        jdrop = j;
        temp = vsig[j];
      }
    }
  }

/* Calculate the step to the new vertex and its sign. */

  temp = gamma * rho * vsig[jdrop];
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dx[i__] = temp * simi[jdrop + i__ * simi_dim1];
  }
  cvmaxp = 0.;
  cvmaxm = 0.;
  i__1 = mp;
  for (k = 1; k <= i__1; ++k) {
    sum = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      sum += a[i__ + k * a_dim1] * dx[i__];
    }
    if (k < mp) {
      temp = datmat[k + np * datmat_dim1];
      d__1 = cvmaxp, d__2 = -sum - temp;
      cvmaxp = max(d__1,d__2);
      d__1 = cvmaxm, d__2 = sum - temp;
      cvmaxm = max(d__1,d__2);
    }
  }
  dxsign = 1.;
  if (parmu * (cvmaxp - cvmaxm) > sum + sum) {
    dxsign = -1.;
  }

/* Update the elements of SIM and SIMI, and set the next X. */

  temp = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dx[i__] = dxsign * dx[i__];
    sim[i__ + jdrop * sim_dim1] = dx[i__];
    temp += simi[jdrop + i__ * simi_dim1] * dx[i__];
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    simi[jdrop + i__ * simi_dim1] /= temp;
  }
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    if (j != jdrop) {
      temp = 0.;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        temp += simi[j + i__ * simi_dim1] * dx[i__];
      }
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        simi[j + i__ * simi_dim1] -= temp * simi[jdrop + i__ * 
            simi_dim1];
      }
    }
    x[j] = sim[j + np * sim_dim1] + dx[j];
  }
  goto L40;

/* Calculate DX=x(*)-x(0). Branch if the length of DX is less than 0.5*RHO. */

L370:
  iz = 1;
  izdota = iz + *n * *n;
  ivmc = izdota + *n;
  isdirn = ivmc + mp;
  idxnew = isdirn + *n;
  ivmd = idxnew + *n;
  trstlp(n, m, &a[a_offset], &con[1], &rho, &dx[1], &ifull, &iact[1], &w[
      iz], &w[izdota], &w[ivmc], &w[isdirn], &w[idxnew], &w[ivmd]);
  if (ifull == 0) {
    temp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      d__1 = dx[i__];
      temp += d__1 * d__1;
    }
    if (temp < rho * .25 * rho) {
      ibrnch = 1;
      goto L550;
    }
  }

/* Predict the change to F and the new maximum constraint violation if the */
/* variables are altered from x(0) to x(0)+DX. */

  resnew = 0.;
  con[mp] = 0.;
  i__1 = mp;
  for (k = 1; k <= i__1; ++k) {
    sum = con[k];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      sum -= a[i__ + k * a_dim1] * dx[i__];
    }
    if (k < mp) {
      resnew = max(resnew,sum);
    }
  }

/* Increase PARMU if necessary and branch back if this change alters the */
/* optimal vertex. Otherwise PREREM and PREREC will be set to the predicted */
/* reductions in the merit function and the maximum constraint violation */
/* respectively. */

  barmu = 0.;
  prerec = datmat[*mpp + np * datmat_dim1] - resnew;
  if (prerec > 0.) {
    barmu = sum / prerec;
  }
  if (parmu < barmu * 1.5) {
    parmu = barmu * 2.;
    if (*iprint >= 2) {
      fprintf(stderr, "cobyla: increase in PARMU to %12.6E\n", parmu);
    }
    phi = datmat[mp + np * datmat_dim1] + parmu * datmat[*mpp + np * 
        datmat_dim1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      temp = datmat[mp + j * datmat_dim1] + parmu * datmat[*mpp + j * 
          datmat_dim1];
      if (temp < phi) {
        goto L140;
      }
      if (temp == phi && parmu == 0.f) {
        if (datmat[*mpp + j * datmat_dim1] < datmat[*mpp + np * 
            datmat_dim1]) {
          goto L140;
        }
      }
    }
  }
  prerem = parmu * prerec - sum;

/* Calculate the constraint and objective functions at x(*). Then find the */
/* actual reduction in the merit function. */

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    x[i__] = sim[i__ + np * sim_dim1] + dx[i__];
  }
  ibrnch = 1;
  goto L40;
L440:
  vmold = datmat[mp + np * datmat_dim1] + parmu * datmat[*mpp + np * 
      datmat_dim1];
  vmnew = f + parmu * resmax;
  trured = vmold - vmnew;
  if (parmu == 0. && f == datmat[mp + np * datmat_dim1]) {
    prerem = prerec;
    trured = datmat[*mpp + np * datmat_dim1] - resmax;
  }

/* Begin the operations that decide whether x(*) should replace one of the */
/* vertices of the current simplex, the change being mandatory if TRURED is */
/* positive. Firstly, JDROP is set to the index of the vertex that is to be */
/* replaced. */

  ratio = 0.;
  if (trured <= 0.f) {
    ratio = 1.f;
  }
  jdrop = 0;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    temp = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      temp += simi[j + i__ * simi_dim1] * dx[i__];
    }
    temp = abs(temp);
    if (temp > ratio) {
      jdrop = j;
      ratio = temp;
    }
    sigbar[j] = temp * vsig[j];
  }

/* Calculate the value of ell. */

  edgmax = delta * rho;
  l = 0;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    if (sigbar[j] >= parsig || sigbar[j] >= vsig[j]) {
      temp = veta[j];
      if (trured > 0.) {
        temp = 0.;
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          d__1 = dx[i__] - sim[i__ + j * sim_dim1];
          temp += d__1 * d__1;
        }
        temp = sqrt(temp);
      }
      if (temp > edgmax) {
        l = j;
        edgmax = temp;
      }
    }
  }
  if (l > 0) {
    jdrop = l;
  }
  if (jdrop == 0) {
    goto L550;
  }

/* Revise the simplex by updating the elements of SIM, SIMI and DATMAT. */

  temp = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sim[i__ + jdrop * sim_dim1] = dx[i__];
    temp += simi[jdrop + i__ * simi_dim1] * dx[i__];
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    simi[jdrop + i__ * simi_dim1] /= temp;
  }
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    if (j != jdrop) {
      temp = 0.;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        temp += simi[j + i__ * simi_dim1] * dx[i__];
      }
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        simi[j + i__ * simi_dim1] -= temp * simi[jdrop + i__ * 
            simi_dim1];
      }
    }
  }
  i__1 = *mpp;
  for (k = 1; k <= i__1; ++k) {
    datmat[k + jdrop * datmat_dim1] = con[k];
  }

/* Branch back for further iterations with the current RHO. */

  if (trured > 0. && trured >= prerem * .1) {
    goto L140;
  }
L550:
  if (iflag == 0) {
    ibrnch = 0;
    goto L140;
  }

/* Otherwise reduce RHO if it is not at its least value and reset PARMU. */

  if (rho > *rhoend) {
    rho *= .5;
    if (rho <= *rhoend * 1.5) {
      rho = *rhoend;
    }
    if (parmu > 0.) {
      denom = 0.;
      i__1 = mp;
      for (k = 1; k <= i__1; ++k) {
        cmin = datmat[k + np * datmat_dim1];
        cmax = cmin;
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          d__1 = cmin, d__2 = datmat[k + i__ * datmat_dim1];
          cmin = min(d__1,d__2);
          d__1 = cmax, d__2 = datmat[k + i__ * datmat_dim1];
          cmax = max(d__1,d__2);
        }
        if (k <= *m && cmin < cmax * .5) {
          temp = max(cmax,0.) - cmin;
          if (denom <= 0.) {
            denom = temp;
          } else {
            denom = min(denom,temp);
          }
        }
      }
      if (denom == 0.) {
        parmu = 0.;
      } else if (cmax - cmin < parmu * denom) {
        parmu = (cmax - cmin) / denom;
      }
    }
    if (*iprint >= 2) {
      fprintf(stderr, "cobyla: reduction in RHO to %12.6E and PARMU =%13.6E\n",
        rho, parmu);
    }
    if (*iprint == 2) {
      fprintf(stderr, "cobyla: NFVALS = %4d, F =%13.6E, MAXCV =%13.6E\n",
        nfvals, datmat[mp + np * datmat_dim1], datmat[*mpp + np * datmat_dim1]);

      fprintf(stderr, "cobyla: X =");
      i__1 = iptem;
      for (i__ = 1; i__ <= i__1; ++i__) {
        if (i__>1) fprintf(stderr, "  ");
        fprintf(stderr, "%13.6E", sim[i__ + np * sim_dim1]);
      }
      if (iptem < *n) {
        i__1 = *n;
        for (i__ = iptemp; i__ <= i__1; ++i__) {
          if (!((i__-1) % 4)) fprintf(stderr, "\ncobyla:  ");
          fprintf(stderr, "%15.6E", x[i__]);
        }
      }
      fprintf(stderr, "\n");
    }
    goto L140;
  }

/* Return the best calculated values of the variables. */

  if (*iprint >= 1) {
    fprintf(stderr, "cobyla: normal return.\n");
  }
  if (ifull == 1) {
    goto L620;
  }
L600:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    x[i__] = sim[i__ + np * sim_dim1];
  }
  f = datmat[mp + np * datmat_dim1];
  resmax = datmat[*mpp + np * datmat_dim1];
L620:
  if (*iprint >= 1) {
    fprintf(stderr, "cobyla: NFVALS = %4d, F =%13.6E, MAXCV =%13.6E\n",
      nfvals, f, resmax);
    i__1 = iptem;
    fprintf(stderr, "cobyla: X =");
    for (i__ = 1; i__ <= i__1; ++i__) {
      if (i__>1) fprintf(stderr, "  ");
      fprintf(stderr, "%13.6E", x[i__]);
    }
    if (iptem < *n) {
      i__1 = *n;
      for (i__ = iptemp; i__ <= i__1; ++i__) {
        if (!((i__-1) % 4)) fprintf(stderr, "\ncobyla:  ");
        fprintf(stderr, "%15.6E", x[i__]);
      }
    }
    fprintf(stderr, "\n");
  }
  *maxfun = nfvals;
  return rc;
} /* cobylb */

/* ------------------------------------------------------------------------- */
int trstlp(int *n, int *m, double *a, 
    double *b, double *rho, double *dx, int *ifull, 
    int *iact, double *z__, double *zdota, double *vmultc,
     double *sdirn, double *dxnew, double *vmultd)
{
  /* System generated locals */
  int a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
  double d__1, d__2;

  /* Local variables */
  double alpha, tempa;
  double beta;
  double optnew, stpful, sum, tot, acca, accb;
  double ratio, vsave, zdotv, zdotw, dd;
  double sd;
  double sp, ss, resold = 0.0, zdvabs, zdwabs, sumabs, resmax, optold;
  double spabs;
  double temp, step;
  int icount;
  int iout, i__, j, k;
  int isave;
  int kk;
  int kl, kp, kw;
  int nact, icon = 0, mcon;
  int nactx = 0;


/* This subroutine calculates an N-component vector DX by applying the */
/* following two stages. In the first stage, DX is set to the shortest */
/* vector that minimizes the greatest violation of the constraints */
/*   A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K=2,3,...,M, */
/* subject to the Euclidean length of DX being at most RHO. If its length is */
/* strictly less than RHO, then we use the resultant freedom in DX to */
/* minimize the objective function */
/*      -A(1,M+1)*DX(1)-A(2,M+1)*DX(2)-...-A(N,M+1)*DX(N) */
/* subject to no increase in any greatest constraint violation. This */
/* notation allows the gradient of the objective function to be regarded as */
/* the gradient of a constraint. Therefore the two stages are distinguished */
/* by MCON .EQ. M and MCON .GT. M respectively. It is possible that a */
/* degeneracy may prevent DX from attaining the target length RHO. Then the */
/* value IFULL=0 would be set, but usually IFULL=1 on return. */

/* In general NACT is the number of constraints in the active set and */
/* IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT */
/* contains a permutation of the remaining constraint indices. Further, Z is */
/* an orthogonal matrix whose first NACT columns can be regarded as the */
/* result of Gram-Schmidt applied to the active constraint gradients. For */
/* J=1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th */
/* column of Z with the gradient of the J-th active constraint. DX is the */
/* current vector of variables and here the residuals of the active */
/* constraints should be zero. Further, the active constraints have */
/* nonnegative Lagrange multipliers that are held at the beginning of */
/* VMULTC. The remainder of this vector holds the residuals of the inactive */
/* constraints at DX, the ordering of the components of VMULTC being in */
/* agreement with the permutation of the indices of the constraints that is */
/* in IACT. All these residuals are nonnegative, which is achieved by the */
/* shift RESMAX that makes the least residual zero. */

/* Initialize Z and some other variables. The value of RESMAX will be */
/* appropriate to DX=0, while ICON will be the index of a most violated */
/* constraint if RESMAX is positive. Usually during the first stage the */
/* vector SDIRN gives a search direction that reduces all the active */
/* constraint violations by one simultaneously. */

  /* Parameter adjustments */
  z_dim1 = *n;
  z_offset = 1 + z_dim1 * 1;
  z__ -= z_offset;
  a_dim1 = *n;
  a_offset = 1 + a_dim1 * 1;
  a -= a_offset;
  --b;
  --dx;
  --iact;
  --zdota;
  --vmultc;
  --sdirn;
  --dxnew;
  --vmultd;

  /* Function Body */
  *ifull = 1;
  mcon = *m;
  nact = 0;
  resmax = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
      z__[i__ + j * z_dim1] = 0.;
    }
    z__[i__ + i__ * z_dim1] = 1.;
    dx[i__] = 0.;
  }
  if (*m >= 1) {
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
      if (b[k] > resmax) {
        resmax = b[k];
        icon = k;
      }
    }
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
      iact[k] = k;
      vmultc[k] = resmax - b[k];
    }
  }
  if (resmax == 0.) {
    goto L480;
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sdirn[i__] = 0.;
  }

/* End the current stage of the calculation if 3 consecutive iterations */
/* have either failed to reduce the best calculated value of the objective */
/* function or to increase the number of active constraints since the best */
/* value was calculated. This strategy prevents cycling, but there is a */
/* remote possibility that it will cause premature termination. */

L60:
  optold = 0.;
  icount = 0;
L70:
  if (mcon == *m) {
    optnew = resmax;
  } else {
    optnew = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      optnew -= dx[i__] * a[i__ + mcon * a_dim1];
    }
  }
  if (icount == 0 || optnew < optold) {
    optold = optnew;
    nactx = nact;
    icount = 3;
  } else if (nact > nactx) {
    nactx = nact;
    icount = 3;
  } else {
    --icount;
    if (icount == 0) {
      goto L490;
    }
  }

/* If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to */
/* the active set. Apply Givens rotations so that the last N-NACT-1 columns */
/* of Z are orthogonal to the gradient of the new constraint, a scalar */
/* product being set to zero if its nonzero value could be due to computer */
/* rounding errors. The array DXNEW is used for working space. */

  if (icon <= nact) {
    goto L260;
  }
  kk = iact[icon];
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dxnew[i__] = a[i__ + kk * a_dim1];
  }
  tot = 0.;
  k = *n;
L100:
  if (k > nact) {
    sp = 0.;
    spabs = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = z__[i__ + k * z_dim1] * dxnew[i__];
      sp += temp;
      spabs += abs(temp);
    }
    acca = spabs + abs(sp) * .1;
    accb = spabs + abs(sp) * .2;
    if (spabs >= acca || acca >= accb) {
      sp = 0.;
    }
    if (tot == 0.) {
      tot = sp;
    } else {
      kp = k + 1;
      temp = sqrt(sp * sp + tot * tot);
      alpha = sp / temp;
      beta = tot / temp;
      tot = temp;
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        temp = alpha * z__[i__ + k * z_dim1] + beta * z__[i__ + kp * 
            z_dim1];
        z__[i__ + kp * z_dim1] = alpha * z__[i__ + kp * z_dim1] - 
            beta * z__[i__ + k * z_dim1];
        z__[i__ + k * z_dim1] = temp;
      }
    }
    --k;
    goto L100;
  }

/* Add the new constraint if this can be done without a deletion from the */
/* active set. */

  if (tot != 0.) {
    ++nact;
    zdota[nact] = tot;
    vmultc[icon] = vmultc[nact];
    vmultc[nact] = 0.;
    goto L210;
  }

/* The next instruction is reached if a deletion has to be made from the */
/* active set in order to make room for the new active constraint, because */
/* the new constraint gradient is a linear combination of the gradients of */
/* the old active constraints. Set the elements of VMULTD to the multipliers */
/* of the linear combination. Further, set IOUT to the index of the */
/* constraint to be deleted, but branch if no suitable index can be found. */

  ratio = -1.;
  k = nact;
L130:
  zdotv = 0.;
  zdvabs = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    temp = z__[i__ + k * z_dim1] * dxnew[i__];
    zdotv += temp;
    zdvabs += abs(temp);
  }
  acca = zdvabs + abs(zdotv) * .1;
  accb = zdvabs + abs(zdotv) * .2;
  if (zdvabs < acca && acca < accb) {
    temp = zdotv / zdota[k];
    if (temp > 0. && iact[k] <= *m) {
      tempa = vmultc[k] / temp;
      if (ratio < 0. || tempa < ratio) {
        ratio = tempa;
        iout = k;
      }
    }
    if (k >= 2) {
      kw = iact[k];
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        dxnew[i__] -= temp * a[i__ + kw * a_dim1];
      }
    }
    vmultd[k] = temp;
  } else {
    vmultd[k] = 0.;
  }
  --k;
  if (k > 0) {
    goto L130;
  }
  if (ratio < 0.) {
    goto L490;
  }

/* Revise the Lagrange multipliers and reorder the active constraints so */
/* that the one to be replaced is at the end of the list. Also calculate the */
/* new value of ZDOTA(NACT) and branch if it is not acceptable. */

  i__1 = nact;
  for (k = 1; k <= i__1; ++k) {
    d__1 = 0., d__2 = vmultc[k] - ratio * vmultd[k];
    vmultc[k] = max(d__1,d__2);
  }
  if (icon < nact) {
    isave = iact[icon];
    vsave = vmultc[icon];
    k = icon;
L170:
    kp = k + 1;
    kw = iact[kp];
    sp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sp += z__[i__ + k * z_dim1] * a[i__ + kw * a_dim1];
    }
    d__1 = zdota[kp];
    temp = sqrt(sp * sp + d__1 * d__1);
    alpha = zdota[kp] / temp;
    beta = sp / temp;
    zdota[kp] = alpha * zdota[k];
    zdota[k] = temp;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = alpha * z__[i__ + kp * z_dim1] + beta * z__[i__ + k * 
          z_dim1];
      z__[i__ + kp * z_dim1] = alpha * z__[i__ + k * z_dim1] - beta * 
          z__[i__ + kp * z_dim1];
      z__[i__ + k * z_dim1] = temp;
    }
    iact[k] = kw;
    vmultc[k] = vmultc[kp];
    k = kp;
    if (k < nact) {
      goto L170;
    }
    iact[k] = isave;
    vmultc[k] = vsave;
  }
  temp = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    temp += z__[i__ + nact * z_dim1] * a[i__ + kk * a_dim1];
  }
  if (temp == 0.) {
    goto L490;
  }
  zdota[nact] = temp;
  vmultc[icon] = 0.;
  vmultc[nact] = ratio;

/* Update IACT and ensure that the objective function continues to be */
/* treated as the last active constraint when MCON>M. */

L210:
  iact[icon] = iact[nact];
  iact[nact] = kk;
  if (mcon > *m && kk != mcon) {
    k = nact - 1;
    sp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sp += z__[i__ + k * z_dim1] * a[i__ + kk * a_dim1];
    }
    d__1 = zdota[nact];
    temp = sqrt(sp * sp + d__1 * d__1);
    alpha = zdota[nact] / temp;
    beta = sp / temp;
    zdota[nact] = alpha * zdota[k];
    zdota[k] = temp;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = alpha * z__[i__ + nact * z_dim1] + beta * z__[i__ + k * 
          z_dim1];
      z__[i__ + nact * z_dim1] = alpha * z__[i__ + k * z_dim1] - beta * 
          z__[i__ + nact * z_dim1];
      z__[i__ + k * z_dim1] = temp;
    }
    iact[nact] = iact[k];
    iact[k] = kk;
    temp = vmultc[k];
    vmultc[k] = vmultc[nact];
    vmultc[nact] = temp;
  }

/* If stage one is in progress, then set SDIRN to the direction of the next */
/* change to the current vector of variables. */

  if (mcon > *m) {
    goto L320;
  }
  kk = iact[nact];
  temp = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    temp += sdirn[i__] * a[i__ + kk * a_dim1];
  }
  temp += -1.;
  temp /= zdota[nact];
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sdirn[i__] -= temp * z__[i__ + nact * z_dim1];
  }
  goto L340;

/* Delete the constraint that has the index IACT(ICON) from the active set. */

L260:
  if (icon < nact) {
    isave = iact[icon];
    vsave = vmultc[icon];
    k = icon;
L270:
    kp = k + 1;
    kk = iact[kp];
    sp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sp += z__[i__ + k * z_dim1] * a[i__ + kk * a_dim1];
    }
    d__1 = zdota[kp];
    temp = sqrt(sp * sp + d__1 * d__1);
    alpha = zdota[kp] / temp;
    beta = sp / temp;
    zdota[kp] = alpha * zdota[k];
    zdota[k] = temp;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = alpha * z__[i__ + kp * z_dim1] + beta * z__[i__ + k * 
          z_dim1];
      z__[i__ + kp * z_dim1] = alpha * z__[i__ + k * z_dim1] - beta * 
          z__[i__ + kp * z_dim1];
      z__[i__ + k * z_dim1] = temp;
    }
    iact[k] = kk;
    vmultc[k] = vmultc[kp];
    k = kp;
    if (k < nact) {
      goto L270;
    }
    iact[k] = isave;
    vmultc[k] = vsave;
  }
  --nact;

/* If stage one is in progress, then set SDIRN to the direction of the next */
/* change to the current vector of variables. */

  if (mcon > *m) {
    goto L320;
  }
  temp = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    temp += sdirn[i__] * z__[i__ + (nact + 1) * z_dim1];
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sdirn[i__] -= temp * z__[i__ + (nact + 1) * z_dim1];
  }
  goto L340;

/* Pick the next search direction of stage two. */

L320:
  temp = 1. / zdota[nact];
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sdirn[i__] = temp * z__[i__ + nact * z_dim1];
  }

/* Calculate the step to the boundary of the trust region or take the step */
/* that reduces RESMAX to zero. The two statements below that include the */
/* factor 1.0E-6 prevent some harmless underflows that occurred in a test */
/* calculation. Further, we skip the step if it could be zero within a */
/* reasonable tolerance for computer rounding errors. */

L340:
  dd = *rho * *rho;
  sd = 0.;
  ss = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    if ((d__1 = dx[i__], abs(d__1)) >= *rho * 1e-6f) {
      d__2 = dx[i__];
      dd -= d__2 * d__2;
    }
    sd += dx[i__] * sdirn[i__];
    d__1 = sdirn[i__];
    ss += d__1 * d__1;
  }
  if (dd <= 0.) {
    goto L490;
  }
  temp = sqrt(ss * dd);
  if (abs(sd) >= temp * 1e-6f) {
    temp = sqrt(ss * dd + sd * sd);
  }
  stpful = dd / (temp + sd);
  step = stpful;
  if (mcon == *m) {
    acca = step + resmax * .1;
    accb = step + resmax * .2;
    if (step >= acca || acca >= accb) {
      goto L480;
    }
    step = min(step,resmax);
  }

/* Set DXNEW to the new variables if STEP is the steplength, and reduce */
/* RESMAX to the corresponding maximum residual if stage one is being done. */
/* Because DXNEW will be changed during the calculation of some Lagrange */
/* multipliers, it will be restored to the following value later. */

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dxnew[i__] = dx[i__] + step * sdirn[i__];
  }
  if (mcon == *m) {
    resold = resmax;
    resmax = 0.;
    i__1 = nact;
    for (k = 1; k <= i__1; ++k) {
      kk = iact[k];
      temp = b[kk];
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        temp -= a[i__ + kk * a_dim1] * dxnew[i__];
      }
      resmax = max(resmax,temp);
    }
  }

/* Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A */
/* device is included to force VMULTD(K)=0.0 if deviations from this value */
/* can be attributed to computer rounding errors. First calculate the new */
/* Lagrange multipliers. */

  k = nact;
L390:
  zdotw = 0.;
  zdwabs = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    temp = z__[i__ + k * z_dim1] * dxnew[i__];
    zdotw += temp;
    zdwabs += abs(temp);
  }
  acca = zdwabs + abs(zdotw) * .1;
  accb = zdwabs + abs(zdotw) * .2;
  if (zdwabs >= acca || acca >= accb) {
    zdotw = 0.;
  }
  vmultd[k] = zdotw / zdota[k];
  if (k >= 2) {
    kk = iact[k];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      dxnew[i__] -= vmultd[k] * a[i__ + kk * a_dim1];
    }
    --k;
    goto L390;
  }
  if (mcon > *m) {
    d__1 = 0., d__2 = vmultd[nact];
    vmultd[nact] = max(d__1,d__2);
  }

/* Complete VMULTC by finding the new constraint residuals. */

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dxnew[i__] = dx[i__] + step * sdirn[i__];
  }
  if (mcon > nact) {
    kl = nact + 1;
    i__1 = mcon;
    for (k = kl; k <= i__1; ++k) {
      kk = iact[k];
      sum = resmax - b[kk];
      sumabs = resmax + (d__1 = b[kk], abs(d__1));
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        temp = a[i__ + kk * a_dim1] * dxnew[i__];
        sum += temp;
        sumabs += abs(temp);
      }
      acca = sumabs + abs(sum) * .1f;
      accb = sumabs + abs(sum) * .2f;
      if (sumabs >= acca || acca >= accb) {
        sum = 0.f;
      }
      vmultd[k] = sum;
    }
  }

/* Calculate the fraction of the step from DX to DXNEW that will be taken. */

  ratio = 1.;
  icon = 0;
  i__1 = mcon;
  for (k = 1; k <= i__1; ++k) {
    if (vmultd[k] < 0.) {
      temp = vmultc[k] / (vmultc[k] - vmultd[k]);
      if (temp < ratio) {
        ratio = temp;
        icon = k;
      }
    }
  }

/* Update DX, VMULTC and RESMAX. */

  temp = 1. - ratio;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dx[i__] = temp * dx[i__] + ratio * dxnew[i__];
  }
  i__1 = mcon;
  for (k = 1; k <= i__1; ++k) {
    d__1 = 0., d__2 = temp * vmultc[k] + ratio * vmultd[k];
    vmultc[k] = max(d__1,d__2);
  }
  if (mcon == *m) {
    resmax = resold + ratio * (resmax - resold);
  }

/* If the full step is not acceptable then begin another iteration. */
/* Otherwise switch to stage two or end the calculation. */

  if (icon > 0) {
    goto L70;
  }
  if (step == stpful) {
    goto L500;
  }
L480:
  mcon = *m + 1;
  icon = mcon;
  iact[mcon] = mcon;
  vmultc[mcon] = 0.;
  goto L60;

/* We employ any freedom that may be available to reduce the objective */
/* function before returning a DX whose length is less than RHO. */

L490:
  if (mcon == *m) {
    goto L480;
  }
  *ifull = 0;
L500:
  return 0;
} /* trstlp */

#endif // GIAC_HAS_STO_38
