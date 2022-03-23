// -*- mode:C++ ; compile-command: "g++-3.4 -I.. -I../include -g -c vecteur.cc -DHAVE_CONFIG_H -DIN_GIAC" -*-
#include "giacPCH.h"

/*
 *  Copyright (C) 2000,7 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#include <cmath>
#include <stdexcept>
#include <map>
#include <iostream>
#include "gen.h"
#include "vecteur.h"
#include "modpoly.h"
#include "unary.h"
#include "symbolic.h"
#include "usual.h"
#include "sym2poly.h"
#include "solve.h"
#include "prog.h"
#include "subst.h"
#include "permu.h"
#include "plot.h"
#include "misc.h"
#include "ti89.h"
#include "giacintl.h"
#ifdef HAVE_LIBGSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_poly.h>
#endif

// for pocketcas compat.
#if defined(HAVE_LIBCLAPACK) && !defined(HAVE_LIBLAPACK)
#define HAVE_LIBLAPACK
// #define integer __CLPK_integer
#endif

#ifdef HAVE_LIBLAPACK
#include <f2c.h>
#include <clapack.h>
#undef abs
#undef min
#endif

#if defined __i386__ && !defined PIC && !defined __APPLE__ && !defined _I386_
// #define _I386_
// commented because it will fail with -O2 optimizations under gcc >= 4.3 
// on Ubuntu 11.04 in Mac VirtualBox
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

#if defined(GIAC_HAS_STO_38) && defined(VISUALC)
  static const int rand_max=2147483647;
#else
  static const int rand_max=RAND_MAX;
#endif

  vecteur makevecteur(const gen & a,const gen & b){
    vecteur v(2);
    v[0]=a;
    v[1]=b;
    return v;
  }

  vecteur makevecteur(const gen & a,const gen & b,const gen & c){
    vecteur v(3);
    v[0]=a;
    v[1]=b;
    v[2]=c;
    return v;
  }

  vecteur makevecteur(const gen & a){
    return vecteur(1,a);
  }

  vecteur makevecteur(const gen & a,const gen & b,const gen & c,const gen & d){
    vecteur v(4);
    v[0]=a;
    v[1]=b;
    v[2]=c;
    v[3]=d;
    return v;
  }

  vecteur makevecteur(const gen & a,const gen & b,const gen & c,const gen & d,const gen & e){
    vecteur v(5);
    v[0]=a;
    v[1]=b;
    v[2]=c;
    v[3]=d;
    v[4]=e;
    return v;
  }

  vecteur makevecteur(const gen & a,const gen & b,const gen & c,const gen & d,const gen & e,const gen & f){
    vecteur v(6);
    v[0]=a;
    v[1]=b;
    v[2]=c;
    v[3]=d;
    v[4]=e;
    v[5]=f;
    return v;
  }

  vecteur makevecteur(const gen & a,const gen & b,const gen & c,const gen & d,const gen & e,const gen & f,const gen & g){
    vecteur v(7);
    v[0]=a;
    v[1]=b;
    v[2]=c;
    v[3]=d;
    v[4]=e;
    v[5]=f;
    v[6]=g;
    return v;
  }

  vecteur makevecteur(const gen & a,const gen & b,const gen & c,const gen & d,const gen & e,const gen & f,const gen & g,const gen & h){
    vecteur v(8);
    v[0]=a;
    v[1]=b;
    v[2]=c;
    v[3]=d;
    v[4]=e;
    v[5]=f;
    v[6]=g;
    v[7]=h;
    return v;
  }

  vecteur makevecteur(const gen & a,const gen & b,const gen & c,const gen & d,const gen & e,const gen & f,const gen & g,const gen & h,const gen & i){
    vecteur v(9);
    v[0]=a;
    v[1]=b;
    v[2]=c;
    v[3]=d;
    v[4]=e;
    v[5]=f;
    v[6]=g;
    v[7]=h;
    v[8]=i;
    return v;
  }

  ref_vecteur * makenewvecteur(const gen & a){
    return new ref_vecteur(vecteur(1,a));
  }

  ref_vecteur * makenewvecteur(const gen & a,const gen & b){
    ref_vecteur *vptr=new ref_vecteur;
    vptr->v.reserve(2);
    vptr->v.push_back(a);
    vptr->v.push_back(b);
    return vptr;
  }

  ref_vecteur * makenewvecteur(const gen & a,const gen & b,const gen & c){
    ref_vecteur * vptr=new ref_vecteur;
    vptr->v.reserve(3);
    vptr->v.push_back(a);
    vptr->v.push_back(b);
    vptr->v.push_back(c);
    return vptr;
  }

  ref_vecteur * makenewvecteur(const gen & a,const gen & b,const gen & c,const gen & d){
    ref_vecteur * vptr=new ref_vecteur;
    vptr->v.reserve(4);
    vptr->v.push_back(a);
    vptr->v.push_back(b);
    vptr->v.push_back(c);
    vptr->v.push_back(d);
    return vptr;
  }

  ref_vecteur * makenewvecteur(const gen & a,const gen & b,const gen & c,const gen & d,const gen & e){
    ref_vecteur * vptr=new ref_vecteur;
    vptr->v.reserve(5);
    vptr->v.push_back(a);
    vptr->v.push_back(b);
    vptr->v.push_back(c);
    vptr->v.push_back(d);
    vptr->v.push_back(e);
    return vptr;
  }

  ref_vecteur * makenewvecteur(const gen & a,const gen & b,const gen & c,const gen & d,const gen & e,const gen & f){
    ref_vecteur * vptr=new ref_vecteur;
    vptr->v.reserve(6);
    vptr->v.push_back(a);
    vptr->v.push_back(b);
    vptr->v.push_back(c);
    vptr->v.push_back(d);
    vptr->v.push_back(e);
    vptr->v.push_back(f);
    return vptr;
  }

  ref_vecteur * makenewvecteur(const gen & a,const gen & b,const gen & c,const gen & d,const gen & e,const gen & f,const gen & g){
    ref_vecteur * vptr=new ref_vecteur;
    vptr->v.reserve(7);
    vptr->v.push_back(a);
    vptr->v.push_back(b);
    vptr->v.push_back(c);
    vptr->v.push_back(d);
    vptr->v.push_back(e);
    vptr->v.push_back(f);
    vptr->v.push_back(g);
    return vptr;
  }

  ref_vecteur * makenewvecteur(const gen & a,const gen & b,const gen & c,const gen & d,const gen & e,const gen & f,const gen & g,const gen & h){
    ref_vecteur * vptr=new ref_vecteur;
    vptr->v.reserve(8);
    vptr->v.push_back(a);
    vptr->v.push_back(b);
    vptr->v.push_back(c);
    vptr->v.push_back(d);
    vptr->v.push_back(e);
    vptr->v.push_back(f);
    vptr->v.push_back(g);
    vptr->v.push_back(h);
    return vptr;
  }

  ref_vecteur * makenewvecteur(const gen & a,const gen & b,const gen & c,const gen & d,const gen & e,const gen & f,const gen & g,const gen & h,const gen & i){
    ref_vecteur * vptr=new ref_vecteur;
    vptr->v.reserve(9);
    vptr->v.push_back(a);
    vptr->v.push_back(b);
    vptr->v.push_back(c);
    vptr->v.push_back(d);
    vptr->v.push_back(e);
    vptr->v.push_back(f);
    vptr->v.push_back(g);
    vptr->v.push_back(h);
    vptr->v.push_back(i);
    return vptr;
  }

  // make a matrix with free rows 
  // (i.e. it is possible to modify the answer in place)
  matrice makefreematrice(const matrice & m){
    matrice res(m);
    int s=m.size();
    for (int i=0;i<s;++i){
      if (m[i].type==_VECT){
	res[i]=makefreematrice(*m[i]._VECTptr);
      }
    }
    return res;
  }

  int alphaposcell(const string & s,int & r){
    int ss=s.size();
    r=0;
    int i=0;
    for (;i<ss;++i){
      if ( (s[i]>='A') && (s[i]<='Z') )
	r=r*26+(s[i]-'A')+1;
      else {
	if ( (s[i]>='a') && (s[i]<='q') )
	  r=r*26+(s[i]-'a')+1;
	else
	  break;
      }
    }
    --r;
    return i;
  }

  bool iscell(const gen & g,int & r,int & c,GIAC_CONTEXT){
    if (g.type!=_IDNT)
      return false;
    const string & s=g._IDNTptr->name();
    int ss=s.size();
    if (ss<2)
      return false;
    int i=alphaposcell(s,r);
    if (!i || (i==ss) )
      return false;
    c=0;
    for (;i<ss;++i){
      if ( (s[i]>='0') && (s[i]<='9') )
	c=c*10+(s[i]-'0');
      else
	break;
    }
    if (xcas_mode(contextptr))
      --c;
    return (i==ss);
  }

  // find all identifiers in g, check if they are of the form
  // Alpha_number, replace them by spread(i,j) if this is the case
  gen spread_convert(const gen & g,int g_row,int g_col,GIAC_CONTEXT){
    // relative cell
    vecteur l(*_lname(g,contextptr)._VECTptr);
    const_iterateur it=l.begin(),itend=l.end();
    vecteur sub_in,sub_out;
    int r,c;
    for (;it!=itend;++it){
      if (iscell(*it,c,r,contextptr)){
	sub_in.push_back(*it);
	sub_out.push_back(symbolic(at_cell,makevecteur(makevecteur(r-g_row),makevecteur(c-g_col))));
      }
    }
    // absolute cell
    l=lop(g,at_dollar);
    itend=l.end();
    for (it=l.begin();it!=itend;++it){
      gen & f=it->_SYMBptr->feuille;
      // cerr << "absolute cell "<< f << endl;
      if ( (f.type!=_VECT) ){
	if (iscell(f,c,r,contextptr)){
	  sub_in.push_back(*it);
	  sub_out.push_back(symbolic(at_cell,makevecteur(makevecteur(r-g_row),c)));
	}
	continue;
      }
      vecteur & v=*f._VECTptr;
      if (v.size()==2){
	gen & a=v.front();
	gen & b=v.back();
	// cerr << "absolute cell "<< a << " " << b <<endl;
	if (b.type!=_INT_)
	  continue;
	if (xcas_mode(contextptr))
	  r=b.val-1;
	else
	  r=b.val;
	if (a.type==_IDNT){
	  const string & chaine=a._IDNTptr->name();
	  int i=alphaposcell(chaine,c);
	  if (i==signed(chaine.size())){
	    sub_in.push_back(*it);
	    sub_out.push_back(symbolic(at_cell,makevecteur(r,makevecteur(c-g_col))));
	  }
	  continue;
	}
	if ( (a.type==_SYMB) && (a._SYMBptr->sommet==at_dollar) && (a._SYMBptr->feuille.type==_IDNT) ){
	  const string & chaine = a._SYMBptr->feuille._IDNTptr->name();
	  int i=alphaposcell(chaine,c);
	  if (i==signed(chaine.size())){
	    sub_in.push_back(*it);
	    sub_out.push_back(symbolic(at_cell,makevecteur(r,c)));
	  }
	}
      }
    }
    if (sub_in.empty())
      return g;
    gen tmp(quotesubst(g,sub_in,sub_out,contextptr));
    tmp.subtype=_SPREAD__SYMB;
    return tmp;
  }


  string printcell(const vecteur & v,GIAC_CONTEXT){
    // cerr << "printcell" << printcell_current_row << " " << printcell_current_col << " " << v << endl;
    string debut,tmp,fin;
    int i;
    // Note: in popular spreadsheet, the column index comes before the row
    // Therefore we translate v.back before v.front
    if (v.back().type==_INT_){
      i=v.back().val;
      debut="$";
    }
    else 
      i=v.back()._VECTptr->front().val+printcell_current_col(contextptr);
    if (i<0)
      return print_INT_(i);
    for(int j=0;;++j){
      tmp=char('A'+i%26-(j!=0))+tmp;
      i=i/26;
      if (!i)
	break;
    }
    debut=debut+tmp;
    if (v.front().type==_INT_){
      i=v.front().val;
      debut=debut+"$";
    }
    else 
      i=v.front()._VECTptr->front().val+printcell_current_row(contextptr);
    if (xcas_mode(contextptr))
      ++i;
    if (i<0)
      return debut+print_INT_(i);
    for (;;){
      fin=char('0'+i%10)+fin;
      i=i/10;
      if (!i)
	break;
    }
    return debut+fin;
  }

  string printascell(const gen & feuille,const char * sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=2) )
      return sommetstr+("("+feuille.print(contextptr)+")");
    return printcell(*feuille._VECTptr,contextptr);
  }
  gen _cell(const gen & args,GIAC_CONTEXT){
    if ( args.type==_STRNG && args.subtype==-1) return  args;
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2) )
      return gensizeerr(contextptr);
    return symbolic(at_cell,args);
  }
  static const char _cell_s []="cell";
  static define_unary_function_eval2 (__cell,&_cell,_cell_s,&printascell);
  define_unary_function_ptr5( at_cell ,alias_at_cell,&__cell,0,true);

  static void lcell(const gen & g,vecteur & res){
    if (g.type==_VECT){
      if (g.subtype==_CELL__VECT){
	if (res.empty())
	  res=*g._VECTptr;
	else { // assumes g._VECTptr has much more elements than res
	  vecteur tmp=res;
	  res=*g._VECTptr;
	  const_iterateur it=tmp.begin(),itend=tmp.end();
	  for (;it!=itend;++it){
	    if (!equalposcomp(res,*it))
	      res.push_back(*it);
	  }
	}
      }
      else {
	const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
	for (;it!=itend;++it)
	  lcell(*it,res);
      }
    }
    if (g.type==_SYMB){
      if (g._SYMBptr->sommet==at_cell || g._SYMBptr->sommet==at_deuxpoints){
	if (!equalposcomp(res,g))
	  res.push_back(g);
      }
      else
	lcell(g._SYMBptr->feuille,res);
    }
  }

  vecteur lcell(const gen & g){
    vecteur res;
    lcell(g,res);
    return res;
  }

  // given g=cell() or its argument at row i, column j 
  // return 0 if not a cell, 1 if a cell, then compute r and c s.t. g refers to (r,c), 
  // return 2 if g is e.g. A1:B4 compute ref of A1 and B4
  int cell2pos(const gen & g,int i,int j,int & r,int & c,int & r2,int & c2){
    if (g.is_symb_of_sommet(at_deuxpoints) && g._SYMBptr->feuille.type==_VECT ){
      vecteur & gf=*g._SYMBptr->feuille._VECTptr;
      if (gf.size()!=2)
	return 0;
      int r1,c1;
      if (cell2pos(gf[0],i,j,r,c,r1,c1)==1 && cell2pos(gf[1],i,j,r2,c2,r1,c1)==1)
	return 2;
      return 0;
    }
    vecteur v;
    if ( (g.type==_SYMB) && (g._SYMBptr->sommet==at_cell))
      v=*g._SYMBptr->feuille._VECTptr;
    else {
      if ( (g.type!=_VECT) || (g._VECTptr->size()!=2) )
	return 0;
      v=*g._VECTptr;
    }
    if (v.front().type==_INT_)
      r=v.front().val;
    else
      r=i+v.front()._VECTptr->front().val;
    if (v.back().type==_INT_)
      c=v.back().val;
    else
      c=j+v.back()._VECTptr->front().val;
    return 1;
  }

  // return cell(r,c) argument at (i,j) with same absolute/relative addressing
  // as g
  gen pos2cell(const gen & g,int i,int j,int r,int c,int r2,int c2){
    if (g.is_symb_of_sommet(at_deuxpoints) && g._SYMBptr->feuille.type==_VECT){
      vecteur & gf=*g._SYMBptr->feuille._VECTptr;
      if (gf.size()!=2)
	return gensizeerr(gettext("pos2cell"));
      return symbolic(at_deuxpoints,makevecteur(pos2cell(gf[0],i,j,r,c,r,c),pos2cell(gf[1],i,j,r2,c2,r2,c2)));
    }
    vecteur v;
    if ( (g.type==_SYMB) && (g._SYMBptr->sommet==at_cell))
      v=*g._SYMBptr->feuille._VECTptr;
    else {
      if ( (g.type!=_VECT) || (g._VECTptr->size()!=2) )
	return gensizeerr(gettext("pos2cell"));
      v=*g._VECTptr;
    }
    vecteur w(2);
    if (v.front().type==_INT_)
      w.front()=r;
    else
      w.front()=vecteur(1,r-i);
    if (v.back().type==_INT_)
      w.back()=c;
    else
      w.back()=vecteur(1,c-j);
    return _cell(w,context0);
  }

  gen freecopy(const gen & g){
    if (g.type!=_VECT)
      return g;
    else
      return *g._VECTptr;
  }

  // insert nrows/ncols of fill in m, e.g. fill= [0,0,2] for a spreadsheet
  // or ["","",2] or 0 for a matrix
  matrice matrice_insert(const matrice & m,int insert_row,int insert_col,int nrows,int ncols,const gen & fill,GIAC_CONTEXT){
    int r,c,cell_r,cell_c;
    int decal_i=0,decal_j;
    mdims(m,r,c);
    matrice res;
    res.reserve(r+nrows);
    // i,j position in the old matrix; i+decal_i,i+decal_j in the new
    for (int i=0;i<r;++i){
      vecteur tmp;
      tmp.reserve(c+ncols);
      if (i==insert_row){ // insert nrows of fill
	for (int j=0;j<nrows;++j){
	  // we must recreate the line each time to have a free line
	  for (int k=0;k<c+ncols;++k)
	    tmp.push_back(freecopy(fill));
	  res.push_back(tmp); 
	  tmp.clear();
	}
	decal_i=nrows;
      }
      decal_j=0;
      for (int j=0;j<c;++j){
	if (j==insert_col){
	  for (int k=0;k<ncols;k++)
	    tmp.push_back(freecopy(fill));
	  decal_j=ncols;
	}
	gen g=m[i][j];
	// find all cells in g
	vecteur sub_in(lcell(g[0])),sub_out;
	if (sub_in.empty()){
	  tmp.push_back(g);
	  continue;
	}
	const_iterateur it=sub_in.begin(),itend=sub_in.end();
	for (;it!=itend;++it){
	  int cell_r2,cell_c2,type=cell2pos(*it,i,j,cell_r,cell_c,cell_r2,cell_c2);
	  if (type){
	    if (cell_r>=insert_row)
	      cell_r += nrows;
	    if (cell_c>=insert_col)
	      cell_c += ncols;
	    if (cell_r2>=insert_row)
	      cell_r2 += nrows;
	    if (cell_c2>=insert_col)
	      cell_c2 += ncols;
	    sub_out.push_back(pos2cell(*it,i+decal_i,j+decal_j,cell_r,cell_c,cell_r2,cell_c2));
	  }
	  else
	    sub_out.push_back(*it);
	}
	g=quotesubst(g,sub_in,sub_out,contextptr);
	if (g.type==_VECT && !g._VECTptr->empty())
	  g._VECTptr->front().subtype=m[i][j][0].subtype;
	tmp.push_back(g);
      } // end for j
      res.push_back(tmp);
    } // end for i
    return res;
  }

  // erase nrows/ncols
  matrice matrice_erase(const matrice & m,int insert_row,int insert_col,int nrows,int ncols,GIAC_CONTEXT){
    int r,c,cell_r,cell_c;
    int decal_i=0,decal_j;
    mdims(m,r,c);
    matrice res;
    if ( (r<=nrows) || (c<=ncols) )
      return res;
    res.reserve(r-nrows);
    for (int i=0;i<r;++i){
      if (i==insert_row){
	i+=nrows;
	if (i>=r)
	  break;
	decal_i=nrows;
      }
      vecteur tmp;
      tmp.reserve(c-ncols);
      decal_j=0;
      for (int j=0;j<c;++j){
	if (j==insert_col){
	  j+=ncols;
	  if (j>=c)
	    break;
	  decal_j=ncols;
	}
	gen g=m[i][j];
	// find all cells in g
	vecteur sub_in(lcell(g)),sub_out;
	if (sub_in.empty()){
	  tmp.push_back(g);
	  continue;
	}
	const_iterateur it=sub_in.begin(),itend=sub_in.end();
	for (;it!=itend;++it){
	  int cell_r2,cell_c2,type=cell2pos(*it,i,j,cell_r,cell_c,cell_r2,cell_c2);
	  if (type){
	    if (cell_r>=insert_row)
	      cell_r -= nrows;
	    if (cell_c>=insert_col)
	      cell_c -= ncols;
	    if (cell_r2>=insert_row)
	      cell_r2 -= nrows;
	    if (cell_c2>=insert_col)
	      cell_c2 -= ncols;
	    sub_out.push_back(pos2cell(*it,i-decal_i,j-decal_j,cell_r,cell_c,cell_r2,cell_c2));
	  }
	  else
	    sub_out.push_back(*it);
	}
	tmp.push_back(quotesubst(g,sub_in,sub_out,contextptr));
      } // end for j
      res.push_back(tmp);
    } // end for i
    return res;
  }

  // extract submatrix
  matrice matrice_extract(const matrice & m,int insert_row,int insert_col,int nrows,int ncols){
    if ( (!nrows) || (!ncols))
      return vecteur(1,vecteur(1,gensizeerr(gettext("matrice_extract"))));
    int mr,mc;
    mdims(m,mr,mc);
    if (mr>insert_row+nrows)
      mr=insert_row+nrows;
    if (mc>insert_col+ncols)
      mc=insert_col+ncols;
    matrice res;
    res.reserve(nrows);
    for (int i=insert_row;i<mr;++i){
      const_iterateur it=m[i]._VECTptr->begin();
      res.push_back(vecteur(it+insert_col,it+mc));
    }
    return res;
  }

  // convert m to a spreadsheet matrix if necessary
  // each cell must be a vector of length 3: v[0] is the formula
  // v[1] is the value and v[2] is 0 (not evaluated), 1 (in eval), 2 (evaled)
  void makespreadsheetmatrice(matrice & m,GIAC_CONTEXT){
    int nr=m.size();
    if (!nr)
      return ;
    int nc=m.front()._VECTptr->size();
    // prepare each cell
    for (int i=0;i<nr;++i){
      gen & g=m[i];
      if (g.type!=_VECT)
	g=vecteur(nc,g);
      vecteur & v=*g._VECTptr;
      for (int j=0;j<nc;++j){
	vecteur w;
	if ((v[j].type==_VECT) && (v[j].subtype==0))
	  w=*v[j]._VECTptr;
	else
	  w=vecteur(2,v[j]);
	int s=w.size();
	if (s>3)
	  w=vecteur(w.begin(),w.begin()+3);
	if (s<1)
	  w.push_back(zero);
	if (s<3)
	  w.push_back(zero);
	if (s<2)
	  w.push_back(w.front());
	/* if (w[2].type!=_INT_)
	   w[2]=0; */
	w[0]=spread_convert(w[0],i,j,contextptr);
	v[j]=w;
      }
    }
    return ;
  }

  matrice extractmatricefromsheet(const matrice & m){
    int I=m.size();
    if (!I)
      return m;
    int J=m.front()._VECTptr->size();
    matrice res(I);
    for (int i=0;i<I;++i){
      vecteur & v=*m[i]._VECTptr;
      vecteur tmp(J);
      for (int j=0;j<J;++j){
	if ( (v[j].type==_VECT) && (v[j]._VECTptr->size()==3) )
	  tmp[j]=(*v[j]._VECTptr)[1];
	else
	  tmp[j]=v[j];
      }
      res[i]=tmp;
    }
    return res;
  }

  static gen evaldeuxpoints(const gen & args,const matrice *mptr,int cr,int cc,int & x,int & y,int & X,int & Y,GIAC_CONTEXT){
    if (args.is_symb_of_sommet(at_deuxpoints))
      return evaldeuxpoints(args._SYMBptr->feuille,mptr,cr,cc,x,y,X,Y,contextptr);
    if (args.type==_VECT && args._VECTptr->size()==2){
      vecteur & w=*args._VECTptr;
      if (!mptr){
	if (w[0].is_symb_of_sommet(at_cell) && w[1].is_symb_of_sommet(at_cell)){
	  if (w[0]._SYMBptr->feuille.type!=_VECT || w[0]._SYMBptr->feuille._VECTptr->size()!=2 || w[1]._SYMBptr->feuille.type!=_VECT || w[1]._SYMBptr->feuille._VECTptr->size()!=2 )
	    return gensizeerr(gettext("Bad cell"));
	  vecteur & w0=*w[0]._SYMBptr->feuille._VECTptr;
	  vecteur & w1=*w[1]._SYMBptr->feuille._VECTptr;
	  // Take absolute types for the returned list
	  int xm,xM,ym,yM;
	  if (w0[0].type==_VECT) 
	    xm=w0[0]._VECTptr->front().val+cr;
	  else 
	    xm=w0[0].val;
	  if (w0[1].type==_VECT) 
	    ym=w0[1]._VECTptr->front().val+cc;
	  else 
	    ym=w0[1].val;
	  // BUG 
	  if (w1[0].type==_VECT) 
	    xM=w1[0]._VECTptr->front().val+cr;
	  else 
	    xM=w1[0].val;
	  if (w1[1].type==_VECT) 
	    yM=w1[1]._VECTptr->front().val+cc;
	  else 
	    yM=w1[1].val;
	  x=min(xm,xM); X=max(xm,xM); y=min(ym,yM); Y=max(ym,yM);
	  return 1;
	}
	return 0;
      } // end if (!mptr)
      int nrows=mptr->size();
      if (X>=nrows)
	X=nrows-1;
      int ncols=nrows?mptr->front()._VECTptr->size():0;
      if (Y>=ncols)
	Y=ncols-1;
      ref_vecteur * resptr=new ref_vecteur;
      resptr->v.reserve((X-x+1)*(Y-y+1));
      ref_vecteur * vptr=0;
      for (int x0=x;x0<=X;++x0){
#ifdef SMARTPTR64
	vptr=((ref_vecteur*)(* (longlong *) &(*mptr)[x0] >> 16));
#else
	vptr=(*mptr)[x0].__VECTptr;
#endif
	for (int y0=y;y0<=Y;++y0){
	  const gen & tmp=(vptr->v)[y0][1];
	  if (tmp.type!=_STRNG || !tmp._STRNGptr->empty())
	    resptr->v.push_back(tmp);
	}
      }
      return resptr;
    }
    return mptr?symbolic(at_deuxpoints,args):zero;
  }

  // find all spread(i,j) that are in m[m_row][m_col], eval them recursively
  static gen spread_eval(matrice & m,int m_row,int m_col,GIAC_CONTEXT){
    control_c();
    if (interrupted){
      *logptr(contextptr) << "Interrupted " << m_row << " " << m_col << endl;
      return undef;
    }
    const gen & g=m[m_row][m_col][0];
    if (g.type!=_SYMB && g.type!=_VECT)
      return protecteval(g,eval_level(contextptr),contextptr);
    int & mr =spread_Row(contextptr);
    mr=m_row;
    int & mc=spread_Col(contextptr);
    mc=m_col;
    // printcell_current_row(contextptr)=m_row; printcell_current_col(contextptr)=m_col;
    vecteur v;
    lcell(g,v);
    if (v.empty()){
      gen temp=g;
      if (temp.type==_SYMB && temp.subtype==_SPREAD__SYMB)
	temp.subtype=0;
      return protecteval(temp,eval_level(contextptr),contextptr);
    }
    vecteur sub_in,sub_out;
    const_iterateur it=v.begin(),itend=v.end();
    int i,j,ms=m.size(),ws,x,y,X,Y;
    for (;it!=itend;++it){
      if (it->_SYMBptr->sommet==at_deuxpoints){
	if (is_one(evaldeuxpoints(*it,0,m_row,m_col,x,y,X,Y,contextptr))){
	  for (i=x;i<ms && i<=X;++i){
	    vecteur & w=*m[i]._VECTptr;
	    ws=w.size();
	    for (j=y;j<ws && j<=Y;++j){
	      vecteur & wj=*w[j]._VECTptr;
	      if (wj.back().val==1)
		return string2gen("Recursive eval",false);
	      if (wj.back().val==0){
		wj.back().val=1;
		wj[1]=spread_eval(m,i,j,contextptr);
		if (interrupted)
		  return undef;
		wj.back().val=2;
	      }
	    }
	  }
	  sub_in.push_back(*it);
	  sub_out.push_back(evaldeuxpoints(*it,&m,m_row,m_col,x,y,X,Y,contextptr));
	}
      } // end at_deuxpoints
      else {
	gen & gi=it->_SYMBptr->feuille._VECTptr->front();
	gen & gj=it->_SYMBptr->feuille._VECTptr->back();
	if (gi.type==_INT_)
	  i=gi.val;
	else
	  i=m_row+gi._VECTptr->front().val;
	if (gj.type==_INT_)
	  j=gj.val;
	else
	  j=m_col+gj._VECTptr->front().val;
	if ( i>=0 && i<ms ){
	  vecteur & w=*m[i]._VECTptr;
	  if ( j>=0 && j<signed(w.size()) ){
	    vecteur & wj=*w[j]._VECTptr;
	    if (wj.back().val==1)
	      return string2gen("Recursive eval",false);
	    if (wj.back().val==0){
	      wj.back().val=1;
	      wj[1]=spread_eval(m,i,j,contextptr);
	      if (interrupted)
		return undef;
	      wj.back().val=2;
	    }
	    sub_in.push_back(*it);
	    sub_out.push_back(wj[1]);
	  }
	}
      } // end at_cell
    }
    // replace evaled cell in g
    // if (sub_in.size()>=1000)
    //  cerr << endl;
    gen temp(quotesubst(g,sub_in,sub_out,contextptr));
    if (temp.type==_SYMB && temp.subtype==_SPREAD__SYMB)
      temp.subtype=0;
    mr=m_row;
    mc=m_col;
    const gen & res=protecteval(temp,eval_level(contextptr),contextptr);
    return res;
  }

  // evaluate a matrix representing a spreadsheet
  // m must be a spreadsheet matrix (see above)
  // lc will contain the list of cell dependances of m
  void spread_eval(matrice & m,GIAC_CONTEXT){
    interrupted=false;
    int nr=m.size();
    if (!nr)
      return;
    int nc=m.front()._VECTptr->size();
    // prepare for evaluation, compute list of cell and set eval flag to 0
    for (int i=0;i<nr;++i){
      vecteur & v=*m[i]._VECTptr;
      for (int j=0;j<nc;++j){
	vecteur & w=*v[j]._VECTptr;
	if (w.front().type<=_POLY){
	  w[1]=w[0];
	  w[2].val=2;
	}
	else {
	  w[2].val=0;
	}
      }
    }
    // eval
    for (int i=0;!interrupted && i<nr;++i){
      vecteur & v=*m[i]._VECTptr;
      for (int j=0;!interrupted && j<nc;++j){
	vecteur & w=*v[j]._VECTptr;
	if (w[2].val==2)
	  continue;
	w[2].val=1;
#ifndef NO_STDEXCEPT
	try {
#endif
	  w[1]=spread_eval(m,i,j,contextptr);
#ifndef NO_STDEXCEPT
	}
	catch (std::runtime_error & e){
	  w[1]=string2gen(e.what(),false);
	}
#endif
	w[2].val=2;
      }
    }
    spread_Row(-1,contextptr);
    spread_Col(-1,contextptr);
    if (interrupted)
      *logptr(contextptr) << "Spreadsheet evaluation interrupted" << endl;
  }

  vecteur mergevecteur(const vecteur & a,const vecteur & b){
    vecteur v(a);
    int as=a.size();
    int bs=b.size();
    v.reserve(as+bs);
    vecteur::const_iterator it=b.begin(),itend=b.end();
    for (;it!=itend;++it)
      v.push_back(*it);
    return v;
  }

  vecteur mergeset(const vecteur & a,const vecteur & b){
    if (a.empty())
      return b;
    vecteur v(a);
    vecteur::const_iterator it=b.begin(),itend=b.end();
    if ( (itend-it)>std::log(double(a.size()))){
      v.reserve(a.size()+(itend-it));
      for (;it!=itend;++it)
	v.push_back(*it);
      sort(v.begin(),v.end(),islesscomplexthanf);
      vecteur res(1,v.front());
      res.reserve(v.size());
      it=v.begin()+1,itend=v.end();
      for (;it!=itend;++it){
	if (*it!=res.back())
	  res.push_back(*it);
      }
      return res;
    }
    for (;it!=itend;++it){
      if (!equalposcomp(v,*it))
	v.push_back(*it);
    }
    return v;
  }

  gen makesuite(const gen & a){
    if ( (a.type==_VECT) && (a.subtype==_SEQ__VECT) )
      return a;
    else 
      return gen(vecteur(1,a),_SEQ__VECT);
  }
  
  gen makesuite_inplace(const gen & a,const gen & b){
    if (a.type!=_VECT || a.subtype!=_VECT || (b.type==_VECT && b.subtype==_SEQ__VECT))
      return makesuite(a,b);
    a._VECTptr->push_back(b);
    return a;
  }

  gen makesuite(const gen & a,const gen & b){
    if ( (a.type==_VECT) && (a.subtype==_SEQ__VECT) ){
      if ( (b.type==_VECT) && (b.subtype==_SEQ__VECT) )
	return gen(mergevecteur(*a._VECTptr,*b._VECTptr),_SEQ__VECT);
      else {
	vecteur va=*a._VECTptr;
	va.push_back(b);
	return gen(va,_SEQ__VECT);
      }
    }
    else {
      if ( (b.type==_VECT) && (b.subtype==_SEQ__VECT) ){
	vecteur vb=*b._VECTptr;
	vb.insert(vb.begin(),a);
	return gen(vb,_SEQ__VECT);
      }
      else
	return gen(makevecteur(a,b),_SEQ__VECT);
    }
  }

  // gluing is done line1 of a with line1 of b and so on
  // look at mergevecteur too
  matrice mergematrice(const matrice & a,const matrice & b){
    if (a.empty())
      return b;
    if (b.empty())
      return a;
    const_iterateur ita=a.begin(),itaend=a.end();
    const_iterateur itb=b.begin(),itbend=b.end();
    matrice res;
    res.reserve(itaend-ita);
    if (itaend-ita!=itbend-itb){
      if (debug_infolevel<1)
	return vecteur(1,vecteur(1,gendimerr(gettext("mergematrice"))));
      if (debug_infolevel<1){
	res.dbgprint();
	std_matrix<gen> M;
	matrice2std_matrix_gen(res,M);
	M.dbgprint();
      }
      return vecteur(1,vecteur(1,gendimerr(gettext("mergematrice"))));
    }
    for (;ita!=itaend;++ita,++itb){
      if (ita->type!=_VECT || itb->type!=_VECT)
	return vecteur(1,vecteur(1,gensizeerr(gettext("mergematrice"))));
      res.push_back(mergevecteur(*ita->_VECTptr,*itb->_VECTptr));
    }
    return res;
  }
  
  static complex<double> horner(const vector< complex<double> > & v, const complex<double> & c){
    vector< complex<double> > :: const_iterator it=v.begin(),itend=v.end();
    complex<double> res(0);
    for (;it!=itend;++it){
      res *= c;
      res += *it;
    }
    // cout << v << "(" << c << ")" << "=" << res << endl;
    return res;
  }

  // find a root of a polynomial with float coeffs
  gen a_root(const vecteur & v,const complex<double> & c0,double eps){
    if (v.empty())
      return gentypeerr(gettext("a_root"));
    vector< complex<double> > v_d,dv_d;
    const_iterateur it=v.begin(),itend=v.end();
    int deg=itend-it-1;
    if (deg==0)
      return gensizeerr(gettext("a_root"));
    if (deg==1)
      return -rdiv(v.back(),v.front());
    if (deg==2){ // use 2nd order equation formula
      return (-v[1]+sqrt(v[1]*v[1]-4*v[0]*v[2],context0))/(2*v[0]); // ok
    }
    v_d.reserve(deg+1);
    dv_d.reserve(deg);
    for (int d=deg;it!=itend;++it,--d){
      gen temp=it->evalf_double(1,context0); // ok
      if (temp.type==_DOUBLE_)
	v_d.push_back(temp._DOUBLE_val);
      else {
	if (temp.type!=_CPLX)
	  return undef;
	v_d.push_back(complex<double>(temp._CPLXptr->_DOUBLE_val,(temp._CPLXptr+1)->_DOUBLE_val));
      }
    }
    // Preconditionning, x->x*lambda
    // a_n x^n + .. + a_0 = a_n*lambda^n x^n + a_[n-1]*lambda^(n-1)*x^(n-1) + 
    // = a_n*lambda^n * ( x^n + a_[n-1]/a_n/lambda * x^(n-1) +
    //                    +  a_[n-2]/a_n/lambda^2 * x^(n-1) + ...)
    // take the largest ratio (a_[n-d]/a_n)^(1/d) for lambda
    double ratio=0.0,tmpratio;
    for (int d=1;d<=deg;++d){
      tmpratio=std::pow(abs(v_d[d]/v_d[0]),1.0/d);
      if (tmpratio>ratio)
	ratio=tmpratio;
    }
    double logratio=std::log(ratio);
    if (debug_infolevel)
      cerr << ratio << endl;
    bool real0=v_d[0].imag()==0;
    // Recompute coefficients
    for (int d=1;d<=deg;++d){
      bool real=real0 && v_d[d].imag()==0;
      v_d[d]=std::exp(std::log(v_d[d]/v_d[0])-d*logratio);
      if (real)
	v_d[d]=v_d[d].real();
    }
    v_d[0]=1;
    for (int d=0;d<deg;++d)
      dv_d.push_back(v_d[d]*(double)(deg-d)) ;
#ifndef __APPLE__
    if (debug_infolevel>1)
      cout << "Aroot init " << c0 << " after renormalization: " << v_d << endl << "Diff " << dv_d << endl;
#endif
    // newton method with prefactor
    complex<double> c(c0),newc,fc,newfc,fprimec,rapport;    
    double prefact=1.0;
    int maxloop=SOLVER_MAX_ITERATE;
    for (double j=1;j<1024;j=2*j,maxloop=(maxloop*3)/2){ // max 10 loop
      double prefactmult=0.5;
      fc=horner(v_d,c);
      for (int i=maxloop; i;--i){
	fprimec=horner(dv_d,c);
	if (fprimec==complex<double>(0,0))
	  break;
	rapport=fc/fprimec;
	if (abs(rapport)>1/eps) // denominator not invertible -> start elsewhere
	  break;
	newc=c-prefact*rapport;
	if (newc==c){
	  if (abs(fc)<eps)
	    return gen(real(newc)*ratio,imag(newc)*ratio);
	  break;
	}
	newfc=horner(v_d,newc);
#ifndef __APPLE__
	if (debug_infolevel>1)
	  cerr << "proot (j=" << j << "i=" << i << "), z'=" << newc << " f(z')=" << newfc << " f(z)=" << fc << " " << prefact << endl;
#endif
	if (abs(rapport)<eps)
	  return gen(real(newc)*ratio,imag(newc)*ratio);
	if (abs(newfc)>abs(fc)){
	  prefact=prefact*prefactmult;
	  // prefactmult = std::max(0.1,prefactmult*prefactmult);
	}
	else { 
	  prefactmult=0.5;
	  c=newc;
	  fc=newfc;
	  if (prefact>0.9)
	    prefact=1;
	  else
	    prefact=prefact*1.1;
	}
      }
      // c=complex<double>(rand()*j/RAND_MAX,rand()*j/RAND_MAX);
      c=complex<double>(rand()*1.0/RAND_MAX,rand()*1.0/RAND_MAX);
    }
    cerr << "proot error "+gen(v).print() << endl;
    return c;
  }

  matrice companion(const vecteur & w){
    vecteur v(w);
    if (!is_one(v.front()))
      v=divvecteur(v,v.front());
    int s=v.size()-1;
    if (s<=0)
      return vecteur(1,gendimerr());
    matrice m;
    m.reserve(s);
    for (int i=0;i<s;++i){
      vecteur w(s);
      w[s-1]=-v[s-i];
      if (i>0)
	w[i-1]=plus_one;
      m.push_back(w);
    }
    return m;
  }

  bool eigenval2(std_matrix<gen> & H,int n2,gen & l1, gen & l2,GIAC_CONTEXT){
    gen a=H[n2-2][n2-2],b=H[n2-2][n2-1],c=H[n2-1][n2-2],d=H[n2-1][n2-1];
    gen delta=a*a-2*a*d+d*d+4*b*c;
    bool save=complex_mode(contextptr);
    complex_mode(true,contextptr);
    delta=sqrt(delta,contextptr);
    complex_mode(save,contextptr);
    l1=(a+d+delta)/2;
    l2=(a+d-delta)/2;
    return is_zero(im(l1,contextptr)) && is_zero(im(l2,contextptr));
  }

  /*
  static gen balance(vecteur &v,GIAC_CONTEXT){
    // Preconditionning, x->x*lambda
    // a_n x^n + .. + a_0 = a_n*lambda^n x^n + a_[n-1]*lambda^(n-1)*x^(n-1) + 
    // = a_n*lambda^n * ( x^n + a_[n-1]/a_n/lambda * x^(n-1) +
    //                    +  a_[n-2]/a_n/lambda^2 * x^(n-1) + ...)
    // take the largest ratio (a_[n-d]/a_n)^(1/d) for lambda
    gen ratio=0,tmpratio;
    int deg=v.size()-1;
    for (int d=1;d<=deg;++d){
      tmpratio=pow(abs(v[d]/v[0]),inv(d,contextptr),contextptr);
      if (is_greater(tmpratio,ratio,contextptr))
	ratio=tmpratio;
    }
    gen logratio=ln(ratio,contextptr);
    if (debug_infolevel)
      cerr << ratio << endl;
    bool real0=is_zero(im(v[0],contextptr));
    // Recompute coefficients
    for (int d=1;d<=deg;++d){
      bool real=real0 && is_zero(im(v[d],contextptr));
      v[d]=exp(log(v[d]/v[0],contextptr)-d*logratio,contextptr);
      if (real)
	v[d]=re(v[d],contextptr);
    }
    return ratio;
  }
  */

  static bool schur_eigenvalues(matrix_double & H1,vecteur & res,double eps,GIAC_CONTEXT){
    int dim=H1.size();
    bool ans=true;
    // read eigenvalues on diagonal of H, using subdiagonal for complex pairs
    for (int i=0;i<dim;++i){
      if (i<dim-1 && std::sqrt(eps)>std::abs(H1[i+1][i])){
	if (dim*eps<std::abs(H1[i+1][i])){
#ifndef GIAC_HAS_STO_38
	  *logptr(contextptr) << "schur row " << i+1 << " " << H1[i+1][i] << endl;
#endif
	  ans=false;
	}
	// subdiagonal element is 0 -> diagonal element is an eigenvalue
	res.push_back(double(H1[i][i]));
	continue;
      }
      if (i==dim-1 && std::sqrt(eps)>std::abs(H1[i][i-1])){
	if (dim*eps<std::abs(H1[i][i-1])){
#ifndef GIAC_HAS_STO_38
	  *logptr(contextptr) << "schur row " << i << " " << H1[i][i-1] << endl;
#endif
	  ans=false;
	}
	// subdiagonal element is 0 -> diagonal element is an eigenvalue
	res.push_back(double(H1[i][i]));
	if (debug_infolevel>2)
	  cerr << "Francis algorithm Success " << res << endl;
	return ans;
      }
      // non-0, next one must be 0
      if (i<dim-2 && dim*eps<std::abs(H1[i+2][i+1])){
#ifndef GIAC_HAS_STO_38
	*logptr(contextptr) << "schur row " << i+2 << " " << H1[i+2][i+1] << endl;
#endif
	ans=false;
	if (std::sqrt(eps)<std::abs(H1[+2][i+1]))
	  continue;
      }
      if (i==dim-1)
	return false;
      giac_double l1,l2;
      if (eigenval2(H1,i+2,l1,l2)){
	res.push_back(double(l1));
	res.push_back(double(l2));
      }
      else {
	res.push_back(gen(double(l1),double(l2)));
	res.push_back(gen(double(l1),-double(l2)));
      }
      ++i;
    }
    return ans;
  }
    
  static bool proot_real1(const vecteur & v,double eps,int rprec,vecteur & res,GIAC_CONTEXT){
    matrice m(companion(v));
    int dim=m.size();
    std_matrix<gen> H,P;
    matrix_double H1,P1;
    matrice2std_matrix_gen(m,H);
    matrice2std_matrix_gen(midn(dim),P);
    if (std_matrix_gen2std_matrix_giac_double(H,H1)){
      std_matrix_gen2std_matrix_giac_double(P,P1);
      if (lapack_schur(H1,P1,false,res))
	return true;
      bool ans=francis_schur(H1,0,dim,P1,SOLVER_MAX_ITERATE,dim*eps,true,false,false);
      return ans && schur_eigenvalues(H1,res,eps,contextptr);
    }
#ifdef HAVE_LIBLAPACK
    vecteur eigenvals;
    if (lapack_schur(H,P,false,eigenvals,contextptr))
      return true;
#endif
    if (francis_schur(H,0,dim,P,SOLVER_MAX_ITERATE,dim*eps,false,false,false,contextptr)){
      eps=dim*dim*eps;
      // read eigenvalues on diagonal of H, using subdiagonal for complex pairs
      for (int i=0;i<dim;++i){
	if (i<dim-1 && is_greater(eps,abs(H[i+1][i],contextptr),contextptr)){
	  // subdiagonal element is 0 -> diagonal element is an eigenvalue
	  res.push_back(H[i][i]);
	  continue;
	}
	if (i==dim-1 && is_greater(eps,abs(H[i][i-1],contextptr),contextptr)){
	  // subdiagonal element is 0 -> diagonal element is an eigenvalue
	  res.push_back(H[i][i]);
	  if (debug_infolevel>2)
	    cerr << "Francis algorithm Success " << res << endl;
	  return true;
	}
	// non-0, next one must be 0
	if (i<dim-2 && !is_greater(eps,abs(H[i+2][i+1],contextptr),contextptr))
	  return false;
	if (i==dim-1)
	  return false;
	gen l1,l2;
	eigenval2(H,i+2,l1,l2,contextptr);
	res.push_back(l1);
	res.push_back(l2);
	++i;
      }
      if (debug_infolevel>2)
	cerr << "Francis algorithm Success " << res << endl;
      return true;
    }
    // old code using GSL
#ifdef HAVE_LIBGSL
    int vsize=v.size();
    int deg2=2*(v.size()-1);
    double *a=new double[vsize];
    for (int j=0;j<vsize;j++){      
      a[vsize-1-j]=evalf_double(v[j],1,contextptr)._DOUBLE_val;
    }
    double *z=new double[deg2];
    gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (vsize);
    int gsl=gsl_poly_complex_solve (a, vsize, w, z);
    gsl_poly_complex_workspace_free (w);
    if (gsl!=GSL_SUCCESS){
      delete [] a; delete [] z;
      return false;
    }
    for (int j=0;j<deg2;j+=2){
      res.push_back(gen(z[j],z[j+1]));
    }
    delete [] a; delete [] z;
    return true;
#else
    return false;
#endif // HAVE_LIBGSL
  }

  static bool proot_real(const vecteur & w,double eps,int rprec,vecteur & res,GIAC_CONTEXT){
    // new code using francis_schur
    gen tmp=evalf(w,1,contextptr);
    if (tmp.type!=_VECT)
      return false;
    vecteur v=*tmp._VECTptr;
    // gen prefact=balance(v,contextptr);
    bool ans=proot_real1(v,eps,rprec,res,contextptr);
    // res=multvecteur(prefact,res);
    return ans;
  }

  static bool improve_root(const vecteur & v,gen & r,int nbits,int rprec){
    int vsize=v.size();
    int deg=vsize-1;
    vecteur cur_v(v);
    double ratiod=0.0,tmpratio;
    for (int d=1;d<=deg;++d){
      tmpratio=std::pow(evalf_double(abs(cur_v[d]/cur_v[0]),1,context0)._DOUBLE_val,1.0/d);
      if (tmpratio>ratiod)
	ratiod=tmpratio;
    }
    gen ratio=accurate_evalf(gen(ratiod),nbits);
    if (ratiod>10 || ratiod<0.1){
      gen logratio=log(ratio,context0);
      if (debug_infolevel)
	cerr << ratio << endl;
      // Recompute coefficients
      for (int d=1;d<=deg;++d){
	cur_v[d]=cur_v[d]/cur_v[0]*exp(-d*logratio,context0);
      }
      cur_v[0]=1;
    }
    else
      ratio=1;
    vecteur dcur_v=derivative(cur_v);
    int j=1;
    gen prefact=accurate_evalf(plus_one,nbits);
    gen oldval,newval,newr,dr,fprimer;
    r=r/ratio;
    oldval=horner(cur_v,r);
    for (;j<SOLVER_MAX_ITERATE*vsize;j++){
      if (!(j%vsize)){
	if (is_zero(im(r,context0),context0))
	  r=r*accurate_evalf(gen(1.,1e-2),nbits);
	// random restart
	else
	  r=accurate_evalf(j/vsize*complex<double>(rand()*1.0/RAND_MAX,rand()*1.0/RAND_MAX),nbits);
	oldval=horner(cur_v,r);
	prefact=accurate_evalf(plus_one,nbits);
      }
      fprimer=horner(dcur_v,r);
      dr=oldval/fprimer;
      newr=r-prefact*dr;
      if (is_positive(-rprec-ln(abs(dr)/abs(r),context0)/std::log(2.0),context0)){
	r=ratio*newr;
	return true;
      }
      newval=horner(cur_v,newr);
      if (is_positive(abs(newval,context0)-abs(oldval,context0),context0)){
	prefact=prefact/2;
      }
      else {
	r=newr;
	oldval=newval;
	prefact=prefact*accurate_evalf(gen(1.1),nbits);
	if (is_positive(prefact-1,context0))
	  prefact=accurate_evalf(plus_one,nbits);
      }
    }
    return false;
  }

  vecteur proot(const vecteur & v,double eps,int rprec){
    int vsize=v.size();
    int deg=vsize-1;
    if (vsize<2)
      return vecteur(0);
    if (vsize==2)
      return vecteur(1,evalf(-v[1]/v[0],1,context0)); // ok
    if (vsize==3){
      gen b2=-v[1]/2;
      gen delta=sqrt(b2*b2-v[0]*v[2],context0); // ok
      return makevecteur(evalf((b2-delta)/v[0],1,context0),evalf((b2+delta)/v[0],1,context0)); // ok
    }
    bool add_conjugate=is_zero(im(v,context0),context0); // ok
    vecteur res,crystalball;
#ifdef HAVE_LIBMPFR
    int nbits = 2*(rprec+vsize);
    vecteur v_accurate(accurate_evalf(v,nbits));
    // GSL call is much faster but not very accurate
#if 0
    if (proot_real(v,eps,rprec,crystalball,context0) && int(crystalball.size())==deg){
      if (rprec<50)
	return crystalball;
    }
#else
    if (!proot_real(v,eps,rprec,crystalball,context0)){
      if (crystalball.size()!=v.size()-1)
	cerr << "Francis algorithm failure for" << v << endl;
      else
	cerr << "Francis algorithm not precise enough for" << v << endl;
    }
    else {
      if (rprec<50 && int(crystalball.size())==deg)
	return crystalball;
    }
#endif
#else // HAVE_LIBMPFR
    int nbits=45;
    rprec = 37;
    vecteur v_accurate(*evalf_double(v,1,context0)._VECTptr);
    if (proot_real(v,eps,rprec,crystalball,context0) && int(crystalball.size())==deg)
      return crystalball;
    // GSL call is much faster but not very accurate
    if (eps<1e-5)
      eps=1e-5;
#endif //HAVE_LIBMPFR
    vecteur dv_accurate(derivative(v_accurate));
    gen r,vr,dr;
    vecteur cur_v(v_accurate),dcur_v(dv_accurate),new_v;
    for (int i=0;;++i,eps*=1.1){
      if (cur_v.size()<2)
	return res;
      // gen scale=linfnorm(cur_v);
      // r=a_root(cur_v,0,scale.evalf_double(1,context0)._DOUBLE_val*eps); // ok
      if (!crystalball.empty()){
	r=crystalball.back();
	crystalball.pop_back();
      }
      else
	r=a_root(*evalf_double(cur_v,1,context0)._VECTptr,0,eps); // ok
      if (debug_infolevel)
	cerr << "Approx float root " << r << endl;
      if (is_undef(r))
	return res;
      r=accurate_evalf(r,nbits);
      int j=1;
      gen prefact=accurate_evalf(plus_one,nbits);
      gen oldval,newval,newr,fprimer;
      oldval=horner(cur_v,r);
      for (;j<SOLVER_MAX_ITERATE*vsize;j++){
	if (!(j%vsize)){
	  if (is_zero(im(r,context0),context0))
	    r=r*accurate_evalf(gen(1.,1e-2),nbits);
	  // random restart
	  else
	    r=accurate_evalf(j/vsize*complex<double>(rand()*1.0/RAND_MAX,rand()*1.0/RAND_MAX),nbits);
	  oldval=horner(cur_v,r);
	  prefact=accurate_evalf(plus_one,nbits);
	}
	fprimer=horner(dcur_v,r);
	dr=oldval/fprimer;
	newr=r-prefact*dr;
	if (is_zero(dr) || is_positive(-rprec-ln(abs(dr)/abs(r),context0)/std::log(2.0),context0)){
	  r=newr;
	  break;
	}
	newval=horner(cur_v,newr);
	if (is_strictly_positive(abs(newval,context0)-abs(oldval,context0),context0)){
	  prefact=prefact/2;
	}
	else {
	  r=newr;
	  oldval=newval;
	  prefact=prefact*accurate_evalf(gen(1.1),nbits);
	  if (is_positive(prefact-1,context0))
	    prefact=accurate_evalf(plus_one,nbits);
	}
      }
      for (j=0;j<vsize;j++){
	dr=horner(v_accurate,r)/horner(dv_accurate,r);
	r=r-dr;
	if (is_zero(dr) || is_positive(-rprec-ln(abs(dr)/abs(r),context0)/std::log(2.0),context0))
	  break;
      }
      if (j==vsize)
	return vecteur(1,gensizeerr(gettext("Proot error : no root found")));
      if (debug_infolevel)
	cerr << "Root found " << evalf_double(r,1,context0) << endl;
      if (add_conjugate && is_greater(abs(im(r,context0),context0),eps,context0) ){ // ok
	res.push_back(rprec<53?evalf_double(conj(r,context0),1,context0):conj(accurate_evalf(r,rprec),context0)); // ok
	if (!crystalball.empty()){
	  gen rcrystal=crystalball.back();
	  if (is_greater(1e-5,abs(rcrystal-res.back(),context0),context0))
	    crystalball.pop_back();
	}	
	vr=horner(cur_v,r,0,new_v);
	horner(new_v,conj(r,context0),0,cur_v); // ok
	cur_v=*(re(cur_v,context0)._VECTptr); // ok
      }
      else {
	if (add_conjugate)
	  r=re(r,context0);
	vr=horner(cur_v,r,0,new_v);
	cur_v=new_v;
      }
      res.push_back(rprec<53?evalf_double(r,1,context0):accurate_evalf(r,rprec)); // ok
      dcur_v=derivative(cur_v);
    } // end i loop
  }

  vecteur proot(const vecteur & v,double eps){
    return proot(v,eps,45);
  }

  vecteur real_proot(const vecteur & v,double eps,GIAC_CONTEXT){
    vecteur w(proot(v,eps));
    if (is_undef(w)) return w;
    vecteur res;
    const_iterateur it=w.begin(),itend=w.end();
    for (;it!=itend;++it){
      if (is_real(*it,contextptr))
	res.push_back(*it);
    }
    return res;
  }

  // eps is defined using the norm of v
  vecteur proot(const vecteur & v){
    double eps=1e-12; 
    return proot(v,eps);
  }

  gen _proot(const gen & v,GIAC_CONTEXT){
    if ( v.type==_STRNG && v.subtype==-1) return  v;
    if (v.type!=_VECT)
      return _proot(makevecteur(v,vx_var),contextptr);
    vecteur & w=*v._VECTptr;
    if (w.size()==2 && w[1].type==_IDNT){
      gen tmp=_e2r(v,contextptr);
      if (is_undef(tmp)) return tmp;
      if (tmp.type==_FRAC)
	tmp=tmp._FRACptr->num;
      if (tmp.type!=_VECT)
	return vecteur(0);
      return _proot(tmp,contextptr);
    }
    for (int i=0;i<w.size();++i){
      if (w[i].type>_REAL && w[i].type!=_FLOAT_ && w[i].type!=_CPLX)
	return gensizeerr(contextptr);
    }
    return proot(w,epsilon(contextptr),int(decimal_digits(contextptr)*3.3));
  }
  gen symb_proot(const gen & e) {
    return symbolic(at_proot,e);
  }
  static const char _proot_s []="proot";
  static define_unary_function_eval (__proot,&giac::_proot,_proot_s);
  define_unary_function_ptr5( at_proot ,alias_at_proot,&__proot,0,true);

  vecteur pcoeff(const vecteur & v){
    vecteur w(1,plus_one),new_w,somme;
    gen a,b;
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (it->type==_CPLX && it+1!=itend && is_zero(*it-conj(*(it+1),context0))){
	a=re(*it,context0);
	b=im(*it,context0);
	b=a*a+b*b;
	a=-2*a;
	w=w*makevecteur(1,a,b);
	++it;
	continue;
      }
      new_w=w;
      new_w.push_back(zero); // new_w=w*x
      mulmodpoly(w,-(*it),w); // w = -w*root
      addmodpoly(new_w,w,somme);
      w=somme;
    }
    return w;
  }
  gen _pcoeff(const gen & v,GIAC_CONTEXT){
    if ( v.type==_STRNG && v.subtype==-1) return  v;
    if (v.type!=_VECT)
      return symb_pcoeff(v);
    return gen(pcoeff(*v._VECTptr),_POLY1__VECT);
  }
  gen symb_pcoeff(const gen & e) {
    return symbolic(at_pcoeff,e);
  }
  static const char _pcoeff_s []="pcoeff";
  static define_unary_function_eval (__pcoeff,&giac::_pcoeff,_pcoeff_s);
  define_unary_function_ptr5( at_pcoeff ,alias_at_pcoeff,&__pcoeff,0,true);

  gen _peval(const gen & e,GIAC_CONTEXT){
    if ( e.type==_STRNG && e.subtype==-1) return  e;
    if (e.type!=_VECT)
      return gentypeerr(contextptr);
    vecteur & args=*e._VECTptr;
    if ( (args.size()==2) && (args.front().type==_VECT) )
      return horner(*(args.front()._VECTptr),args.back());
    if ( (args.size()!=3) || (args[1].type!=_VECT) || (args[2].type!=_VECT) )
      return gentypeerr(contextptr);
    gen pol(args.front());
    vecteur vars(*args[1]._VECTptr);
    vecteur vals(*args[2]._VECTptr);
    if (vars.size()!=vals.size())
      return gendimerr(contextptr);
    for (int i=0;i<signed(vars.size());++i){
      if (vars[i].type!=_IDNT)
	return gensizeerr(contextptr);
    }
    // convert to internal form: 
    // now put vars at the beginning of the list of variables
    vecteur lv(vars);
    lvar(e,lv);
    vecteur lv1(lv.begin()+vars.size(),lv.end());
    pol=sym2r(pol,lv,contextptr);
    gen polnum,polden;
    fxnd(pol,polnum,polden);
    for (int i=0;i<signed(vals.size());++i){
      if (debug_infolevel)
	cerr << "// Peval conversion of var " << i << " " << clock() << endl;
      vals[i]=e2r(vals[i],lv1,contextptr);
    }
    if (debug_infolevel)
      cerr << "// Peval conversion to internal form completed " << clock() << endl;
    if (polnum.type==_POLY)
      polnum=peval(*polnum._POLYptr,vals,0);
    if (polden.type==_POLY)
      polden=peval(*polden._POLYptr,vals,0);
    pol=rdiv(polnum,polden);
    return r2sym(pol,lv1,contextptr);
  }
  gen symb_peval(const gen & arg1,const gen & arg2) {
    return symbolic(at_peval,makevecteur(arg1,arg2));
  }
  static const char _peval_s []="peval";
  static define_unary_function_eval (__peval,&giac::_peval,_peval_s);
  define_unary_function_ptr5( at_peval ,alias_at_peval,&__peval,0,true);
  
  int vrows(const vecteur & a){
    return a.size();
  }

  // addvecteur is different from addmodpoly if a and b have != sizes
  // because it always start adding at the beginning of a and b
  void addvecteur(const vecteur & a,const vecteur & b,vecteur & res){
    if (&b==&res && &b!=&a){
      addvecteur(b,a,res);
      return ;
    }
    vecteur::const_iterator itb=b.begin(), itbend=b.end();
    if (&a==&res){ // in-place addition
      vecteur::iterator ita=res.begin(), itaend=res.end();
      for (;(ita!=itaend)&&(itb!=itbend);++ita,++itb){
	*ita=*ita+*itb;
      }
      for (;itb!=itbend;++itb)
	res.push_back(*itb);
      return;
    }
    vecteur::const_iterator ita=a.begin(), itaend=a.end();
    res.clear();
    res.reserve(max(itbend-itb,itaend-ita));
    for (;(ita!=itaend)&&(itb!=itbend);++ita,++itb){
      res.push_back(*ita+*itb);
    }
    for (;ita!=itaend;++ita)
      res.push_back(*ita);
    for (;itb!=itbend;++itb)
      res.push_back(*itb);
  }

  // subvecteur is different from submodpoly if a and b have != sizes
  // because it always start substr. at the beginning of a and b
  void subvecteur(const vecteur & a,const vecteur & b,vecteur & res){
    if (&b==&res){
      vecteur::const_iterator ita=a.begin(), itaend=a.end();
      vecteur::iterator itb=res.begin(), itbend=res.end();
      for (;(ita!=itaend)&&(itb!=itbend);++ita,++itb){
	*itb=*ita-*itb;
      }
      for (;ita!=itaend;++ita)
	res.push_back(*ita);
      return;
    }
    vecteur::const_iterator itb=b.begin(), itbend=b.end();
    if (&a==&res){ // in-place substract
      vecteur::iterator ita=res.begin(), itaend=res.end();
      for (;(ita!=itaend)&&(itb!=itbend);++ita,++itb){
	*ita=*ita+*itb;
      }
      for (;itb!=itbend;++itb)
	res.push_back(-*itb);
      return;
    }
    vecteur::const_iterator ita=a.begin(), itaend=a.end();
    res.clear();
    res.reserve(max(itbend-itb,itaend-ita));
    for (;(ita!=itaend)&&(itb!=itbend);++ita,++itb){
      res.push_back(*ita-*itb);
    }
    for (;ita!=itaend;++ita)
      res.push_back(*ita);
    for (;itb!=itbend;++itb)
      res.push_back(-*itb);
  }

  vecteur addvecteur(const vecteur & a,const vecteur & b){
    vecteur res;
    addvecteur(a,b,res);
    return res;
  }

  vecteur subvecteur(const vecteur & a,const vecteur & b){
    vecteur res;
    subvecteur(a,b,res);
    return res;
  }

  vecteur negvecteur(const vecteur & v){
    vecteur w;
    negmodpoly(v,w);
    return w;
  }

  gen dotvecteur(const vecteur & a,const vecteur & b){
    vecteur::const_iterator ita=a.begin(), itaend=a.end();
    vecteur::const_iterator itb=b.begin(), itbend=b.end();
    gen res,tmp;
    for (;(ita!=itaend)&&(itb!=itbend);++ita,++itb){
      type_operator_times((*ita),(*itb),tmp);
      res += tmp;
    }
    return res;
  }

  gen dotvecteur(const gen & g1,const gen & g2){
    gen a=remove_at_pnt(g1);
    gen b=remove_at_pnt(g2);
    if (a.type!=_VECT || b.type!=_VECT)
      return gensizeerr(gettext("dotvector"));
    if (a.subtype==_VECTOR__VECT)
      return dotvecteur(vector2vecteur(*a._VECTptr),b);
    if (b.subtype==_VECTOR__VECT)
      return dotvecteur(a,vector2vecteur(*b._VECTptr));
    return dotvecteur(*a._VECTptr,*b._VECTptr);
  }

  void multvecteur(const gen & a,const vecteur & b,vecteur & res){
    if (b.empty()){
      res.clear();
      return;
    }
    if (b.front().type==_VECT && ckmatrix(b)){
      vecteur temp;
      const_iterateur it=b.begin(),itend=b.end();
      res.clear();
      res.reserve(itend-it);
      for (;it!=itend;++it){
	if (it->type==_VECT){
	  multvecteur(a,*it->_VECTptr,temp);
	  res.push_back(temp);
	}
	else
	  res.push_back(a*(*it));
      }
      return;
    }
    if (is_zero(a,context0)){
      const_iterateur it=b.begin(),itend=b.end();
      res.clear();
      res.reserve(itend-it);
      for (;it!=itend;++it)
	res.push_back((*it)*zero);
    }
    else {
      environment * env=new environment;
      mulmodpoly(b,a,env,res);
      delete env;
    }
  }

  vecteur multvecteur(const gen & a,const vecteur & b){
    vecteur res;
    multvecteur(a,b,res);
    return res;
  }

  void divvecteur(const vecteur & b,const gen & a,vecteur & res){
    if (b.empty()){
      res.clear();
      return;
    }
    if (&b==&res){
      iterateur it=res.begin(),itend=res.end();
      for (;it!=itend;++it){
	if (it->type==_VECT){
	  vecteur temp;
	  divvecteur(*it->_VECTptr,a,*it->_VECTptr);
	}
	else
	  *it=rdiv(*it,a);
      }
      return;
    }
    if (b.front().type==_VECT && ckmatrix(b)){
      const_iterateur it=b.begin(),itend=b.end();
      res.clear();
      res.reserve(itend-it);
      for (;it!=itend;++it){
	if (it->type==_VECT){
	  vecteur temp;
	  divvecteur(*it->_VECTptr,a,temp);
	  res.push_back(temp);
	}
	else
	  res.push_back(rdiv(*it,a));
      }
      return;
    }
    environment * env=new environment;
    divmodpoly(b,a,res);
    delete env;
  }

  vecteur divvecteur(const vecteur & b,const gen & a){
    vecteur res;
    divvecteur(b,a,res);
    return res;
  }

  void multmatvecteur(const matrice & a,const vecteur & b,vecteur & res){
    vecteur::const_iterator ita=a.begin(), itaend=a.end();
    res.clear();
    res.reserve(itaend-ita);
    for (;ita!=itaend;++ita)
      res.push_back(dotvecteur(*(ita->_VECTptr),b));
  }

  vecteur multmatvecteur(const matrice & a,const vecteur & b){
    vecteur res;
    multmatvecteur(a,b,res);
    return res;
  }

  void multvecteurmat(const vecteur & a,const matrice & b,vecteur & res){
    matrice btran;
    mtran(b,btran);
    multmatvecteur(btran,a,res);
  }

  vecteur multvecteurmat(const vecteur & a,const matrice & b){
    vecteur res;
    multvecteurmat(a,b,res);
    return res;
  }

  gen ckmultmatvecteur(const vecteur & a,const vecteur & b){
    if (ckmatrix(a)){
      if (ckmatrix(b)){
	matrice res;
	mmultck(a,b,res);
	return _simplifier(res,context0);
      }
      // matrice * vecteur
      vecteur res;
      if (a.front()._VECTptr->size()!=b.size())
	return gendimerr(gettext("dotvecteur"));
      multmatvecteur(a,b,res);
      return _simplifier(res,context0);
    }
    if (ckmatrix(b)){
      vecteur res;
      multvecteurmat(a,b,res);
      return _simplifier(res,context0);
    }
    if (xcas_mode(context0)==3)
      return apply(a,b,prod);
    return dotvecteur(a,b);
  }

  // *********************
  // ***   Matrices    ***
  // *********************

  bool ckmatrix(const matrice & a,bool allow_embedded_vect){
    vecteur::const_iterator it=a.begin(),itend=a.end();
    if (itend==it)
      return false;
    int s=-1;
    int cur_s;
    for (;it!=itend;++it){
      if (it->type!=_VECT)
	return false;
      cur_s=it->_VECTptr->size();
      if (!cur_s)
	return false;
      if (s<0)
	s = cur_s;
      else {
	if (s!=cur_s)
	  return false;
	if (s && it->_VECTptr->front().type==_VECT && !allow_embedded_vect)
	  return false;
      }
    }
    return true;
  }

  bool ckmatrix(const matrice & a){
    return ckmatrix(a,false);
  }

  bool ckmatrix(const gen & a,bool allow_embedded_vect){
    if (a.type!=_VECT)
      return false;
    return ckmatrix(*a._VECTptr,allow_embedded_vect);
  }

  bool ckmatrix(const gen & a){
    return ckmatrix(a,false);
  }

  bool is_squarematrix(const matrice & a){
    if (!ckmatrix(a))
      return false;
    return a.size()==a.front()._VECTptr->size();
  }

  bool is_squarematrix(const gen & a){
    if (!ckmatrix(a))
      return false;
    return a._VECTptr->size()==a._VECTptr->front()._VECTptr->size();
  }

  bool is_fully_numeric(const vecteur & v){
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (!is_fully_numeric(*it))
	return false;
    }
    return true;
  }

  bool is_fully_numeric(const gen & a){
    switch (a.type){
    case _DOUBLE_: case _FLOAT_:
      return true;
    case _CPLX:
      return is_fully_numeric(*a._CPLXptr) && is_fully_numeric(*(a._CPLXptr+1));
    case _VECT:
      return is_fully_numeric(*a._VECTptr);
    case _IDNT:
      return *a._IDNTptr==_IDNT_pi;
    default:
      return false;
    }
  }

  int mrows(const matrice & a){
    return a.size();
  }

  int mcols(const matrice & a){
    return a.begin()->_VECTptr->size();
  }

  void mdims(const matrice &m,int & r,int & c){
    r=m.size();
    c=0;
    if (r){
      const gen & g=m.front();
      if (g.type==_VECT)
	c=g._VECTptr->size();
    }
  }

  void mtran(const matrice & a,matrice & res,int ncolres){
    vecteur::const_iterator it=a.begin(),itend=a.end();
    int n=itend-it; // nrows of a = ncols of res if ncolres was 0
    res.clear();
    if (!n)
      return;
    if (!ncolres)
      ncolres=n;
    int c=it->_VECTptr->size(); // ncols of a = rows of res
    res.reserve(c);
    // find begin of each row
#ifdef VISUALC
    vecteur::const_iterator * itr=new vecteur::const_iterator[ncolres];
#else
    vecteur::const_iterator itr[ncolres];
#endif
    vecteur::const_iterator * itrend= itr+ncolres;
    vecteur::const_iterator * itrcur;
    int i;
    for (i=0;(i<n) && (it!=itend);++it,++i)
      itr[i]=it->_VECTptr->begin();
    for (;(i<ncolres) ;++i)
#ifdef VISUALC
      * (int *) &itr[i]=0;
#else
      itr[i]=(vecteur::const_iterator) NULL;
#endif
    vecteur cur_row; 
    // make current row of res with currents elements of itr[]
    for (int j=0;j<c;++j){
      cur_row.clear();
      cur_row.reserve(ncolres);
      for (itrcur=itr;itrcur!=itrend;++itrcur){
	if
#ifdef VISUALC
	  (* (int *)itrcur!=0)
#else
	  (*itrcur!=(vecteur::const_iterator)NULL)
#endif
	    {
	      cur_row.push_back(**itrcur);
	      ++(*itrcur);
	    }
	else
	  cur_row.push_back(0);
      }
      res.push_back(cur_row);
    }
#ifdef VISUALC
    delete [] itr;
#endif
  }

  matrice mtran(const matrice & a){
    matrice res;
    mtran(a,res);
    return res;
  }

  gen _tran(const gen & a,GIAC_CONTEXT){
    if ( a.type==_STRNG && a.subtype==-1) return  a;
    vecteur v;
    if (!ckmatrix(a)){
      if (a.type==_VECT && !a._VECTptr->empty())
	v=vecteur(1,a);
      else
	return symb_tran(a);
    }
    else
      v=*a._VECTptr;
    matrice res;
    mtran(v,res);
    return res;
  }
  static const char _tran_s []="tran";
  static define_unary_function_eval (__tran,&giac::_tran,_tran_s);
  define_unary_function_ptr5( at_tran ,alias_at_tran,&__tran,0,true);
  
  bool matrice2std_matrix_double(const matrice & m,matrix_double & M){
    int n=m.size(),c;
    gen g;
    M.resize(n);
    for (int i=0;i<n;++i){
      const vecteur & mi=*m[i]._VECTptr;
      c=mi.size();
      std::vector<giac_double> & v =M[i];
      v.clear();
      v.reserve(c);
      const_iterateur it=mi.begin(),itend=mi.end();
      for (;it!=itend;++it){
	g=evalf(*it,1,context0);
	if (g.type==_FLOAT_){
	  v.push_back(get_double(g._FLOAT_val));
	  continue;
	}
	if (g.type!=_DOUBLE_)
	  return false;
	v.push_back(g._DOUBLE_val);
      }
    }
    return true;
  }

  giac_double dotvecteur(const std::vector<giac_double> & a,const std::vector<giac_double> & b){
    std::vector<giac_double>::const_iterator ita=a.begin(), itaend=a.end();
    std::vector<giac_double>::const_iterator itb=b.begin(), itbend=b.end();
    giac_double res=0;
    for (;(ita!=itaend)&&(itb!=itbend);++ita,++itb){
      res += (*ita)*(*itb);
    }
    return res;
  }

  bool matrice2std_matrix_complex_double(const matrice & m,matrix_complex_double & M){
    int n=m.size(),c;
    gen g;
    M.resize(n);
    for (int i=0;i<n;++i){
      const vecteur & mi=*m[i]._VECTptr;
      c=mi.size();
      std::vector<complex_double> & v =M[i];
      v.clear();
      v.reserve(c);
      const_iterateur it=mi.begin(),itend=mi.end();
      for (;it!=itend;++it){
	g=evalf_double(*it,1,context0);
	if (g.type==_CPLX && g._CPLXptr->type==_DOUBLE_ && (g._CPLXptr+1)->type==_DOUBLE_){
	  v.push_back(complex_double(g._CPLXptr->_DOUBLE_val,(g._CPLXptr+1)->_DOUBLE_val));
	  continue;
	}
	if (g.type!=_DOUBLE_)
	  return false;
	v.push_back(g._DOUBLE_val);
      }
    }
    return true;
  }

  complex_double dotvecteur(const std::vector<complex_double> & a,const std::vector<complex_double> & b){
    std::vector<complex_double>::const_iterator ita=a.begin(), itaend=a.end();
    std::vector<complex_double>::const_iterator itb=b.begin(), itbend=b.end();
    complex_double res=0;
    for (;(ita!=itaend)&&(itb!=itbend);++ita,++itb){
      res += (*ita)*(*itb);
    }
    return res;
  }

  bool mmult_float(const matrice & a,const matrice & btran,matrice & res){
    matrix_double ad,btrand;
    if (matrice2std_matrix_double(a,ad) && matrice2std_matrix_double(btran,btrand)){
      matrix_double::const_iterator ita=ad.begin(),itaend=ad.end();
      matrix_double::const_iterator itbbeg=btrand.begin(),itb,itbend=btrand.end();
      int resrows=mrows(a);
      int rescols=mrows(btran);
      res.clear();
      res.reserve(resrows);
      for (;ita!=itaend;++ita){
	res.push_back(new ref_vecteur(rescols));
	vecteur & cur_row=*res.back()._VECTptr;
	for (itb=itbbeg;itb!=itbend;++itb)
	  cur_row[itb-itbbeg]=double(dotvecteur(*ita,*itb));
      }
      return true;
    }
    matrix_complex_double zad,zbtrand;
    if (!matrice2std_matrix_complex_double(a,zad) || !matrice2std_matrix_complex_double(btran,zbtrand))
      return false;
    matrix_complex_double::const_iterator ita=zad.begin(),itaend=zad.end();
    matrix_complex_double::const_iterator itbbeg=zbtrand.begin(),itb,itbend=zbtrand.end();
    int resrows=mrows(a);
    int rescols=mrows(btran);
    res.clear();
    res.reserve(resrows);
    for (;ita!=itaend;++ita){
      res.push_back(new ref_vecteur(rescols));
      vecteur & cur_row=*res.back()._VECTptr;
      for (itb=itbbeg;itb!=itbend;++itb)
	cur_row[itb-itbbeg]=dotvecteur(*ita,*itb);
    }
    return true;
  }

  void mmult(const matrice & a,const matrice & b,matrice & res){
    matrice btran;
    mtran(b,btran);
    // now make the (dotvecteur) product of row i of a with rows of btran to get
    // row i of res
#ifndef BCD
    if ( (has_num_coeff(a) || has_num_coeff(b)) && mmult_float(a,btran,res))
      return;
#endif
    vecteur::const_iterator ita=a.begin(),itaend=a.end();
    vecteur::const_iterator itbbeg=btran.begin(),itb;//itbend=btran.end(),
    int resrows=mrows(a);
    int rescols=mrows(btran);
    res.clear();
    res.reserve(resrows);
    /* old code replaced to enhance product of sparse matrices
    vecteur cur_row;
    for (;ita!=itaend;++ita){
      cur_row.clear();
      cur_row.reserve(rescols);
      for (itb=itbbeg;itb!=itbend;++itb)
	cur_row.push_back(dotvecteur(*(ita->_VECTptr),*(itb->_VECTptr)));
      res.push_back(cur_row);
    }
    */
    int s=btran.size();
    gen tmp;
    const_iterateur it,itend;
    vector<const_iterateur> itbb(s);
    iterateur itc;
    for (;ita!=itaend;++ita){
      vecteur c(rescols,zero);
      it=ita->_VECTptr->begin();
      itend=ita->_VECTptr->end();
      itb=itbbeg;
      for (int i=0;i<s;++i,++itb)
	itbb[i]=itb->_VECTptr->begin();
      for (;it!=itend;++it){
	const gen & acur=*it;
	if (is_zero(acur,context0)){
	  int p=1;
	  ++it;
	  for (; (it!=itend) && is_zero(*it,context0);++it,++p){
	  }
	  if (it==itend)
	    break;
	  else
	    --it;
	  for (int i=0;i<s;++i)
	    itbb[i]+=p;
	}
	else {
	  itc=c.begin();
	  gen tmp;
	  for (int i=0;i<s;++itc,++(itbb[i]),++i){
	    type_operator_times(acur, *(itbb[i]),tmp);
	    *itc += tmp;
	  }
	}
      }
      res.push_back(c);
    }
  }

  matrice mmult(const matrice & a,const matrice & b){
    matrice res;
    mmult(a,b,res);
    return res;
  }

  bool mmultck(const matrice & a, const matrice & b,matrice & res){
    if (mcols(a)!=mrows(b))
      return false; 
    mmult(a,b,res);
    return true;
  }

  matrice mmultck(const matrice & a, const matrice & b){
    matrice res;
    if (!mmultck(a,b,res))
      return vecteur(1,vecteur(1,gendimerr(gettext("mmultck"))));
    return res;
  }

  gen mtrace(const matrice & a){
    gen res(0);
    vecteur::const_iterator it=a.begin(),itend=a.end();
    for (int i=0;it!=itend;++it,++i)
      res = res + (*it)[i];
    return res;
  }
  
  gen ckmtrace(const gen & a,GIAC_CONTEXT){
    if (!is_squarematrix(a))
      return symbolic(at_trace,a); // gendimerr(contextptr); required to keep trace for geometry
    return mtrace(*a._VECTptr);
  }
  static const char _trace_s []="trace";
  static define_unary_function_eval (__trace,&giac::ckmtrace,_trace_s);
  define_unary_function_ptr5( at_trace ,alias_at_trace,&__trace,0,true);

  gen common_deno(const vecteur & v){
    const_iterateur it=v.begin(),itend=v.end();
    gen lcm_deno(1);
    for (;it!=itend;++it){
      if (it->type==_FRAC)
	lcm_deno=rdiv(lcm_deno,gcd(lcm_deno,it->_FRACptr->den))*(it->_FRACptr->den);
    }
    return lcm_deno;
  }

  static gen common_num(const vecteur & v){
    const_iterateur it=v.begin(),itend=v.end();
    gen gcd_num(0);
    for (;it!=itend;++it){
      if (it->type!=_FRAC)
	gcd_num=gcd(gcd_num,*it);
    }
    return gcd_num;
  }

  static inline gen trim(const gen & a,const gen & b,double eps){
    if (eps && a.type==_DOUBLE_ && b.type==_DOUBLE_ &&
	fabs(a._DOUBLE_val)<eps*fabs(b._DOUBLE_val)) 
      return 0;
    else
      return a;
  }

  static gen exact_div(const gen & a,const gen & b){
    if (a.type==_POLY && b.type==_POLY){
      polynome *quoptr=new polynome, rem;
      if (!divrem1(*a._POLYptr,*b._POLYptr,*quoptr,rem,2)) 
	cerr << "bad quo("+a.print()+","+b.print()+")" << endl;
      gen res= *quoptr;
      // if (!is_zero(a-b*res))
      //	cerr << "Bad division" << endl;
      return res;
      polynome quo;
      if (!a._POLYptr->Texactquotient(*b._POLYptr,quo))
	cerr << "bad quo("+a.print()+","+b.print()+")" << endl;
      return quo;
    }
    return rdiv(a,b);
  }

  // v=(c1*v1+c2*v2)/c
  // Set cstart to 0, or to c+1 for lu decomposition
  void linear_combination(const gen & c1,const vecteur & v1,const gen & c2,const vecteur & v2,const gen & c,vecteur & v,double eps,int cstart){
    const_iterateur it1=v1.begin()+cstart,it1end=v1.end(),it2=v2.begin()+cstart;
    iterateur jt1=v.begin()+cstart;
#ifdef DEBUG_SUPPORT
    if (it1end-it1!=v2.end()-it2)
      setdimerr();
#endif
    if (it2==jt1)
      linear_combination(c2,v2,c1,v1,c,v,eps,cstart);
    else {
      if (it1==jt1){
	if (is_one(c)){
	  for (;jt1!=it1end;++jt1,++it2){
	    *jt1=trim(c1*(*jt1)+c2*(*it2),c1,eps);
	  }
	}
	else {
	  for (;jt1!=it1end;++jt1,++it2){
	    *jt1=trim(exact_div(c1*(*jt1)+c2*(*it2),c),c1,eps);
	  }
	}
      }
      else {
	v.clear();
	v.reserve(it1end-it1);
	if (is_one(c)){
	  for (;it1!=it1end;++it1,++it2)
	    v.push_back(trim(c1*(*it1)+c2*(*it2),c1,eps));
	}
	else {
	  for (;it1!=it1end;++it1,++it2)
	    v.push_back(trim(exact_div(c1*(*it1)+c2*(*it2),c),c1,eps));
	}
      }
    }
  }

  // v1=v1+c2*v2 smod modulo
  void modlinear_combination(vecteur & v1,const gen & c2,const vecteur & v2,const gen & modulo,int cstart,int cend){
    if (!is_exactly_zero(c2)){
      iterateur it1=v1.begin()+cstart,it1end=v1.end();
      if (cend && cend>=cstart && cend<it1end-v1.begin())
	it1end=v1.begin()+cend;
      const_iterateur it2=v2.begin()+cstart;
      for (;it1!=it1end;++it1,++it2)
	*it1=smod((*it1)+c2*(*it2),modulo);
    }
  }

  // v1=v1+c2*v2 
  void double_linear_combination(std::vector<giac_double> & v1,double c2,const std::vector<giac_double> & v2,int cstart,int cend){
    if (c2){
      std::vector<giac_double>::iterator it1=v1.begin()+cstart,it1end=v1.end();
      if (cend && cend>=cstart && cend<it1end-v1.begin())
	it1end=v1.begin()+cend;
      std::vector<giac_double>::const_iterator it2=v2.begin()+cstart;
      for (;it1!=it1end;++it1,++it2)
	*it1 += c2*(*it2);
    }
  }

#ifdef _I386_
  // a->a+b*c mod m
  inline void mod(int & a,int b,int c,int m){
    if (c){
      asm volatile("testl %%ebx,%%ebx\n\t" /* sign bit=1 if negative */
		   "jns .Lok%=\n\t"
		   "addl %%edi,%%ebx\n" /* a+=m*/
		   ".Lok%=:\t"
		   "imull %%ecx; \n\t" /* b*c in edx:eax */
		   "addl %%ebx,%%eax; \n\t" /* b*c+a */
		   "adcl $0x0,%%edx; \n\t" /* b*c+a carry */
		   "idivl %%edi; \n\t"
		   :"=d"(a)
		   :"a"(b),"b"(a),"c"(c),"D"(m)
		   );
    }
  }

  // a->a+b*c mod m
  inline int smod(int a,int b,int c,int m){
    if (c){
      if (a<0) a+=m;
      asm volatile("imull %%ecx; \n\t" /* b*c in edx:eax */
		   "addl %%ebx,%%eax; \n\t" /* b*c+a */
		   "adcl $0x0,%%edx; \n\t" /* b*c+a carry */
		   "idivl %%edi; \n\t"
		   :"=d"(a)
		   :"a"(b),"b"(a),"c"(c),"D"(m)
		   );
    }
    return a;
  }
#else
  // a->a+b*c mod m
  inline void mod(int & a,int b,int c,int m){
    a = (a + longlong(b)*c)%m;
  }

  // a->a+b*c mod m
  inline int smod(int a,int b,int c,int m){
    return (a + longlong(b)*c)%m;
  }

#endif

  int dotvecteur(const vecteur & a,const vecteur & b,int modulo){
    vecteur::const_iterator ita=a.begin(), itaend=a.end();
    vecteur::const_iterator itb=b.begin();
    int res=0;
    for (;ita!=itaend;++ita,++itb){
#ifdef _I386_
      mod(res,ita->val,itb->val,modulo);
#else
      res = (res + ita->val*itb->val) % modulo; 
#endif
    }
    return res;
  }

  void multmatvecteur(const matrice & a,const vecteur & b,vecteur & res,int modulo){
    vecteur::const_iterator ita=a.begin(), itaend=a.end();
    res.clear();
    res.reserve(itaend-ita);
    for (;ita!=itaend;++ita)
      res.push_back(dotvecteur(*(ita->_VECTptr),b,modulo));
  }

  // v1=v1+c2*v2 smod modulo
  void modlinear_combination(vector<int> & v1,int c2,const vector<int> & v2,int modulo,int cstart,int cend){
    if (c2){
      vector<int>::iterator it1=v1.begin()+cstart,it1end=v1.end();
      if (cend && cend>=cstart && cend<it1end-v1.begin())
	it1end=v1.begin()+cend;
      vector<int>::const_iterator it2=v2.begin()+cstart;
      for (;it1!=it1end;++it1,++it2)
#ifdef _I386_
	// *it1=( (*it1) + (longlong) c2*(*it2)) % modulo ; // replace smod
	mod(*it1,c2,*it2,modulo);
#else
	*it1=( (*it1) + c2*(*it2)) % modulo ; // replace smod
#endif
    }
  }

  void matrice2std_matrix_gen(const matrice & m,std_matrix<gen> & M){
    int n=m.size();
    M.clear();
    M.reserve(n);
    for (int i=0;i<n;++i)
      M.push_back(*m[i]._VECTptr);
  }

  void std_matrix_gen2matrice(const std_matrix<gen> & M,matrice & m){
    int n=M.size();
    m.clear();
    m.reserve(n);
    for (int i=0;i<n;++i)
      m.push_back(M[i]);
  }

  bool vecteur2index(const vecteur & v,index_t & i){
    i.clear();
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (it->type!=_INT_)
	return false;
      i.push_back(it->val);
    }
    return true;
  }

  static void print_debug_info(const gen & pivot){
    if ( (pivot.type==_POLY) && !pivot._POLYptr->coord.empty())
      cerr << "poly(" << sum_degree(pivot._POLYptr->coord.front().index) << "," << pivot._POLYptr->coord.size() << ") ";
    else
      cerr << pivot << " ";
  }

  bool is_integer_vecteur(const vecteur & m){
    const_iterateur it=m.begin(),itend=m.end();
    for (;it!=itend;++it)
      if (!is_integer(*it)) return false;
    return true;
  }

  bool is_integer_matrice(const matrice & m){
    const_iterateur it=m.begin(),itend=m.end();
    for (;it!=itend;++it)
      if (it->type!=_VECT || !is_integer_vecteur(*it->_VECTptr)) return false;
    return true;
  }

  gen modproduct(const vecteur & v, const gen & modulo){
    const_iterateur it=v.begin(),itend=v.end();
    gen res(1);
    for (;it!=itend;++it){
      res = smod(res * (*it),modulo);
    }
    return res;
  }

  static gen untrunc1(const gen & g){
    return g.type==_POLY?g._POLYptr->untrunc1():g;
  }

  vecteur fracmod(const vecteur & v,const gen & modulo){
    const_iterateur it=v.begin(),itend=v.end();
    vecteur res;
    res.reserve(itend-it);
    for (;it!=itend;++it){
      if (it->type==_VECT)
	res.push_back(fracmod(*it->_VECTptr,modulo));
      else
	res.push_back(fracmod(*it,modulo));
    }
    return res;
  }
  // row reduction from line l and column c to line lmax and column cmax
  // lmax and cmax are not included
  // line are numbered starting from 0
  // if fullreduction is false, reduction occurs under the diagonal only
  // if dont_swap_below !=0, for line numers < dont_swap_below
  // the pivot is searched in the line instead of the column
  // hence no line swap occur
  // convert_internal=false if we do not want conversion to rational fractions
  // algorithm=0 Gauss-Jordan, 1 guess, 2 Bareiss, 3 modular, 4 p-adic, 5 interp
  // rref_or_det_or_lu = 0 for rref, 1 for det, 2 for lu, 
  // 3 for lu without pemutation
  bool mrref(const matrice & a, matrice & res, vecteur & pivots, gen & det,int l, int lmax, int c,int cmax,
	     bool fullreduction,int dont_swap_below,bool convert_internal,int algorithm,int rref_or_det_or_lu,GIAC_CONTEXT){
    if (!ckmatrix(a))
      return false;
    int modular=(algorithm==RREF_MODULAR || algorithm==RREF_PADIC);
    unsigned as=a.size(),a0s=a.front()._VECTptr->size();
    if (algorithm==RREF_GUESS && rref_or_det_or_lu==0 && as>10 && as==a0s-1 && int(as)==lmax && int(a0s)==cmax)
      modular=2;
    if (algorithm==RREF_GUESS && rref_or_det_or_lu<0){
      modular=1;
      rref_or_det_or_lu=-rref_or_det_or_lu;
    }
    if (rref_or_det_or_lu==2 || rref_or_det_or_lu == 3){ // LU decomposition
      algorithm=RREF_GAUSS_JORDAN;
      dont_swap_below=0;
      convert_internal=false;
      fullreduction=false;
    }
    vector<int> permutation(lmax);
    for (int i=0;i<lmax;++i)
      permutation[i]=i;
    // The following works also for non-det operations but it seems
    // to be slower than Bareiss
    if ( ( (algorithm==RREF_GUESS &&rref_or_det_or_lu==1) || modular ) && is_integer_matrice(a) && as<=a0s){
      // Modular algorithm for matrix integer reduction
      // Find Hadamard bound
      gen h2=4*square_hadamard_bound(a);
#ifdef _I386_
      gen p(536870923),det_mod_p;
#else
      gen p(36007),det_mod_p;
#endif
      gen pi_p;
      int done=0;
      bool failure=false;
      gen factdet(1); // find a divisor of the determinant
      // by solving a random linear system having a as matrix
      // using a p-adic method 
#ifdef _I386_
      double p0=3037000500./std::sqrt(double(as))/5.; // so that p0^2*rows(a)<2^63
#else
      double p0=46340./std::sqrt(double(as))/5.; // so that p0^2*rows(a)<2^31
#endif
      if (modular==2){ // rref is like linsolve
	p=nextprime(int(p0));
	matrice A(mtran(a));
	vecteur b=*A.back()._VECTptr,x;
	A.pop_back();
	A=mtran(A);
	int done=padic_linsolve(A,b,x,p,det);
	if (done>0){
	  res=midn(as);
	  res.push_back(x);
	  res=mtran(res);
	  return true;
	}
	failure=true;
      }
      if (as>50 && as==a0s){
	p=nextprime(int(p0));
	vecteur b(vranm(as,8,contextptr)),resb;
	// reconstruct at most 5 components of res for lcm
	if ( (done=padic_linsolve(a,b,resb,p,det,5)) ){ 
	  if (done==-1){
	    det=0;
	    return true;
	  }
	  lcmdeno(resb,factdet,contextptr);
	  cerr << "lif=" << factdet << endl;
	  h2=iquo(h2,factdet*factdet)+1;
	  det=smod(det*invmod(factdet,p),p);
	  pi_p=p;
	}
#ifdef _I386_
	p=536870923;
#endif
      }
      if (!failure){
	double proba=1.0;
	if (!done){
	  pi_p=p;
	  if (!modrref(a,res,pivots,det,l,lmax,c,cmax,
		       false,dont_swap_below,p,1 /* det */))
	    return false;
	}
	// First find det to avoid bad primes
	for (;is_strictly_greater(h2,pi_p*pi_p,contextptr);){
	  for (;;){
	    p=nextprime(p+1);
	    if (is_one(gcd(factdet,p))) // keep p prime with factdet
	      break;
	  }
	  if (as>10 && debug_infolevel)
	    cerr << "Modrref, % done " << evalf_double(_evalf(gen(makevecteur(200*ln(pi_p,contextptr)/ln(h2,contextptr),20),_SEQ__VECT),contextptr),1,contextptr)<< ", prime " << p << ", det/lif=" << det << endl;
	  if (!modrref(a,res,pivots,det_mod_p,l,lmax,c,cmax,
		       false,dont_swap_below,p,1 /* det */))
	    return false;
	  det_mod_p=smod(det_mod_p*invmod(factdet,p),p);
	  gen old_det=det;
	  det=ichinrem(det,det_mod_p,pi_p,p);
	  if (old_det==det)
	    proba=proba/evalf_double(p,1,contextptr)._DOUBLE_val;
	  else
	    proba=1.0;
	  pi_p=pi_p*p;
	  if (proba<proba_epsilon(contextptr))
	    break;
	}
	det=smod(det,pi_p)*factdet;
	if (rref_or_det_or_lu==1)
	  return true;
	if (is_zero(det,contextptr))
	  failure=true;
      }
      if (!failure){
	// Improve: currently permutation should always be the idn for lu
	// instead of by det (det works for rref)
	if (rref_or_det_or_lu==2){
	  rref_or_det_or_lu=3;
	  h2=h2*h2; // need to square for LU decomp (rational reconstruction)
	}
	// Now do the reduction again, avoiding bad primes
	p=36007;
	gen q;
	while (is_zero(irem(det,p,q),contextptr))
	  p=nextprime(p+1);
	pi_p=p;
	gen det1;
	if (!modrref(a,res,pivots,det1,l,lmax,c,cmax,
		     fullreduction,dont_swap_below,p,rref_or_det_or_lu))
	  return false;
	// Multiply res by product of pivots in order to have the det
	// as initial non-zero element of each line after the reduction
	if (rref_or_det_or_lu==0)
	  res=smod(multvecteur(det,res),p);
	for (;is_strictly_greater(h2,pi_p*pi_p,contextptr);){
	  p=nextprime(p+1);
	  while (is_zero(irem(det,p,q),contextptr))
	    p=nextprime(p+1);
	  matrice res_mod_p,pivots_mod_p;
	  if (!modrref(a,res_mod_p,pivots_mod_p,det_mod_p,l,lmax,c,cmax,
		       fullreduction,dont_swap_below,p,rref_or_det_or_lu))
	    return false;
	  if (rref_or_det_or_lu==3){
	    if (is_zero(det_mod_p,contextptr))
	      continue;
	  }
	  else
	    res_mod_p=smod(multvecteur(det,res_mod_p),p);
	  res=*ichinrem(gen(res),gen(res_mod_p),pi_p,p)._VECTptr;
	  pivots=*ichinrem(gen(pivots),gen(pivots_mod_p),pi_p,p)._VECTptr;
	  pi_p=pi_p*p;
	}
	res=smod(res,pi_p);
	pivots=smod(pivots,pi_p);
	if (rref_or_det_or_lu==3) // rational reconstruction
	  res=fracmod(res,pi_p);
	if (rref_or_det_or_lu==2 || rref_or_det_or_lu == 3){
	  vecteur P;
	  vector_int2vecteur(permutation,P);
	  pivots.push_back(P);
	}
	return true;
      } // end if !failure
    } // end modular/padic algorithm
    gen tmp=a.front();
    if (tmp.type==_VECT && !tmp._VECTptr->empty()){
      tmp=tmp._VECTptr->front();
      if (tmp.type==_MOD){
	gen modulo=*(tmp._MODptr+1);
	if (!modrref(*unmod(a)._VECTptr,res,pivots,det,l,lmax,c,cmax,
		     fullreduction,dont_swap_below,modulo,rref_or_det_or_lu))
	  return false;
	res=*makemod(res,modulo)._VECTptr;
	pivots=*makemod(pivots,modulo)._VECTptr;
	det=makemod(det,modulo);
	return true;
      }
    }
    int linit=l;//,previous_l=l;
    vecteur lv;
    bool num_mat=has_num_coeff(a);
    if (num_mat){
      res=*evalf_VECT(a,0,1,contextptr)._VECTptr;
      if (algorithm==RREF_GUESS)
	algorithm=RREF_LAGRANGE;
      num_mat=is_fully_numeric(res);
#if 1 // ndef BCD
      matrix_double N;
      if (num_mat && matrice2std_matrix_double(res,N)){
	// specialization for double
	double ddet;
	vector<int> maxrankcols;
	doublerref(N,pivots,permutation,maxrankcols,ddet,l,lmax,c,cmax,fullreduction,dont_swap_below,rref_or_det_or_lu);
	if (rref_or_det_or_lu!=1){
	  std_matrix<gen> RES;
	  std_matrix_giac_double2std_matrix_gen(N,RES);
	  std_matrix_gen2matrice(RES,res);
	}
	det=ddet;
	if (rref_or_det_or_lu==2 || rref_or_det_or_lu == 3){
	  vecteur P;
	  vector_int2vecteur(permutation,P);
	  pivots.push_back(P);
	}
	return true;
      }
#endif
    }
    else
      res=a;
    if (convert_internal){
      // convert a to internal form
      lv=alg_lvar(res);
      if (!lv.empty() && lv.front().type==_VECT && lv.front()._VECTptr->size()>1){
	vecteur lw=*tsimplify(lv.front(),contextptr)._VECTptr;
	if (lvar(lw).size()<lw.front()._VECTptr->size()){
	  res=*subst(gen(res),lv.front(),lw,false,contextptr)._VECTptr;
	  lv=alg_lvar(res);
	}
      }
      res = *(e2r(res,lv,contextptr)._VECTptr);
    }
    int lvs=lv.size();
    // cout << res << endl;
    gen lcm_deno,gcd_num;
    gen detnum = plus_one;
    gen detden = plus_one;
    if (algorithm!=RREF_GAUSS_JORDAN){
      // remove common denominator of each line (fraction-free elim)
      iterateur it=res.begin(),itend=res.end();
      for (;it!=itend;++it){
	if (num_mat){ // divide each line by max coeff in abs value
	  lcm_deno=linfnorm(*it,contextptr);
	  detnum=lcm_deno*detnum;
	  multvecteur(inv(lcm_deno,contextptr),*it->_VECTptr,*it->_VECTptr);
	}
	else { // non num mat
	  lcm_deno=common_deno(*it->_VECTptr);
	  iterateur jt=it->_VECTptr->begin(),jtend=it->_VECTptr->end();
	  for (;jt!=jtend;++jt){
	    if (jt->type==_FRAC){
	      gen nm(jt->_FRACptr->num);
	      gen dn(jt->_FRACptr->den);
	      // *jt -> lcmdeno* (nm/dn) = nm * tmp/dn
	      gen tmp(lcm_deno);
	      simplify(tmp,dn);
	      if (dn.type<=_CPLX){
		*jt=nm*tmp/dn;
		continue;
	      }
	      if (dn.type==_POLY){
		*jt=nm*tmp/dn._POLYptr->coord.front().value;
		continue;
	      }
	      return false; // settypeerr();
	    }
	    else
	      *jt=(*jt) * lcm_deno;
	  }
	  detden = detden * lcm_deno;
	  gcd_num=common_num(*it->_VECTptr);
	  if (!is_zero(gcd_num,contextptr))
	    *it=rdiv(*it,gcd_num);
	  detnum=detnum*gcd_num;
	} // end else (non num mat)
      } // end for (;it!=itend;)
      // check if res is integer or polynomial
      if (lvs==1 && lv.front().type==_VECT && lv.front()._VECTptr->empty() && (rref_or_det_or_lu==1 || modular ) && is_integer_matrice(res) && as<=a0s){
	matrice res1;
	if (!mrref(res,res1,pivots,det,l,lmax,c,cmax,fullreduction,dont_swap_below,false,algorithm,rref_or_det_or_lu,contextptr))
	  return false;
	res=res1;
	det=detnum*det/detden;
	if (rref_or_det_or_lu==2 || rref_or_det_or_lu == 3){
	  vecteur P;
	  vector_int2vecteur(permutation,P);
	  pivots.push_back(P);
	}
	return true;
      }
      if (rref_or_det_or_lu==1 && as==a0s && as>4 && algorithm==RREF_GUESS && convert_internal && lvs==1 && lv.front().type==_VECT){
	// guess if Bareiss or Lagrange interpolation is faster
	// Bareiss depends on the total degree, Lagrange on partial degrees
	// gather line/columns statistics
	int polydim=lv.front()._VECTptr->size();
	index_t col_totaldeg(as);
	vector< index_t > col_partialdeg(as,index_t(polydim));
	int maxtotaldeg=0,summaxtotaldeg=0;
	if (polydim){
	  index_t summaxdeg(polydim);
	  for (unsigned int i=0;i<as;++i){
	    index_t maxdeg(polydim);
	    for (unsigned int j=0;j<as;++j){
	      const gen & tmp = (*res[i]._VECTptr)[j];
	      if (tmp.type==_POLY){
		const index_t & degij=tmp._POLYptr->degree();
		maxdeg=index_lcm(degij,maxdeg);
		col_partialdeg[j]=index_lcm(degij,col_partialdeg[j]);
		int totaldeg=tmp._POLYptr->total_degree();
		if (maxtotaldeg<totaldeg)
		  maxtotaldeg=totaldeg;
		if (col_totaldeg[j]<totaldeg)
		  col_totaldeg[j]=totaldeg;
	      }
	    }
	    summaxtotaldeg += maxtotaldeg;
	    summaxdeg=summaxdeg+maxdeg;
	  }
	  maxtotaldeg=std::min(summaxtotaldeg,int(total_degree(col_totaldeg)));
	  index_t col_sumpartialdeg(polydim);
	  for (unsigned int j=0;j<as;++j)
	    col_sumpartialdeg = col_sumpartialdeg+col_partialdeg[j];
	  for (int i=0;i<polydim;++i){
	    summaxdeg[i]=std::min(summaxdeg[i],col_sumpartialdeg[i]);
	  }
	  if (debug_infolevel)
	    cerr << "Total degree " << maxtotaldeg << ", partial degrees " << summaxdeg << endl;
	  // Now modify algorithm to RREF_LAGRANGE if it's faster
	  double lagrange_time=std::pow(double(as),2)*(as*10+160); 
	  // coeffs of as*.+. are guess
	  for (int j=0;j<polydim;j++){
	    lagrange_time *= (summaxdeg[j]+1);
	  }
	  double bareiss_time=0;
	  // time is almost proportionnal to sum( comb(maxtotaldeg*j/as+polydim,polydim)^2, j=1..as-1)
	  for (unsigned int j=1;j<as;++j){
	    int tmpdeg=int(double(maxtotaldeg*j)/as+.5);
	    double tmp = evalf_double(comb(tmpdeg+polydim,polydim),1,contextptr)._DOUBLE_val;
	    tmp = tmp*tmp*std::log(tmp)*(as-j);
	    bareiss_time += tmp;
	  }
	  bareiss_time *= as; // take account of the size of the coefficients
	  if (debug_infolevel)
	    cerr << "lagrange " << lagrange_time << " bareiss " << bareiss_time << endl;
	  if (lagrange_time<bareiss_time){
	    algorithm=RREF_LAGRANGE;
	  }
	} // end if (polydim)
      }
      if ( algorithm==RREF_LAGRANGE && rref_or_det_or_lu==1 && as==a0s){
	vecteur lva=lvar(a);
	if ( (!convert_internal && lva.empty()) || (lvs==1 && lv.front()==lva) ){
	  // find degrees wrt main variable
	  int polydim=0;
	  int totaldeg=0;
	  vector<int> maxdegj(as);
	  for (unsigned int i=0;i<as;++i){
	    int maxdegi=0;
	    for (unsigned int j=0;j<a0s;++j){
	      gen & tmp = (*res[i]._VECTptr)[j];
	      if (tmp.type==_POLY){
		polydim=tmp._POLYptr->dim;
		const int & curdeg=tmp._POLYptr->lexsorted_degree();
		if (curdeg>maxdegi)
		  maxdegi=tmp._POLYptr->lexsorted_degree();
		if (curdeg>maxdegj[j])
		  maxdegj[j]=curdeg;
		tmp=polynome2poly1(tmp,1);
	      }
	    }
	    totaldeg+=maxdegi;
	  }
	  if (polydim){
	    totaldeg=std::min(totaldeg,total_degree(maxdegj));
	    proba_epsilon(contextptr) /= totaldeg;
	    vecteur X(totaldeg+1),Y(totaldeg+1);
	    for (int x=0;x<=totaldeg;++x){
	      X[x]=x;
	      vecteur resx;
	      resx.reserve(totaldeg+1);
	      for (unsigned int i=0;i<as;++i){
		vecteur resxi;
		resxi.reserve(totaldeg+1);
		for (unsigned int j=0;j<a0s;++j){
		  const gen & tmp = (*res[i]._VECTptr)[j];
		  resxi.push_back(horner(tmp,x));
		}
		resx.push_back(resxi);
	      }
	      matrice res1;
	      if (!mrref(resx,res1,pivots,det,l,lmax,c,cmax,fullreduction,dont_swap_below,false,algorithm,1,contextptr))
		return false;
	      Y[x]=det;
	    } // end for x
	    proba_epsilon(contextptr) *= totaldeg;
	    // Lagrange interpolation
	    vecteur L=divided_differences(X,Y);
	    det=untrunc1(L[totaldeg]);
	    monomial<gen> mtmp(1,1,polydim);
	    gen xpoly=polynome(mtmp);
	    for (int i=totaldeg-1;i>=0;--i){
	      det = det*(xpoly-untrunc1(X[i]))+untrunc1(L[i]);
	    }
	    det=det*detnum/detden;
	    if (convert_internal)
	      det=r2sym(det,lva,contextptr);
	    if (rref_or_det_or_lu==2 || rref_or_det_or_lu == 3){
	      vecteur P;
	      vector_int2vecteur(permutation,P);
	      pivots.push_back(P);
	    }
	    return true;
	  } // end if polydim
	}
      }
    }

    std_matrix<gen> M;
    matrice2std_matrix_gen(res,M);
    gen bareiss (1),pivot,temp;
    // vecteur vtemp;
    int pivotline,pivotcol;
    pivots.clear();
    pivots.reserve(cmax-c);
    for (;(l<lmax) && (c<cmax);){
      if ( (!fullreduction) && (l==lmax-1) )
	break;
      if (debug_infolevel)
	cerr <<  "// mrref line " << l << ":" << clock() <<endl;
      pivot=M[l][c];
      if (debug_infolevel ){
	cerr << "// ";
	print_debug_info(pivot);
      }
      pivotline=l;
      pivotcol=c;
      if (l<dont_swap_below){ // scan current line for the best pivot available
	for (int ctemp=c+1;ctemp<cmax;++ctemp){
	  temp=M[l][ctemp];
	  if (debug_infolevel)
	    print_debug_info(temp);
	  if (!is_zero(temp,contextptr) && temp.islesscomplexthan(pivot)){
	    pivot=temp;
	    pivotcol=ctemp;
	  }
	}	
      }
      else {      // scan M current column for the best pivot available
	if (rref_or_det_or_lu == 3){ // LU without line permutation
	  if (is_zero(pivot,contextptr)){
	    det = 0;
	    vecteur P;
	    vector_int2vecteur(permutation,P);
	    pivots.push_back(P);
	    return true;
	  }
	}
	else {
	  for (int ltemp=l+1;ltemp<lmax;++ltemp){
	    temp=M[ltemp][c];
	    if (debug_infolevel)
	      print_debug_info(temp);
	    if (!is_zero(temp,contextptr) && temp.islesscomplexthan(pivot)){
	      pivot=temp;
	      pivotline=ltemp;
	    }
	  }
	}
      }
      if (debug_infolevel)
	cerr << endl;
      //cout << M << endl << pivot << endl;
      if (!is_zero(pivot,contextptr)){
	// exchange lines if needed
	if (l!=pivotline){
	  swap(M[l],M[pivotline]);
	  swap(permutation[l],permutation[pivotline]);
	  // temp = M[l];
	  // M[l] = M[pivotline];
	  // M[pivotline] = temp;
	  detnum = -detnum;
	}
	// make the reduction
	if (fullreduction){
	  for (int ltemp=linit;ltemp<lmax;++ltemp){
	    if (debug_infolevel>=2)
	      cerr << "// " << l << "," << ltemp << " "<< endl;
	    if (ltemp!=l){
	      if (algorithm!=RREF_GAUSS_JORDAN) // M[ltemp] = rdiv( pivot * M[ltemp] - M[ltemp][pivotcol]* M[l], bareiss);
		linear_combination(pivot,M[ltemp],-M[ltemp][pivotcol],M[l],bareiss,M[ltemp],1e-12,0);
	      else // M[ltemp]=M[ltemp]-rdiv(M[ltemp][pivotcol],pivot)*M[l];
		linear_combination(plus_one,M[ltemp],-rdiv(M[ltemp][pivotcol],pivot),M[l],plus_one,M[ltemp],1e-12,0);
	    }
	  }
	}
	else { // subdiagonal reduction
	  for (int ltemp=l+1;ltemp<lmax;++ltemp){
	    if (debug_infolevel>=2)
	      cerr << "// " << l << "," << ltemp << " "<< endl;
	    if (algorithm!=RREF_GAUSS_JORDAN)
	      linear_combination(pivot,M[ltemp],-M[ltemp][pivotcol],M[l],bareiss,M[ltemp],1e-12,(c+1)*(rref_or_det_or_lu>0));
	    else {
	      gen coeff=M[ltemp][pivotcol]/pivot;
	      linear_combination(plus_one,M[ltemp],-coeff,M[l],plus_one,M[ltemp],1e-12,(c+1)*(rref_or_det_or_lu>0));
	      if (rref_or_det_or_lu==2 || rref_or_det_or_lu == 3)
		M[ltemp][pivotcol]=coeff;
	    }
	  }
	  if (rref_or_det_or_lu==1 && algorithm!=RREF_GAUSS_JORDAN) {
	    if (debug_infolevel)
	      cerr << "//mrref clear line " << l << endl;
	    // clear pivot line to save memory
	    M[l].clear();
	  }
	} // end else
	// cout << M << endl;
	// increment column number if swap was allowed
	if (l>=dont_swap_below)
	  ++c;
	// increment line number since reduction has been done
	++l;	  
	// multiply det
	// set new bareiss for next reduction round
	if (algorithm!=RREF_GAUSS_JORDAN)
	  bareiss=pivot;
	// save pivot for annulation test purposes
	if (rref_or_det_or_lu!=1){
	  if (convert_internal)
	    pivots.push_back(r2sym(pivot,lv,contextptr));
	  else
	    pivots.push_back(pivot);
	  if (debug_infolevel)
	    cerr << pivots.back() << endl;
	}
      }
      else { // if pivot is 0 increment either the line or the col
	if (rref_or_det_or_lu==1){
	  det=0;
	  return true;
	}
	if (l>=dont_swap_below)
	  c++;
	else
	  l++;
      }
    } // end for reduction loop
    if (debug_infolevel)
      cerr << "// mrref reduction end:" << clock() << endl;
    if (algorithm!=RREF_GAUSS_JORDAN){
      int last=min(lmax,cmax);
      det=M[last-1][last-1];
      if ( (debug_infolevel) && (det.type==_POLY) )
	cerr << "// polynomial size " << det._POLYptr->coord.size() << endl;
      if (rref_or_det_or_lu==1) // erase last line of the matrix
	M[lmax-1].clear();
      det=rdiv(det*detnum,detden);
      if (convert_internal)
	det=r2sym(det,lv,contextptr);
      // cerr << det << endl;
    }
    else {
      // adjust determinant by multiplication by all diagonal coeffs
      for (int i=linit;i<lmax;++i)
	detnum = detnum * M[i][i];
      det = rdiv(detnum,detden);
      if (convert_internal)
	det = r2sym(det,lv,contextptr);
    }
    std_matrix_gen2matrice(M,res);
    if (convert_internal)
      res = *(r2sym (res,lv,contextptr)._VECTptr);
    if (rref_or_det_or_lu==2 || rref_or_det_or_lu == 3){
      vecteur P;
      vector_int2vecteur(permutation,P);
      pivots.push_back(P);
    }
    if (debug_infolevel)
      cerr << "// mrref end:" << clock() << " " << M << endl;
    return true;
  }

  // convert a to vector< vector<int> > with modular reduction (if modulo!=0)
  void vect_vecteur_2_vect_vector_int(const std_matrix<gen> & M,int modulo,vector< vector<int> > & N){
    int Msize=M.size();
    N.clear();
    N.reserve(Msize);
    for (int k=0;k<Msize;k++){
      const vecteur & v = M[k];
      const_iterateur it=v.begin(),itend=v.end();
      vector<int> vi(itend-it);
      vector<int>::iterator jt=vi.begin();
      for (;it!=itend;++jt,++it){
	if (!modulo)
	  *jt=it->val;
	else
	  *jt=smod(*it,modulo).val;
      }
      N.push_back(vi);
    }
  }

  void vect_vector_int_2_vect_vecteur(const vector< vector<int> > & N,std_matrix<gen> & M){
    // Back convert N to M
    int Msize=N.size();
    M = std_matrix<gen>(Msize);
    for (int k=0;k<Msize;k++){
      const vector<int> & v = N[k];
      vector<int>::const_iterator it=v.begin(),itend=v.end();
      vecteur vi(itend-it);
      iterateur jt=vi.begin();
      for (;it!=itend;++jt,++it){
	*jt=*it;
      }
      M[k]=vi;
    }
  }

  //transforme un vecteur en vector<int>  
  void vecteur2vector_int(const vecteur & v,int modulo,vector<int> & res){
    vecteur::const_iterator it=v.begin(),itend=v.end();
    res.clear();
    res.reserve(itend-it);
    for (;it!=itend;++it){
      if (it->type==_MOD)
	res.push_back(it->_MODptr->val);
      else {
	if (modulo) 
	  res.push_back(smod((*it),modulo).val); 
	else
	  res.push_back(it->val); 
      }
    }
  } 

  bool vecteur2vectvector_int(const vecteur & v,int modulo,vector< vector<int> > & res){
    vecteur::const_iterator it=v.begin(),itend=v.end();
    res.clear();
    res.reserve(itend-it);
    vector<int> w;
    for (;it!=itend;++it){
      if (it->type!=_VECT)
	return false;
      vecteur2vector_int(*it->_VECTptr,modulo,w);
      res.push_back(w);
    }
    return true;
  }

  void vector_int2vecteur(const vector<int> & v,vecteur & res){
    //transforme un vector<int> en vecteur 
    vector<int>::const_iterator it=v.begin(),itend=v.end();
    res.clear();
    res.reserve(itend-it);
    for (;it!=itend;++it)
      res.push_back(*it);
  } 

  void vectvector_int2vecteur(const vector< vector<int> > & v,vecteur & res){
    //transforme un vector< vector<int> > en vecteur  
    int s=v.size();
    res.clear();
    res.reserve(s);
    vecteur w;
    for (int i=0;i<s;++i){
      vector_int2vecteur(v[i],w);
      res.push_back(w);
    }
  }

  int dotvector_int(const vector<int> & v,const vector<int> & w,int modulo){
    vector<int>::const_iterator it=v.begin(),itend=v.end(),it1,jt=w.begin();
    unsigned n=itend-it;
    if ( ((longlong(modulo)*modulo)/RAND_MAX)*n>RAND_MAX){
      int res=0;
      for (;it!=itend;++jt,++it){
#ifdef _I386_
	mod(res,*it,*jt,modulo);
#else
	res = (res + longlong(*it)*(*jt))% modulo;
#endif
      }
      return smod(res,modulo) ;
    }
    longlong res=0;
    it1 = it + ((n>>2) <<2);
    for (;it!=it1;){
      res += (longlong (*it))*(*jt);
      ++jt; ++it;
      res += (longlong (*it))*(*jt);
      ++jt; ++it;
      res += (longlong (*it))*(*jt);
      ++jt; ++it;
      res += (longlong (*it))*(*jt);
      ++jt; ++it;
    }
    for (;it!=itend;++jt,++it){
      res += (longlong (*it))*(*jt);
    }
    return smod(res,modulo) ;
  }

  bool multvectvector_int_vector_int(const vector< vector<int> > & M,const vector<int> & v,int modulo,vector<int> & Mv){
    unsigned n=M.size();
    Mv.clear();
    if (!n)
      return true;
    if (M.front().size()!=v.size())
      return false; 
    Mv.reserve(n);
    vector< vector<int> >::const_iterator it=M.begin(),itend=M.end();
    for (;it!=itend;++it){
      Mv.push_back(dotvector_int(*it,v,modulo));
    }
    return true;
  }

  void tran_vect_vector_int(const vector< vector<int> > & N,vector< vector<int> > & tN){
    tN.clear();
    unsigned r=N.size();
    if (!r)
      return;
    unsigned c=N.front().size();
    tN.reserve(c);
    for (unsigned int i=0;i<c;++i){
      vector<int> current;
      current.reserve(r);
      for (unsigned int j=0;j<r;++j){
	current.push_back(N[j][i]);
      }
      tN.push_back(current);
    }
  }

  void apply_permutation(const vector<int> & permutation,const vector<int> &x,vector<int> & y){
    unsigned n=x.size();
    y.clear();
    y.reserve(n);
    for (unsigned int i=0;i<n;++i)
      y.push_back(x[permutation[i]]);
  }

  /*
  vector<int> perminv(const vector<int> & p);
  // solve LU x= b (permutation P)
  void smallsolvelu(const vector< vector<int> > & LU,const vector<int> & P,const vector<int> & b,vector<int> & x,int modulo){
    unsigned n=P.size();
    vector<int> bp(n),y(n);
    apply_permutation(P,b,bp);
    // solve U y=bp
    for (int i=n-1;i>=0;--i){
      // y[i]=LU[i,i]^(-1)*(bp[i]-sum(j>i)LU[i,j]*y[j])
      int res=0;
      const vector<int> & li=LU[i];
      for (int j=i+1;j<n;++j)
	mod(res,li[j],y[j],modulo);
      y[i]=(invmod(li[i],modulo)*longlong(bp[i]-res))%modulo;
    }
    // solve L bp = y
    for (int i=0;i<n;++i){
      // bp[i]=(y[i]-sum(j<i)LU[i,j]*bp[j])
      int res=0;
      const vector<int> & li=LU[i];
      for (int j=0;j<i;++j)
	mod(res,li[j],y[j],modulo);
      y[i]=longlong(bp[i]-res)%modulo;
    }
    // reorder bp
    apply_permutation(perminv(P),bp,x);
  }
  */

  // if dont_swap_below !=0, for line numers < dont_swap_below
  // the pivot is searched in the line instead of the column
  // hence no line swap occur
  // rref_or_det_or_lu = 0 for rref, 1 for det, 2 for lu, 
  // 3 for lu without permutation
  // fullreduction=0 or 1, use 2 if the right part of a is idn
  void smallmodrref(vector< vector<int> > & N,vecteur & pivots,vector<int> & permutation,vector<int> & maxrankcols,longlong & idet,int l, int lmax, int c,int cmax,int fullreduction,int dont_swap_below,int modulo,int rref_or_det_or_lu){
    bool use_cstart=!c;
    bool inverting=fullreduction==2;
    int linit=l;//,previous_l=l;
    // Reduction
    int pivot,temp;
    idet=1;
    // vecteur vtemp;
    int pivotline,pivotcol;
    pivots.clear();
    pivots.reserve(cmax-c);
    permutation.clear();
    maxrankcols.clear();
    bool noswap=true;
    for (int i=0;i<lmax;++i)
      permutation.push_back(i);
    for (;(l<lmax) && (c<cmax);){
      pivot=N[l][c];
      if (rref_or_det_or_lu==3 && !pivot){
	idet=0;
	return;
      }
      if ( rref_or_det_or_lu==1 && l==lmax-1 ){
	idet = (idet * pivot) % modulo ;
	break;
      }
      pivotline=l;
      pivotcol=c;
      if (!pivot){ // scan current line
	noswap=false;
	if (l<dont_swap_below){ 
	  for (int ctemp=c+1;ctemp<cmax;++ctemp){
	    temp=N[l][ctemp];
	    if (temp){
	      pivot=smod(temp,modulo);
	      pivotcol=ctemp;
	      break;
	    }
	  }
	}
	else {      // scan N current column for the best pivot available
	  for (int ltemp=l+1;ltemp<lmax;++ltemp){
	    temp=N[ltemp][c];
	    if (debug_infolevel)
	      print_debug_info(temp);
	    if (temp){
	      pivot=smod(temp,modulo);
	      pivotline=ltemp;
	      break;
	    }
	  }
	}
      } // end if is_zero(pivot), true pivot found on line or column
      if (pivot){
	maxrankcols.push_back(c);
	if (l!=pivotline){
	  swap(N[l],N[pivotline]);
	  swap(permutation[l],permutation[pivotline]);
	  pivotline=l;
	  idet = -idet;
	}
	// save pivot for annulation test purposes
	if (rref_or_det_or_lu!=1)
	  pivots.push_back(pivot);
	// invert pivot 
	temp=invmod(pivot,modulo);
	// multiply det
	idet = (idet * pivot) % modulo ;
	if (fullreduction || rref_or_det_or_lu<2){ // not LU decomp
	  vector<int>::iterator it=N[pivotline].begin(),itend=N[pivotline].end();
	  for (;it!=itend;++it){
	    *it=(longlong(temp) * *it)%modulo;
	    if (*it>0){
	      if (2* *it>modulo)
		*it -= modulo;
	    }
	    else {
	      if (2* *it<-modulo)
		*it += modulo;
	    }
	  }
	}
	// make the reduction
	if (fullreduction){
	  for (int ltemp=linit;ltemp<lmax;++ltemp){
	    if (ltemp!=l)
	      modlinear_combination(N[ltemp],-N[ltemp][pivotcol],N[l],modulo,(use_cstart?c:0),(inverting && noswap)?(c+1+lmax):0);
	  }
	}
	else {
	  for (int ltemp=l+1;ltemp<lmax;++ltemp){
	    if (rref_or_det_or_lu>=2) // LU decomp
	      N[ltemp][pivotcol]= (N[ltemp][pivotcol]*longlong(temp)) % modulo;
	    modlinear_combination(N[ltemp],-N[ltemp][pivotcol],N[l],modulo,(rref_or_det_or_lu>0)?(c+1):(use_cstart?c:0));
	  }
	} // end else
	  // increment column number if swap was allowed
	if (l>=dont_swap_below)
	  ++c;
	// increment line number since reduction has been done
	++l;	  
      } // end if (!is_zero(pivot)
      else { // if pivot is 0 increment either the line or the col
	idet = 0;
	if (rref_or_det_or_lu==1)
	  return;
	if (l>=dont_swap_below)
	  c++;
	else
	  l++;
      }
    } // end for reduction loop
    if (rref_or_det_or_lu!=1){
      for (int i=0;i<lmax;i++){
	for (int j=0;j<cmax;j++){
	  N[i][j]=smod(N[i][j],modulo);
	}
      }
    }
  }

  // if dont_swap_below !=0, for line numers < dont_swap_below
  // the pivot is searched in the line instead of the column
  // hence no line swap occur
  // rref_or_det_or_lu = 0 for rref, 1 for det, 2 for lu, 
  // 3 for lu without permutation
  // fullreduction=0 or 1, use 2 if the right part of a is idn
  void doublerref(matrix_double & N,vecteur & pivots,vector<int> & permutation,vector<int> & maxrankcols,double & idet,int l, int lmax, int c,int cmax,int fullreduction,int dont_swap_below,int rref_or_det_or_lu){
    bool use_cstart=!c;
    bool inverting=fullreduction==2;
    int linit=l;//,previous_l=l;
    // Reduction
    double pivot,temp;
    idet=1;
    // vecteur vtemp;
    int pivotline,pivotcol;
    pivots.clear();
    pivots.reserve(cmax-c);
    permutation.clear();
    maxrankcols.clear();
    bool noswap=true;
    for (int i=0;i<lmax;++i)
      permutation.push_back(i);
    for (;(l<lmax) && (c<cmax);){
      pivot=N[l][c];
      if (rref_or_det_or_lu==3 && !pivot){
	idet=0;
	return;
      }
      if ( rref_or_det_or_lu==1 && l==lmax-1 ){
	idet = (idet * pivot);
	break;
      }
      pivotline=l;
      pivotcol=c;
      noswap=false;
      if (l<dont_swap_below){ 
	for (int ctemp=c+1;ctemp<cmax;++ctemp){
	  temp=N[l][ctemp];
	  if (std::abs(temp)>std::abs(pivot)){
	    pivot=temp;
	    pivotcol=ctemp;
	  }
	}
      }
      else {      // scan N current column for the best pivot available
	for (int ltemp=l+1;ltemp<lmax;++ltemp){
	  temp=N[ltemp][c];
	  if (debug_infolevel)
	    print_debug_info(temp);
	  if (std::abs(temp)>std::abs(pivot)){
	    pivot=temp;
	    pivotline=ltemp;
	  }
	}
      }
      if (pivot){
	maxrankcols.push_back(c);
	if (l!=pivotline){
	  swap(N[l],N[pivotline]);
	  swap(permutation[l],permutation[pivotline]);
	  pivotline=l;
	  idet = -idet;
	}
	// save pivot for annulation test purposes
	if (rref_or_det_or_lu!=1)
	  pivots.push_back(pivot);
	// invert pivot 
	temp=1./pivot;
	// multiply det
	idet = idet * pivot ;
	if (fullreduction || rref_or_det_or_lu<2){ // not LU decomp
	  std::vector<giac_double>::iterator it=N[pivotline].begin(),itend=N[pivotline].end();
	  for (;it!=itend;++it){
	    *it /= pivot;
	  }
	}
	// make the reduction
	if (fullreduction){
	  for (int ltemp=linit;ltemp<lmax;++ltemp){
	    if (ltemp!=l)
	      double_linear_combination(N[ltemp],-N[ltemp][pivotcol],N[l],(use_cstart?c:0),(inverting && noswap)?(c+1+lmax):0);
	  }
	}
	else {
	  for (int ltemp=l+1;ltemp<lmax;++ltemp){
	    if (rref_or_det_or_lu>=2) // LU decomp
	      N[ltemp][pivotcol]= N[ltemp][pivotcol]*temp;
	    double_linear_combination(N[ltemp],-N[ltemp][pivotcol],N[l],(rref_or_det_or_lu>0)?(c+1):(use_cstart?c:0),0);
	  }
	} // end else
	  // increment column number if swap was allowed
	if (l>=dont_swap_below)
	  ++c;
	// increment line number since reduction has been done
	++l;	  
      } // end if (!is_zero(pivot)
      else { // if pivot is 0 increment either the line or the col
	idet = 0;
	if (rref_or_det_or_lu==1)
	  return;
	if (l>=dont_swap_below)
	  c++;
	else
	  l++;
      }
    } // end for reduction loop
  }

  // if dont_swap_below !=0, for line numers < dont_swap_below
  // the pivot is searched in the line instead of the column
  // hence no line swap occur
  // rref_or_det_or_lu = 0 for rref, 1 for det, 2 for lu, 
  // 3 for lu without permutation
  // fullreduction=0 or 1, use 2 if the right part of a is idn
  bool modrref(const matrice & a, matrice & res, vecteur & pivots, gen & det,int l, int lmax, int c,int cmax,int fullreduction,int dont_swap_below,const gen & modulo,int rref_or_det_or_lu){
    if (modulo.type==_INT_ && 
#ifndef _I386_
	modulo.val<46340 &&
#endif
	is_integer_matrice(a) ){ // Small mod reduction
      vector< vector<int> > N;
      if (!vecteur2vectvector_int(a,modulo.val,N))
	return false;
      longlong idet=1;
      vector<int> permutation,maxrankcol;
      smallmodrref(N,pivots,permutation,maxrankcol,idet,l,lmax,c,cmax,fullreduction,dont_swap_below,modulo.val,rref_or_det_or_lu);
      det = smod(longlong(idet),modulo.val);
      if (rref_or_det_or_lu!=1)
	vectvector_int2vecteur(N,res);
      if (rref_or_det_or_lu==2){
	vecteur P;
	vector_int2vecteur(permutation,P);
	pivots.push_back(P);
      }
      return true;
    }
    // bool use_cstart=!c;
    // bool inverting=fullreduction==2;
    det = 1;
    int linit=l;//,previous_l=l;
    vecteur lv;
    // Large mod reduction (coeff do not fit in an int)
    res=a;
    // cout << res << endl;
    std_matrix<gen> M;
    matrice2std_matrix_gen(res,M);
    gen pivot,temp;
    // vecteur vtemp;
    int pivotline,pivotcol;
    pivots.clear();
    pivots.reserve(cmax-c);
    for (;(l<lmax) && (c<cmax);){
      if ( (!fullreduction) && (l==lmax-1) ){
	det = smod(det*M[l][c],modulo);
	break;
      }
      pivot=M[l][c];
      pivotline=l;
      pivotcol=c;
      if (is_exactly_zero(pivot)){ // scan current line
	if (rref_or_det_or_lu==3){
	  det=0;
	  return true;
	}
	if (l<dont_swap_below){ 
	  for (int ctemp=c+1;ctemp<cmax;++ctemp){
	    temp=M[l][ctemp];
	    if (!is_exactly_zero(temp)){
	      pivot=temp;
	      pivotcol=ctemp;
	      break;
	    }
	  }
	}
	else {      // scan M current column for the best pivot available
	  for (int ltemp=l+1;ltemp<lmax;++ltemp){
	    temp=M[ltemp][c];
	    if (debug_infolevel)
	      print_debug_info(temp);
	    if (!is_exactly_zero(temp)){
	      pivot=temp;
	      pivotline=ltemp;
	      break;
	    }
	  }
	}
      } // end if is_zero(pivot), true pivot found on line or column
      if (!is_exactly_zero(pivot)){
	if (l!=pivotline){
	  swap(M[l],M[pivotline]);
	  det = -det;
	}
	det = smod(det*pivot,modulo);
	// save pivot for annulation test purposes
	if (rref_or_det_or_lu!=1)
	  pivots.push_back(pivot);
	// invert pivot 
	temp=invmod(pivot,modulo);
	if (fullreduction || rref_or_det_or_lu<2){
	  iterateur it=M[pivotline].begin(),itend=M[pivotline].end();
	  for (;it!=itend;++it)
	    *it=smod(temp * *it,modulo);
	}
	// make the reduction
	if (fullreduction){
	  for (int ltemp=linit;ltemp<lmax;++ltemp){
	    if (ltemp!=l)
	      modlinear_combination(M[ltemp],-M[ltemp][pivotcol],M[l],modulo,0);
	  }
	}
	else {
	  for (int ltemp=l+1;ltemp<lmax;++ltemp){
	    if (rref_or_det_or_lu>=2)
	      M[ltemp][pivotcol]=smod(M[ltemp][pivotcol]*temp,modulo);
	    modlinear_combination(M[ltemp],-M[ltemp][pivotcol],M[l],modulo,(c+1)*(rref_or_det_or_lu>0));
	  }
	} // end else
	// increment column number if swap was allowed
	if (l>=dont_swap_below)
	  ++c;
	// increment line number since reduction has been done
	++l;	  
      } // end if (!is_zero(pivot)
      else { // if pivot is 0 increment either the line or the col
	det = 0;
	if (rref_or_det_or_lu==1)
	  return true;
	if (l>=dont_swap_below)
	  c++;
	else
	  l++;
      }
    } // end for reduction loop
    std_matrix_gen2matrice(M,res);
    return true;
  }

  bool mrref(const matrice & a, matrice & res, vecteur & pivots, gen & det,GIAC_CONTEXT){
    return mrref(a,res,pivots,det,0,a.size(),0,a.front()._VECTptr->size(),
	  true,0,true,1,0,
	  contextptr);
  }

  bool modrref(const matrice & a, matrice & res, vecteur & pivots, gen & det,const gen& modulo){
    return modrref(a,res,pivots,det,0,a.size(),0,a.front()._VECTptr->size(),
	    true /* full reduction */,0 /* dont_swap_below*/,modulo,0 /* rref */);
  }

  // add identity matrix, modifies arref in place
  void add_identity(matrice & arref){
    int s=arref.size();
    vecteur v;
    for (int i=0;i<s;++i){
      v = *arref[i]._VECTptr;
      v.reserve(2*s);
      for (int j=0;j<s;++j)
	v.push_back(i==j);
      arref[i] = v;
    }
  }

  // add identity matrix, modifies arref in place
  void add_identity(vector< vector<int> > & arref){
    int s=arref.size();
    for (int i=0;i<s;++i){
      vector<int> & v= arref[i];
      v.reserve(2*s);
      for (int j=0;j<s;++j)
	v.push_back(i==j);
    }
  }

  bool remove_identity(matrice & res){
    int s=res.size();
    // "shrink" res
    for (int i=0;i<s;++i){
      vecteur & v = *res[i]._VECTptr;
      if (is_zero(v[i],context0))
	return false;
      vecteur w(v.begin()+s,v.end());
      res[i] = divvecteur(w,v[i]);
    }
    return true;
  }

  bool remove_identity(vector< vector<int> > & res,int modulo){
    int s=res.size();
    // "shrink" res
    for (int i=0;i<s;++i){
      vector<int> & v = res[i];
      if (!v[i])
	return false;
      longlong inv=invmod(v[i],modulo);
      v = vector<int>(v.begin()+s,v.end());
      for (int j=0;j<s;++j){
	longlong tmp=v[j]*inv;
	v[j] = tmp % modulo;
      }
    }
    return true;
  }

  bool modinv(const matrice & a,matrice & res,const gen & modulo,gen & det_mod_p){
    matrice arref = a;
    add_identity(arref);
    int s=a.size();
    vecteur pivots;
    if (!modrref(arref,res,pivots,det_mod_p,0,s,0,2*s,
		 2/* full reduction*/,0/*dont_swap_below*/,modulo,0/* rref */))
      return false;
    return remove_identity(res);
  }

  bool smallmodinv(const vector< vector<int> > & a,vector< vector<int> > & res,int modulo,longlong & det_mod_p){
    res = a;
    add_identity(res);
    int s=a.size();
    vecteur pivots;
    vector<int> permutation,rankcols;
    smallmodrref(res,pivots,permutation,rankcols,det_mod_p,0,s,0,2*s,
	    2/* full reduction*/,0/*dont_swap_below*/,modulo,0/* rref */);
    return remove_identity(res,modulo);
  }

  // works if |v|^2,|w|^2<2^31
  static gen dotvecteur_int(const vecteur & a,const vecteur & b,bool smallint){
    vecteur::const_iterator ita=a.begin(), itaend=a.end();
    vecteur::const_iterator itb=b.begin(), itbend=b.end();
    if (smallint) {
      longlong res=0;
      for (;(ita!=itaend)&&(itb!=itbend);++ita,++itb){
	res += longlong (ita->val)*(itb->val);
      }
      return res;
    }
    ref_mpz_t * e = new ref_mpz_t;
    mpz_set_ui(e->z,0);
    gen tmp;
    for (;(ita!=itaend)&&(itb!=itbend);++ita,++itb){
      type_operator_times(*ita,*itb,tmp);
      if (tmp.type==_INT_){
	if (tmp.val<0)
	  mpz_sub_ui(e->z,e->z,-tmp.val);
	else
	  mpz_add_ui(e->z,e->z,tmp.val);
      }
      else 
	mpz_add(e->z,e->z,*tmp._ZINTptr);
    }
    return e;
  }

  static vecteur multmatvecteur_int(const matrice & a,const vecteur & b,bool smallint){
    vecteur::const_iterator ita=a.begin(), itaend=a.end();
    vecteur res;
    res.reserve(itaend-ita);
    for (;ita!=itaend;++ita)
      res.push_back(dotvecteur_int(*(ita->_VECTptr),b,smallint));
    return res;
  }

  // solve a*x=b where a and b have integer coeffs
  // using a p-adic algorithm, n is the precision required
  // c is the inverse of a mod p 
  // FIXME make it work on Z[i]
  vecteur padic_linsolve(const matrice & a,const vecteur & b,const matrice & c,unsigned n,const gen & p){
    vecteur res,y(b),x; // initialize y_0=b
    int smallint = p.type==_INT_ && 
      ((ulonglong) p.val*p.val < ((ulonglong) 1 << 63)/a.size() );
    if (smallint && ck_is_greater(p,linfnorm(a,context0),context0)) // ok
      smallint=2;
    for (unsigned i=0;i<n;++i){
      x=smod(multmatvecteur_int(c,smod(y,p),smallint!=0),p); // x_{n}=c*y_n mod p
      y=divvecteur(subvecteur(y,multmatvecteur_int(a,x,smallint==2)),p); // y_{n+1}=(y_n-Ax_n)/p
      // should use below on Z[i]
      // x=smod(multmatvecteur(c,y),p); // x_{n+1}=c*y_n mod p
      // y=divvecteur(subvecteur(y,multmatvecteur(a,x)),p); // y_{n+1}=(y_n-Ax_n)/p
      res.push_back(x);
    }
    // res=[x_0,...,x_n-1]
    reverse(res.begin(),res.end());
    // res=[x_n-1,...,x_0]
    res=*horner(res,p)._VECTptr;
    return res;
  }  
  
  bool iszero(const vector<int> & p){
    vector<int>::const_iterator it=p.begin(),itend=p.end();
    for (;it!=itend;++it){
      if (*it)
	return false;
    }
    return true;    
  }
  
  static void extract(const vector< vector<int> > & source,const vector<int> & lines,int rang,vector< vector<int> > & target,vector< vector<int> > & excluded){
    int s=lines.size(),i;
    target.clear();
    target.reserve(rang);
    for (i=0;i<rang;++i){
      const vector<int> & si=source[lines[i]];
      if (iszero(si))
	break;
      target.push_back(si);
    }
    for (;i<s;++i)
      excluded.push_back(source[lines[i]]);
  }

  static void extract1(const matrice & source,const vector<int> & lines,int rang,matrice & target,vector<vecteur> & excluded){
    int s=lines.size(),i;
    target.clear();
    target.reserve(rang);
    for (i=0;i<rang;++i){
      const vecteur & si=*source[lines[i]]._VECTptr;
      if (is_zero(si,context0))
	break;
      target.push_back(si);
    }
    for (;i<s;++i)
      excluded.push_back(*source[lines[i]]._VECTptr);
  }

  static void extract2(const matrice & source,const vector<int> & lines,int rang,matrice & target,vector<vecteur> & k_excluded,vector<int> & k_excluded_col){
    int s=source.size(),i,j=0;
    target.clear();
    target.reserve(rang);
    for (i=0;i<rang;++i,++j){
      int li=lines[i];
      for (;j<li;++j){
	k_excluded.push_back(*source[j]._VECTptr);
	k_excluded_col.push_back(j);
      }
      const vecteur & si=*source[li]._VECTptr;
      target.push_back(si);
    }
    for (;j<s;++j){
      k_excluded.push_back(*source[j]._VECTptr);
      k_excluded_col.push_back(j);
    }
  }

  static gen init_modulo(int n,double logbound){
#ifdef _I386_
    double pinit= double(longlong(1) << 62);
    pinit /=n ;
    pinit = std::sqrt(pinit);
    // so that pinit^2*n<2^63 
    pinit -= 3*logbound; // keep enough primes satisfying p^2*n<2^63
    if (pinit> 1<<30)
      pinit /=2;
    return nextprime(int(pinit)); 
#else
    return 36007;
#endif
  }

  // a is a matrix with integer coeffs
  // find p such that a mod p has the same rank
  // rankline and rankcols are the lines/cols used for the submatrix
  // asub of max rank, ainv is the inverse of asub mod p
  // return -1 or the rank
  int padic_linsolve_prepare(const matrice & a,gen & p,vector<int> & ranklines, vector<int> & rankcols,matrice & asub,matrice & ainv,vecteur & compat,vecteur & kernel){
    if (!is_integer_matrice(a))
      return -1;
    vector< vector<int> > N;
    if (!vecteur2vectvector_int(a,p.val,N))
      return -1;
    int nrows=N.size();
    int ncols=N.front().size();
    int n0=max(ncols,nrows);
    gen h2=4*square_hadamard_bound(a);
    gen B=evalf_double(linfnorm(a,context0),0,context0);
    double Bd=B._DOUBLE_val;
    double logbound=n0*(std::log10(double(n0))/2+std::log10(Bd));
    if (is_exactly_zero(p) || p.type!=_INT_)
      p=init_modulo(n0,logbound);
    for (;;){ // break loop as soon as a good p is found
      matrice atmp;
      vector< vector<int> > Ntmp,Ninv,Nsub,Nmaxrank,Nexcluded;
      vector<vecteur> excluded,k_excluded;
      Nsub=N;
      int nstep=int(min(nrows,ncols)*std::log(evalf_double(h2,1,context0)._DOUBLE_val)/2/std::log(double(p.val)))+1;
      vecteur pivots;
      longlong idet;
      smallmodrref(Nsub,pivots,ranklines,rankcols,idet,0,nrows,0,ncols,0 /* fullreduction*/,0 /* dont_swap_below */,p.val,0 /* rref_or_det_or_lu */);
      int rang=rankcols.size();
      /* extract maxrank submatrix */
      extract(N,ranklines,rang,Ntmp,Nexcluded);
      tran_vect_vector_int(Ntmp,Nsub);
      extract1(a,ranklines,rang,atmp,excluded);
      mtran(atmp,asub);
      // next call does not change excluded 
      extract(Nsub,rankcols,rang,Ntmp,Nexcluded);
      tran_vect_vector_int(Ntmp,Nsub);  
      vector<int> k_excluded_col;
      extract2(asub,rankcols,rang,atmp,k_excluded,k_excluded_col);
      mtran(atmp,asub);  
      // now asub contains the invertible matrix
      // we must truncate excluded, then rewrite each line of excluded
      // as a linear combination of lines of Np
      int es=excluded.size();
      for (int i=0;i<es;++i){
	vecteur & v=excluded[i];
	vecteur tmp=v;
	v.clear();
	for (int j=0;j<rang;++j)
	  v.push_back(tmp[rankcols[j]]);
      }
      /* invert Nsub in Ninv */
      longlong det;
      if (!smallmodinv(Nsub,Ninv,p.val,det))
	break; // should not happen!
      /* find compatibility with the lines of a that are not in Nmaxrank */
      vectvector_int2vecteur(Ninv,ainv);
      tran_vect_vector_int(Ninv,Ntmp);  
      matrice c;
      vectvector_int2vecteur(Ntmp,c);
      int i;
      gen pn=pow(p,nstep,context0);
      for (i=0;i<es;++i){
	/* p-adic lift of each compatibility equation */
	vecteur current,cond(nrows);
	cond[ranklines[rang+i]]=-1;
	current=padic_linsolve(atmp,excluded[i],c,nstep,p);
	int cs=current.size();
	for (int j=0;j<cs;++j)
	  current[j]=fracmod(current[j],pn);
	/* rewrite using ranklines */
	for (int j=0;j<rang;++j)
	  cond[ranklines[j]]=current[j];
	/* check that it is correct with the original matrix */
	if (!is_zero(multvecteurmat(cond,a),context0))
	  break;
	compat.push_back(cond);
      }
      if (i==es){
	/* compute kernel using k_excluded */
	es=k_excluded.size();
	for (i=0;i<es;++i){
	  /* p-adic lift of each kernel element basis */
	  vecteur current,cond(ncols);
	  cond[k_excluded_col[i]]=-1;
	  current=padic_linsolve(asub,k_excluded[i],ainv,nstep,p);
	  int cs=current.size();
	  for (int j=0;j<cs;++j)
	    current[j]=fracmod(current[j],pn);
	  /* rewrite using ranklines */
	  for (int j=0;j<rang;++j)
	    cond[rankcols[j]]=current[j];
	  kernel.push_back(cond);
	}
	return rang;
      }
      /* bad modulo */
      p=nextprime(p+1);
    }
    return -1;
  }

  // solve
  bool padic_linsolve_solve(const matrice & a,const gen & p,const vector<int> & ranklines,const vector<int> & rankcols,const matrice & asub,const matrice & ainv,const vecteur & compat,const vecteur & b,vecteur & sol){
    // first check that b verifies compat
    int es=compat.size();
    for (int i=0;i<es;++i){
      if (!is_exactly_zero(dotvecteur(compat[i],b))){
	return false;
      }
    }
    /* padic solve asub*x=part of b (ranklines) */
    int rang=asub.size();
    vecteur newb(rang);
    for (int i=0;i<rang;++i)
      newb[i]=b[ranklines[i]];
    gen h2=4*square_hadamard_bound(asub)*l2norm2(newb);
    int nstep=int(rang*std::log(evalf_double(h2,1,context0)._DOUBLE_val)/2/std::log(double(p.val)))+1;
    vecteur res=padic_linsolve(asub,newb,ainv,nstep,p);
    gen pn=pow(p,nstep,context0);
    int ress=res.size();
    for (int i=0;i<ress;++i)
      res[i]=fracmod(res[i],pn);
    /* find x (using rankcols) */
    int xs=a.front()._VECTptr->size();
    sol=vecteur(xs);
    for (int i=0;i<rang;++i){
      sol[rankcols[i]]=res[i];
    }
    return true;
  }

  gen _padic_linsolve(const gen & g,GIAC_CONTEXT){
    if ( g.type==_STRNG && g.subtype==-1) return  g;
    if (g.type!=_VECT || g._VECTptr->empty())
      return gensizeerr(contextptr);
    if (g.subtype==_SEQ__VECT && g._VECTptr->size()==2){
      gen a=g._VECTptr->front();
      gen b=g._VECTptr->back();
      if (!ckmatrix(a) || b.type!=_VECT)
	return gensizeerr(contextptr);
      if (a._VECTptr->front()._VECTptr->size()!=b._VECTptr->size())
	return gendimerr(contextptr);
      matrice & A=*a._VECTptr;
      gen p;
      matrice asub,ainv;
      vecteur compat,kernel;
      vector<int> ranklines,rankcols;
      if (!padic_linsolve_prepare(A,p,ranklines,rankcols,asub,ainv,compat,kernel))
	return gensizeerr(gettext("Unable to find a modulus to solve"));
      vecteur & B=*b._VECTptr;
      vecteur res;
      if (padic_linsolve_solve(A,p,ranklines,rankcols,asub,ainv,compat,B,res))
	return makevecteur(res,kernel);
      else
	return gensizeerr(gettext("Incompatible system"));
    }
    return gensizeerr(contextptr);
  }
  static const char _padic_linsolve_s []="padic_linsolve";
  static define_unary_function_eval (__padic_linsolve,&giac::_padic_linsolve,_padic_linsolve_s);
  define_unary_function_ptr5( at_padic_linsolve ,alias_at_padic_linsolve,&__padic_linsolve,0,true);

  // solve a*x=b where a and b have integer coeffs using a p-adic algorithm
  // lcmdeno of the answer may be used to give an estimate of the 
  // least divisor element of a if b is random
  // returns 0 if no invertible found, -1 if det==0, 1 otherwise
  int padic_linsolve(const matrice & a,const vecteur & b,vecteur & res,gen & p,gen & det_mod_p,unsigned reconstruct,int maxtry){
    // first find p such that a mod p is invertible
    // find a bound on the num/den of x
    // let c=(a mod p)^(-1)
    matrice c;
    matrice ab(a);
    ab.push_back(b);
    if (is_exactly_zero(p))
      p=36007;
    gen h2=4*square_hadamard_bound(ab),pip(1);
    cerr << "Modinv begin " << clock() << endl;
    for (int tryinv=0;;++tryinv){
      if (modinv(a,c,p,det_mod_p))
	break;
      pip=pip*p;
      if (tryinv>maxtry)
	return 0;
      if (is_strictly_greater(pip*pip,h2,context0)) // ok
	return -1;
      p=nextprime(p+1);
    }
    cerr << "Modinv end " << clock() << endl;
    unsigned n=1;
    gen pn=p;
    while (is_strictly_greater(h2,pn,context0)){ // ok
      ++n;
      pn = pn * p;
    }
    vecteur resp=padic_linsolve(a,b,c,n,p);
    cerr << "Padic end " << clock() << endl;
    // rational reconstruction
    unsigned s=resp.size();
    if (reconstruct)
      s=std::min(s,reconstruct);
    res.clear();
    res.reserve(s);
    for (unsigned j=0;j<s;++j){
      res.push_back(fracmod(resp[j],pn));
    }
    return 1;    
  }

  matrice mrref(const matrice & a,GIAC_CONTEXT){
    if (a.empty())
      return vecteur(vecteur(1,gendimerr(contextptr)));
    gen det;
    vecteur pivots;
    matrice res;
    if (!mrref(a,res,pivots,det,0,a.size(),0,a.front()._VECTptr->size(),
	  true,0,true,1,0,
	       contextptr))
      return vecteur(1,vecteur(1,gendimerr(contextptr)));
    return res;
  }

  bool read_reduction_options(const gen & a_orig,matrice & a,bool & convert_internal,int & algorithm,bool & minor_det,bool & keep_pivot,int & last_col){
    convert_internal=true;
    algorithm=RREF_GUESS;
    minor_det=false;
    keep_pivot=false;
    last_col=-1;
    if (ckmatrix(a_orig)){
      a=*a_orig._VECTptr;
    }
    else { // rref with options
      if (a_orig.type!=_VECT)
	return false;
      vecteur & v=*a_orig._VECTptr;
      int s=v.size();
      if (!s || !ckmatrix(v[0]))
	return false;
      a=*v[0]._VECTptr;
      for (int i=1;i<s;++i){
	if (v[i]==at_lagrange)
	  algorithm=RREF_LAGRANGE;
	if (v[i].type==_INT_){
	  if (v[i].subtype==_INT_SOLVER){
	    switch (v[i].val){
	    case _RATIONAL_DET:
	      convert_internal=false;
	      algorithm=RREF_GAUSS_JORDAN;
	      break;
	    case _BAREISS:
	      algorithm=RREF_BAREISS;
	      break;
	    case _KEEP_PIVOT:
	      keep_pivot=true;
	      break;
	    case _MINOR_DET:
	      minor_det=true;
	    }
	  }
	  else
	    last_col=v[i].val;
	}
      }
    }
    return true;
  }
  gen _rref(const gen & a_orig,GIAC_CONTEXT) {
    if ( a_orig.type==_STRNG && a_orig.subtype==-1) return  a_orig;
    matrice a;
    bool convert_internal,minor_det,keep_pivot;
    int algorithm,last_col;
    if (!read_reduction_options(a_orig,a,convert_internal,algorithm,minor_det,keep_pivot,last_col))
      return gensizeerr(contextptr);
    if (minor_det)
      return gensizeerr(gettext("minor_det option applies only to det"));
    gen det;
    vecteur pivots;
    matrice res;
    int ncols=a.front()._VECTptr->size();
    if (last_col>=0)
      ncols=min(ncols,last_col);
    if (!mrref(a,res,pivots,det,0,a.size(),0,ncols,
	  true,0,convert_internal,algorithm,0,
	       contextptr))
      return gendimerr(contextptr);
    if (!keep_pivot){
      bool reducelast = int(a.size())!=ncols-1;
      mdividebypivot(res,reducelast);
    }
    return ratnormal(res);
  }
  static const char _rref_s []="rref";
  static define_unary_function_eval (__rref,&giac::_rref,_rref_s);
  define_unary_function_ptr5( at_rref ,alias_at_rref,&__rref,0,true);

  // returns 0 if all elements are 0
  gen first_non_zero(const vecteur & v,bool uselast){
    vecteur::const_iterator it=v.begin(),itend=v.end();
    if (it!=itend && !uselast)
      --itend;
    for (;it!=itend;++it){
      if (!is_zero(*it,context0))
	return *it;
    }
    return 0;
  }

  void mdividebypivot(matrice & a,bool uselast){
    vecteur::const_iterator ita=a.begin(),itaend=a.end();
    gen pivot;
    for (;ita!=itaend;++ita){
      pivot=first_non_zero(*(ita->_VECTptr),uselast);
      if (!is_zero(pivot,context0))
	divvecteur(*(ita->_VECTptr),pivot,*(ita->_VECTptr));
    }
  }

  void midn(int n,matrice & res){
    if (n<=0 || longlong(n)*n>LIST_SIZE_LIMIT){
      res= vecteur(1,vecteur(1,gendimerr(gettext("idn"))));
      return ;
    }
    res.clear();
    res.reserve(n);
    vecteur v;
    for (int i=0;i<n;++i){
      v.clear();
      v.reserve(n);
      for (int j=0;j<n;++j)
	v.push_back(i==j);
      res.push_back(v);
    }
  }

  matrice midn(int n){
    matrice res;
    midn(n,res);
    return res;
  }

  gen _idn(const gen & e,GIAC_CONTEXT) {
    if ( e.type==_STRNG && e.subtype==-1) return  e;
    matrice res;
    if (e.type==_INT_)
      midn(e.val,res);
    else {
      if (e.type==_DOUBLE_)
	midn(int(e._DOUBLE_val),res);
      else {
	if ((e.type==_VECT) && is_squarematrix(*e._VECTptr))
	  midn(e._VECTptr->size(),res);
	else
	  return gensizeerr(contextptr);
      }
    }
    return res;
  }
  static const char _idn_s []="idn";
  static define_unary_function_eval (__idn,&giac::_idn,_idn_s);
  define_unary_function_ptr5( at_idn ,alias_at_idn,&__idn,0,true);

  vecteur vranm(int n,const gen & f,GIAC_CONTEXT){
    n=max(1,n);
    if (n>LIST_SIZE_LIMIT)
      setstabilityerr();
    vecteur res;
    for (int i=0;i<n;++i){
      if (is_zero(f,contextptr))
	res.push_back((int) (2*randrange*giac_rand(contextptr)/(rand_max+1.0)-randrange));
      else {
	if (f.type==_INT_)
	  res.push_back(_rand(f,contextptr));
	else {
	  if (f.is_symb_of_sommet(at_interval) && f._SYMBptr->feuille.type==_VECT){
	    res.push_back(_rand(f._SYMBptr->feuille,contextptr));
	  }
	  else
	    if (f.is_symb_of_sommet(at_program))
	      res.push_back(f(vecteur(0),contextptr));
	    else 
	      res.push_back(eval(f,eval_level(contextptr),contextptr));
	}
      }
    }
    return res;
  }

  matrice mranm(int n,int m,const gen & f,GIAC_CONTEXT){
    n=max(1,n);
    m=max(1,m);
    if (longlong(n)*m>LIST_SIZE_LIMIT)
      setstabilityerr();
    matrice res;
    res.reserve(n);
    for (int i=0;i<n;++i)
      res.push_back(vranm(m,f,contextptr));
    return res;
  }

  gen _ranm(const gen & e,GIAC_CONTEXT){
    if ( e.type==_STRNG && e.subtype==-1) return  e;
    int n=0,m=0;
    switch (e.type){
    case _INT_:
      return vranm(e.val,zero,contextptr);
    case _DOUBLE_:
      return vranm(int(e._DOUBLE_val),zero,contextptr);
    case _VECT:
      if (e._VECTptr->size()==1)
	return _ranm(e._VECTptr->front(),contextptr);
      if (e._VECTptr->size()>=2){
	if (e._VECTptr->front().type==_INT_)
	  n=e._VECTptr->front().val;
	else {
	  if (e._VECTptr->front().type==_DOUBLE_)
	    n=int(e._VECTptr->front()._DOUBLE_val);
	  else
	    return gensizeerr(contextptr);
	}
	if ((*e._VECTptr)[1].type==_INT_)
	  m=(*e._VECTptr)[1].val;
	else {
	  if ((*e._VECTptr)[1].type==_DOUBLE_)
	    m=int((*e._VECTptr)[1]._DOUBLE_val);
	  else
	    return gensizeerr(contextptr);
	}
	if (e._VECTptr->size()==3)
	  return mranm(n,m,e._VECTptr->back(),contextptr);
	return mranm(n,m,0,contextptr);
      }
    default:
      return gensizeerr(contextptr);
    }
    return undef;
  }
  static const char _ranm_s []="ranm";
  static define_unary_function_eval (__ranm,&giac::_ranm,_ranm_s);
  define_unary_function_ptr5( at_ranm ,alias_at_ranm,&__ranm,0,true);

  gen _randvector(const gen & e,GIAC_CONTEXT){
    if ( e.type==_STRNG && e.subtype==-1) return  e;
    int n=0;
    switch (e.type){
    case _INT_:
      return vranm(e.val,zero,contextptr);
    case _DOUBLE_:
      return vranm(int(e._DOUBLE_val),zero,contextptr);
    case _VECT:
      if (e._VECTptr->size()==1)
	return _randvector(e._VECTptr->front(),contextptr);
      if (e._VECTptr->size()==2){
	if (e._VECTptr->front().type==_INT_)
	  n=e._VECTptr->front().val;
	else {
	  if (e._VECTptr->front().type==_DOUBLE_)
	    n=int(e._VECTptr->front()._DOUBLE_val);
	  else
	    return gensizeerr(contextptr);
	}
	return vranm(n,e._VECTptr->back(),contextptr);
      }
    default:
      return gensizeerr(contextptr);
    }
    return undef;
  }
  static const char _randvector_s []="randvector";
  static define_unary_function_eval (__randvector,&giac::_randvector,_randvector_s);
  define_unary_function_ptr5( at_randvector ,alias_at_randvector,&__randvector,0,true);

#ifdef HAVE_LIBLAPACK
  bool matrix2zlapack(const std_matrix<gen> & m,doublef2c_complex * A,GIAC_CONTEXT){
    std_matrix<gen>::const_iterator it=m.begin(),itend=m.end();
    gen g;
    int rows=itend-it;
    for (int i = 0; it!=itend; ++i,++it){
      const_iterateur jt=it->begin(),jtend=it->end();
      for (int j = 0; jt!=jtend;++j, ++jt){
	g=evalf_double(*jt,1,contextptr);
	if (g.type==_DOUBLE_){
	  A[i + j * rows].r=g._DOUBLE_val;
	  A[i + j * rows].i=0;
	  continue;
	}
	if (g.type==_CPLX && g._CPLXptr->type==_DOUBLE_ && (*g._CPLXptr+1).type==_DOUBLE_){
	  A[i + j * rows].r = g._CPLXptr->_DOUBLE_val;
	  A[i + j * rows].i = (g._CPLXptr+1)->_DOUBLE_val;
	}
	else 
	  return false;
      }
    }
    return true;
  }

  bool matrice2zlapack(const matrice & m,doublef2c_complex * A,GIAC_CONTEXT){
    const_iterateur it=m.begin(),itend=m.end();
    gen g;
    int rows=itend-it;
    for (int i = 0; it!=itend; ++i,++it){
      if (it->type!=_VECT)
	return false;
      const_iterateur jt=it->_VECTptr->begin(),jtend=it->_VECTptr->end();
      for (int j = 0; jt!=jtend;++j, ++jt){
	g=evalf_double(*jt,1,contextptr);
	if (g.type==_DOUBLE_){
	  A[i + j * rows].r=g._DOUBLE_val;
	  A[i + j * rows].i=0;
	  continue;
	}
	if (g.type==_CPLX && g._CPLXptr->type==_DOUBLE_ && (*g._CPLXptr+1).type==_DOUBLE_){
	  A[i + j * rows].r = g._CPLXptr->_DOUBLE_val;
	  A[i + j * rows].i = (g._CPLXptr+1)->_DOUBLE_val;
	}
	else 
	  return false;
      }
    }
    return true;
  }

  void zlapack2matrix(doublef2c_complex * A,unsigned rows,unsigned cols,std_matrix<gen> & R){
    R.resize(rows);
    for (int i=0;i<rows;++i){
      vecteur r(cols);
      for (int j=0;j<cols;++j)
	r[j] = gen(A[i + j * rows].r,A[i + j * rows].i);
      R[i]=r;
    }
  }

  void zlapack2matrice(doublef2c_complex * A,unsigned rows,unsigned cols,matrice & R){
    R.resize(rows);
    for (int i=0;i<rows;++i){
      vecteur r(cols);
      for (int j=0;j<cols;++j)
	r[j] = gen(A[i + j * rows].r,A[i + j * rows].i);
      R[i]=r;
    }
  }
#endif

  bool minv(const matrice & a,matrice & res,bool convert_internal,int algorithm,GIAC_CONTEXT){
#ifdef HAVE_LIBLAPACK
    if (is_squarematrix(a) && is_fully_numeric(a)){
      integer N,LDA,INFO,LWORK;
      int n=a.size();
      LDA=n; N=n; LWORK=N*N;
      integer * IPIV=new integer[n];
      if (is_zero(im(a,contextptr))){
	double * A = new double[N*N];
	matrice2lapack(a,A,contextptr);
	dgetrf_( &N, &N, A, &LDA, IPIV, &INFO );
	if (INFO){
	  delete [] IPIV;
	  delete [] A;
	  return false;
	}
	double * WORK=new double [LWORK];
	/* DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO ) */
	dgetri_(&N,A,&LDA,IPIV,WORK,&LWORK,&INFO);
	delete [] IPIV; delete [] WORK;
	if (INFO){
	  delete [] A; 
	  return false;
	}
	lapack2matrice(A,N,N,res);
	delete [] A; 
	return true;
      }
      doublef2c_complex * A = new doublef2c_complex[N*N];
      matrice2zlapack(a,A,contextptr);
      zgetrf_( &N, &N, A, &LDA, IPIV, &INFO );
      if (INFO){
	delete [] IPIV;
	delete [] A;
	return false;
      }
      doublef2c_complex * WORK=new doublef2c_complex [LWORK];
      /* ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO ) */
      zgetri_(&N,A,&LDA,IPIV,WORK,&LWORK,&INFO);
      delete [] IPIV; delete [] WORK;
      if (INFO){
	delete [] A; 
	return false;
      }
      zlapack2matrice(A,N,N,res);
      delete [] A; 
      return true;
    }
#endif
    matrice arref = a;
    add_identity(arref);
    int s=a.size();
    gen det;
    vecteur pivots;
    if (!mrref(arref,res,pivots,det,0,s,0,2*s,
	  true,0,convert_internal,algorithm,0,
	       contextptr))
      return false;
    if (!remove_identity(res))
      return false;
    return true;
  }

  matrice minv(const matrice & a,GIAC_CONTEXT){
    matrice res;
    if (!minv(a,res,/*convert_internal */true,/* algorithm */ 1,contextptr))
      return vecteur(1,vecteur(1,gensizeerr(gettext("Not invertible"))));
    return res;
  }

  // determinant by expanding wrt last column
  gen det_minor(const matrice & a,bool convert_internal,GIAC_CONTEXT){
    int n=a.size();
    if (n==1)
      return a.front()._VECTptr->front();
    std_matrix<gen> A;
    vecteur lv;
    if (convert_internal){
      lv=alg_lvar(a);
      matrice2std_matrix_gen(*(e2r(a,lv,contextptr)._VECTptr),A);
    }
    else
      matrice2std_matrix_gen(a,A);
    index_t index(n),mineur_index(n);
    map< index_t, gen > tab_mineurs,old_tab;
    // int s=int(std::exp(lgamma(n+1)-2*lgamma(n/2+1)))+1;
    // init: compute 2*2 determinants lines i,j columns 1,2
    gen res;
    for (int i=0;i<n;++i){
      index[0]=i;
      for (int j=i+1;j<n;++j){
	index[1]=j;
	res=A[i][0]*A[j][1]-A[i][1]*A[j][0];
	tab_mineurs[index]=res;
      }
    }
    // compute all possibles i*i det with columns 0..i sing (i-1)*(i-1) det
    for (int i=2;i<n;++i){
      if (debug_infolevel)
	cerr << "// Computing " << i+1 << "*" << i+1 << "minors " << clock() << endl;
      swap(old_tab,tab_mineurs);
      tab_mineurs.clear();
      // initialize index
      for (int j=0;j<=i;++j)
	index[j]=j;
      // computation loop
      for (;;){
	res=zero;
	for (int j=0;j<=i;++j){
	  // make mineur without line index[j]
	  for (int k=0;k<=i;++k){
	    if (k==j)
	      continue;
	    if (k>j)
	      mineur_index[k-1]=index[k];
	    else
	      mineur_index[k]=index[k];
	  }
	  if ((i+j)%2)
	    res = res-A[index[j]][i]*old_tab[mineur_index];
	  else
	    res = res+A[index[j]][i]*old_tab[mineur_index];
	}
	tab_mineurs[index]=res;
	// increment index and test for breaking loop
	int j=i;
	for (;j>=0;--j){
	  ++index[j];
	  if (index[j]!=n+j-i)
	    break;
	}
	if (j<0)
	  break;
	for (;j<i;++j)
	  index[j+1]=index[j]+1;
      }
    }
    if (debug_infolevel)
      cerr << "// Computation done " << clock() << endl;
    if (convert_internal)
      return r2sym(res,lv,contextptr);
    else
      return res;
  }

  gen _det_minor(const gen & a,GIAC_CONTEXT){
    if ( a.type==_STRNG && a.subtype==-1) return  a;
    if (!is_squarematrix(a))
      return symbolic(at_det_minor,a);
    return det_minor(*a._VECTptr,true,contextptr);
  }
  static const char _det_minor_s []="det_minor";
  static define_unary_function_eval (__det_minor,(const gen_op_context)giac::_det_minor,_det_minor_s);
  define_unary_function_ptr5( at_det_minor ,alias_at_det_minor,&__det_minor,0,true);

  gen mdet(const matrice & a,GIAC_CONTEXT){
    if (!is_squarematrix(a))
      return gendimerr(contextptr);
    vecteur pivots;
    matrice res;
    gen determinant;
    int s=a.size();
    if (!mrref(a,res,pivots,determinant,0,s,0,s,
	  false,0,true,1/* guess algorithm */,1/* determinant */,
	       contextptr))
      return gendimerr(contextptr);
    return determinant;
  }

  gen _det(const gen & a_orig,GIAC_CONTEXT){
    if ( a_orig.type==_STRNG && a_orig.subtype==-1) return  a_orig;
    matrice a;
    bool convert_internal,minor_det,keep_pivot;
    int algorithm,last_col;
    if (!read_reduction_options(a_orig,a,convert_internal,algorithm,minor_det,keep_pivot,last_col))
      return gensizeerr(contextptr);
    if (keep_pivot)
      return gensizeerr(gettext("Option keep_pivot not applicable"));
    if (minor_det)
      return det_minor(a,convert_internal,contextptr);
    vecteur pivots;
    matrice res;
    gen determinant;
    int s=a.size();
    if (!mrref(a,res,pivots,determinant,0,s,0,s,
	  false,0,convert_internal,algorithm,1/* det */,
	       contextptr))
      return gendimerr(contextptr);
    return determinant;
  }
  static const char _det_s []="det";
  static define_unary_function_eval (__det,&giac::_det,_det_s);
  define_unary_function_ptr5( at_det ,alias_at_det,&__det,0,true);

  // Find minimal poly by trying with 3 random vectors
  bool probabilistic_pmin(const matrice & m,vecteur & w,bool check,GIAC_CONTEXT){
    int n=m.size();
    modpoly p;
    for (int i=0;i<3;++i){
      vecteur v(vranm(n,0,0));
      // /* Old algorithm
      matrice temp(1,v);
      for (int j=0;j<n;++j){
	v=multmatvecteur(m,v);
	temp.push_back(v);
      }
      temp=mtran(temp);
      temp=mker(temp,contextptr);
      if (temp.empty() || is_undef(temp))
	return false; // setsizeerr();
      w=-*temp.front()._VECTptr;
      reverse(w.begin(),w.end());
      w=trim(w,0);
      // */
      /*
      // New algorithm using A^(2n-1)v and Pade
      vecteur temp(1,v[0]);
      for (int j=1;j<2*n;++j){
	v=multmatvecteur(m,v);
	temp.push_back(v[0]);
      }
      w=reverse_rsolve(temp,false);
      // End new algorithù
      */
      if (signed(w.size())!=n+1 && !p.empty())
	w=lcm(w,p,0);
      p=w;
      if (signed(w.size())==n+1){
	w=w/w.front();
	return true;
      }
    }
    if (!check)
      return false;
    gen res=horner(w,m);
    return is_zero(res,contextptr);
  }

  // Reduction to Hessenberg form, see e.g. Cohen algorithm 2.2.9
  // (with C array indices)
  // integer modulo case
  void mhessenberg(vector< vector<int> > & H,vector< vector<int> > & P,int modulo){
    int t,u,n=H.size();
    vecteur vtemp;
    for (int m=0;m<n-2;++m){
      if (debug_infolevel>=2)
	cerr << "// hessenberg reduction line " << m << endl;
      // check for a non zero coeff in the column m below ligne m+1
      int i=m+1;
      for (;i<n;++i){
	t=H[i][m];
	if (t)
	  break;
      }
      if (i==n) //not found
	continue;
      t=invmod(t,modulo);
      // permutation of lines m+1 and i and columns m+1 and i
      if (i>m+1){
	vector<int>::iterator Hib=H[i].begin(),Hie=H[i].end(),Hmp1b=H[m+1].begin();
	for (;Hib!=Hie;++Hmp1b,++Hib){
	  swapint(*Hib,*Hmp1b);
	  // tmp=H[i][j]; H[i][j]=H[m+1][j]; H[m+1][j]=tmp;
	}
	Hib=P[i].begin(); Hie=P[i].end();Hmp1b=P[m+1].begin();
	for (;Hib!=Hie;++Hmp1b,++Hib){
	  swapint(*Hib,*Hmp1b);
	  // tmp=H[i][j]; H[i][j]=H[m+1][j]; H[m+1][j]=tmp;
	}
	for (int j=0;j<n;++j){
	  swapint(H[j][i],H[j][m+1]);
	  // tmp=H[j][i]; H[j][i]=H[j][m+1]; H[j][m+1]=tmp;
	}
      }
      // now coeff at line m+1 column m is H[m+1][m]=t!=0
      // creation of zeros in column m+1, lines i=m+2 and below
      vector<int> & Hmp1=H[m+1];
      vector<int> & Pmp1=P[m+1];
      for (i=m+2;i<n;++i){
	// line operation
	vector<int> & Hi=H[i];
	u=( (longlong) t*Hi[m]) % modulo;
	if (debug_infolevel>=2)
	  cerr << "// i=" << i << " " << u <<endl;
	modlinear_combination(Hi,-u,Hmp1,modulo,0); // H[i]=H[i]-u*H[m+1];
	vector<int> & Pi=P[i];
	modlinear_combination(Pi,-u,Pmp1,modulo,0); // P[i]=P[i]-u*P[m+1];
	// column operation
	for (int j=0;j<n;++j){
	  vector<int> & Hj=H[j];
#ifdef _I386_
	  mod(Hj[m+1],u,Hj[i],modulo);
#else
	  Hj[m+1]=(Hj[m+1]+u*Hj[i])%modulo;
#endif
	}
      }
    }
  }

  // Hessenberg reduction, P is not orthogonal
  // P^(-1)*H*P = original
  void hessenberg(std_matrix<gen> & H,std_matrix<gen> & P,GIAC_CONTEXT){
    int n=H.size();
    gen t,tabs,u,tmp;
    vecteur vtemp;
    for (int m=0;m<n-2;++m){
      if (debug_infolevel>=2)
	cerr << "// hessenberg reduction line " << m << endl;
      // check for a non zero coeff in the column m below ligne m+1
      int i=m+1;
      gen pivot=0;
      int pivotline=0;
      for (;i<n;++i){
	t=H[i][m];
	tabs=abs(t,contextptr);
	if (is_strictly_greater(tabs,pivot,contextptr)){
	  pivotline=i;
	  pivot=tabs;
	}
      }
      if (is_zero(pivot)) //not found
	continue;
      i=pivotline;
      t=H[i][m];
      // permutation of lines m+1 and i and columns m+1 and i
      /*
      if (i>m+1){
	for (int j=0;j<n;++j){
	  tmp=H[i][j];
	  H[i][j]=H[m+1][j];
	  H[m+1][j]=tmp;
	}
	for (int j=0;j<n;++j){
	  tmp=H[j][i];
	  H[j][i]=H[j][m+1];
	  H[j][m+1]=tmp;
	}
      }
      */
      if (i>m+1){
	swap(H[i],H[m+1]);
	swap(P[i],P[m+1]);
	for (int j=0;j<n;++j){
	  vecteur & Hj=H[j];
	  swapgen(Hj[i],Hj[m+1]);
	}
      }
      // now coeff at line m+1 column m is H[m+1][m]=t!=0
      // creation of zeros in column m+1, lines i=m+2 and below
      for (i=m+2;i<n;++i){
	// line operation
	u=rdiv(H[i][m],t);
	if (debug_infolevel>=2)
	  cerr << "// i=" << i << " " << u <<endl;
	linear_combination(plus_one,H[i],-u,H[m+1],plus_one,vtemp,1e-12,0); // H[i]=H[i]-u*H[m+1];
	swap(H[i],vtemp);
	linear_combination(plus_one,P[i],-u,P[m+1],plus_one,vtemp,1e-12,0); // H[i]=H[i]-u*H[m+1];
	swap(P[i],vtemp);
	// column operation
	for (int j=0;j<n;++j){
	  vecteur & Hj=H[j];
	  tmp=Hj[m+1]+u*Hj[i];
	  Hj[m+1]=tmp;
	}
      }
    }
  }

  // QR reduction, P is orthogonal and should be initialized to identity
  // trn(P)*H=original
  void qr_ortho(std_matrix<gen> & H,std_matrix<gen> & P,GIAC_CONTEXT){
    int n=H.size();
    gen t,tn,tc,tabs,u,un,uc,tmp1,tmp2,norme;
    vecteur v1,v2;
    for (int m=0;m<n-1;++m){
      if (debug_infolevel>=2)
	cerr << "// hessenberg reduction line " << m << endl;
      // check for a non zero coeff in the column m below ligne m
      int i=m;
      gen pivot=0;
      int pivotline=0;
      for (;i<n;++i){
	t=H[i][m];
	tabs=abs(t,contextptr);
	if (is_strictly_greater(tabs,pivot,contextptr)){
	  pivotline=i;
	  pivot=tabs;
	}
      }
      if (is_zero(pivot)) //not found
	continue;
      i=pivotline;
      // exchange lines 
      if (i>m){
	swap(H[i],H[m]);
	swap(P[i],P[m]);
      }
      // now coeff at line m column m is H[m][m]=t!=0
      // creation of zeros in lines i=m+1 and below
      for (i=m+1;i<n;++i){
	// line operation
	t=H[m][m];
	if (is_zero(t)){
	  swap(H[i],H[m]);
	  swap(P[i],P[m]);
	  t=H[m][m];
	}
	u=H[i][m];
	if (is_zero(u))
	  continue;
	uc=conj(u,contextptr);
	tc=conj(t,contextptr);
	norme=sqrt(u*uc+t*tc,contextptr);
	un=u/norme; tn=t/norme; uc=conj(un,contextptr);	tc=conj(tn,contextptr); 
	if (debug_infolevel>=2)
	  cerr << "// i=" << i << " " << u <<endl;
	// H[m]=un*H[i]+tn*H[m] and H[i]=tn*H[i]-un*H[m];
	linear_combination(uc,H[i],tc,H[m],plus_one,v1,1e-12,0); 
	linear_combination(tn,H[i],-un,H[m],plus_one,v2,1e-12,0); 
	swap(H[m],v1);
	swap(H[i],v2);
	linear_combination(uc,P[i],tc,P[m],plus_one,v1,1e-12,0); 
	linear_combination(tn,P[i],-un,P[m],plus_one,v2,1e-12,0); 
	swap(P[m],v1);
	swap(P[i],v2);
      }
    }
  }

  // IMPROVE: don't do operations with 0
  void qr_rq(std_matrix<gen> & H,std_matrix<gen> & P,const gen & shift,int n,int & nitershift0,GIAC_CONTEXT){
    gen t,tn,tc,tabs,uabs,t2,u,un,uc,tmp1,tmp2,norme;
    int n_orig=H.size();
    vecteur v1,v2,TN(n_orig),UN(n_orig);
    if (is_zero(shift)){
      nitershift0++;
    }
    else{
      for (int i=0;i<n_orig;++i){
	H[i][i] -= shift;
      }
    }
    // H -> H-shift*identity
    for (int m=0;m<n-1;++m){
      // reduce coeff line m+1, col m
      t=H[m][m];
      // if (is_zero(t)) *logptr(contextptr) << "qr iteration: 0 on diagonal");
      int i=m+1;
      u=H[i][m];
      // now coeff at line m+1 column m is H[m+1][m]=t
      // creation of zeros in column m+1, lines i=m+2 and below
      // normalization of t and u such that t is real positive
      tabs=abs(t,contextptr);
      uabs=abs(u,contextptr);
      if (is_strictly_greater(uabs/tabs,1,contextptr))
	t2=uabs/u;
      else
	t2=tabs/t;
      t=t*t2; 
      u=u*t2;
      // compute unitary matrix coefficients
      tc=conj(t,contextptr);
      uc=conj(u,contextptr);
      norme=sqrt(re(u*uc+t*tc,contextptr),contextptr);
      un=u/norme; tn=t/norme; uc=conj(un,contextptr);	tc=conj(tn,contextptr); 
      // line operation
      // H[m]=uc*H[i]+tc*H[m] and H[i]=tn*H[i]-un*H[m];
      linear_combination(uc,H[i],tc,H[m],plus_one,v1,1e-12,0); 
      linear_combination(tn,H[i],-un,H[m],plus_one,v2,1e-12,0); 
      swap(H[m],v1);
      swap(H[i],v2);
      linear_combination(uc,P[i],tc,P[m],plus_one,v1,1e-12,0); 
      linear_combination(tn,P[i],-un,P[m],plus_one,v2,1e-12,0); 
      swap(P[m],v1);
      swap(P[i],v2);
      TN[m]=tn;
      UN[m]=un;
    } // end for m
    for (int m=0;m<n-1;++m){
      tn=TN[m];
      un=UN[m];
      tc=conj(tn,contextptr);
      uc=conj(un,contextptr);
      // column operation
      // int nmax=n_orig>m+3?m+3:n_orig;
      for (int j=0;j<n_orig;++j){
	vecteur & Hj=H[j];
	gen & Hjm=Hj[m];
	gen & Hjm1=Hj[m+1];
	tmp1=tn*Hjm+un*Hjm1;
	tmp2=-uc*Hjm+tc*Hjm1;
	Hjm=tmp1;
	Hjm1=tmp2;
      }
    }
    if (!is_zero(shift)){
      for (int i=0;i<n_orig;++i){
	H[i][i] += shift;
      }
    }
  }

  void re(std_matrix<gen> & H,int n,GIAC_CONTEXT){
    for (int i=0;i<n;i++){
      for (int j=0;j<n;j++){
	H[i][j]=re(H[i][j],contextptr);
      }
    }
  }

  bool convert(const vecteur & v,vector<giac_double> & v1){
    int n=v.size();
    v1.clear();
    v1.reserve(n);
    for (int i=0;i<n;++i){
      if (v[i].type==_INT_){
	v1.push_back(v[i].val);
	continue;
      }
      if (v[i].type==_FLOAT_){
	v1.push_back(get_double(v[i]._FLOAT_val));
	continue;
      }
      if (v[i].type!=_DOUBLE_)
	return false;
      v1.push_back(v[i]._DOUBLE_val);
    }
    return true;
  }

  bool std_matrix_gen2std_matrix_giac_double(const std_matrix<gen> & H,matrix_double & H1){
#ifdef VISUALC
    // return false;
#endif
    int n=H.size();
    H1.resize(n);
    for (int i=0;i<n;++i){
      if (!convert(H[i],H1[i]))
	return false;
    }
    return true;
  }

  bool convert(const vector<giac_double> & v,vecteur & v1){
    int n=v.size();
    v1.resize(n);
    for (int i=0;i<n;++i){
      v1[i]=double(v[i]);
    }
    return true;
  }

  bool std_matrix_giac_double2std_matrix_gen(const matrix_double & H,std_matrix<gen> & H1){
    int n=H.size();
    H1.resize(n);
    for (int i=0;i<n;++i){
      if (!convert(H[i],H1[i]))
	return false;
    }
    return true;
  }

  void hessenberg_ortho(std_matrix<gen> & H,std_matrix<gen> & P,GIAC_CONTEXT){
    hessenberg_ortho(H,P,-1,-1,true,0,contextptr);
  }

  // v=(c1*v1+c2*v2), begin at cstart
  void linear_combination(const gen & c1,const vecteur & v1,const gen & c2,const vecteur & v2,vecteur & v,int cstart,double eps){
    eps=0;
    if (cstart<0)
      cstart=0;
    const_iterateur it1=v1.begin()+cstart,it1end=v1.end(),it2=v2.begin()+cstart;
    iterateur jt1=v.begin()+cstart;
#ifdef DEBUG_SUPPORT
    if (it1end-it1!=v2.end()-it2)
      setdimerr();
#endif
    if (it2==jt1){
      linear_combination(c2,v2,c1,v1,v,cstart,eps);
      return;
    }
    if (it1==jt1){
      for (;jt1!=it1end;++jt1,++it2){
	*jt1=trim(c1*(*jt1)+c2*(*it2),c1,eps);
      }
      return;
    }
    if (v.size()==it1end-it1){
      jt1=v.begin();
      for (int i=0;i<cstart;++i,++jt1)
	*jt1=0;
      for (;it1!=it1end;++it1,++it2,++jt1)
	*jt1=trim(c1*(*it1)+c2*(*it2),c1,eps);
      return;
    }
    v.clear();
    v.reserve(it1end-it1);
    for (int i=0;i<cstart;++i)
      v.push_back(0);
    for (;it1!=it1end;++it1,++it2)
      v.push_back(trim(c1*(*it1)+c2*(*it2),c1,eps));
  }

  matrice H0;

  void dbg_schur(const std_matrix<gen> & H,const std_matrix<gen> & P){
    matrice Hg,Pg;
    std_matrix_gen2matrice(H,Hg);
    std_matrix_gen2matrice(P,Pg);
    matrice res=mmult(mtran(Pg),Hg);
    res=mmult(res,Pg);
    gen t=subvecteur(res,H0);
    gen t1=_max(_abs(t,context0),context0);
    if (t1._DOUBLE_val>1e-5)
      cerr << "Error" << endl;
  }

  // Hessenberg reduction, P is orthogonal and should be initialized to identity
  // trn(P)*H*P=original
  // already_zero is either <=0 or an integer such that H[i][j]==0 if i>j+already_zero
  // (already_zero==1 if H is hessenberg, ==3 for Francis algorithm)
  void hessenberg_ortho(std_matrix<gen> & H,std_matrix<gen> & P,int firstrow,int n,bool compute_P,int already_zero,GIAC_CONTEXT){
    int nH=H.size();
    if (n<0 || n>nH) 
      n=nH;
    if (firstrow<0 || firstrow>n)
      firstrow=0;
    gen t,tn,tc,tabs,u,un,uc,tmp1,tmp2,norme;
    vecteur v1(nH),v2(nH),TN(n,1),UN(n);
    for (int m=firstrow;m<n-2;++m){
      if (debug_infolevel>=2)
	cerr << "// hessenberg reduction line " << m << endl;
      // check for a non zero coeff in the column m below ligne m+1
      int i=m+1;
      gen pivot=0;
      int pivotline=0;
      for (;i<n;++i){
	t=H[i][m];
	tabs=abs(t,contextptr);
	if (is_strictly_greater(tabs,pivot,contextptr)){
	  pivotline=i;
	  pivot=tabs;
	}
      }
      if (is_zero(pivot)) //not found
	continue;
      i=pivotline;
      // exchange line and columns
      if (i>m+1){
	swap(H[i],H[m+1]);
	if (compute_P)
	  swap(P[i],P[m+1]);
	for (int j=0;j<n;++j){
	  vecteur & Hj=H[j];
	  swapgen(Hj[i],Hj[m+1]);
	}
      }
      // now coeff at line m+1 column m is H[m+1][m]=t!=0
      // creation of zeros in column m+1, lines i=m+2 and below
      // if (firstrow==100) dbg_schur(H,P);
      int nprime=n;
      for (i=m+2;i<nprime;++i){
	// line operation
	t=H[m+1][m];
	u=H[i][m];
	// cerr << t << " " << u << endl;
	uc=conj(u,contextptr);
	tc=conj(t,contextptr);
	norme=sqrt(u*uc+t*tc,contextptr);
	un=u/norme; tn=t/norme; uc=conj(un,contextptr);	tc=conj(tn,contextptr); 
	if (is_zero(un)){
	  UN[i]=0;
	  continue;
	}
	if (debug_infolevel>=2)
	  cerr << "// i=" << i << " " << u <<endl;
	// H[m+1]=tc*H[m+1]+uc*H[i] and H[i]=tn*H[i]-un*H[m+1];
	linear_combination(uc,H[i],tc,H[m+1],v1,0,1e-12); 
	linear_combination(tn,H[i],-un,H[m+1],v2,0,1e-12); 
	swap(H[m+1],v1);
	swap(H[i],v2);
	if (compute_P){
	  linear_combination(uc,P[i],tc,P[m+1],v1,0,1e-12); 
	  linear_combination(tn,P[i],-un,P[m+1],v2,0,1e-12); 
	  swap(P[m+1],v1);
	  swap(P[i],v2);
	}
	TN[i]=tn;
	UN[i]=un;
      }
      for (i=m+2;i<nprime;++i){      
	if (is_zero( (un=UN[i]) ) )
	  continue;
	tn=TN[i];
	tc=conj(tn,contextptr);
	uc=conj(un,contextptr);
	// column operation
	for (int j=0;j<nH;++j){
	  vecteur & Hj=H[j];
	  tmp1=tn*Hj[m+1]+un*Hj[i];
	  tmp2=-uc*Hj[m+1]+tc*Hj[i];
	  Hj[m+1]=tmp1;
	  Hj[i]=tmp2;
	}
      }
      // if (firstrow==100) dbg_schur(H,P);
    }
  }

  void tri_linear_combination(const gen & c1,const vecteur & x1,const gen & c2,const vecteur & x2,const gen & c3,const vecteur & x3,vecteur & y){
    const_iterateur it1=x1.begin(),it2=x2.begin(),it3=x3.begin(),it3end=x3.end();
    iterateur jt=y.begin();
    for (;it3!=it3end;++jt,++it1,++it2,++it3){
      *jt=c1*(*it1)+c2*(*it2)+c3*(*it3);
    }
  }

  // Francis algorithm on submatrix rows and columns n1..n2-1
  // Invariant: trn(P)*H*P=orig matrix
  bool francis_schur(std_matrix<gen> & H,int n1,int n2,std_matrix<gen> & P,int maxiter,double eps,bool is_hessenberg,bool complex_schur,bool compute_P,GIAC_CONTEXT){
    vecteur eigenv;
    if (n1==0 && eps>1e-15 && lapack_schur(H,P,compute_P,eigenv,contextptr))
      return true;
#ifdef VISUALC // tested on 100 rand poly of degree 10, should work...
    // return false; 
#endif
    int n_orig=H.size(),nitershift0=0;
    if (!is_hessenberg){
      std_matrix_gen2matrice(H,H0);
      hessenberg_ortho(H,P,0,n_orig,compute_P,0,contextptr); // insure Hessenberg form (on the whole matrix)
    }
    if (n2-n1<=1)
      return true; // nothing to do
    if (n2-n1==2){ // 2x2 submatrix, we know how to diagonalize
      gen l1,l2;
      if (eigenval2(H,n2,l1,l2,contextptr) || complex_schur){
	// [[a,b],[c,d]] -> [b,l1-a] or [l1-d,c] as first eigenvector
	gen a=H[n2-2][n2-2],b=H[n2-2][n2-1],c=H[n2-1][n2-2],d=H[n2-1][n2-1];
	gen l1a=l1-a,l1d=l1-d,x,y,xy;
	if (is_greater(abs(l1a,contextptr),abs(l1d,contextptr),contextptr)){
	  x=b; y=l1a;
	}
	else {
	  x=l1d; y=c;
	}
	xy=sqrt(x*conj(x,contextptr)+y*conj(y,contextptr),contextptr);
	// normalize eigenvector
	x = x/xy; y = y/xy;
	if (is_greater(eps,abs(y,contextptr),contextptr))
	  return true;
	if (fabs(evalf_double(re(x-1,contextptr),1,contextptr)._DOUBLE_val)<0.5)
	  x=-y*y/(x+1);
	else
	  x-=1;
	xy=sqrt(x*conj(x,contextptr)+y*conj(y,contextptr),contextptr);
	x = x/xy; y=y/xy;
	// compute reflection matrix let n=[[x],[y]] trn(n)=conj([[x,y]])
	// Q=idn(2)-2*n*trn(n);
	// i.e. [[ 1-2x*conj(x), -2x*conj(y)   ]   -2x*conj(z)  ],
	//       [ -2*y*conj(x), 1-2*y*conj(y) ]]
	// apply Q on H and P: line operations on H and P
	gen c11=1-2*x*conj(x,contextptr), c12=-2*x*conj(y,contextptr),
	  c21=conj(c12,contextptr), c22=1-2*y*conj(y,contextptr),tmp1,tmp2;
	vecteur v1(n_orig),v2(n_orig);
	linear_combination(c11,H[n1],c12,H[n1+1],v1,0,eps);
	linear_combination(c21,H[n1],c22,H[n1+1],v2,0,eps);
	swap(H[n1],v1);
	swap(H[n1+1],v2);
	if (compute_P){
	  linear_combination(c11,P[n1],c12,P[n1+1],v1,0,eps);
	  linear_combination(c21,P[n1],c22,P[n1+1],v2,0,eps);
	  swap(P[n1],v1);
	  swap(P[n1+1],v2);
	}
	// now columns operations on H (not on P)
	for (int j=0;j<n_orig;++j){
	  vecteur & Hj=H[j];
	  gen & Hjm1=Hj[n1];
	  gen & Hjm2=Hj[n1+1];
	  tmp1=Hjm1*c11+Hjm2*c21;
	  tmp2=Hjm1*c12+Hjm2*c22;
	  Hjm1=tmp1;
	  Hjm2=tmp2;
	}
      }
      return true;
    }
    vecteur v1(n_orig),v2(n_orig),v3(n_orig);
    gen tmp1,tmp2,tmp3;
    for (int niter=0;n2-n1>2 && niter<maxiter;niter++){
      // make 0 below subdiagonal
      for (int i=2;i<n_orig;i++){
	vecteur & Hi=H[i];
	for (int j=0;j<i-1;j++){
	  Hi[j]=0;
	}
      }
      if (debug_infolevel>=2){
	cerr << "// qr iteration number " << niter << " " << endl;
	H.dbgprint();
      }
      // check if one subdiagonal element is sufficiently small, if so 
      // we can increase n1 or decrease n2 or split
      for (int i=n1;i<n2-1;++i){
	gen ratio=abs(H[i+1][i]/H[i][i],contextptr);
	ratio=evalf_double(ratio,1,contextptr);
	if (ratio.type==_DOUBLE_ && fabs(ratio._DOUBLE_val)<eps){
	  if (debug_infolevel>2)
	    cerr << "Francis split " << n1 << " " << i+1 << " " << n2 << endl;
	  // submatrices n1..i and i+1..n2-1
	  if (!francis_schur(H,n1,i+1,P,maxiter,eps,true,complex_schur,compute_P,contextptr))
	    return false;
	  return francis_schur(H,i+1,n2,P,maxiter,eps,true,complex_schur,compute_P,contextptr);
	}
      }
      // now H is proper hessenberg (indices n1 to n2-1)
      // find eigenvalues l1 and l2 of last 2x2 matrix, they will be taken as shfits
      gen s,p; // s=l1+l2, p=l1*l2
      s=H[n2-2][n2-2]+H[n2-1][n2-1];
      p=H[n2-2][n2-2]*H[n2-1][n2-1]-H[n2-1][n2-2]*H[n2-2][n2-1];
      // compute (H-l2)(H-l1)=(H-s)*H+p on n1-th basis vector (if n1==0, on [1,0,...,0])
      gen ha=H[n1][n1],hb=H[n1][n1+1],
	hd=H[n1+1][n1],he=H[n1+1][n1+1],
	hh=H[n1+2][n1+1];
      gen x=hb*hd+ha*(ha-s)+p,y=hd*(he-s+ha),z=hd*hh;
      // normalize, substract [1,0,0] and normalize again
      gen xyz=sqrt(x*conj(x,contextptr)+y*conj(y,contextptr)+z*conj(z,contextptr),contextptr);
      // if x/xyz is near 1, improve precision:
      // x/xyz-1 = ((x/xyz)^2-1)/(x/xyz+1)=-((y/xyz)^2+(z/xyz)^2)/(x/xyz+1)
      x=x/xyz; y=y/xyz; z=z/xyz;
      if (fabs(evalf_double(re(x,contextptr)-1,1,contextptr)._DOUBLE_val)<0.5)
	x=-(y*y+z*z)/(x+1);
      else
	x-=1;
      xyz=sqrt(x*conj(x,contextptr)+y*conj(y,contextptr)+z*conj(z,contextptr),contextptr);
      x=x/xyz; y=y/xyz; z=z/xyz;
      // compute reflection matrix let n=[[x],[y],[z]] trn(n)=conj([[x,y,z]])
      // Q=idn(3)-2*n*trn(n);
      // i.e. [[ 1-2x*conj(x), -2x*conj(y),   -2x*conj(z)  ],
      //       [ -2*y*conj(x), 1-2*y*conj(y), -2*y*conj(z) ],
      //       [ -2*z*conj(x), -2*z*conj(y),  1-2*z*conj(z)]]
      // apply Q on H and P: line operations on H and P
      gen c11=1-2*x*conj(x,contextptr),c12=-2*x*conj(y,contextptr),c13=-2*x*conj(z,contextptr);
      gen c21=-2*y*conj(x,contextptr),c22=1-2*y*conj(y,contextptr),c23=-2*y*conj(z,contextptr);
      gen c31=-2*z*conj(x,contextptr),c32=-2*z*conj(y,contextptr),c33=1-2*z*conj(z,contextptr);
      // cerr << "[[" << c11 <<"," << c12 << "," << c13 << "],[" <<  c21 <<"," << c22 << "," << c23 << "],[" << c31 <<"," << c32 << "," << c33 << "]]" << endl;
      tri_linear_combination(c11,H[n1],c12,H[n1+1],c13,H[n1+2],v1);
      tri_linear_combination(c21,H[n1],c22,H[n1+1],c23,H[n1+2],v2);
      tri_linear_combination(c31,H[n1],c32,H[n1+1],c33,H[n1+2],v3);
      swap(H[n1],v1);
      swap(H[n1+1],v2);
      swap(H[n1+2],v3);
      // cerr << H << endl;
      if (compute_P){
	tri_linear_combination(c11,P[n1],c12,P[n1+1],c13,P[n1+2],v1);
	tri_linear_combination(c21,P[n1],c22,P[n1+1],c23,P[n1+2],v2);
	tri_linear_combination(c31,P[n1],c32,P[n1+1],c33,P[n1+2],v3);
	swap(P[n1],v1);
	swap(P[n1+1],v2);
	swap(P[n1+2],v3);
      }
      // now columns operations on H (not on P)
      for (int j=0;j<n_orig;++j){
	vecteur & Hj=H[j];
	gen & Hjm1=Hj[n1];
	gen & Hjm2=Hj[n1+1];
	gen & Hjm3=Hj[n1+2];
	tmp1=Hjm1*c11+Hjm2*c21+Hjm3*c31;
	tmp2=Hjm1*c12+Hjm2*c22+Hjm3*c32;
	tmp3=Hjm1*c13+Hjm2*c23+Hjm3*c33;
	Hjm1=tmp1;
	Hjm2=tmp2;
	Hjm3=tmp3;
      }
      if (n1==100)
	dbg_schur(H,P);
      // cerr << H << endl;
      // chase the bulge: Hessenberg reduction on 2 subdiagonals
      hessenberg_ortho(H,P,n1,n2,compute_P,3,contextptr); // <- improve
    } // end for loop on niter
    return false;
  }

  // trn(P)*H*P=orig matrix
  // IMPROVE: real iterations
  void hessenberg_schur(std_matrix<gen> & H,std_matrix<gen> & P,int maxiter,double eps,GIAC_CONTEXT){
    int n_orig=H.size(),n=n_orig,nitershift0=0;
    bool real=true,is_double=true; // true does not work because the conjugate iteration might be bad conditionned
#ifdef VISUALC
    // is_double=false;
#endif
    for (int i=0;real && i<n;i++){
      for (int j=0;j<n;j++){
	if (is_double){
	  if (H[i][j].type==_DOUBLE_)
	    continue;
	  else
	    is_double=false;
	}
	if (!is_zero(im(H[i][j],contextptr))){
	  real=false;
	  break;
	}
      }
    }
    if (real && is_double){
      matrix_double H1,P1;
      std_matrix_gen2std_matrix_giac_double(H,H1);
      std_matrix_gen2std_matrix_giac_double(P,P1);
      if (francis_schur(H1,0,n_orig,P1,maxiter,std::sqrt(double(n_orig))*eps,false,false,true)){
	std_matrix_giac_double2std_matrix_gen(P1,P);
	std_matrix_giac_double2std_matrix_gen(H1,H);
	return;
      }
    }
    else {
      if (francis_schur(H,0,n_orig,P,maxiter,std::sqrt(double(n_orig))*eps,false,!real,true,contextptr)){
	return ;
      }
    }
    hessenberg_ortho(H,P,contextptr); // insure
    // make 0 below subdiagonal
    for (int i=2;i<n_orig;i++){
      for (int j=0;j<i-1;j++){
	H[i][j]=0;
      }
    }
    gen shift=0,ratio,oldratio=0,maxi,absmaxi,tmp,abstmp;
    vecteur SHIFT;
    for (int niter=0;n>1 && niter<maxiter;niter++){
      if (debug_infolevel>=2)
	cerr << "// qr iteration number " << niter << endl;
      shift=0;
      gen test=abs(H[n-1][n-2],contextptr);
      ratio=test/abs(H[n-1][n-1],contextptr);
      bool Small=is_strictly_greater(0.01,ratio,contextptr);
      if (Small)
	shift=H[n-1][n-1];
      else {
	if (n==2 || is_strictly_greater(0.01,(ratio=abs(H[n-2][n-3]/H[n-2][n-2],contextptr)),contextptr)){
	  // define shift according to the smallest eigenvalues 
	  // of the last 2x2 submatrix bloc
	  gen a=H[n-2][n-2],b=H[n-2][n-1],c=H[n-1][n-2],d=H[n-1][n-1];
	  gen delta=a*a-2*a*d+d*d+4*b*c;
	  if (real && n==2 && is_strictly_positive(-delta,0))
	    break;
	  delta=sqrt(delta,contextptr);
	  gen l1=(a+d+delta)/2,l2=(a+d-delta)/2;
	  if (is_strictly_greater(abs(l1,contextptr),abs(l2,contextptr),contextptr))
	    shift=l2;
	  else
	    shift=l1;
	}
	else {
	  if (niter>=maxiter/4)
	    shift=H[n-1][n-1]/2;
	}
      }
      oldratio=ratio;
      qr_rq(H,P,shift,n,nitershift0,contextptr);
      if (real && !is_zero(im(shift,contextptr)))
	SHIFT.push_back(shift);
      test=abs(H[n-1][n-2],contextptr);
      ratio=test/abs(H[n-1][n-1],contextptr);
      if (is_strictly_greater(gen(eps)/oldratio,1,contextptr) && is_strictly_greater(gen(eps)/ratio,1,contextptr)){
	// eigenvalue has been found
	niter=0;
	oldratio=0;
	n--;
	if (real && !SHIFT.empty()){
	  // int ni=SHIFT.size();
	  for(int i=0;i<(int)SHIFT.size();++i)
	    qr_rq(H,P,conj(SHIFT[i],contextptr),n,nitershift0,contextptr);
	  for (int i=0;i<n-1;i++){
	    vecteur & Pi=P[i];
	    maxi=Pi.front();
	    absmaxi=abs(maxi,contextptr);
	    for (int j=1;j<n-1;j++){
	      tmp=Pi[j];
	      abstmp=abs(tmp,contextptr);
	      if (abstmp>absmaxi){
		absmaxi=abstmp;
		maxi=tmp;
	      }
	    }
	    tmp=absmaxi/maxi;
	    multvecteur(tmp,Pi,Pi);
	    multvecteur(tmp,H[i],H[i]);
	    tmp=maxi/absmaxi;
	    for (int j=0;j<n_orig;j++){
	      gen & Hji= H[j][i];
	      Hji = tmp*Hji;
	    }
	  }
	  re(H,n-1,contextptr); re(P,n-1,contextptr);
	  n--;
	  SHIFT.clear();
	} // end if (real)
      } // end eigenvalue detected
    } // end loop on n for 0 subdiagonal elements
  }

  // Reduction to Hessenberg form, see e.g. Cohen algorithm 2.2.9
  // (with C array indices)
  // general case
  // if modulo==-1 Schur reduction up to precision eps and maxiterations maxiter
  // if modulo<-1, using orthogonal/unitary matrices
  bool mhessenberg(const matrice & M,matrice & h,matrice & p,int modulo,int maxiter,double eps,GIAC_CONTEXT){
    int n=M.size();
    if (!n || n!=mcols(M))
      return false; // setdimerr();
    bool modularize=!modulo && M[0][0].type==_MOD && (M[0][0]._MODptr+1)->type==_INT_;
    if (modularize)
      modulo=(M[0][0]._MODptr+1)->val;
    if (modulo>0){
      vector< vector<int> > H;
      if (!vecteur2vectvector_int(M,modulo,H))
	return false;
      vector< vector<int> > P;
      if (!vecteur2vectvector_int(midn(n),modulo,P))
	return false;
      mhessenberg(H,P,modulo);
      vectvector_int2vecteur(H,h);
      vectvector_int2vecteur(P,p);
      if (modularize){
	h=*makemod(h,modulo)._VECTptr;
	p=*makemod(p,modulo)._VECTptr;
      }
      return true;
    }
    std_matrix<gen> H,P;
    matrice2std_matrix_gen((modulo<0?*evalf(gen(M),1,contextptr)._VECTptr:M),H);
    p=midn(n);
    matrice2std_matrix_gen(p,P);
    if (modulo==-1)
      hessenberg_schur(H,P,maxiter,eps,contextptr);
    else {
      if (modulo<0)
	hessenberg_ortho(H,P,contextptr);
      else
	hessenberg(H,P,contextptr);
    }
    // store result
    std_matrix_gen2matrice(H,h);
    std_matrix_gen2matrice(P,p);
    return true;
  }
  gen _hessenberg(const gen & g0,GIAC_CONTEXT){
    if ( g0.type==_STRNG && g0.subtype==-1) return  g0;
    gen g(g0);
    int modulo=0;
    double eps=epsilon(contextptr);
    int maxiter=500;
    if (g.type==_VECT && g._VECTptr->size()>=2 && g.subtype==_SEQ__VECT){
      vecteur & v = *g._VECTptr;
      gen v1=v[1];
      if (v1.type==_INT_)
	modulo=v1.val;
      else {
	v1=evalf_double(v1,1,contextptr);
	if (v1.type==_DOUBLE_){
	  modulo=-1;
	  eps=v1._DOUBLE_val;
	  if (v.size()>2 && v[2].type==_INT_)
	    maxiter=v[2].val;
	}
      }
      g=v.front();
    }
    if (!is_squarematrix(g))
      return symbolic(at_hessenberg,g);
    matrice m(*g._VECTptr),h,p;
    if (!mhessenberg(m,h,p,modulo,maxiter,eps,contextptr))
      return gensizeerr(contextptr);
    if (modulo<0)
      return makevecteur(_trn(p,contextptr),h); // p,h such that p*h*p^-1=orig
    else
      return makevecteur(inv(p,contextptr),h); // p,h such that p*h*p^-1=orig
  }
  static const char _hessenberg_s []="hessenberg";
  static define_unary_function_eval (__hessenberg,&giac::_hessenberg,_hessenberg_s);
  define_unary_function_ptr5( at_hessenberg ,alias_at_hessenberg,&__hessenberg,0,true);

  dense_POLY1 mpcar_hessenberg(const matrice & A,int modulo,GIAC_CONTEXT){
    int n=A.size();
    if (modulo || is_integer_matrice(A)){
      if (modulo){ // try Krylov pmin
	vector< vector<int> > N,temp(n+1),ttemp;
	if (debug_infolevel)
	  cerr << "Charpoly mod " << modulo << " A*v" << clock() << endl;
	if (!vecteur2vectvector_int(A,modulo,N))
	  return vecteur(1,gendimerr(contextptr));
	vector<int> & t0=temp[0];
	t0.reserve(n);
	for (int i=0;i<n;++i)
	  t0.push_back(rand()%modulo);
	for (int j=0;j<n;++j){
	  if (!multvectvector_int_vector_int(N,temp[j],modulo,temp[j+1]))
	    return vecteur(1,gendimerr(contextptr));
	}
	if (debug_infolevel)
	  cerr << "Charpoly mod " << modulo << " tran " << clock() << endl;
	tran_vect_vector_int(temp,ttemp);
	vecteur pivots;
	longlong det;
	vector<int> permutation,maxrankcol;
	if (debug_infolevel)
	  cerr << "Charpoly mod " << modulo << " rref " << clock() << endl;
	smallmodrref(ttemp,pivots,permutation,maxrankcol,det,0,n,0,n+1,false/* LU decomp */,0,modulo,2/* LU */);
	if (debug_infolevel)
	  cerr << "Charpoly mod " << modulo << " det=" << det << " " << clock() << endl;
	// if det==0 we will use Hessenberg
	// If rank==n-1 we could extract the min polynomial and find charpoly using the trace
	if (
	    // false 
	    det
	    ){ 
	  // U*charpol=last column
	  for (int i=n-1;i>=0;--i){
	    // charpol[i]=LU[i,i]^(-1)*(bp[i]-sum(j>i)LU[i,j]*charpol[j])
	    int res=0;
	    vector<int> & li=ttemp[i];
	    for (int j=i+1;j<n;++j)
	      mod(res,li[j],ttemp[j][n],modulo);
	    li[n]=(invmod(li[i],modulo)*longlong(li[n]-res))%modulo;
	  }
	  // the last column is the min poly
	  modpoly charpol(n+1);
	  for (int i=0;i<n;++i)
	    charpol[n-i]=smod(-ttemp[i][n],modulo);
	  charpol[0]=1;
	  return charpol;
	}
	else
	  if (debug_infolevel)
	    cerr << "Singular, back to Hessenberg " << endl;
      }
      else {
	gen B=evalf_double(linfnorm(A,contextptr),0,contextptr);
	double Bd=B._DOUBLE_val;
	if (!Bd){
	  modpoly charpol(n+1);
	  charpol[0]=1;
	  return charpol;
	}
	// max value of any coeff in the charpoly
	// max eigenval is <= sqrt(n)||A|| hence bound is in n (log(B)+log(n)/2)
	// we must add combinatorial (n k)<2^n
	double logbound=n*(std::log10(double(n))/2+std::log10(Bd)+std::log10(2.0));
	double proba=proba_epsilon(contextptr),currentprob=1;
	gen currentp(init_modulo(n,logbound));
	gen pip(currentp);
	double pipd=std::log10(pip.val/2+1.0);
	modpoly charpol=*makemod(mpcar_hessenberg(A,currentp.val,contextptr),0)._VECTptr;
	if (is_undef(charpol)) return charpol;
	for (;pipd<logbound && currentprob>proba;){
	  currentp=nextprime(currentp.val+2);
	  modpoly currentcharpol=*makemod(mpcar_hessenberg(A,currentp.val,contextptr),0)._VECTptr;
	  if (is_undef(currentcharpol)) return currentcharpol;
	  modpoly newcharpol=ichinrem(charpol,currentcharpol,pip,currentp);
	  if (newcharpol==charpol)
	    currentprob=currentprob/currentp.val;
	  else {
	    charpol=newcharpol;
	    currentprob=1.0;
	  }
	  pip=pip*currentp;
	  pipd += std::log10(double(currentp.val));
	}
	if (debug_infolevel && pipd<logbound)
	  cerr << "Probabilistic answer" << endl;
	return charpol;
      }
    } // end if (is_integer_matrix)
    matrice H,P;
    if (!mhessenberg(A,H,P,modulo,500,1e-10,contextptr))
      return vecteur(1,gensizeerr(contextptr));
    if (modulo)
      H=*makemod(H,modulo)._VECTptr;
    dense_POLY1 p0(1,plus_one),pX(2,plus_one);
    vector< dense_POLY1 > p(1,p0);
    for (int m=1;m<=n;++m){
      pX[1]=-H[m-1][m-1];
      p0=pX*p0;
      gen t(plus_one);
      for (int i=1;i<m;++i){
	t=t*H[m-i][m-i-1];
	p0=p0-t*H[m-i-1][m-1]*p[m-i-1];
      }
      p.push_back(p0);
    }
    return p0;
  }
  gen _pcar_hessenberg(const gen & g,GIAC_CONTEXT){
    if ( g.type==_STRNG && g.subtype==-1) return  g;
    if (!is_squarematrix(g)){
      if (g.type==_VECT && g._VECTptr->size()==2){
	gen m=g._VECTptr->front(),x=g._VECTptr->back();
	if (is_squarematrix(m))
	  return symb_horner(mpcar_hessenberg(*m._VECTptr,0,contextptr),x);
      }
      return symbolic(at_pcar_hessenberg,g);
    }
    matrice m(*g._VECTptr);
    return mpcar_hessenberg(m,0,contextptr);
  }
  static const char _pcar_hessenberg_s []="pcar_hessenberg";
  static define_unary_function_eval (__pcar_hessenberg,&giac::_pcar_hessenberg,_pcar_hessenberg_s);
  define_unary_function_ptr5( at_pcar_hessenberg ,alias_at_pcar_hessenberg,&__pcar_hessenberg,0,true);


  // Fadeev algorithm to compute the char poly of a matrix
  // B is a vector of matrices
  // the returned value is the vector of coeff of the char poly
  // see modpoly.h for polynomial operations on vecteur
  dense_POLY1 mpcar(const matrice & a,vecteur & Bv,bool compute_Bv,bool convert_internal,GIAC_CONTEXT){
    int n=a.size();
    if (n && a[0]._VECTptr->front().type==_MOD){
      vecteur P(mpcar_hessenberg(a,0,contextptr));
      // do Horner to compute Bv
      if (compute_Bv){
	horner(P,a,0,Bv);
	Bv[0]=midn(n);
      }
      return P;
    }
    matrice A,Bi,Ai,I,lv;
    if (convert_internal){
      // convert a to internal form
      lv=alg_lvar(a);
      A = *(e2r(a,lv,contextptr)._VECTptr);
    }
    else
      A=a;
    midn(n,I);
    Bi=I; // B0=Id
    Bv.push_back(Bi); 
    vecteur P;
    gen pk;
    P.push_back(1); // p0= 1
    for (int i=1;i<=n;++i){
      mmult(A,Bi,Ai); // Ai = A*Bi
      pk = rdiv(-mtrace(Ai),i); 
      P.push_back(convert_internal?r2e(pk,lv,contextptr):pk);
      addvecteur( Ai,multvecteur(pk,I),Bi); // Bi = Ai+pk*I
      // cout << i << ":" << Bi << endl;
      if (i!=n)
	Bv.push_back(convert_internal?r2e(Bi,lv,contextptr):Bi);
    }
    return P;
  }

  dense_POLY1 mpcar(const matrice & a,vecteur & Bv,bool compute_Bv,GIAC_CONTEXT){
    return mpcar(a,Bv,compute_Bv,false,contextptr);
  }

  gen _lagrange(const gen & g,GIAC_CONTEXT);

  static gen pcar_interp(const matrice & a,gen & g,GIAC_CONTEXT){
    int m=a.size();
    vecteur v1,v2,I(midn(m));
    for (int j=0;j<=m;++j){
      v1.push_back(j);
      v2.push_back(mdet(addvecteur(a,multvecteur(-j,I)),contextptr));
    }
    return _lagrange(makevecteur(v1,v2,g),contextptr);
  }

  gen _pcar(const gen & a,GIAC_CONTEXT){
    if ( a.type==_STRNG && a.subtype==-1) return  a;
    vecteur Bv;
    matrice M;
    gen b(undef);
    if (!is_squarematrix(a)){
      if (a.type!=_VECT)
	return symb_pcar(a);
      vecteur v=*a._VECTptr;
      int s=v.size();
      if (s<2 || !is_squarematrix(v.front()))
	return gensizeerr(contextptr);
      matrice &m=*v.front()._VECTptr;
      if (v.back().type==_INT_ && v.back().val==_FADEEV){
	vecteur res=mpcar(m,Bv,false,true,contextptr);
	return s==2?res:symb_horner(res,v[1]);
      }
      if (v.back()==at_pmin && probabilistic_pmin(m,Bv,false,contextptr))
	return s==2?Bv:symb_horner(Bv,v[1]); 
      if (v.back()==at_lagrange)
	return pcar_interp(m,s==2?vx_var:v[1],contextptr);
      if (v.back()==at_hessenberg || v.back()==at_pcar_hessenberg){
	Bv=mpcar_hessenberg(m,0,contextptr);
	return s==2?Bv:symb_horner(Bv,v[1]);
      }
      b=v[1];
      M=m;
    }
    else
      M=*a._VECTptr;
    // int n=M.size();
    // search for the best algorithm
    if (is_integer_matrice(M)){
      vecteur res=mpcar_hessenberg(M,0,contextptr);
      if (is_undef(b))
	return gen(res,_POLY1__VECT);
      return symb_horner(res,b);	
    }
    else {
      vecteur res;
      res=mpcar(M,Bv,false,true,contextptr);
      if (is_undef(b))
	return res;
      return symb_horner(res,b);
    }
  }
  static const char _pcar_s []="pcar";
  static define_unary_function_eval (__pcar,&giac::_pcar,_pcar_s);
  define_unary_function_ptr5( at_pcar ,alias_at_pcar,&__pcar,0,true);

  static vecteur polymat2matpoly(const vecteur & v){
    if (v.empty() || v.front().type!=_VECT)
      return vecteur(1,gensizeerr(gettext("polymat2matpoly")));
    int l,c,s=v.size();
    mdims(*v.front()._VECTptr,l,c);
    vecteur mat;
    mat.reserve(l);
    for (int i=0;i<l;++i){
      vecteur ligne;
      ligne.reserve(c);
      for (int j=0;j<c;++j){
	bool trim=true;
	vecteur res;
	res.reserve(s);
	for (int k=0;k<s;++k){
	  gen g=v[k];
	  gen tmp=g[i][j];
	  if (trim && is_zero(tmp,context0))
	    continue;
	  trim=false;
	  res.push_back(tmp);
	}
	ligne.push_back(gen(res,_POLY__VECT));
      }
      mat.push_back(ligne);
    }
    return mat;
  }

  static vecteur polymat2mat(const vecteur & v){
    if (v.empty()) 
      return v;
    if (v.front().type!=_VECT)
      return vecteur(1,gensizeerr(gettext("polymat2mat")));
    int l,c,s=v.size();
    vecteur w(v);
    for (int i=0;i<s;++i)
      w[i]=mtran(*v[i]._VECTptr);
    mdims(*v.front()._VECTptr,l,c);
    vecteur mat;
    mat.reserve(l*s);
    for (int k=0;k<s;++k){
      gen & g=w[k];
      for (int i=0;i<l;++i){
	mat.push_back(g[i]);
      }
    }
    return mat;
  }

  // dot product of a[0..a.size()-1] and b[pos..pos+a.size()-1]
  gen generalized_dotvecteur(const vecteur & a,const vecteur & b,int pos){
    vecteur::const_iterator ita=a.begin(), itaend=a.end();
    vecteur::const_iterator itb=b.begin()+pos;
    gen res;
    for (;(ita!=itaend);++ita,++itb){
      res = res + (*ita)*(*itb);
    }
    return res;
  }
  
  vecteur generalized_multmatvecteur(const matrice & a,const vecteur & b){
    vecteur::const_iterator ita=a.begin(), itaend=a.end();
    int s=b.size();
    int n=itaend-ita; // number of vectors stored in b=s/n
    vecteur res;
    res.reserve(s);
    for (int i=0;i<s;i+=n){
      for (ita=a.begin();ita!=itaend;++ita){
	res.push_back(generalized_dotvecteur(*(ita->_VECTptr),b,i));
      }
    }
    return res;
  }

  // [almost] rational jordan block
  matrice rat_jordan_block(const vecteur & v,int n,bool pseudo){
    if (n<1)
      return vecteur(1,gendimerr(gettext("rat_jordan_block")));
    int s=v.size()-1;
    // Size of the matrix is s*n
    vecteur ligne(s*n,zero);
    std_matrix<gen> M(s*n,ligne);
    for (int i=0;i<n;++i){
      // Fill the block-diagonal part with companion block
      for (int j=0;j<s;++j){
	M[i*s+j][i*s+s-1]=-v[s-j];
	if (j>0)
	  M[i*s+j][i*s+j-1]=plus_one;
      }
      // Fill the upper diagonal with idn or a single 1
      if (i!=n-1){
	if (pseudo)
	  M[i*s][i*s+s+s-1]=1;
	else {
	  for (int j=0;j<s;++j){
	    M[i*s+j][i*s+s+j]=1;
	  }
	}
      }
    }
    matrice res;
    std_matrix_gen2matrice(M,res);
    return res;
  }

  gen _rat_jordan_block(const gen &args,GIAC_CONTEXT){
    if (args.type==_STRNG && args.subtype==-1) return args;
    if (args.type==_VECT && args._VECTptr->size()==3){
      vecteur & v=*args._VECTptr;
      gen Px=_e2r(makevecteur(v[0],v[1]),contextptr);
      if (Px.type==_VECT && v[2].type==_INT_){
	int n=v[2].val;
	return rat_jordan_block(*Px._VECTptr,absint(n),n<0);
      }
    }
    return gensizeerr(contextptr);
  }
  static const char _rat_jordan_block_s []="rat_jordan_block";
  static define_unary_function_eval (__rat_jordan_block,&giac::_rat_jordan_block,_rat_jordan_block_s);
  define_unary_function_ptr5( at_rat_jordan_block ,alias_at_rat_jordan_block,&__rat_jordan_block,0,true);

  matrice pseudo_rat_to_rat(const vecteur & v,int n){
    if (n<1)
      return vecteur(1,gendimerr(gettext("pseudo_rat_ro_rat")));
    matrice A(rat_jordan_block(v,n,true));
    if (is_undef(A)) return A;
    // lines of A are initial v
    vecteur q(v);
    int d=q.size()-1; // degree of the polynomial
    matrice res(midn(n*d));
    reverse(q.begin(),q.end());
    for (int j=1;j<n;++j){
      // compute Q(A) v_{j,0}
      vecteur QAvj0(n*d);
      for (int l=1;l<=d;++l){
	int mmax=min(l,j);
	for (int m=1;m<=mmax;++m){
	  QAvj0=addvecteur(QAvj0,multvecteur(q[l]*comb((unsigned long) l,(unsigned long)m),*res[(j-m)*d+(l-m)]._VECTptr));
	}
      }
      // shift
      vecteur vj0=mergevecteur(vecteur(d),vecteur(QAvj0.begin(),QAvj0.begin()+(n-1)*d));
      // replace in res
      res[j*d]=vj0;
      // compute images by A, ..., A^[d-1]
      for (int l=1;l<d;++l){
	vj0=multmatvecteur(A,vj0);
	vecteur tmp(vj0);
	int mmax=min(l,j);
	for (int m=1;m<=mmax;++m)
	  tmp=subvecteur(tmp,multvecteur(comb((unsigned long) l,(unsigned long) m),*res[(j-m)*d+(l-m)]._VECTptr));
	res[j*d+l]=tmp;
      }
    }
    return res;
  }

  // input trn(p)*d*p=original matrix, d upper triangular
  // output p*d*inv(p)=original matrix, d diagonal
  bool schur_eigenvectors(matrice &p,matrice & d,double eps,GIAC_CONTEXT){
    bool ans=true;
    int dim=p.size();
    matrice m(midn(dim)); 
    // columns of m are the vector of the basis of the Schur decomposition
    // in terms of the eigenvector
    for (int k=1;k<dim;++k){
      // compute column k of m
      for (int j=0;j<k;++j){
	gen tmp=0;
	for (int i=0;i<k;++i){
	  tmp += d[i][k]*m[j][i];
	}
	tmp = tmp/(d[j][j]-d[k][k]);
	(*m[j]._VECTptr)[k]=tmp;
      }
    }
    m=minv(m,contextptr);
    p=mmult(*_trn(p,contextptr)._VECTptr,m);
    // set d to its diagonal
    for (int i=0;i<dim;++i){
      vecteur & di= *d[i]._VECTptr;
      for (int j=0;j<dim;++j){
	if (j==i) continue;
#ifndef GIAC_HAS_STO_38
	if (ans && j==i-1 && is_greater(abs(di[j]/di[j+1],contextptr),eps,contextptr)){
	  *logptr(contextptr) << "Schur row " << j << " " << d[i] << endl;
	  ans=false;
	}
#endif
	di[j]=0;
      }
    }
    return ans;
  }

  void matrice_double2lapack(const matrix_double & H,double * A){
    matrix_double::const_iterator it=H.begin(),itend=H.end();
    int rows=itend-it;
    for (int i = 0; it!=itend; ++i,++it){
      const vector<giac_double> & v =*it;
      vector<giac_double>::const_iterator jt=v.begin(),jtend=v.end();
      for (int j = 0; jt!=jtend;++j, ++jt){
	A[i + j * rows] = *jt;
      }
    }
  }

  void lapack2matrice_double(double * A, int rows,int cols,matrix_double & R){
    R.resize(rows);
    for (int i=0;i<rows;++i){
      vector<double> & r=R[i];
      r.resize(cols);
      for (int j=0;j<cols;++j)
	r[j] = A[i + j * rows];
    }
  }

  void transpose_square_matrix(matrix_double & R){
    int n=R.size();
    for (int i=0;i<n;++i){
      for (int j=0;j<i;++j)
	swapdouble(R[i][j],R[j][i]);
    }
  }

#ifdef HAVE_LIBLAPACK
  bool lapack_schur(matrix_double & H,matrix_double & P,bool compute_P,vecteur & eigenvalues){
    if (!CAN_USE_LAPACK)
      return false;
    /* int dgees_(char *jobvs, char *sort, L_fp select, integer *n, 
       doublereal *a, integer *lda, integer *sdim, doublereal *wr, 
       doublereal *wi, doublereal *vs, integer *ldvs, doublereal *work, 
       integer *lwork, logical *bwork, integer *info)
       jobvs="n" or "v" (compute Schur vectors)
       sort="n" or "s" (sort eigenvals), 
       select(real,imag) should return true if eigenval is selected
       n order of matrix A (==size(A))
       a input/output matrix (Schur form on output)
       lda leading dimension of A >= max(1,n)
       sdim (output) =0 if sort=="n" or number of eigenval selected if =="s" 
       wr and wi (output) contains real and imaginary part of eigenvals
       vs (output) if jobs=="v: contains the orthogonal matrix
       ldvs leading dimension of vs (>=n if jobvs=="v")
       lwork >=3*n, dimension of the work array (should be larger for good results)
       bwork array of dimension n
       info ==0 success, <0 bad arg value, >0 runtime error
    */
    integer n=H.size(),sdim,ldvs=n,lwork=max(20,n)*n,info;
    doublef2c_real * Hlapack=new doublef2c_real[n*n], * Plapack=new doublef2c_real[n*n];
    doublef2c_real * Wreal=new doublef2c_real[n], * Wimag=new doublef2c_real[n], * work=new doublef2c_real[lwork];
    logical * bwork=new logical[n];
    char ch[2]={0,0};
    ch[0]=(compute_P?'v':'n');
    matrice_double2lapack(H,Hlapack);
    char ch2[]="n";
    dgees_(ch,ch2,0,&n,
	   Hlapack,&n,&sdim,Wreal,
	   Wimag,Plapack,&ldvs,work,&lwork,
	   bwork,&info);
    lapack2matrice_double(Hlapack,n,n,H);
    if (compute_P){
      lapack2matrice_double(Plapack,n,n,P);
      transpose_square_matrix(P);
    }
    delete [] Hlapack; delete [] Plapack;
    eigenvalues.resize(n);
    for (int i=0;i<n;++i)
      eigenvalues[i]=gen(Wreal[i],Wimag[i]);
    delete [] Wreal; delete [] Wimag; delete [] work;
    delete [] bwork;
    return info==0;
  }
  
  bool lapack_schur(std_matrix<gen> & H,std_matrix<gen> & P,bool compute_P,vecteur & eigenvalues,GIAC_CONTEXT){
    if (!CAN_USE_LAPACK)
      return false;
  /* int zgees_(char *jobvs, char *sort, L_fp select, integer *n, 
     doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w, 
     doublecomplex *vs, integer *ldvs, doublecomplex *work, integer *lwork, 
     doublereal *rwork, logical *bwork, integer *info) 
       jobvs="n" or "v" (compute Schur vectors)
       sort="n" or "s" (sort eigenvals), 
       select(real,imag) should return true if eigenval is selected
       n order of matrix A (==size(A))
       a input/output matrix (Schur form on output)
       lda leading dimension of A >= max(1,n)
       sdim (output) =0 if sort=="n" or number of eigenval selected if =="s" 
       w (output) contains eigenvals
       vs (output) if jobs=="v: contains the orthogonal matrix
       ldvs leading dimension of vs (>=n if jobvs=="v")
       lwork >=3*n, dimension of the work array (should be larger for good results)
       bwork array of dimension n
       info ==0 success, <0 bad arg value, >0 runtime error
    */
    integer n=H.size(),sdim,ldvs=n,lwork=max(20,n)*n,info;
    doublef2c_complex * Hlapack=new doublef2c_complex[n*n], * Plapack=new doublef2c_complex[n*n];
    doublef2c_complex * W=new doublef2c_complex[n], * work=new doublef2c_complex[lwork];
    doublef2c_real * rwork=new doublef2c_real[lwork];
    logical * bwork=new logical[n];
    matrix2zlapack(H,Hlapack,contextptr);
    char ch[2]={0,0};
    ch[0]=(compute_P?'v':'n');
    char ch2[]="n";
    zgees_(ch,ch2,0,&n,
	   Hlapack,&n,&sdim,W,
	   Plapack,&ldvs,work,&lwork,
	   rwork,bwork,&info);
    zlapack2matrix(Hlapack,n,n,H);
    if (compute_P){
      zlapack2matrix(Plapack,n,n,P);
      P=P.transconjugate();
    }
    delete [] Hlapack; delete [] Plapack;
    eigenvalues.resize(n);
    for (int i=0;i<n;++i)
      eigenvalues[i]=gen(W[i].r,W[i].i);
    delete [] W; delete [] work; delete [] rwork;
    delete [] bwork;
    return info==0;
  }
#else

  bool lapack_schur(matrix_double & H,matrix_double & P,bool compute_P,vecteur & eigenvalues){
    return false;
  }

  bool lapack_schur(std_matrix<gen> & H,std_matrix<gen> & P,bool compute_P,vecteur & eigenvalues,GIAC_CONTEXT){
    return false;
  }
#endif // HAVE_LIBLAPACK

  // if jordan is false, errors for non diagonalizable matrices
  // if jordan is true, d is a matrix, not a vector
  bool egv(const matrice & m,matrice & p,vecteur & d, GIAC_CONTEXT,bool jordan,bool rational_jordan_form,bool eigenvalues_only){
    bool numeric_matrix=is_fully_numeric(m);
    bool sym=(m==mtran(*conj(m,contextptr)._VECTptr));
    double eps=epsilon(contextptr);
    // check for symmetric numeric matrix
    if (numeric_matrix){
#ifdef HAVE_LIBLAPACK
      if ( !is_zero(im(m,contextptr),contextptr) ){ 
	// complex matrix, try hermitian
	if (m==conj(mtran(m),contextptr)){
	  // call to zheev
	  // ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,INFO)
	  char JOBZ=eigenvalues_only?'N':'V';
	  char UPLO='U';
	  integer N=m.size();
	  doublef2c_complex * A=new doublef2c_complex[N*N];
	  matrice2zlapack(m,A,contextptr);
	  integer LDA=N;
	  double * W = new double[N];
	  integer LWORK=max(10,N)*N,INFO;
	  doublef2c_complex * WORK=new doublef2c_complex[LWORK];
	  double * RWORK=new double[LWORK];
	  zheev_( &JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, RWORK,&INFO );
	  delete [] WORK; delete [] RWORK;
	  if (INFO){ 
	    delete [] A; delete [] W;  return false;
	  }
	  d.resize(N);
	  for (int i=0;i<N;++i){
	    if (jordan){
	      vecteur v(N);
	      v[i]=W[i];
	      d[i]=v;
	    }
	    else
	      d[i]=W[i];
	  }
	  zlapack2matrice(A,N,N,p);
	  delete [] A; delete [] W;
	  return true;
	}
	// call to zgeev
	// ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,WORK, LWORK, RWORK, INFO )
	char JOBVL=eigenvalues_only?'N':'V';
	char JOBVR=eigenvalues_only?'N':'V';
	integer N=m.size();
	doublef2c_complex * A=new doublef2c_complex[N*N];
	matrice2zlapack(m,A,contextptr);
	integer LDA=N;
	doublef2c_complex * W=new doublef2c_complex[N];
	doublef2c_complex * VL=new doublef2c_complex[N*N];
	doublef2c_complex * VR=new doublef2c_complex[N*N];
	integer LDVL=N,LDVR=N,LWORK=max(10,N)*N,INFO;
	double * RWORK=new double[LWORK];
	doublef2c_complex * WORK=new doublef2c_complex[LWORK];
	zgeev_(&JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR,&LDVR, WORK, &LWORK, RWORK,&INFO );
	delete [] WORK; delete [] RWORK; delete [] A; delete [] VL; 
	if (INFO){
	  delete [] W; delete [] VR; return false;
	}
	zlapack2matrice(VR,N,N,p);
	d.resize(N);
	for (int i=0;i<N;++i){
	  if (jordan){
	    vecteur v(N);
	    v[i]=gen(W[i].r,W[i].i);
	    d[i]=v;
	  }
	  else
	    d[i]=gen(W[i].r,W[i].i);
	}
	delete [] W; delete [] VR; 
	return true;
      }
      if (sym){ // call to dsyev
	// DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
	char JOBZ=eigenvalues_only?'N':'V';
	char UPLO='U';
	integer N=m.size();
	double * A=new double[N*N];
	matrice2lapack(m,A,contextptr);
	integer LDA=N;
	double * W = new double[N];
	integer LWORK=max(10,N)*N,INFO;
	double * WORK=new double[LWORK];
	dsyev_( &JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO );
	delete [] WORK;
	if (INFO){ 
	  delete [] A; delete [] W;  return false;
	}
	d.resize(N);
	for (int i=0;i<N;++i){
	  if (jordan){
	    vecteur v(N);
	    v[i]=W[i];
	    d[i]=v;
	  }
	  else
	    d[i]=W[i];
	}
	lapack2matrice(A,N,N,p);
	delete [] A; delete [] W;
	return true;
      }
      // call to dgeev
      // DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,LDVR, WORK, LWORK, INFO )
      char JOBVL=eigenvalues_only?'N':'V';
      char JOBVR=eigenvalues_only?'N':'V';
      integer N=m.size();
      double * A=new double[N*N];
      matrice2lapack(m,A,contextptr);
      integer LDA=N;
      double * WR=new double[N];
      double * WI=new double[N];
      double * VL=new double[N*N];
      double * VR=new double[N*N];
      integer LDVL=N,LDVR=N,LWORK=max(10,N)*N,INFO=0;
      double * WORK=new double[LWORK];
      dgeev_(&JOBVL, &JOBVR, &N, A, &LDA, WR, WI, VL, &LDVL, VR,&LDVR, WORK, &LWORK, &INFO );
      delete [] WORK; delete [] A; delete [] VL;
      if (INFO){
	delete [] WR; delete [] WI; delete [] VR; return false;
      }
      lapack2matrice(VR,N,N,p);
      p=mtran(p);
      d.resize(N);
      for (int i=0;i<N;++i){
	if (jordan){
	  vecteur v(N);
	  v[i]=gen(WR[i],WI[i]);
	  d[i]=v;
	  if (WI[i]!=0 && i!=N-1){
	    v[i]=0;
	    v[i+1]=gen(WR[i],-WI[i]);
	    d[i+1]=v;
	  }
	}
	else {
	  d[i]=gen(WR[i],WI[i]);
	  if (WI[i]!=0 && i!=N-1)
	    d[i+1]=gen(WR[i],-WI[i]);
	}
	if (WI[i]==0 || i==N-1)
	  continue;
	gen tmp=p[i]+cst_i*p[i+1];
	p[i]=tmp;
	p[i+1]=conj(tmp,contextptr);
	++i;
      }
      p=mtran(p);
      delete [] WR; delete [] WI; delete [] VR;
      return true;
#endif // HAVE_LIBLAPACK
#ifdef HAVE_LIBGSL
      if (sym){
	gsl_matrix * a=matrice2gsl_matrix(m,contextptr);
	int s=a->size1;
	gsl_matrix * eigenvectors= gsl_matrix_alloc(s,s);
	gsl_vector * eigenvalues =gsl_vector_alloc(s);
	gsl_eigen_symmv_workspace * w=gsl_eigen_symmv_alloc(s);
	gsl_eigen_symmv (a, eigenvalues,eigenvectors,w);
	gsl_eigen_symmv_free(w);
	p=gsl_matrix2matrice(eigenvectors);
	d=gsl_vector2vecteur(eigenvalues);
	if (jordan){
	  for (int i=0;i<s;++i){
	    vecteur tmp(s);
	    tmp[i]=d[i];
	    d[i]=tmp;
	  }
	}
	gsl_matrix_free(eigenvectors);
	gsl_vector_free(eigenvalues);
	return true;
      } // end sym. matrix
#endif // HAVE_LIBGSL
      std_matrix<gen> H,P;
      matrice2std_matrix_gen(m,H);
      int dim(H.size());
      matrix_double H1,P1;
      if (!eigenvalues_only){
	matrice p(midn(dim));
	matrice2std_matrix_gen(p,P);
	std_matrix_gen2std_matrix_giac_double(P,P1);
      }
      std_matrix_gen2std_matrix_giac_double(H,H1);
      bool ans=francis_schur(H1,0,dim,P1,SOLVER_MAX_ITERATE,dim*eps,false,false,!eigenvalues_only);
      if (eigenvalues_only){
	vecteur res;
	ans = ans && schur_eigenvalues(H1,res,eps,contextptr);
	gen tmp=_diag(res,contextptr);
	if (tmp.type!=_VECT)
	  return false;
	d=*tmp._VECTptr;
	return ans;
      }
      else {
	std_matrix_giac_double2std_matrix_gen(H1,H);
	std_matrix_giac_double2std_matrix_gen(P1,P);
	// finish Schur with complex entries
	ans=francis_schur(H,0,dim,P,SOLVER_MAX_ITERATE,dim*eps,true,true,true,contextptr);
	std_matrix_gen2matrice(P,p);
	std_matrix_gen2matrice(H,d);
	return ans && schur_eigenvectors(p,d,eps,contextptr);
      }
    } // end if (numeric_matrix)
    int taille=m.size();
    vecteur lv(lvar(m));
    numeric_matrix=has_num_coeff(m) && is_fully_numeric(evalf(m,1,contextptr));
    matrice mr=*(e2r(numeric_matrix?exact(m,contextptr):m,lv,contextptr)._VECTptr); // convert to internal form
    // vecteur lv;
    // matrice mr = m;
    matrice m_adj;
    vecteur p_car;
    p_car=mpcar(mr,m_adj,true,contextptr);
    p_car=common_deno(p_car)*p_car; // remove denominators
    // factorizes p_car
    factorization f;
    polynome p_content(lv.size()+1);
    if (!factor(poly12polynome(p_car,1),p_content,f,false,rational_jordan_form?false:withsqrt(contextptr),complex_mode(contextptr),1))
      return false;
    factorization::const_iterator f_it=f.begin(),f_itend=f.end();
    int total_char_found=0;
    for (;f_it!=f_itend;++f_it){
      // find roots of it->fact
      // works currently only for 1st order factors
      // vecteur v=solve(f_it->fact);
      vecteur v;
      vecteur w=polynome2poly1(f_it->fact,1);
      int s=w.size();
      if (s<2)
	continue;
      if (s==2)
	v.push_back(rdiv(-w.back(),w.front()));
      if (is_undef(v))
	return false;
      gen x;
      vecteur cur_m_adj(m_adj),cur_lv(lv),new_m_adj,char_m;
      if (s>=3 && rational_jordan_form){
	int mult=f_it->mult;
	int qdeg=s-1;
	int n=mult*qdeg; // number of vectors to find
	// Divide cur_m_adj by w f_it->mult times
	// Collect the remainders matrices in C
	vecteur C,quo,rem;
	int char_line=0,char_found=0,cycle_size=mult; 
	for (int i=0;i<mult;++i){
	  DivRem(cur_m_adj,w,0,quo,rem);
	  // rem is a polynomial made of matrices
	  // we convert it to a matrix (explode the polys)
	  if (rem.empty()){
	    --cycle_size;
	  }
	  else {
	    C=mergematrice(C,polymat2mat(rem));
	    if (is_undef(C)) return false;
	  }
	  cur_m_adj=quo;
	}
	// char_line is the line where the reduction begins
	vecteur Ccopy(C),pivots;
	gen det;
	for (;char_found<n;){
	  // Reduce
	  if (!mrref(Ccopy,C,pivots,det,0,Ccopy.size(),0,taille,
		true,char_line,true,1,0,
		     contextptr))
	    return false;
	  // Extract a non-0 line at char_line
	  vecteur line=*C[char_line]._VECTptr;
	  if (is_zero(vecteur(line.begin(),line.begin()+taille),contextptr)){
	    // Keep lines 0 to char_line-1, remove last taille columns
	    Ccopy=mtran(vecteur(C.begin(),C.begin()+char_line));
	    if (signed(Ccopy.size())<taille)
	      return false; // setdimerr();
	    vecteur debut(Ccopy.begin(),Ccopy.end()-taille);
	    debut=mtran(debut);
	    // Cut first taille columns of the remainder of the matrix
	    Ccopy=mtran(vecteur(C.begin()+char_line,C.end()));
	    if (signed(Ccopy.size())<taille)
	      return false; // setdimerr();
	    vecteur fin(Ccopy.begin()+taille,Ccopy.end());
	    fin=mtran(fin);
	    Ccopy=mergevecteur(debut,fin);
	    --cycle_size;
	    continue;
	  }
	  Ccopy=vecteur(C.begin(),C.begin()+char_line);
	  // make a bloc with line and A, A^2, ..., A^[qdeg-1]*line
	  // and put them into Ccopy and in ptmp
	  vecteur ptmp;
	  for (int i=0;i<qdeg;++i){
	    Ccopy.push_back(line);
	    ptmp.push_back(line);
	    line=generalized_multmatvecteur(mr,line);
	  }
	  // finish Ccopy by copying the remaining lines of C
	  const_iterateur ittmp=C.begin()+char_line+1,ittmpend=C.end();
	  for (;ittmp!=ittmpend;++ittmp)
	    Ccopy.push_back(*ittmp);
	  // update d (with a ratjord bloc) 
	  int taille_bloc=qdeg*cycle_size;
	  matrice tmp=mtran(rat_jordan_block(w,cycle_size,false));
	  tmp=mergematrice(vecteur(qdeg*cycle_size,vecteur(total_char_found)),tmp);
	  tmp=mergematrice(tmp,vecteur(qdeg*cycle_size,vecteur(taille-total_char_found-taille_bloc)));
	  if (is_undef(tmp)) return false;
	  d=mergevecteur(d,tmp);
	  // update p with ptmp 
	  matrice padd;
	  for (int j=0;j<cycle_size;++j){
	    for (int i=0;i<qdeg;++i){
	      vecteur & ptmpi=*ptmp[i]._VECTptr;
	      padd.push_back(vecteur(ptmpi.begin()+taille*j,ptmpi.begin()+taille*(j+1)));
	    }  
	  }
	  matrice AA(pseudo_rat_to_rat(w,cycle_size));
	  if (is_undef(AA)) return false;
	  padd=mmult(AA,padd);
	  p=mergevecteur(p,padd);
	  char_found += taille_bloc;
	  total_char_found += taille_bloc;
	  char_line += cycle_size;
	}
	continue;
      } // end if s>=3 and rational_jordan_form
      if (s>=3){ // recompute cur_m_adj using new extensions
	cur_m_adj=*r2sym(m_adj,lv,contextptr)._VECTptr;
	identificateur tmpx(" x");
	v=solve(horner(r2sym(w,lv,contextptr),tmpx),tmpx,complex_mode(contextptr),contextptr); 
	// compute new lv and update v and m_adj accordingly
	cur_lv=alg_lvar(v);
	alg_lvar(cur_m_adj,cur_lv);
	cur_m_adj=*(e2r(cur_m_adj,cur_lv,contextptr)._VECTptr);
	v=*(e2r(v,cur_lv,contextptr)._VECTptr);
      }
      const_iterateur it=v.begin(),itend=v.end();
      gen cur_m;
      for (;it!=itend;++it){
	vecteur cur_m_adjx(cur_m_adj);
	char_m.clear();
	int n=f_it->mult;
	x=r2sym(*it,cur_lv,contextptr);
	// compute Taylor expansion of m_adj at roots of it->fact
	// at order n-1
	for (;;){
	  --n;
	  if (n){
	    cur_m=horner(cur_m_adjx,*it,0,new_m_adj);
	    if (char_m.empty())
	      char_m=mtran(*cur_m._VECTptr);
	    else
	      char_m=mergematrice(char_m,mtran(*cur_m._VECTptr));
	    if (is_undef(char_m) || (!jordan && !is_zero(cur_m,contextptr)) ){
#ifndef NO_STDEXCEPT
	      throw(std::runtime_error("Not diagonalizable at eigenvalue "+x.print()));
#endif
	      return false;
	    }
	    cur_m_adjx=new_m_adj;
	  }
	  else {
	    cur_m=horner(cur_m_adjx,*it);
	    char_m=mergematrice(char_m,mtran(*cur_m._VECTptr));
	    if (is_undef(char_m)) return false;
	    break;
	  }
	}
	n=f_it->mult;
	if (n==1){ 
	  char_m=mtran(*cur_m._VECTptr);
	  iterateur ct=char_m.begin(),ctend=char_m.end();
	  for (;ct!=ctend;++ct){
	    if (!is_zero(*ct,contextptr))
	      break;
	  }
	  if (ct==ctend)
	    return false; // setsizeerr(gettext("egv/jordan bug"));
	  // FIXME take 1st non-0 col as eigenvector
	  *ct=*ct/lgcd(*ct->_VECTptr);
	  gen eigenvector=r2sym(*ct,cur_lv,contextptr);
	  if (is_fully_numeric(eigenvector) || numeric_matrix)
	    eigenvector=_normalize(eigenvector,contextptr);
	  p.push_back(eigenvector);
	  if (jordan){
	    vecteur vegv(taille,zero);
	    if (total_char_found>taille)
	      return false; // setsizeerr(gettext("Bug in egv/jordan"));
	    vegv[total_char_found]=x;
	    d.push_back(vegv);
	  }
	  else
	    d.push_back(x);
	  ++total_char_found;
	  continue;
	}
	if (jordan){
	  // back to external form
	  char_m=*r2sym(char_m,cur_lv,contextptr)._VECTptr;
	  int egv_found=0;
	  int char_found=0;
	  vecteur char_m_copy(char_m),pivots;
	  gen det;
	  for (;char_found<n;){ 
	    if (!mrref(char_m_copy,char_m,pivots,det,0,taille,0,taille,
		  true,egv_found,true,1,0,
		       contextptr))
	      return false;
	    if (sym)
	      char_m=gramschmidt(char_m,false,contextptr);
	    char_m_copy.clear();
	    // extract non-0 lines starting from line number egv_found
	    vecteur vegv;
	    int j=0;
	    for (;j<egv_found;++j)
	      char_m_copy.push_back(vecteur(char_m[j]._VECTptr->begin(),char_m[j]._VECTptr->end()-taille));
	    for (;j<taille;++j){
	      vegv=vecteur( char_m[j]._VECTptr->begin(),char_m[j]._VECTptr->begin()+taille);
	      if (is_zero(vegv,contextptr) || (numeric_matrix && evalf(abs(vegv,contextptr),1,contextptr)._DOUBLE_val<10*taille*epsilon(contextptr)) ) 
		break;
	      // cycle found! 
	      // update char_m_copy with all the cycle except first vector
	      char_m_copy.push_back(vecteur(char_m[j]._VECTptr->begin(),char_m[j]._VECTptr->end()-taille));
	      // Store cycle
	      const_iterateur c_it=char_m[j]._VECTptr->begin(),c_itend=char_m[j]._VECTptr->end();
	      for (;c_it!=c_itend;c_it+=taille){
		p.push_back(vecteur(c_it,c_it+taille)); // char vector
		// update d
		vegv=vecteur(taille,zero);
		if (total_char_found>=taille)
		  return false; // setsizeerr(gettext("Bug in egv/jordan"));
		if (c_it==char_m[j]._VECTptr->begin()){
		  vegv[total_char_found]=x;
		  ++egv_found;
		}
		else {
		  vegv[total_char_found-1]=1;
		  vegv[total_char_found]=x;
		}
		++char_found;
		++total_char_found;
		d.push_back(vegv);
	      }
	    }
	    for (;j<taille;++j){
	      char_m_copy.push_back(vecteur(char_m[j]._VECTptr->begin()+taille,char_m[j]._VECTptr->end()));
	    }
	  }
	} // end if (jordan)
	else {
	  d=mergevecteur(d,vecteur(n,x));
	  // back to external form
	  cur_m=r2sym(cur_m,cur_lv,contextptr);
	  // column reduction
	  matrice m_egv=mrref(mtran(*cur_m._VECTptr),contextptr);
	  if (sym){
	    // orthonormalize basis
	    m_egv=gramschmidt(matrice(m_egv.begin(),m_egv.begin()+f_it->mult),false,contextptr);
	  }
	  // non zero rows of cur_m are eigenvectors
	  const_iterateur m_it=m_egv.begin(),m_itend=m_egv.end();
	  for (; m_it!=m_itend;++m_it){
	    if (!is_zero(*m_it,contextptr))
	      p.push_back(*m_it);
	  }
	}
      }
    } // end for factorization
    p=mtran(p);
    if (jordan)
      d=mtran(d);
    return true;
  }
  matrice megv(const matrice & e,GIAC_CONTEXT){
    matrice m;
    vecteur d;
    if (!egv(e,m,d,contextptr))
      *logptr(contextptr) << gettext("Low accuracy") << endl;
    return m;
  }

  gen symb_egv(const gen & a){
    return symbolic(at_egv,a);
  }
  gen _egv(const gen & a,GIAC_CONTEXT){
    if ( a.type==_STRNG && a.subtype==-1) return  a;
    if (!is_squarematrix(a))
      return symb_egv(a);
    return megv(*a._VECTptr,contextptr);
  }
  static const char _egv_s []="egv";
  static define_unary_function_eval (__egv,&giac::_egv,_egv_s);
  define_unary_function_ptr5( at_egv ,alias_at_egv,&__egv,0,true);


  vecteur megvl(const matrice & e,GIAC_CONTEXT){
    matrice m;
    vecteur d;
    if (!egv(e,m,d,contextptr,true,false,true))
      *logptr(contextptr) << gettext("Low accuracy") << endl;
    return d;
  }
  gen symb_egvl(const gen & a){
    return symbolic(at_egvl,a);
  }
  gen _egvl(const gen & a,GIAC_CONTEXT){
    if ( a.type==_STRNG && a.subtype==-1) return  a;
    if (!is_squarematrix(a))
      return symb_pcar(a);
    return megvl(*a._VECTptr,contextptr);
  }
  static const char _egvl_s []="egvl";
  static define_unary_function_eval (__egvl,&giac::_egvl,_egvl_s);
  define_unary_function_ptr5( at_egvl ,alias_at_egvl,&__egvl,0,true);

  vecteur mjordan(const matrice & e,bool rational_jordan,GIAC_CONTEXT){
    matrice m;
    vecteur d;
    if (!egv(e,m,d,contextptr,true,rational_jordan))
      *logptr(contextptr) << "Low accuracy" << endl;
    return makevecteur(m,d);
  }
  gen symb_jordan(const gen & a){
    return symbolic(at_jordan,a);
  }
  gen jordan(const gen & a,bool rational_jordan,GIAC_CONTEXT){
    if (a.type==_VECT && a.subtype==_SEQ__VECT && a._VECTptr->size()==2 && is_squarematrix(a._VECTptr->front()) ){
      vecteur v(mjordan(*a._VECTptr->front()._VECTptr,rational_jordan,contextptr));
      if (is_undef(v))
	return v;
      gen tmpsto=sto(v[0],a._VECTptr->back(),contextptr);
      if (is_undef(tmpsto)) return tmpsto;
      return v[1];
    }
    if (!is_squarematrix(a))
      return symb_jordan(a);
    vecteur v(mjordan(*a._VECTptr,rational_jordan,contextptr));
    if (is_undef(v))
      return v;
    if (xcas_mode(contextptr)==1)
      return v[1];
    else
      return gen(v,_SEQ__VECT);
  }

  gen _jordan(const gen & a,GIAC_CONTEXT){
    if ( a.type==_STRNG && a.subtype==-1) return  a;
    bool mode=complex_mode(contextptr);
    complex_mode(true,contextptr);
    gen res=jordan(a,false,contextptr);
    complex_mode(mode,contextptr);
    return res;
  }
  static const char _jordan_s []="jordan";
  static define_unary_function_eval (__jordan,&giac::_jordan,_jordan_s);
  define_unary_function_ptr5( at_jordan ,alias_at_jordan,&__jordan,0,true);

  gen _rat_jordan(const gen & a,GIAC_CONTEXT){
    if ( a.type==_STRNG && a.subtype==-1) return  a;
    return jordan(a,true,contextptr);
  }
  static const char _rat_jordan_s []="rat_jordan";
  static define_unary_function_eval (__rat_jordan,&giac::_rat_jordan,_rat_jordan_s);
  define_unary_function_ptr5( at_rat_jordan ,alias_at_rat_jordan,&__rat_jordan,0,true);

  matrice diagonal_apply(const gen & g,const gen & x,const matrice & m,GIAC_CONTEXT){
    if (!is_squarematrix(m))
      return vecteur(1,gensizeerr(contextptr));
    int n=m.size();
    matrice res;
    for (int i=0;i<n;++i){
      vecteur v=*m[i]._VECTptr;
      v[i]=subst(g,x,v[i],false,contextptr);
      res.push_back(v);
    }
    return res;
  }

  matrice analytic_apply(const gen &ux,const gen & x,const matrice & m,GIAC_CONTEXT){
    if (!is_squarematrix(m))
      return vecteur(1,gensizeerr(contextptr));
    int n=m.size();
    matrice p,d,N,v(n),D;
    bool cplx=complex_mode(contextptr),sqrtb=withsqrt(contextptr);
    complex_mode(true,contextptr);
    withsqrt(true,contextptr);
    if (!egv(m,p,d,contextptr,true))
      return vecteur(1,gensizeerr(contextptr));
    complex_mode(cplx,contextptr);
    withsqrt(sqrtb,contextptr);
    if (int(p.size())!=n)
      return vecteur(1,gensizeerr(gettext("Unable to find all eigenvalues")));
    // search for distance of 1st non-zero non-diagonal element
    int dist=0;
    for (int i=0;i<n;++i){
      for (int j=0;j<n;++j){
	const gen & g=d[i][j];
	if (!is_zero(g,contextptr) && i!=j)
	  dist=max(dist,n-absint(i-j));
	if (i==j)
	  v[j]=g;
	else
	  v[j]=zero;
      }
      D.push_back(v);
    }
    identificateur y(" y");
    if (!dist) {// u(d) should be replaced with applying u to elements of d
      d=diagonal_apply(ux,x,d,contextptr); 
      if (is_undef(d)) return d;
      return mmult(mmult(p,d),minv(p,contextptr));
    }
    N=subvecteur(d,D);
    vecteur pol;
    if (!taylor(ux,x,y,dist,pol,contextptr)) 
      return vecteur(1,gensizeerr(ux.print()+" is not analytic"));
    if (is_undef(pol.back()))
      pol.pop_back();
    reverse(pol.begin(),pol.end());
    // subst y with D (i.e. diagonal element by diagonal element)
    int pols=pol.size();
    for (int i=0;i<pols;++i){
      if (is_undef( (pol[i]=diagonal_apply(pol[i],y,D,contextptr)) ))
	return gen2vecteur(pol[i]);
    }
    gen res=horner(pol,N);
    if (res.type!=_VECT)
      return vecteur(1,gensizeerr(contextptr));
    d=mmult(p,*res._VECTptr);
    d=mmult(d,minv(p,contextptr));
    return d;
  }

  matrice analytic_apply(const unary_function_ptr *u,const matrice & m,GIAC_CONTEXT){
    identificateur x(" x");
    gen ux=(*u)(x,contextptr);
    return analytic_apply(ux,x,m,contextptr);
  }

  // return a vector which elements are the basis of the ker of a
  bool mker(const matrice & a,vecteur & v,GIAC_CONTEXT){
    v.clear();
    gen det;
    vecteur pivots;
    matrice res;
    if (!mrref(a,res,pivots,det,0,a.size(),0,a.front()._VECTptr->size(),
	  true,0,true,1,0,
	       contextptr))
      return false;
    mdividebypivot(res);
    // put zero lines in res at their proper place, so that
    // non zero pivot are on the diagonal
    int s=res.size(),c=res.front()._VECTptr->size();
    matrice newres;
    newres.reserve(s);
    matrice::const_iterator it=res.begin(),itend=res.end();
    int i;
    for (i=0;(i<c) && (it!=itend);++i){
      if (is_zero(((*(it->_VECTptr))[i]),contextptr)){
	newres.push_back(vecteur(c,zero));
      }
      else {
	newres.push_back(*it);
	++it;
      }
    }
    for (;i<c;++i)
      newres.push_back(vecteur(c,zero));
    // now tranpose newres & resize, keep the ith line if it's ith coeff is 0
    // replace 0 by -1 to get an element of the basis
    matrice restran;
    mtran(newres,restran,res.front()._VECTptr->size());
    it=restran.begin();
    itend=restran.end();
    bool modular=!pivots.empty() && pivots.front().type==_MOD;
    for (int i=0;it!=itend;++it,++i){
      if (is_zero((*(it->_VECTptr))[i],contextptr)){
	(*(it->_VECTptr))[i]=modular?makemod(-1,*(pivots.front()._MODptr+1)):-1;
	v.push_back(*it);
      }
    }
    return true;
  }

  vecteur mker(const matrice & a,GIAC_CONTEXT){
    vecteur v;
    if (!mker(a,v,contextptr))
      return vecteur(1,gendimerr(contextptr));
    return v;
  }
  gen _ker(const gen & a,GIAC_CONTEXT){
    if ( a.type==_STRNG && a.subtype==-1) return  a;
    if (!ckmatrix(a))
      return symb_ker(a);
    vecteur v;
    if (!mker(*a._VECTptr,v,contextptr))
      return vecteur(1,gendimerr(contextptr));
    return v;    
  }
  static const char _ker_s []="ker";
  static define_unary_function_eval (__ker,&giac::_ker,_ker_s);
  define_unary_function_ptr5( at_ker ,alias_at_ker,&__ker,0,true);

  bool mimage(const matrice & a, vecteur & v,GIAC_CONTEXT){
    matrice atran;
    mtran(a,atran);
    v.clear();
    gen det;
    vecteur pivots;
    matrice res;
    if (!mrref(atran,res,pivots,det,0,atran.size(),0,atran.front()._VECTptr->size(),
	  true,0,true,1,0,
	       contextptr))
      return false;
    matrice::const_iterator it=res.begin(),itend=res.end();
    for (int i=0;it!=itend;++it,++i){
      if (!is_zero(*(it),contextptr))
	v.push_back(*it);
    }
    return true;
  }

  vecteur mimage(const matrice & a,GIAC_CONTEXT){
    vecteur v;
    if (!mimage(a,v,contextptr))
      return vecteur(1,gendimerr(contextptr));
    return v;
  }

  gen _image(const gen & a,GIAC_CONTEXT){
    if ( a.type==_STRNG && a.subtype==-1) return  a;
    if (!ckmatrix(a))
      return symb_image(a);
    vecteur v;
    if (!mimage(*a._VECTptr,v,contextptr))
      return gensizeerr(contextptr);
    return v;    
  }
  static const char _image_s []="image";
  static define_unary_function_eval (__image,&giac::_image,_image_s);
  define_unary_function_ptr5( at_image ,alias_at_image,&__image,0,true);

  vecteur cross(const vecteur & v_orig,const vecteur & w_orig){
    vecteur v(v_orig),w(w_orig);
    int s1=v.size(),s2=w.size();
    if (s1==2){
      v.push_back(0);
      ++s1;
    }
    if (s2==2){
      w.push_back(0);
      ++s2;
    }
    if (s1!=3 || s2!=3)
      return vecteur(1,gendimerr(gettext("cross")));
    vecteur res;
    res.push_back(v[1]*w[2]-v[2]*w[1]);
    res.push_back(v[2]*w[0]-v[0]*w[2]);
    res.push_back(v[0]*w[1]-v[1]*w[0]);
    return res;
  }
  gen symb_cross(const gen & arg1,const gen & arg2){
    return symbolic(at_cross,makevecteur(arg1,arg2));
  }
  gen symb_cross(const gen & args){
    return symbolic(at_cross,args);
  }
  gen cross(const gen & a,const gen & b){
    gen g1=remove_at_pnt(a);
    gen g2=remove_at_pnt(b);
    if (is_undef(g1) || g1.type!=_VECT || is_undef(g2) || g2.type!=_VECT)
      return gensizeerr(gettext("cross"));
    if (g1.subtype==_VECTOR__VECT)
      return cross(vector2vecteur(*g1._VECTptr),g2);
    if (g2.subtype==_VECTOR__VECT)
      return cross(g1,vector2vecteur(*g2._VECTptr));
    return cross(*g1._VECTptr,*g2._VECTptr);
  }
  gen _cross(const gen &args,GIAC_CONTEXT){
    if (args.type==_STRNG && args.subtype==-1) return args;
    if (args.type!=_VECT)
      return symb_cross(args);
    if (args._VECTptr->size()!=2)
      return gendimerr(contextptr);
    if (ckmatrix(args._VECTptr->front()) || ckmatrix(args._VECTptr->back()))
      return gensizeerr(contextptr);
    return cross(args._VECTptr->front(),args._VECTptr->back());
  }
  static const char _cross_s []="cross";
  string texprintascross(const gen & feuille,const char * sommetstr,GIAC_CONTEXT){
    return texprintsommetasoperator(feuille," \\wedge ",contextptr);
  }
  static define_unary_function_eval4 (__cross,&giac::_cross,_cross_s,0,texprintascross);
  define_unary_function_ptr5( at_cross ,alias_at_cross,&__cross,0,true);

  static string printassize(const gen & feuille,const char * sommetstr,GIAC_CONTEXT){
    string res(sommetstr);
    if (xcas_mode(contextptr)>0)
      res="nops";
    return res+"("+feuille.print(contextptr)+")";
  }
  gen symb_size(const gen & args){
    return symbolic(at_size,args);
  }
  gen _size(const gen &args,GIAC_CONTEXT){
    if (args.type==_STRNG && args.subtype==-1) return args;
    if (args.type==_STRNG)
      return (int) args._STRNGptr->size();
    if (args.type==_SYMB){
      if (args._SYMBptr->feuille.type==_VECT)
	return (int) args._SYMBptr->feuille._VECTptr->size();
      else
	return 1;
    }
    if (args.type==_POLY)
      return args._POLYptr->coord.size();
    if (args.type!=_VECT)
      return 1;
    return (int) args._VECTptr->size();
  }
  static const char _size_s []="size";
  static define_unary_function_eval2 (__size,&giac::_size,_size_s,&printassize);
  define_unary_function_ptr5( at_size ,alias_at_size,&__size,0,true);

#ifdef HAVE_LIBGSL  

  int vecteur2gsl_vector(const_iterateur it,const_iterateur itend,gsl_vector * w,GIAC_CONTEXT){
#ifdef DEBUG_SUPPORT
    if (itend-it!=signed(w->size))
      setsizeerr(gettext("vecteur.cc vecteur2gsl_vector"));
#endif
    gen g;
    int res=GSL_SUCCESS;
    for (int i = 0; it!=itend; ++i,++it){
      g=it->evalf(1,contextptr);
      if (g.type==_DOUBLE_)
	gsl_vector_set (w, i, g._DOUBLE_val);
      else {
	gsl_vector_set (w, i, nan());	
	res=!GSL_SUCCESS;
      }
    }
    return res;
  }

  int vecteur2gsl_vector(const vecteur & v,gsl_vector * w,GIAC_CONTEXT){
    const_iterateur it=v.begin(),itend=v.end();
    return vecteur2gsl_vector(it,itend,w,contextptr);
  }
  // this function allocate all space needed for the gsl_vector
  gsl_vector * vecteur2gsl_vector(const vecteur & v,GIAC_CONTEXT){
    const_iterateur it=v.begin(),itend=v.end();
    gsl_vector * w = gsl_vector_alloc (itend-it);
    vecteur2gsl_vector(it,itend,w,contextptr);
    return w;
  }

  // this function does not deallocate the gsl vector
  // call gsl_vector_free(v) for this
  vecteur gsl_vector2vecteur(const gsl_vector * v){
    vecteur res;
    int s=v->size;
    res.reserve(s);
    for (int i=0;i<s;++i)
      res.push_back(gsl_vector_get(v,i));
    return res;
  }

  int matrice2gsl_matrix(const matrice & m,gsl_matrix * w,GIAC_CONTEXT){
    int s1=w->size1,s2=w->size2;
#ifdef DEBUG_SUPPORT
    ckmatrix(m);
    if (mrows(m)!=s1 || mcols(m)!=s2)
      setdimerr();
#endif
    gen g;
    const_iterateur it=m.begin(),itend=m.end();
    int res=GSL_SUCCESS;
    for (int i = 0; it!=itend; ++i,++it){
      if (it->type!=_VECT)
	res=!GSL_SUCCESS;
      vecteur & v =*it->_VECTptr;
      const_iterateur jt=v.begin(),jtend=v.end();
      for (int j=0;jt!=jtend;++j,++jt){
	g=evalf(*jt,1,contextptr);
	if (g.type==_DOUBLE_)
	  gsl_matrix_set(w,i,j,g._DOUBLE_val);
	else {
	  res=!GSL_SUCCESS;
	  gsl_matrix_set(w,i,j,nan());	  
	}
      }
    }
    return res;
  }

  // this function allocate all space needed for the gsl_matrix
  gsl_matrix * matrice2gsl_matrix(const matrice & m,GIAC_CONTEXT){
    int n1=mrows(m),n2=mcols(m);
    gsl_matrix * w = gsl_matrix_alloc (n1,n2);
    matrice2gsl_matrix(m,w,contextptr);
    return w;
  }
  
  // this function does not deallocate the gsl vector
  // call gsl_matrix_free(v) for this
  matrice gsl_matrix2matrice(const gsl_matrix * v){
    matrice res;
    int s1=v->size1,s2=v->size2;
    res.reserve(s1);
    for (int i=0;i<s1;++i){
      vecteur tmp;
      tmp.reserve(s2);
      for (int j=0;j<s2;++j){
	tmp.push_back(gsl_matrix_get(v,i,j));
      }
      res.push_back(tmp);
    }
    return res;
  }

  vecteur gsl_permutation2vecteur(const gsl_permutation * p,GIAC_CONTEXT){
    int s=p->size;
    vecteur res(s);
    for (int i=0;i<s;++i)
      res[i]=(int)gsl_permutation_get(p,i)+(xcas_mode(contextptr)?1:0);
    return res;
  }
#endif // HAVE_LIBGSL

  bool mlu(const matrice & a0,vecteur & P,matrice & L,matrice & U,GIAC_CONTEXT){
    matrice a(a0);
    bool modular=false;
    if (!is_squarematrix(a)){
      if (a.front().type==_VECT && !a.front()._VECTptr->empty() && (a.back()==at_irem || a.back()==at_ichinrem)){
	modular=true;
	a=*a.front()._VECTptr;
      }
      if (!is_squarematrix(a))
	return false; // setsizeerr(gettext("Expecting a square matrix"));
    }
    gen det;
    vecteur pivots;
    matrice res;
    int s=a.size();
    if (!mrref(a,res,pivots,det,0,s,0,s,
	  false,0,false,(modular?3:0) /* algorithm */,2 /* lu */,
	       contextptr))
      return false;
    if (pivots.empty())
      return false;
    gen tmp=pivots.back();
    if (tmp.type!=_VECT)
      return false; // setsizeerr();
    P=*tmp._VECTptr;
    // Make L and U from res
    L.reserve(s); U.reserve(s);
    for (int i=0;i<s;++i){
      vecteur & v=*res[i]._VECTptr;
      vecteur wl(s);
      for (int j=0;j<i;++j){ // L part
	wl[j]=v[j];
      }
      wl[i]=1;
      L.push_back(wl);
      vecteur wu(s);
      for (int j=i;j<s;++j){ // U part
	wu[j]=v[j];
      }
      U.push_back(wu);
    }
    return true;
  }

  // in: l= r rows, c cols
  // out: l= r,r and u=r rows, c cols
  void splitlu(matrice & l,matrice & u){
    u=l;
    int r,c;
    mdims(l,r,c);
    for (int i=0;i<r;++i){
      vecteur li=*l[i]._VECTptr;
      li.resize(r);
      vecteur & ui=*u[i]._VECTptr;
      for (int j=0;j<i;++j){
	ui[j]=0;
      }
      li[i]=1;
      for (int j=i+1;j<r;++j){
	li[j]=0;
      }
      l[i]=li;
    }
  }

  gen lu(const gen &args,GIAC_CONTEXT){
    matrice L,U,P;
#ifdef HAVE_LIBLAPACK
    if (ckmatrix(args) && is_fully_numeric(args)){
      
      /* DGETRF( M, N, A, LDA, IPIV, INFO ), ZGETRF( M, N, A, LDA, IPIV, INFO )
       *  Purpose
       *  =======
       *
       *  DGETRF/ZGETRF computes an LU factorization of a general M-by-N matrix A
       *  using partial pivoting with row interchanges.
       *
       *  The factorization has the form
       *     A = P * L * U
       *  where P is a permutation matrix, L is lower triangular with unit
       *  diagonal elements (lower trapezoidal if m > n), and U is upper
       *  triangular (upper trapezoidal if m < n).
       *
       *  This is the right-looking Level 3 BLAS version of the algorithm.
       *
       *  Arguments
       *  =========
       *
       *  M       (input) INTEGER
       *          The number of rows of the matrix A.  M >= 0.
       *
       *  N       (input) INTEGER
       *          The number of columns of the matrix A.  N >= 0.
       *
       *  A       (input/output) DOUBLE PRECISION/COMPLEX array, dimension (LDA,N)
       *          On entry, the M-by-N matrix to be factored.
       *          On exit, the factors L and U from the factorization
       *          A = P*L*U; the unit diagonal elements of L are not stored.
       *
       *  LDA     (input) INTEGER
       *          The leading dimension of the array A.  LDA >= max(1,M).
       *
       *  IPIV    (output) INTEGER array, dimension (min(M,N))
       *          The pivot indices; for 1 <= i <= min(M,N), row i of the
       *          matrix was interchanged with row IPIV(i).
       */
      integer M,N,LDA,INFO;
      int m,n;
      mdims(*args._VECTptr,m,n);
      int mn=min(m,n);
      M=m; N=n;
      LDA=M;
      integer * IPIV=new integer[mn];
      for (int i=0;i<mn;++i)
	P.push_back(i+1);
      if (is_zero(im(args,contextptr))){
	double * A = new double[M*N];
	matrice2lapack(*args._VECTptr,A,contextptr);
	dgetrf_( &M, &N, A, &LDA, IPIV, &INFO );
	if (INFO){
	  delete [] IPIV;
	  delete [] A;
	  return gensizeerr(gettext("LAPACK LU error"));
	}
	lapack2matrice(A,M,N,L);
	// get U from upper part of L and clear
	splitlu(L,U);
	// get P
	for (int i=1;i<mn;++i){
	  if (IPIV[i-1]!=i)
	    swapgen(P[i-1],P[IPIV[i-1]-1]);
	}
	if (!xcas_mode(contextptr)){
	  for (int i=0;i<mn;++i)
	    P[i] -= 1;
	}
	delete [] IPIV;
	delete [] A;
	return gen(makevecteur(P,L,U),_SEQ__VECT);
      }
      doublef2c_complex * A = new doublef2c_complex[M*N];
      matrice2zlapack(*args._VECTptr,A,contextptr);
      zgetrf_( &M, &N, A, &LDA, IPIV, &INFO );
      if (INFO){
	delete [] IPIV;
	delete [] A;
	return gensizeerr(gettext("LAPACK LU error"));
      }
      zlapack2matrice(A,M,N,L);
      // get U from upper part of L and clear
      splitlu(L,U);
      // get P
      for (int i=1;i<mn;++i){
	if (IPIV[i-1]!=i)
	  swapgen(P[i-1],P[IPIV[i-1]-1]);
      }
      if (!xcas_mode(contextptr)){
	for (int i=0;i<mn;++i)
	  P[i] -= 1;
      }
      delete [] IPIV;
      delete [] A;
      return gen(makevecteur(P,L,U),_SEQ__VECT);
    }
#endif
#ifdef HAVE_LIBGSL
    bool gsl_lu = is_fully_numeric(args) && is_zero(im(args,contextptr),contextptr);
    if (gsl_lu){
      if (!is_squarematrix(args))
	return gensizeerr(gettext("Expecting a square matrix"));
      gsl_matrix * m=matrice2gsl_matrix(*args._VECTptr,contextptr);
      int s1=m->size1;
      gsl_permutation * p=gsl_permutation_alloc (s1);
      int sign;
      gsl_linalg_LU_decomp (m, p,&sign);
      P=gsl_permutation2vecteur(p,contextptr);
      L.reserve(s1);
      U.reserve(s1);
      // get L and U
      for (int i=0;i<s1;++i){
	vecteur l(s1),u(s1);
	for (int j=0;j<i;++j){
	  l[j]=gsl_matrix_get(m,i,j);
	}
	l[i]=1.0;
	for (int j=i;j<s1;++j){
	  u[j]=gsl_matrix_get(m,i,j);
	}
	L.push_back(l);
	U.push_back(u);
      }
      gsl_permutation_free(p);
      gsl_matrix_free(m);
      return gen(makevecteur(P,L,U),_SEQ__VECT);
    }
#endif // HAVE_LIBGSL
    if (args.type!=_VECT)
      return gentypeerr(contextptr);
    // Giac LU decomposition
    if (!mlu(*args._VECTptr,P,L,U,contextptr))
      return gendimerr(contextptr);
    if (xcas_mode(contextptr)){
      int s=P.size();
      for (int i=0;i<s;++i){
	P[i]=P[i]+1;
      }
    }
    return gen(makevecteur(P,L,U),_SEQ__VECT);
  }
  static const char _lu_s []="lu";
  static define_unary_function_eval (__lu,&giac::lu,_lu_s);
  define_unary_function_ptr5( at_lu ,alias_at_lu,&__lu,0,true);

  bool matrice2lapack(const matrice & m,double * A,GIAC_CONTEXT){
    const_iterateur it=m.begin(),itend=m.end();
    gen g;
    int rows=itend-it;
    for (int i = 0; it!=itend; ++i,++it){
      if (it->type!=_VECT)
	return false;
      vecteur & v =*it->_VECTptr;
      const_iterateur jt=v.begin(),jtend=v.end();
      for (int j = 0; jt!=jtend;++j, ++jt){
	g=evalf_double(*jt,1,contextptr);
	if (g.type==_DOUBLE_)
	  A[i + j * rows] = g._DOUBLE_val;
	else 
	  return false;
      }
    }
    return true;
  }

  void lapack2matrice(double * A,unsigned rows,unsigned cols,matrice & R){
    R.reserve(rows);
    for (int i=0;i<rows;++i){
      vecteur r(cols);
      for (int j=0;j<cols;++j)
	r[j] = A[i + j * rows];
      R.push_back(r);
    }
  }

  gen qr(const gen &args_orig,GIAC_CONTEXT){
    gen args;
    int method=0; // use -1 to check built-in qr
    if ( (args_orig.type==_VECT) && (args_orig._VECTptr->size()==2) && (args_orig._VECTptr->back().type==_INT_)){
      args=args_orig._VECTptr->front();
      method=args_orig._VECTptr->back().val;
    }
    else
      args=args_orig;
    if (!ckmatrix(args))
      return symbolic(at_qr,args);
    if (!is_zero(im(args,contextptr),contextptr))
      return gensizeerr(gettext("Complex entry!"));
    if (method<0 || !is_fully_numeric(evalf_double(args,1,contextptr))){
      matrice r;
      if (is_fully_numeric(args)){ 
	// qr decomposition using rotations, numerically stable
	// but not suited to exact computations
	matrice h=*args._VECTptr,p(midn(h.size()));
	std_matrix<gen> H,P;
	matrice2std_matrix_gen(h,H);
	matrice2std_matrix_gen(p,P);
	qr_ortho(H,P,contextptr);
	std_matrix_gen2matrice(H,h);
	std_matrix_gen2matrice(P,p);
	return makevecteur(_trn(p,contextptr),h,midn(H[0].size()));
      }
      // qr decomposition using GramSchmidt (not numerically stable)
      matrice res(gramschmidt(*_trn(args,contextptr)._VECTptr,r,method==-1,contextptr));
      return gen(makevecteur(_trn(res,contextptr),r,midn(r.size())),_SEQ__VECT);
    }
#ifdef HAVE_LIBLAPACK
    if (!CAN_USE_LAPACK
	|| dgeqrf_ == NULL
	|| dorgqr_ == NULL)
      return gensizeerr(gettext("LAPACK not available"));
    
    const matrice &m = *args._VECTptr;
    integer rows = mrows(m), cols = mcols(m);
    if (rows < cols)
      return gendimerr(contextptr);
    integer info;
    double *A = new double[rows * cols];
    if (!matrice2lapack(m,A,contextptr))
      return gensizeerr(gettext("Lapack conversion error"));
    double *tau = new double[cols];
    integer lwork = -1;
    double worktmp;
    if (!is_zero(im(args,contextptr))){
      // complex QR decomposition, currently disabled above
    }
    // first call to determine optimum work vector size lwork
    dgeqrf_(&rows, &cols, A, &rows, tau, &worktmp, &lwork, &info);
    if (info){
      delete [] A;
      delete [] tau;
      return gensizeerr(gettext("LAPACK error ") + (-info));
    }
    
    lwork = (int)worktmp;
    double *work = new double[lwork];
    
    // second call, computes the QR-decomposition
    dgeqrf_(&rows, &cols, A, &rows, tau, work, &lwork, &info);
    if (info){
      delete [] A;
      delete [] tau;
      delete [] work;
      return gensizeerr(gettext("LAPACK error ") + (-info));
    }
    
    // load R from the upper right part of A
    matrice R;
    R.reserve(rows);
    for (int i=0;i<rows;++i){
      vecteur r(cols);
      for (int j=i;j<cols;++j)
	r[j] = A[i + j * rows];
      R.push_back(r);
    }
    
    // compute Q and store it in A
    dorgqr_(&rows, &cols, &cols, A, &rows, tau, work, &lwork, &info);
    delete [] tau;
    delete [] work;
    if (info){
      delete [] A;
      return gensizeerr(gettext("LAPACK error ") + (-info));
    }
    
    // load Q from A
    matrice Q;
    Q.reserve(rows);
    for (int i=0;i<rows;++i){
      vecteur q(cols);
      for (int j=0;j<cols;++j)
	q[j] = A[i + j * rows];
      Q.push_back(q);
    }
    
    delete [] A;
    
    // I prefer to get Q and R returned. Your mileage may vary.
    return gen(makevecteur(Q,R),_SEQ__VECT);
#endif // HAVE_LIBLAPACK
#ifdef HAVE_LIBGSL
    {
    gsl_matrix * m=matrice2gsl_matrix(*args._VECTptr,contextptr);
    int s1=m->size1,s2=m->size2;
    gsl_vector * tau=gsl_vector_alloc(min(s1,s2));
    gsl_linalg_QR_decomp (m,tau);
    matrice R;
    R.reserve(s1);
    // get R
    for (int i=0;i<s1;++i){
      vecteur r(s2);
      for (int j=i;j<s2;++j){
	r[j]=gsl_matrix_get(m,i,j);
      }
      R.push_back(r);
    }
    // get the list of tau_i,v_i
    vecteur Q;
    for (int i=0;i<signed(tau->size);++i){
      vecteur tmp(m->size2);
      tmp[i]=1.0;
      for (int j=i+1;j<signed(m->size2);++j)
	tmp[j]=gsl_matrix_get(m,j,i);
      Q.push_back(makevecteur(gsl_vector_get(tau,i),tmp));
    }
    gsl_vector_free(tau);
    gsl_matrix_free(m);
    // return gen(makevecteur(Q,R),_SEQ__VECT);
    return R;
    }
#endif // HAVE_LIBGSL

    return symbolic(at_qr,args);
  }
  static const char _qr_s []="qr";
  static define_unary_function_eval (__qr,&giac::qr,_qr_s);
  define_unary_function_ptr5( at_qr ,alias_at_qr,&__qr,0,true);

  matrice thrownulllines(const matrice & res){
    int i=res.size()-1;
    for (;i>=0;--i){
      if (!is_zero(res[i],context0))
	break;
    }
    return vecteur(res.begin(),res.begin()+i+1);
  }
  gen _basis(const gen &args,GIAC_CONTEXT){
    if (args.type==_STRNG && args.subtype==-1) return args;
    if (!ckmatrix(args))
      return symbolic(at_basis,args);
    matrice res=mrref(*args._VECTptr,contextptr);
    return gen(thrownulllines(res),_SET__VECT);
  }
  static const char _basis_s []="basis";
  static define_unary_function_eval (__basis,&giac::_basis,_basis_s);
  define_unary_function_ptr5( at_basis ,alias_at_basis,&__basis,0,true);

  // Sylvester matrix, in lines line0=v1 0...0, line1=0 v1 0...0, etc.
  matrice sylvester(const vecteur & v1,const vecteur & v2){
    int m=v1.size()-1;
    int n=v2.size()-1;
    if (m<0 || n<0)
      return vecteur(0);
    matrice res(m+n);
    for (int i=0;i<n;++i){
      vecteur w(m+n);
      for (int j=0;j<=m;++j)
	w[i+j]=v1[j];
      res[i]=w;
    }
    for (int i=0;i<m;++i){
      vecteur w(m+n);
      for (int j=0;j<=n;++j)
	w[i+j]=v2[j];
      res[n+i]=w;
    }
    return res;
  }

  gen _sylvester(const gen &args,GIAC_CONTEXT){
    if (args.type==_STRNG && args.subtype==-1) return args;
    if (args.type!=_VECT || args._VECTptr->size()<2)
      return gensizeerr(contextptr);
    vecteur & v = *args._VECTptr;
    gen x(vx_var);
    if (v.size()>2)
      x=v[2];
    gen p1(_e2r(makevecteur(v[0],x),contextptr));
    gen p2(_e2r(makevecteur(v[1],x),contextptr));
    if (p1.type!=_VECT || p2.type!=_VECT)
      return gensizeerr(contextptr);
    vecteur & v1 =*p1._VECTptr;
    vecteur & v2 =*p2._VECTptr;
    return sylvester(v1,v2);
  }
  static const char _sylvester_s []="sylvester";
  static define_unary_function_eval (__sylvester,&giac::_sylvester,_sylvester_s);
  define_unary_function_ptr5( at_sylvester ,alias_at_sylvester,&__sylvester,0,true);

  gen _ibasis(const gen &args,GIAC_CONTEXT){
    if (args.type==_STRNG && args.subtype==-1) return args;
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2) )
      return symbolic(at_basis,args);
    gen g=args._VECTptr->front(),h=args._VECTptr->back();
    if (!ckmatrix(g) || !ckmatrix(h))
      return gensizeerr(contextptr);
    vecteur & v1=*g._VECTptr;
    vecteur & v2=*h._VECTptr;
    if (v1.empty() || v2.empty())
      return vecteur(0);
    vecteur v=mker(mtran(mergevecteur(v1,v2)),contextptr);
    if (is_undef(v)) return v;
    // if v is not empty compute each corresponding vector of the basis
    int s=v1.size();
    int l=v1.front()._VECTptr->size();
    matrice res;
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      vecteur tmp(l);
      vecteur & i=*it->_VECTptr;
      for (int j=0;j<s;++j)
	tmp=addvecteur(tmp,multvecteur(i[j],*v1[j]._VECTptr));
      res.push_back(tmp);
    }
    return gen(thrownulllines(mrref(res,contextptr)),_SET__VECT);
  }
  static const char _ibasis_s []="ibasis";
  static define_unary_function_eval (__ibasis,&giac::_ibasis,_ibasis_s);
  define_unary_function_ptr5( at_ibasis ,alias_at_ibasis,&__ibasis,0,true);

  gen _svd(const gen &args_orig,GIAC_CONTEXT){
    if (args_orig.type==_STRNG && args_orig.subtype==-1) return args_orig;
    gen args;
    int method=0; // use -1 to check built-in svd, -2 for svl (singular values only) 
    if ( (args_orig.type==_VECT) && (args_orig._VECTptr->size()==2) && (args_orig._VECTptr->back().type==_INT_)){
      args=args_orig._VECTptr->front();
      method=args_orig._VECTptr->back().val;
    }
    else
      args=args_orig;
    if (!ckmatrix(args))
      return symbolic(at_svd,args);
    if (!is_zero(im(args,contextptr),contextptr))
      return gensizeerr(gettext("Complex entry!"));
    if (method>=0 && is_fully_numeric(evalf_double(args,1,contextptr))){
#ifdef HAVE_LIBLAPACK
      if (!CAN_USE_LAPACK
	  || dgeqrf_ == NULL
	  || dorgqr_ == NULL)
	return gensizeerr(gettext("LAPACK not available"));
      
      const matrice &m = *args._VECTptr;
      gen g;
      const_iterateur it=m.begin(),itend=m.end();
      integer rows = mrows(m), cols = mcols(m);
      if (rows < cols)
	return gendimerr(contextptr);
      double *A = new double[rows * cols];
      matrice2lapack(m,A,contextptr);
      integer info;
      char jobU = 'A', jobVT = 'A';
      double *S = new double[cols];
      double *U = new double[rows * rows];
      double *VT = new double[cols * cols];
      integer lwork = -1;
      double worktmp;
      
      // first call to determine optimum work vector size lwork
      dgesvd_(&jobU, &jobVT, &rows, &cols, A, &rows, S, U, &rows, VT, &cols, &worktmp, &lwork, &info);
      if (info){
	delete [] A;
	delete [] S;
	delete [] U;
	delete [] VT;
	return gensizeerr(gettext("LAPACK error ") + (-info));
      }
      
      lwork = (int)worktmp;
      double *work = new double[lwork];
      
      // second call, computes the SVD
      dgesvd_(&jobU, &jobVT, &rows, &cols, A, &rows, S, U, &rows, VT, &cols, work, &lwork, &info);
      if (info){
	delete [] A;
	delete [] S;
	delete [] U;
	delete [] VT;
	delete [] work;
	setsizeerr(gettext("LAPACK error ") + (-info));
      }
      
      // load S
      vecteur s(cols);
      for (int j=0;j<cols;++j)
	s[j] = S[j];
      
      // load U
      matrice mU;
      lapack2matrice(U,rows,rows,mU);
      
      // load VT
      matrice mVT;
      lapack2matrice(VT,cols,cols,mVT);
      
      delete [] A;
      delete [] S;
      delete [] U;
      delete [] VT;
      delete [] work;
      
      return gen(makevecteur(mU,s,mVT),_SEQ__VECT);
#endif // HAVE_LIBLAPACK
#ifdef HAVE_LIBGSL
      {   
	gsl_matrix * u=matrice2gsl_matrix(*args._VECTptr,contextptr);
	int s1=u->size1,s2=u->size2;
	gsl_vector * work=gsl_vector_alloc (s1);
	gsl_matrix * v=gsl_matrix_alloc(s2,s2);
	gsl_vector * s=gsl_vector_alloc(s1);
	gsl_matrix * x=gsl_matrix_alloc(s1,s1);
	switch(method){
	case _GOLUB_REINSCH_MOD_DECOMP:
	  gsl_linalg_SV_decomp_mod(u,x,v,s,work);
	  break;
	case _JACOBI_DECOMP:
	  gsl_linalg_SV_decomp_jacobi(u,v,s);
	  break;
	default:
	  gsl_linalg_SV_decomp (u, v,s,work);
	  break;
	}
	gsl_vector_free(work);
	gsl_matrix_free(x);
	matrice U(gsl_matrix2matrice(u)),S(gsl_vector2vecteur(s)),V(gsl_matrix2matrice(v)); // A=U*S*tran(V)
	gsl_matrix_free(u);
	gsl_matrix_free(v);
	gsl_vector_free(s);
	return gen(makevecteur(U,S,V),_SEQ__VECT); // M=U*diag(S)*tran(V)
      }
#endif // HAVE_LIBGSL
    }    
    // non numeric code
    matrice tm,mtm,p,d,pd,invpd,v,tv;
    mtran(*args._VECTptr,tm);
    mmult(*args._VECTptr,tm,mtm);
    if (!egv(mtm,p,d,contextptr,true,false))
      return gensizeerr(contextptr);
    // should reorder eigenvalue (decreasing order)
    int s=d.size();
    vecteur svl;
    gen eps=epsilon(contextptr);
    for (int i=0;i<s;++i){
      vecteur vi=*d[i]._VECTptr;
      gen & di=vi[i];
      di=sqrt(di,contextptr);
      svl.push_back(di);
      d[i]=vi;
      if (method==-1 && is_greater(eps,di,contextptr))
	return gensizeerr(gettext("0 as singular value, not implemented"));
    }
    if (method==-2)
      return svl;
    // now m*tran(m)=p*d*tran(p)
    mmult(p,d,pd);
    invpd=minv(pd,contextptr);
    mmult(invpd,*args._VECTptr,v);
    mtran(v,tv);
    return gen(makevecteur(p,tv,svl),_SEQ__VECT); // M=p*diag(svl)*tran(tv)
  }
  static const char _svd_s []="svd";
  static define_unary_function_eval (__svd,&giac::_svd,_svd_s);
  define_unary_function_ptr5( at_svd ,alias_at_svd,&__svd,0,true);

  gen _cholesky(const gen &_args,GIAC_CONTEXT){
    if (_args.type==_STRNG && _args.subtype==-1) return _args;
    if (!is_squarematrix(_args))
      return gensizeerr(contextptr);
    gen args;
    if (_args==_tran(_args,contextptr))
      args=_args;
    else
      args=(_args+_tran(_args,contextptr))/2;
#ifdef HAVE_LIBGSL
    if (is_fully_numeric(args) && is_zero(im(args,contextptr),contextptr)){
      gsl_matrix * m=matrice2gsl_matrix(*args._VECTptr,contextptr);
      int s1=m->size1;
      int i=gsl_linalg_cholesky_decomp (m);
      if (i==GSL_EDOM)
	return gensizeerr(gettext("Non positive definite"));
      // clear upper part
      for (i=0;i<s1;++i){
	for (int j=i+1;j<s1;++j)
	  gsl_matrix_set(m,i,j,0.0);
      }
      matrice LL(gsl_matrix2matrice(m));
      gsl_matrix_free(m);
      return LL;
    }
#endif // HAVE_LIBGSL
    matrice &A=*args._VECTptr;
    int n=A.size(),j,k,l;
    std_matrix<gen> C(n,vecteur(n));
    for (j=0;j<n;j++) {
      gen s;
      for (l=j;l<n;l++) {
	s=0;
	for (k=0;k<j;k++) {
	  if (is_zero(C[k][k],contextptr)) 
	    return gensizeerr(gettext("Not invertible matrice"));
	  //if (is_strictly_positive(-C[k][k])) setsizeerr(gettext("Not a positive define matrice"));
	  s=s+C[l][k]*C[j][k]/C[k][k];
	}
	C[l][j]=ratnormal(A[l][j]-s);
      }
    }
    for (k=0;k<n;k++) {
      gen c=normal(inv(sqrt(C[k][k],contextptr),contextptr),contextptr);
      for (j=k;j<n;j++) {
	C[j][k]=C[j][k]*c;
      }
    }
    matrice Cmat;
    std_matrix_gen2matrice(C,Cmat);
    return Cmat;
/*
    matrice & A = *args._VECTptr;
    int n=A.size(),j,k,l;
    // Use LU decomposition without line permutation
    matrice LU,pivots;
    gen det;
    mrref(A,LU,pivots,det,0,n,0,n,false,0,false,false,3,contextptr);
    if (is_zero(det)) return gensizeerr(contextptr"Not a positive defined matrix");
    matrice D,L;
    for (int i=0;i<n;++i){
      vecteur v(n);
      v[i]=sqrt(LU[i][i]);
      D.push_back(v);
      vecteur w(n);
      w[i]=1;
      for (j=0;j<i;j++)
	w[j]=LU[i][j];
      L.push_back(w);
    }
    return ckmultmatvecteur(L,D);
*/
    /*
    std_matrix<gen> C(n,vecteur(n));
    for (j=0;j<n;++j){
      gen s;
      for (k=0;k<j;++k){
	s=s+pow(C[j][k],2);
      }
      gen c2=A[j][j]-s;
      if (is_strictly_positive(-c2,contextptr))
	return gensizeerr(contextptr"Not a positive defined matrix");
      gen c=normal(sqrt(c2,contextptr),contextptr);
      C[j][j]=c;
      for (l=j+1;l<n;++l){
	s=0;
	for (k=0;k<j;++k)
	  s=s+C[l][k]*C[j][k];
	C[l][j]=normal((A[l][j]-s)/c,contextptr);
      }
    }
    matrice Cmat;
    std_matrix_gen2matrice(C,Cmat);
    return Cmat;
    */
  }
  static const char _cholesky_s []="cholesky";
  static define_unary_function_eval (__cholesky,&giac::_cholesky,_cholesky_s);
  define_unary_function_ptr5( at_cholesky ,alias_at_cholesky,&__cholesky,0,true);

  gen l2norm(const vecteur & v,GIAC_CONTEXT){
    const_iterateur it=v.begin(),itend=v.end();
    gen res;
    for (;it!=itend;++it)
      res = res + (*it)*conj(*it,contextptr);
    return sqrt(res,contextptr);
  }

  matrice gramschmidt(const matrice & m,matrice & r,bool normalize,GIAC_CONTEXT){
    r.clear();
    vecteur v(m);
    int s=v.size();
    if (!s)
      return v;
    vecteur sc(1,dotvecteur(*conj(v[0],contextptr)._VECTptr,*v[0]._VECTptr));
    vecteur rcol0(s);
    rcol0[0]=1;
    r.push_back(rcol0);
    for (int i=1;i<s;++i){
      gen cl,coeff;
      vecteur rcol(s);
      rcol[i]=1;
      for (int j=0;j<i;++j){
	coeff=rdiv(dotvecteur(*conj(v[j],contextptr)._VECTptr,*v[i]._VECTptr),sc[j]);
	cl=cl+coeff*v[j];
	rcol[j]=coeff;
      }
      v[i]=v[i]-cl;
      sc.push_back(dotvecteur(*conj(v[i],contextptr)._VECTptr,*v[i]._VECTptr));
      r.push_back(rcol);
      if (is_zero(sc.back(),contextptr))
	break;
    }
    r=mtran(r);
    if (normalize){
      gen coeff;
      for (int i=0;i<s;++i){
	if (is_zero(sc[i],contextptr))
	  break;
	coeff=sc[i]=sqrt(sc[i],contextptr);
	v[i]=rdiv(v[i],coeff);
      }
      for (int i=0;i<s;++i){
	if (is_zero(sc[i],contextptr))
	  break;
	r[i]=sc[i]*r[i];
      }
    }
    return v;
  }

  matrice gramschmidt(const matrice & m,bool normalize,GIAC_CONTEXT){
    matrice r;
    return gramschmidt(m,r,normalize,contextptr);
  }

  // lll decomposition of M, returns S such that S=A*M=L*O
  // L is lower and O is orthogonal
  matrice lll(const matrice & M,matrice & L,matrice & O,matrice &A,GIAC_CONTEXT){
    if (!ckmatrix(M))
      return vecteur(1,gensizeerr(contextptr));
    matrice res(M);
    int n=res.size();
    if (!n)
      return res;
    int c=res[0]._VECTptr->size();
    if (c<n)
      return vecteur(1,gendimerr(contextptr));
    A=midn(c);
    A=vecteur(A.begin(),A.begin()+n);
    int k=0;
    for (;k<n;){
      if (!k){ // push first vector
	vecteur tmp(c);
	tmp[0]=1;
	L.push_back(tmp);
	O.push_back(res.front());
	++k;
	continue;
      }
      // Find new vector in L,O
      vecteur tmp(c);
      gen Otmp(res[k]);
      for (int j=0;j<k;++j){
	// tmp[j]=dotvecteur(res[j],res[k])/dotvecteur(res[j],res[j]);
	tmp[j]=dotvecteur(conj(O[j],contextptr),Otmp)/dotvecteur(conj(O[j],contextptr),O[j]);
	Otmp=subvecteur(*Otmp._VECTptr,multvecteur(tmp[j],*O[j]._VECTptr));
      }
      tmp[k]=1;
      L.push_back(tmp);
      O.push_back(Otmp);
      // Compare norm of O[k] and O[k-1]
      for (int j=k-1;j>=0;--j){
	gen alpha=dotvecteur(conj(O[j],contextptr),res[k])/dotvecteur(conj(O[j],contextptr),O[j]);
	alpha=_round(alpha,contextptr);
	res[k]=subvecteur(*res[k]._VECTptr,multvecteur(alpha,*res[j]._VECTptr));
	A[k]=subvecteur(*A[k]._VECTptr,multvecteur(alpha,*A[j]._VECTptr));
	L[k]=subvecteur(*L[k]._VECTptr,multvecteur(alpha,*L[j]._VECTptr));
      }
      gen lastalpha=dotvecteur(conj(O[k-1],contextptr),res[k])/dotvecteur(conj(O[k-1],contextptr),O[k-1]);
      if (ck_is_greater(dotvecteur(conj(O[k],contextptr),O[k]),(gen(3)/4-lastalpha*lastalpha)*dotvecteur(conj(O[k-1],contextptr),O[k-1]),contextptr)){
	// Ok, continue the reduction
	++k;
      }
      else {
	swapgen(res[k],res[k-1]);
	swapgen(A[k],A[k-1]);
	--k;
	L.pop_back();
	L.pop_back();
	O.pop_back();
	O.pop_back();
      }
    }
    return res;
  }
  matrice lll(const matrice & m,GIAC_CONTEXT){
    matrice L,O,A;
    return lll(m,L,O,A,contextptr);
  }
  gen _lll(const gen & g,GIAC_CONTEXT){
    if ( g.type==_STRNG && g.subtype==-1) return  g;
    if (g.type!=_VECT)
      return gensizeerr(contextptr);
    matrice L,O,A;
    matrice S=lll(*g._VECTptr,L,O,A,contextptr);
    return gen(makevecteur(S,A,L,O),_SEQ__VECT);
  }
  static const char _lll_s []="lll";
  static define_unary_function_eval (__lll,&_lll,_lll_s);
  define_unary_function_ptr5( at_lll ,alias_at_lll,&__lll,0,true);

  // Utilities for Hermite and Smith normal forms
  static gen rem(const gen & p,const gen & q,environment * env){
    if (!env)
      return smod(p,q);
    if (p.type!=_VECT)
      return zero;
    if (q.type!=_VECT)
      return q;
    return operator_mod(*p._VECTptr,*q._VECTptr,env);
  }

  static gen quo(const gen & p,const gen & q,environment * env){
    if (!env)
      return (p-smod(p,q))/q;
    if (p.type!=_VECT)
      return zero;
    if (q.type!=_VECT)
      return q;
    return operator_div(*p._VECTptr,*q._VECTptr,env);
  }

  static void egcd(const gen & a,const gen & b,gen & u,gen & v,gen & d,environment * env){
    if (!env){
      egcd(a,b,u,v,d);
      return ;
    }
    if (a.type!=_VECT){
      d=a;
      u=plus_one;
      v=zero;
      return;
    }
    if (b.type!=_VECT){
      d=b;
      v=plus_one;
      u=zero;
      return;
    }
    modpoly U,V,D;
    egcd(*a._VECTptr,*b._VECTptr,env,U,V,D);
    u=U; v=V; d=D;
  }

  // degree + 1 for poly, abs for integer, 1 otherwise
  static gen smith_deg(const gen & a,environment * env,GIAC_CONTEXT){
    if (!env)
      return abs(a,contextptr); 
    if (a.type!=_VECT)
      return is_zero(a,contextptr)?zero:plus_one;
    return int(a._VECTptr->size());
  }


  // If Aorig has integer coefficients, hermite
  // finds U and A such that A=U*Aorig with U invertible in Z and A
  // is upper triangular, with non zero coeff || <= |pivot|/2
  bool hermite(const std_matrix<gen> & Aorig,std_matrix<gen> & U,std_matrix<gen> & A,environment * env,GIAC_CONTEXT){
    A=Aorig;
    int n=A.size();
    if (!n) return false;
    int m=A.front().size();
    matrice2std_matrix_gen(midn(n),U);
    gen u,v,d;
    vecteur B1(n),B2(m);
    int i0=0;
    for (int j=0;j<m ;j++ ){
      // Find non zero entry of smallest abs value in column j
      int k=-1;
      gen min_val=plus_inf,tmp,q;
      for (int i=i0;i<n;++i){
	tmp=smith_deg(A[i][j],env,contextptr);
	if (!is_zero(tmp,contextptr) && is_strictly_greater(min_val,tmp,contextptr)){
	  k=i;
	  min_val=tmp;
	}
      }
      if (k>=0 && !is_zero(min_val,contextptr)){
	if (i0!=k){ // Exchange lines i0 and k in A and U
	  swap(A[i0],A[k]);
	  swap(U[i0],U[k]);
	}
	for (int i=n-1;i>=0;--i){
	  if (i==i0 || is_zero(A[i][j],contextptr) )
	    continue;
	  if (i<i0){
	    // Above diag do: L_i <- L_i - q*L_j
	    q=quo(A[i][j],A[i0][j],env);
	    linear_combination(plus_one,U[i],-q,U[i0],plus_one,U[i],0.0,0);
	    linear_combination(plus_one,A[i],-q,A[i0],plus_one,A[i],0.0,0);
	  }
	  else {
	    // Below diag: we use Bezout u*a+v*b=d where a=coeff, b="pivot"
	    // L_i0 <- v*L_i0 + u*L_i
	    // L_i <- (-a * L_i0 + b * L_i)/d
	    // This transformation is Z-invertible since det=(U*a+b*v)/d=1
	    // it will cancel the leading coeff of L_i
	    // We should use the smallest possible |u| and |v|
	    gen a = A[i][j];
	    gen b = A[i0][j];
	    egcd(a,b,u,v,d,env);
	    linear_combination(v,U[i0],u,U[i],plus_one,B1,0.0,0);
	    linear_combination(-a,U[i0],b,U[i],d,U[i],0.0,0);
	    U[i0]=B1;
	    linear_combination(v,A[i0],u,A[i],plus_one,B2,0.0,0);
	    linear_combination(-a,A[i0],b,A[i],d,A[i],0.0,0);
	    A[i0]=B2;	    
	  }
	} // end for (column reduced)
	// cerr << A << endl;
	if (!env && is_strictly_positive(-A[i0][i0],contextptr)){ 
	  A[i0]=-A[i0];
	  U[i0]=-U[i0];
	}
	++i0;
      }
    }
    return true;
  }

  // fonction ihermite
  // Forme normale de Hermite pour une matrice a coeff entiers
  // effectue la reduction sous forme echelonnee (de type Gauss)
  // d'une matrice d'entiers en utilisant uniquement des operations
  // de lignes inversibles dans les entiers, en d'autres termes si A0
  // est la matrice originale, on calcule une matrice U inversible dans Z
  // et une matrice A triangulaire superieure telles que
  //   A = U*A0
  // De plus les coefficients au-dessus de la diagonale de A sont en module
  // inferieurs au pivot de la colonne /2 .
  // exemple
  // A0:=[[9,-36,30], [-36,192,-180], [30,-180,180]];
  // U,A:=ihermite(A0);
  // U*A0-A (renvoie 0)
  // det(U) = 1 donc on passe aussi de A a A0 uniquement avec des
  // manipulations de ligne a coeffs entiers
  // Application: calcul d'une Z-base d'un noyau
  // Soit M la matrice dont on cherche le noyau
  // U,A:=ihermite(transpose(M)) -> A=U*transpose(M)
  // -> transpose(A)=M*transpose(U)
  // les colonnes nulles de transpose(A) correspondent aux colonnes 
  // de transpose(U) dans Ker(M) -> les lignes nulles de A aux lignes de U
  // dans le noyau. 
  // Exemple: M:=[[1,2,3],[4,5,6],[7,8,9]]
  // U,A:=ihermite(M) renvoie
  // [[-3,1,0],[4,-1,0],[-1,2,-1]],[[1,-1,-3],[0,3,6],[0,0,0]]
  // A[2]==0 donc base de Ker(M) composee de U[2], on a bien
  // M*U[2]==0

  bool ihermite(const matrice & Aorig, matrice & U,matrice & A,GIAC_CONTEXT){
    std_matrix<gen> aorig,u,a;
    matrice2std_matrix_gen(Aorig,aorig);
    if (!hermite(aorig,u,a,0,contextptr))
      return false;
    std_matrix_gen2matrice(u,U);
    std_matrix_gen2matrice(a,A);
    return true;
  }

  gen _ihermite(const gen & g,GIAC_CONTEXT){
    if ( g.type==_STRNG && g.subtype==-1) return  g;
    if (g.type!=_VECT)
      return gensizeerr(contextptr);
    if (!is_integer_matrice(*g._VECTptr))
      return gensizeerr(gettext("Integer matrix expected"));
    matrice U,A;
    if (!ihermite(*g._VECTptr,U,A,contextptr))
      return gensizeerr(contextptr);
    return gen(makevecteur(U,A),_SEQ__VECT);
  }
  static const char _ihermite_s []="ihermite";
  static define_unary_function_eval (__ihermite,&_ihermite,_ihermite_s);
  define_unary_function_ptr5( at_ihermite ,alias_at_ihermite,&__ihermite,0,true);

  // A=U*Aorig*V, U and V Z-invertible, A diagonal, A[i,i] divides A[i+1,i+1]
  bool smith(const std_matrix<gen> & Aorig,std_matrix<gen> & U,std_matrix<gen> & A,std_matrix<gen> & V,environment * env,GIAC_CONTEXT){
    A=Aorig;
    int n=A.size();
    if (!n) return false; // setsizeerr();
    int m=A.front().size();
    matrice2std_matrix_gen(midn(n),U);
    matrice2std_matrix_gen(midn(m),V);
    // FIXME: possible improvement if only A is computed
    // do ihermite, compute det, 
    // and make computations below mod 2*det
    // It is also possible at increment step to divide by the pivot
    // the remaining coeffs of the matrix (and multiply back later)
    gen u,v,d;
    vecteur B1(n),B2(m);
    int i0=0,j0=0; // row below i0 and col below j0 done
    for (;j0<m && i0<n; ){
      bool increment=true;
      if (j0<m){
	// Find non zero entry of smallest abs value in column j0
	int k=-1;
	gen min_val=plus_inf,tmp,q;
	for (int i=i0;i<n;++i){
	  tmp=smith_deg(A[i][j0],env,contextptr);
	  if (!is_zero(tmp,contextptr) && is_strictly_greater(min_val,tmp,contextptr)){
	    k=i;
	    min_val=tmp;
	  }
	}
	if (k>=0 && !is_zero(min_val,contextptr)){
	  if (i0!=k){ // Exchange lines i0 and k in A and U
	    swap(A[i0],A[k]);
	    swap(U[i0],U[k]);
	  }
	  for (int i=n-1;i>i0;--i){
	    if (is_zero(A[i][j0],contextptr) )
	      continue;
	    increment=false;
	    // we use Bezout u*a+v*b=d where a=coeff, b="pivot"
	    // L_i0 <- v*L_i0 + u*L_i
	    // L_i <- (-a * L_i0 + b * L_i)/d
	    // This transformation is Z-invertible since det=(U*a+b*v)/d=1
	    // it will cancel the leading coeff of L_i
	    // We should use the smallest possible |u| and |v|
	    gen a = A[i][j0];
	    gen b = A[i0][j0];
	    egcd(a,b,u,v,d,env);
	    linear_combination(v,U[i0],u,U[i],plus_one,B1,0.0,0);
	    linear_combination(-a,U[i0],b,U[i],d,U[i],0.0,0);
	    U[i0]=B1;
	    linear_combination(v,A[i0],u,A[i],plus_one,B2,0.0,0);
	    linear_combination(-a,A[i0],b,A[i],d,A[i],0.0,0);
	    A[i0]=B2;	    
	  } // end for (row reduced)
	  if (!env && is_strictly_positive(-A[i0][j0],contextptr)){
	    A[i0]=-A[i0];
	    U[i0]=-U[i0];
	  }
	} // end if k>=0 && !is_zero(min_val)
      } // end if (j0<m)
      if (i0<n){
	// Column reduction
	A=A.transpose();
	// Find non zero entry of smallest abs value in transposed col i0
	int k=-1;
	gen min_val=plus_inf,tmp,q;
	for (int i=j0;i<m;++i){
	  tmp=smith_deg(A[i][i0],env,contextptr);
	  if (!is_zero(tmp,contextptr) && is_strictly_greater(min_val,tmp,contextptr)){
	    k=i;
	    min_val=tmp;
	  }
	}
	if (k>=0 && !is_zero(min_val,contextptr)){
	  if (j0!=k){ // Exchange transposed rows j0 and k in A and V
	    swap(A[j0],A[k]);
	    swap(V[j0],V[k]);
	  }
	  for (int i=m-1;i>j0;--i){
	    if (is_zero(A[i][i0],contextptr) )
	      continue;
	    increment=false;
	    // we use Bezout u*a+v*b=d where a=coeff, b="pivot"
	    // L_j0 <- v*L_j0 + u*L_i
	    // L_i <- (-a * L_j0 + b * L_i)/d
	    // This transformation is Z-invertible since det=(U*a+b*v)/d=1
	    // it will cancel the leading coeff of L_i
	    // We should use the smallest possible |u| and |v|
	    gen a = A[i][i0];
	    gen b = A[j0][i0];
	    egcd(a,b,u,v,d,env);
	    linear_combination(v,V[j0],u,V[i],plus_one,B2,0.0,0);
	    linear_combination(-a,V[j0],b,V[i],d,V[i],0.0,0);
	    V[j0]=B2;
	    linear_combination(v,A[j0],u,A[i],plus_one,B1,0.0,0);
	    linear_combination(-a,A[j0],b,A[i],d,A[i],0.0,0);
	    A[j0]=B1;	    
	  } // end for (row reduced)
	  if (!env && is_strictly_positive(-A[j0][i0],contextptr)){
	    A[j0]=-A[j0];
	    V[j0]=-V[j0];
	  }
	} // end if (k>=0 && !is_zero(min_val) )
	// End column reduction
	A=A.transpose();
      } // end if (i0<n)
      // Now check that all remaining elements are divisible by A[i0][j0]
      // otherwise replace A[i0] by A[i0]+A[i]
      if (i0<n && j0<m){
	gen pivot=A[i0][j0];
	int i=i0+1;
	for (;i<n;++i){
	  int j=j0+1;
	  for (;j<m;++j){
	    if (!is_zero(rem(A[i][j],pivot,env),contextptr))
	      break;
	  }
	  if (j!=m)
	    break;
	}
	if (i!=n){
	  increment=false;
	  A[i0]=addvecteur(A[i0],A[i]);
	  U[i0]=addvecteur(U[i0],U[i]);
	}
      }
      if (increment){
	++i0;
	++j0;
      }
    } // end for (;j0<m && i0<n;)
    V=V.transpose();
    return true;
  }

  bool ismith(const matrice & Aorig, matrice & U,matrice & A,matrice & V,GIAC_CONTEXT){
    std_matrix<gen> aorig,u,a,v;
    matrice2std_matrix_gen(Aorig,aorig);
    if (!smith(aorig,u,a,v,0,contextptr))
      return false;
    std_matrix_gen2matrice(u,U);
    std_matrix_gen2matrice(a,A);
    std_matrix_gen2matrice(v,V);
    return true;
  }

  gen _ismith(const gen & g,GIAC_CONTEXT){
    if ( g.type==_STRNG && g.subtype==-1) return  g;
    if (g.type!=_VECT)
      return gensizeerr(contextptr);
    if (!is_integer_matrice(*g._VECTptr))
      return gensizeerr(gettext("Integer matrix expected"));
    matrice U,A,V;
    if (!ismith(*g._VECTptr,U,A,V,contextptr))
      return gensizeerr(contextptr);
    return gen(makevecteur(U,A,V),_SEQ__VECT);
  }
  static const char _ismith_s []="ismith";
  static define_unary_function_eval (__ismith,&_ismith,_ismith_s);
  define_unary_function_ptr5( at_ismith ,alias_at_ismith,&__ismith,0,true);

  //   ismith, calcule la forme normale de Smith d'une
  //   matrice, A0 a coefficients entiers
  //   U,A,V := ismith(A0);
  //   calcule U,V Z-inversibles et A=U*Aorig*V, A est diagonale avec
  //   A[i,i] divise A[i+1,i+1]
  //   Les A[i,i] s'appellent diviseurs elementaires et permettent entre
  //   autre de trouver la structure des groupes abeliens de type fini

  // FIXME: Hermite and Smith normal form, same code except for smod/iquo/egcd
  // For polynomials use egcd(a,b,env,u,v,d)

  // Read a CSV file (comma separated) with separator, newline, end of file
  // decsep = decimal separator (, -> .)
  matrice csv2gen(istream & i,char sep,char nl,char decsep,char eof,GIAC_CONTEXT){
    vecteur res,line;
    size_t nrows=0,ncols=0;
    char c;
    string s;
    for (;i;){
      c=i.get();
      if (i.eof() || c==eof)
	break;
      if (c=='%')
	c=' ';
      if (c==sep || c==nl){
	// remove spaces at beginning of s
	while (!s.empty() && s[0]==' ')
	  s=s.substr(1,s.size()-1);
	// if sep==' ' remove spaces in i
	if (sep==' '){
	  char c2;
	  for (;;){
	    c2=i.get();
	    if (i.eof() || c2!=' '){
	      i.putback(c2);
	      break;
	    }
	  }
	}
	// if 1st char is = or digit parse, else string
	int ss=s.size();
	if (s.empty())
	  line.push_back(string2gen(s,false));
	else {
	  if (ss>2 && s[0]=='"' && s[1]=='=' && s[ss-1]=='"'){
	    s=s.substr(1,ss-2);
	    ss -= 2;
	  }
	  if (s[0]=='=' || s[0]=='-'){
	    line.push_back(gen(s,contextptr));
	  }
	  else {
	    if (s[0]>='0' && s[0]<='9'){
	      line.push_back(gen(s,contextptr));
	    }
	    else
	      line.push_back(string2gen(s,s[0]=='"'));
	  }
	}
	s="";
	if (c==nl){
	  res.push_back(line);
	  ncols=max(ncols,line.size());
	  line.clear();
	  nrows++;
	  continue;
	}
      } // end if c==sep || nl
      else  {
	if (c==decsep)
	  s += '.';
	else
	  s += c;
      }
    } // end reading stream
    // now make a matrix from res
    for (unsigned j=0;j<nrows;j++){
      res[j]=mergevecteur(*res[j]._VECTptr,vecteur(ncols-res[j]._VECTptr->size(),0));
    }
    return res;
  }
  gen _csv2gen(const gen & g,GIAC_CONTEXT){
    if ( g.type==_STRNG && g.subtype==-1) return  g;
    char sep(';'),nl('\n'),eof(0),decsep(',');
    gen tmp,gs;
    if (g.type==_VECT && !g._VECTptr->empty()){
      gs=g._VECTptr->front();
      int s=g._VECTptr->size();
      if (s>1){
	tmp=g[1];
	if (tmp.type==_STRNG && !tmp._STRNGptr->empty())
	  sep=(*tmp._STRNGptr)[0];
      }
      if (s>2){
	tmp=g[2];
	if (tmp.type==_STRNG && !tmp._STRNGptr->empty())
	  nl=(*tmp._STRNGptr)[0];
      }
      if (s>3){
	tmp=g[1];
	if (tmp.type==_STRNG && !tmp._STRNGptr->empty())
	  decsep=(*tmp._STRNGptr)[0];
      }
      if (s>4){
	tmp=g[4];
	if (tmp.type==_STRNG && !tmp._STRNGptr->empty())
	  eof=(*tmp._STRNGptr)[0];
      }
    }
    else
      gs=g;
    if (gs.type!=_STRNG)
      return gensizeerr(gettext("Expecting file name to convert"));
    string file=*gs._STRNGptr;
    ifstream i(file.c_str());
    return csv2gen(i,sep,nl,decsep,eof,contextptr);
  }
  static const char _csv2gen_s []="csv2gen";
  static define_unary_function_eval (__csv2gen,&_csv2gen,_csv2gen_s);
  define_unary_function_ptr5( at_csv2gen ,alias_at_csv2gen,&__csv2gen,0,true);


  matrice matpow(const matrice & m,const gen & n,GIAC_CONTEXT){
    identificateur x("x");
    gen ux=symbolic(at_pow,gen(makevecteur(x,n),_SEQ__VECT));
    return analytic_apply(ux,x,m,contextptr);
  }

      // FIXME: pow should not always call egv stuff
  gen _matpow(const gen & a,GIAC_CONTEXT){
    if ( a.type==_STRNG && a.subtype==-1) return  a;
    if (a.type==_VECT && a._VECTptr->size()==2 && ckmatrix(a._VECTptr->front()))
      return matpow(*a._VECTptr->front()._VECTptr,a._VECTptr->back(),contextptr);
    return gensizeerr(contextptr);
  }
  static const char _matpow_s []="matpow";
  static define_unary_function_eval (__matpow,&giac::_matpow,_matpow_s);
  define_unary_function_ptr5( at_matpow ,alias_at_matpow,&__matpow,0,true);

  // EIGENVALUES for double coeff
  bool eigenval2(matrix_double & H,int n2,giac_double & l1, giac_double & l2){
    giac_double a=H[n2-2][n2-2],b=H[n2-2][n2-1],c=H[n2-1][n2-2],d=H[n2-1][n2-1];
    giac_double delta=a*a-2*a*d+d*d+4*b*c;
    if (delta<0){
      l1=(a+d)/2;
      l2=std::sqrt(-delta)/2;
      return false;
    }
    delta=std::sqrt(delta);
    l1=(a+d+delta)/2;
    l2=(a+d-delta)/2;
    return true;
  }

  static inline giac_double trim(giac_double a,giac_double b,giac_double eps){
    if (eps && std::abs(a)<eps*std::abs(b)) 
      return 0;
    else
      return a;
  }

  // v=(c1*v1+c2*v2), begin at cstart
  void linear_combination(giac_double c1,const vector<giac_double> & v1,giac_double c2,const vector<giac_double> & v2,vector<giac_double> & v,int cstart,double eps){
    eps=0;
    if (cstart<0)
      cstart=0;
    vector<giac_double>::const_iterator it1=v1.begin()+cstart,it1end=v1.end(),it2=v2.begin()+cstart;
    vector<giac_double>::iterator jt1=v.begin()+cstart;
#ifdef DEBUG_SUPPORT
    if (it1end-it1!=v2.end()-it2)
      setdimerr();
#endif
    if (it2==jt1){
      linear_combination(c2,v2,c1,v1,v,cstart,eps);
      return;
    }
    if (it1==jt1){
      for (;jt1!=it1end;++jt1,++it2){
	*jt1=c1*(*jt1)+c2*(*it2); // trim(c1*(*jt1)+c2*(*it2),c1,eps);
      }
      return;
    }
    if (v.size()==it1end-it1){
      jt1=v.begin();
      for (int i=0;i<cstart;++i,++jt1)
	*jt1=0;
      for (;it1!=it1end;++it1,++it2,++jt1)
	*jt1=c1*(*it1)+c2*(*it2); // trim(c1*(*it1)+c2*(*it2),c1,eps);
      return;
    }
    v.clear();
    v.reserve(it1end-it1);
    for (int i=0;i<cstart;++i)
      v.push_back(0);
    for (;it1!=it1end;++it1,++it2)
      v.push_back(c1*(*it1)+c2*(*it2)); // trim(c1*(*it1)+c2*(*it2),c1,eps);
  }

  void exchange(matrix_double & H,matrix_double &P,bool compute_P,int i,int m1){
    int n=H.size();
    H[i].swap(H[m1]);
    if (compute_P)
      P[i].swap(P[m1]);
    for (int j=0;j<n;++j){
      vector<giac_double> & Hj=H[j];
      swap<giac_double>(Hj[i],Hj[m1]);
    }
  }

  // Hessenberg reduction, P is orthogonal and should be initialized to identity
  // trn(P)*H*P=original
  // already_zero is either <=0 or an integer such that H[i][j]==0 if i>j+already_zero
  // (already_zero==1 if H is hessenberg, ==3 for Francis algorithm)
  void hessenberg_ortho(matrix_double & H,matrix_double & P,int firstrow,int n,bool compute_P,int already_zero){
    double eps=1e-10;
    matrix_double::iterator Hbegin=H.begin();
    vector<giac_double>::iterator Hiterator,Hjmptr,Hjiptr;
    int nH=H.size();
    if (n<0 || n>nH) 
      n=nH;
    if (firstrow<0 || firstrow>n)
      firstrow=0;
    giac_double t,tn,tc,tabs,u,un,uc,tmp1,tmp2,norme;
    vector<giac_double> v1(nH),v2(nH),TN(n,1),UN(n);
    for (int m=firstrow;m<n-2;++m){
      if (debug_infolevel>=2)
	cerr << "// hessenberg reduction line " << m << endl;
      // check for a non zero coeff in the column m below ligne m+1
      int i=m+1;
      giac_double pivot=0;
      int pivotline=0;
      for (;i<n;++i){
	t=H[i][m];
	tabs=std::abs(t);
	if (tabs>pivot){
	  pivotline=i;
	  pivot=tabs;
	}
      }
      if (std::abs(pivot)<eps) //no pivot found
	continue;
      // if (std::abs(pivot/H[m][m])<eps) //not found
      //	continue;
      // if (std::abs(pivot)<1e-8) cerr << pivot << " " << H[m][m] << endl;
      i=pivotline;
      // exchange line and columns
      if (i>m+1)
	exchange(H,P,compute_P,i,m+1);
      // now coeff at line m+1 column m is H[m+1][m]=t!=0
      // creation of zeros in column m+1, lines i=m+2 and below
      for (i=m+2;i<n;++i){
	// line operation
	t=H[m+1][m];
	u=H[i][m];
	uc=u;
	tc=t;
	norme=std::sqrt(u*uc+t*tc);
	un=u/norme; tn=t/norme; uc=un; tc=tn;
	if (std::abs(un)<eps){
	  UN[i]=0;
	  continue;
	}
	if (debug_infolevel>=2)
	  cerr << "// i=" << i << " " << u <<endl;
	// H[m+1]=un*H[i]+tn*H[m+1] and H[i]=tn*H[i]-un*H[m+1];
#if 0
	if (m){
	  for (int k=0;k<m;++k){
	    if (std::abs(H[i][k])>eps || std::abs(H[m+1][k])>eps)
	      cerr << H[i] << " " << H[m+1] << " " << v1 << " " << v2 << endl;
	  }
	}
#endif
	linear_combination(uc,H[i],tc,H[m+1],v1,0,1e-12); // m or 0
	linear_combination(tn,H[i],-un,H[m+1],v2,0,1e-12); 
	H[m+1].swap(v1);
	H[i].swap(v2);
	if (compute_P){
	  linear_combination(uc,P[i],tc,P[m+1],v1,0,1e-12); 
	  linear_combination(tn,P[i],-un,P[m+1],v2,0,1e-12); 
	  P[m+1].swap(v1);
	  P[i].swap(v2);
	}
	TN[i]=tn;
	UN[i]=un;
      }
      for (i=m+2;i<n;++i){      
	if (!(un=UN[i]))
	  continue;
	tn=TN[i];
	tc=tn;
	uc=un;
	// column operation
	for (int j=0;j<nH;++j){
#if 1
	  vector<giac_double> & Hj=H[j];
	  giac_double & Hjm=Hj[m+1];
	  giac_double & Hji=Hj[i];
	  tmp1=tn*Hjm+un*Hji;
	  Hji=-uc*Hjm+tc*Hji;
	  Hjm=tmp1;
#else
	  Hiterator=(Hbegin+j)->begin();
	  Hjmptr=Hiterator+m+1;
	  Hjiptr=Hiterator+i;
	  tmp1=tn*(*Hjmptr)+un*(*Hjiptr);
	  *Hjiptr=-uc*(*Hjmptr)+tc*(*Hjiptr);
	  *Hjmptr=tmp1;	  
#endif
	}
      }
    }
    // make 0 below subdiagonal (i<nH all matrix, i<n only relevant lines/column)
    for (int i=2;i<n;i++){
      vector<giac_double>::iterator it=H[i].begin(),itend=it+i-1; // or min(i-1,n);
      for (;it!=itend;++it){
	if (debug_infolevel>2 && std::abs(*it>1e-10))
	  cerr << "Precision " << i << " " << *it << endl;
	*it=0;
      }
    }
  }

  void tri_linear_combination(const giac_double & c1,const vector<giac_double> & x1,const giac_double & c2,const vector<giac_double> & x2,const giac_double & c3,const vector<giac_double> & x3,vector<giac_double> & y){
    vector<giac_double>::const_iterator it1=x1.begin(),it2=x2.begin(),it3=x3.begin(),it3end=x3.end();
    vector<giac_double>::iterator jt=y.begin();
    for (;it3!=it3end;++jt,++it1,++it2,++it3){
      *jt=c1*(*it1)+c2*(*it2)+c3*(*it3);
    }
  }

  ostream & operator << (ostream & os,const vector<giac_double> & m){
    int s=m.size();
    for (int i=0;i<s;++i)
      os << m[i] << " ";
    return os;
  }

  ostream & operator << (ostream & os,const matrix_double & m){
    int s=m.size();
    for (int i=0;i<s;++i)
      os << m[i] << endl;
    return os;
  }

  void matrix_double::dbgprint() const { std::cout << *this << std::endl; }

  void francis_iterate1(matrix_double & H,int n1,int n2,matrix_double & P,int maxiter,double eps,bool is_hessenberg,bool complex_schur,bool compute_P,vector<giac_double> & v1,vector<giac_double> & v2,giac_double l1){
    int n_orig=H.size();
    // [[a,b],[c,d]] -> [b,l1-a] or [l1-d,c] as first eigenvector
    giac_double a=H[n2-2][n2-2],b=H[n2-2][n2-1],c=H[n2-1][n2-2],d=H[n2-1][n2-1];
    giac_double l1a=l1-a,l1d=l1-d,x,y,xy;
    if (std::abs(l1a)>std::abs(l1d)){
      x=b; y=l1a;
    }
    else {
      x=l1d; y=c;
    }
    xy=std::sqrt(x*x+y*y);
    // normalize eigenvector
    x = x/xy; y = y/xy;
    if (eps>std::abs(y))
      return;
    if (std::abs(x-1)<0.5)
      x=-y*y/(x+1);
    else
      x-=1;
    xy=std::sqrt(x*x+y*y);
    x = x/xy; y=y/xy;
    // compute reflection matrix let n=[[x],[y]] trn(n)=conj([[x,y]])
    // Q=idn(2)-2*n*trn(n);
    // i.e. [[ 1-2x*conj(x), -2x*conj(y)   ]   -2x*conj(z)  ],
    //       [ -2*y*conj(x), 1-2*y*conj(y) ]]
    // apply Q on H and P: line operations on H and P
    giac_double c11=1-2*x*x, c12=-2*x*y,
      c21=c12, c22=1-2*y*y,tmp1,tmp2;
    linear_combination(c11,H[n1],c12,H[n1+1],v1,0,eps);
    linear_combination(c21,H[n1],c22,H[n1+1],v2,0,eps);
    H[n1].swap(v1);
    H[n1+1].swap(v2);
    if (compute_P){
      linear_combination(c11,P[n1],c12,P[n1+1],v1,0,eps);
      linear_combination(c21,P[n1],c22,P[n1+1],v2,0,eps);
      P[n1].swap(v1);
      P[n1+1].swap(v2);
    }
    // now columns operations on H (not on P)
    for (int j=0;j<n_orig;++j){
      vector<giac_double> & Hj=H[j];
      giac_double & Hjm1=Hj[n1];
      giac_double & Hjm2=Hj[n1+1];
      tmp1=Hjm1*c11+Hjm2*c21;
      tmp2=Hjm1*c12+Hjm2*c22;
      Hjm1=tmp1;
      Hjm2=tmp2;
    }
    hessenberg_ortho(H,P,n1,n2,compute_P,2); 
  }

  void francis_iterate2(matrix_double & H,int n1,int n2,matrix_double & P,int maxiter,double eps,bool is_hessenberg,bool complex_schur,bool compute_P,vector<giac_double> & v1,vector<giac_double> & v2,vector<giac_double> & v3){
    giac_double tmp1,tmp2,tmp3;
    int n_orig(H.size());
    // now H is proper hessenberg (indices n1 to n2-1)
    // find eigenvalues l1 and l2 of last 2x2 matrix, they will be taken as shfits
    giac_double s,p; // s=l1+l2, p=l1*l2
    s=H[n2-2][n2-2]+H[n2-1][n2-1];
    p=H[n2-2][n2-2]*H[n2-1][n2-1]-H[n2-1][n2-2]*H[n2-2][n2-1];
    // compute (H-l2)(H-l1)=(H-s)*H+p on n1-th basis vector (if n1==0, on [1,0,...,0])
    giac_double ha=H[n1][n1],hb=H[n1][n1+1],
      hd=H[n1+1][n1],he=H[n1+1][n1+1],
      hh=H[n1+2][n1+1];
    giac_double x=hb*hd+ha*(ha-s)+p,y=hd*(he-s+ha),z=hd*hh;
    // normalize, substract [1,0,0] and normalize again
    giac_double xyz=std::sqrt(x*x+y*y+z*z);
    // if x/xyz is near 1, improve precision:
    // x/xyz-1 = ((x/xyz)^2-1)/(x/xyz+1)=-((y/xyz)^2+(z/xyz)^2)/(x/xyz+1)
    x=x/xyz; y=y/xyz; z=z/xyz;
    if (std::abs(x-1)<0.5)
      x=-(y*y+z*z)/(x+1);
    else
      x-=1;
    xyz=std::sqrt(x*x+y*y+z*z);
    x=x/xyz; y=y/xyz; z=z/xyz;
    // compute reflection matrix let n=[[x],[y],[z]] trn(n)=conj([[x,y,z]])
    // Q=idn(3)-2*n*trn(n);
    // i.e. [[ 1-2x*conj(x), -2x*conj(y),   -2x*conj(z)  ],
    //       [ -2*y*conj(x), 1-2*y*conj(y), -2*y*conj(z) ],
    //       [ -2*z*conj(x), -2*z*conj(y),  1-2*z*conj(z)]]
    // apply Q on H and P: line operations on H and P
    giac_double c11=1-2*x*x,c12=-2*x*y,c13=-2*x*z;
    giac_double c21=c12,c22=1-2*y*y,c23=-2*y*z;
    giac_double c31=c13,c32=c23,c33=1-2*z*z;
    // cerr << "[[" << c11 <<"," << c12 << "," << c13 << "],[" <<  c21 <<"," << c22 << "," << c23 << "],[" << c31 <<"," << c32 << "," << c33 << "]]" << endl;
    tri_linear_combination(c11,H[n1],c12,H[n1+1],c13,H[n1+2],v1);
    tri_linear_combination(c21,H[n1],c22,H[n1+1],c23,H[n1+2],v2);
    tri_linear_combination(c31,H[n1],c32,H[n1+1],c33,H[n1+2],v3);
    H[n1].swap(v1);
    H[n1+1].swap(v2);
    H[n1+2].swap(v3);
    // cerr << H << endl;
    if (compute_P){
      tri_linear_combination(c11,P[n1],c12,P[n1+1],c13,P[n1+2],v1);
      tri_linear_combination(c21,P[n1],c22,P[n1+1],c23,P[n1+2],v2);
      tri_linear_combination(c31,P[n1],c32,P[n1+1],c33,P[n1+2],v3);
      P[n1].swap(v1);
      P[n1+1].swap(v2);
      P[n1+2].swap(v3);
    }
    // now columns operations on H (not on P)
    for (int j=0;j<n_orig;++j){
      vector<giac_double> & Hj=H[j];
      giac_double & Hjm1=Hj[n1];
      giac_double & Hjm2=Hj[n1+1];
      giac_double & Hjm3=Hj[n1+2];
      tmp1=Hjm1*c11+Hjm2*c21+Hjm3*c31;
      tmp2=Hjm1*c12+Hjm2*c22+Hjm3*c32;
      tmp3=Hjm1*c13+Hjm2*c23+Hjm3*c33;
      Hjm1=tmp1;
      Hjm2=tmp2;
      Hjm3=tmp3;
    }
    // cerr << H << endl;
    // chase the bulge: Hessenberg reduction on 2 subdiagonals
    hessenberg_ortho(H,P,n1,n2,compute_P,3); // <- improve
  }


  // Francis algorithm on submatrix rows and columns n1..n2-1
  // Invariant: trn(P)*H*P=orig matrix, complex_schur not used for giac_double coeffs
  bool francis_schur(matrix_double & H,int n1,int n2,matrix_double & P,int maxiter,double eps,bool is_hessenberg,bool complex_schur,bool compute_P){
    vecteur eigenv;
    if (n1==0 && !complex_schur && lapack_schur(H,P,compute_P,eigenv))
      return true;
#ifdef VISUALC // tested on 100 rand poly of degree 10, should work...
    // return false; 
#endif
    int n_orig=H.size(),nitershift0=0;
    if (!is_hessenberg){
      hessenberg_ortho(H,P,0,n_orig,compute_P,0); // insure Hessenberg form (on the whole matrix)
    }
    if (n2-n1<=1)
      return true; // nothing to do
    vector<giac_double> v1(n_orig),v2(n_orig),v3(n_orig);
    if (n2-n1==2){ // 2x2 submatrix, we know how to diagonalize
      giac_double l1,l2;
      if (eigenval2(H,n2,l1,l2)){
	francis_iterate1(H,n1,n2,P,maxiter,eps,true,complex_schur,compute_P,v1,v2,l1);
      }
      return true;
    }
    for (int niter=0;n2-n1>2 && niter<maxiter;niter++){
      if (debug_infolevel>=2){
	cerr << "// qr iteration number " << niter << " " << endl;
	H.dbgprint();
      }
      // check if one subdiagonal element is sufficiently small, if so 
      // we can increase n1 or decrease n2 or split
      giac_double ratio,coeff=1;
      if (niter>maxiter-3)
	coeff=100;
      for (int i=n1;i<n2-1;++i){
	ratio=std::abs(H[i+1][i]/H[i][i]);
	if (ratio<coeff*eps){ 
	  // do a final iteration
	  // francis_iterate2(H,n1,n2,P,maxiter,eps,true,complex_schur,compute_P,v1,v2,v3);
	  // submatrices n1..i and i+1..n2-1
	  if (debug_infolevel>2)
	    cerr << "Francis split double " << n1 << " " << i+1 << " " << n2 << endl;
	  if (!francis_schur(H,n1,i+1,P,maxiter,eps,true,complex_schur,compute_P)){
	    francis_schur(H,i+1,n2,P,maxiter,eps,true,complex_schur,compute_P);
	    return false;
	  }
	  return francis_schur(H,i+1,n2,P,maxiter,eps,true,complex_schur,compute_P);
	}
	if (i==n1 && ratio<std::sqrt(eps)){
	  // exchange lines/columns n1/n2-1
	  // exchange(H,P,compute_P,n1,n2-1);
	  // break;
	}
	// IMPROVE: in that case we should iterate_n using the eigenvalues of the
	// submatrix i+1..n2-1
      }
      if (ratio<std::sqrt(eps)){
	francis_iterate1(H,n1,n2,P,maxiter,eps,true,complex_schur,compute_P,v1,v2,H[n2-1][n2-1]);
      }
      francis_iterate2(H,n1,n2,P,maxiter,eps,true,complex_schur,compute_P,v1,v2,v3);
    } // end for loop on niter
    return false;
  }


#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
