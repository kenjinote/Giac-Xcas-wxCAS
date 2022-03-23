// -*- compile-command: "g++ -I.. -I. -g -c Editeur.cc " -*-
#ifndef _EDITEUR_H
#define _EDITEUR_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef IN_GIAC
#include <giac/first.h>
#else
#include "first.h"
#endif
/*
 *  Copyright (C) 2002 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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

#ifndef IN_GIAC
#include <giac/gen.h>
#else
#include "gen.h"
#endif
#include "Input.h"
#include "Graph.h"
#ifdef HAVE_LIBFLTK
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Value_Input.H>
#include <FL/fl_draw.H>
#include <FL/fl_message.H>
#endif

#ifndef NO_NAMESPACE_XCAS
namespace xcas {
#endif // ndef NO_NAMESPACE_XCAS

#ifdef HAVE_LIBFLTK

  extern Fl_Menu_Item Editeur_menu[] ;
  extern Fl_Text_Display::Style_Table_Entry styletable[];
  extern int styletable_n;

  extern void (*alt_ctrl_cb)(int);
  extern int alt_ctrl;

  class Save_Focus_Button;

  class Xcas_Text_Editor:public Fl_Text_Editor {
  public:
    Fl_Text_Buffer * stylebuf;
    Xcas_Text_Editor(int X, int Y, int W, int H, Fl_Text_Buffer * b,const char* l = 0);
    virtual ~Xcas_Text_Editor() { if (Xcas_input_focus==this) Xcas_input_focus=0;}
    virtual int handle(int event);
    virtual void draw();
    int indent(int pos); // indent current line, return new cursor position
    void indent(); // indent all
    void match();
  };
  
  class Editeur :public Fl_Group {
  public:
    Log_Output * log ;
    Xcas_Text_Editor * editor;
    Fl_Output * output;
    Fl_Menu_Bar * menubar;
    Fl_Button * button,*save_button,*nxt_button;
    Save_Focus_Button * exec_button;
    Fl_Value_Input * linenumber;
    giac::context * contextptr;
    std::string extension;
    std::string search;
    Editeur(int x,int y,int w,int h,const char * l=0);
    // ~Editeur(){ delete editor->buffer(); delete editor; delete output; }
    void position(int & i,int & j) const ; // returns begin and end selection or pos
    std::string value() const ;
    bool eval() ; // parse and eval
  };

  std::string editeur_load(Fl_Text_Editor * e);

  // A button that does not take focus, useful for script exec
  class Save_Focus_Button:public Fl_Button {
  public:
    static Fl_Widget * widget;
    virtual int handle(int event);
    Save_Focus_Button(int X,int Y,int W,int H,const char * l=0):Fl_Button(X,Y,W,H,l){};
  };

  class No_Focus_Button:public Fl_Button {
  public:
    virtual int handle(int event);
    No_Focus_Button(int X,int Y,int W,int H,const char * l=0):Fl_Button(X,Y,W,H,l){};
  };

  // Exported callback for Xcas1.cc for tortue loading
  void cb_Editeur_Exec_All(Fl_Widget * m , void*) ;
  // exported callback for hist.fl
    
  // void Xcas_input_1arg(Save_Focus_Button * m , void*) ;
  void Xcas_input_1arg(No_Focus_Button * m , void*) ;

  void in_Xcas_input_1arg(Fl_Widget * widget,const std::string & chaine ,bool eval=false);
  // void Xcas_input_arg(Save_Focus_Button * m , const std::string & chaine) ;
  void Xcas_input_arg(No_Focus_Button * m , const std::string & chaine);

  void in_Xcas_input_char(Fl_Widget * widget,const std::string & chaine,char keysym);
  // void Xcas_input_char(Save_Focus_Button * m , void*) ;
  void Xcas_input_char(No_Focus_Button * m , void*) ;

  // void Xcas_binary_op(Save_Focus_Button * m , void*) ;
  void Xcas_binary_op(No_Focus_Button * m , void*) ;
  // void Xcas_input_0arg(Save_Focus_Button * m , const std::string & chaine) ;
  void Xcas_input_0arg(No_Focus_Button * m , const std::string & chaine) ;

#endif // FLTK

#ifndef NO_NAMESPACE_XCAS
} // namespace xcas
#endif // ndef NO_NAMESPACE_XCAS

#endif // _EDITEUR_H
