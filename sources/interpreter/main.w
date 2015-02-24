% Copyright (C) 2006-2012 Marc van Leeuwen
% This file is part of the Atlas of Lie Groups and Representations (the Atlas)

% This program is made available under the terms stated in the GNU
% General Public License (GPL), see http://www.gnu.org/licences/licence.html

% The Atlas is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

% The Atlas is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with the Atlas; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


\def\point{\item{$\bullet$}}

@* Introduction to the {\tt realex} interpreter.
%
While the Atlas software initially consisted of a single executable
program \.{atlas} that essentially allows running each of the functionalities
provided by the software library separately, a second executable program
called \.{realex} has been developed since the year~$2006$, which program also
gives access to the functionalities provided by the software library, but in a
manner that provides the user with means to capture the results from previous
computations and use them in subsequent computations. It started out as an
interface giving the user just assignable variables and an expression language
for calling library functions, but has grown out to a complete programming
language, which at the time of writing is still being extended in many
directions. The name of the program might stand for Redesigned
Expression-based Atlas of Lie groups EXecutable, until a better explanation
for the name comes along.

The sources specific to the interpreter are situated in the
current (\.{sources/interpreter/}) subdirectory, and almost all its names are
placed in the |atlas::interpreter| namespace. This part of the program has
been split into several modules, only some of which are centred around
specific data structures; in general these modules are more interdependent
than those of the main Atlas library, and the subdivision is more arbitrary
and subject to change. We list here the source files for these modules, and
their dependencies notably on the level of their header files.

\point The file \.{buffer.w} defines the classes |BufferedInput| providing an
interface to input streams, and |Hash_table| storing and providing a
translation from identifiers to small integers |id_type|.

\point The file \.{parsetree.w} defines the |expr| structure representing parsed
expressions, and many related types, and node constructing functions for use
by the parser. Its header file includes \.{buffer.h}, and declares some
pointer types to types that will be defined in \.{types.h}.

\point The file \.{parser.y} is the source for \.{bison}-generated the parser
file \.{parser.tab.c}. The header file \.{parser.tab.h} includes nothing, but
contains definitions depending on \.{parse\_types.h} having been included.

\point The file \.{lexer.w} defines the lexical analyser class
|Lexical_analyser| and the readline completion function |id_completion_func|.
Its header file includes \.{buffer.h}, \.{parse\_types.h} and \.{parser.tab.h}.

\point The file \.{types.w} defines the main base classes for the evaluator,
|type_expr| for representing (realex) types, |value_base| for dynamic values,
|context| for dynamic evaluation contexts, |expression_base| for ``compiled''
expressions, |program_error| for exceptions, and numerous types related to
these. Its header file includes \.{parsetree.h} (which is needed for the error
classes only).

\point The file \.{built-in-types.w} defines types primitive to realex which
encapsulate types of the Atlas library, and their interface functions. Its
header file includes \.{types.h} and many headers from the Atlas library.

\point The file \.{global.w} other primitive types less related to Atlas, like
integers, rationals, matrices. Also some global aspects of the interpreter
like operating the global identifier tables. Its
header file includes \.{types.h} and some headers from the Atlas library.

\point The file \.{evaluator.w} Defines the realex type-checker and the
evaluator proper. It defines many classes derived from |expression_base|. Its
header file includes \.{lexer.h} and \.{global.h}.

\point The file \.{main.w}, the current module, brings everything together and
defining the main program. It has no header file.


@* Main program. This file defines a small main program to test the parser
under development. It is written in \Cpp, but is it mainly concerned with
interfacing to the parser that is generated by~\.{bison}, and which is
therefore written in~\Cee. However, we now compile the generated parser file
\.{parser.tab.c} using a \Cpp\ compiler, which means we can write our code
here as an ordinary \Cpp\ program, including use of namespaces.

Since depending on the readline libraries still gives difficulties on some
platforms, we arrange for the possibility of compiling this program in the
absence of that library. This means that we should refrain from any reference
to its header files, and so the corresponding \&{\#include} statements cannot
be given in the usual way, which would cause their inclusion unconditionally.
Like for the \.{atlas} program, the compile time flag |NREADLINE|, if defined
by setting \.{-DNREADLINE} as a flag to the compiler, will prevent any
dependency on the readline library.

@d realex_version "0.8.5"
 // numbering from 0.5 (on 27/11/2010); last change December 1, 2014

@c

@< Conditionally include the header files for the readline library @>

@< Declaration of interface to the parser @>@;
namespace { @< Local static data @>@; }@;
@< Definitions of local functions @>@;
@< Main program @>

@ Since the file \.{parser.y} declares \.{\%pure-parser} and \.{\%locations},
the prototype of the lexical analyser (wrapper) function |yylex| is the one
below. Curiously, the program~\.{bison} does not write this prototype to
\.{parser.tab.h}, but it does write the definitions of the types |YYSTYPE| and
|YYLTYPE| there; these require that \.{parse\_types.h} be included first. We
also declare ``{\tt\%parse-param \char`\{} |int* verbosity, expr_p*
parsed_expr@;| {\tt\char`\}}'' in~\.{parser.y}, so that the parser itself,
|yyparse|, takes an integer pointer as parameter, which it uses to signal
special requests from the user (such as verbose output but also termination or
output redirection), and a pointer to an expression, in which it writes the
result of parsing.

The definitions below used to start with |extern "C"|, but no longer do so
since the parser is now compiled as a \Cpp\ program.

@h "parse_types.h"
@h "parser.tab.h"

@< Declaration of interface to the parser @>=

int yylex (YYSTYPE *, YYLTYPE *);
@/int yyparse( atlas::interpreter::expr_p* parsed_expr, int* verbosity );

@ Here is an array that declares the keywords that the lexical scanner is to
recognise, terminated by a null pointer. Currently the lexical analyser adds
the offset of the keyword in this list to |QUIT|, so the recognition depends
on the fact that |"quit"| is the first keyword, and that they are listed below
in the same order as in the \.{\%token} declarations in \.{parser.y}.

@< Local static data @>=

const char* keywords[] =
 {"quit"
 ,"set","let","in","begin","end"
 ,"if","then","else","elif","fi"
 ,"and","or","not"
 ,"while","do","od","next","for","from","downto"
 ,"true","false"
 ,"quiet","verbose"
 ,"whattype","showall","forget"
 ,nullptr};

@ Here are the wrapper function for the lexical analyser and the error
reporting function, which are necessary because the parser cannot directly
call a class method. The prototypes are imposed, in particular the second and
third arguments to |yyerror| are those passed to |yyparse|, even though they
are not used in |yyerror|. In |yyerror| we close any open include files, as
continuing to execute their commands is undesirable.

@< Definitions of local functions @>=

int yylex(YYSTYPE *valp, YYLTYPE *locp)
{@; return atlas::interpreter::lex->get_token(valp,locp); }
@)

void yyerror (YYLTYPE* locp, atlas::interpreter::expr_p* ,int* ,char const *s)
{ atlas::interpreter::main_input_buffer->show_range@|
  (std::cerr,
   locp->first_line, locp->first_column,
   locp->last_line,  locp->last_column);
  std::cerr << s << std::endl;
  atlas::interpreter::main_input_buffer->close_includes();
}

@ Here are some header files which need to be included for this main program.
As we discussed above, the inclusion of header files for the readline
libraries is made dependent on the flag |NREADLINE|. In case the flag is set,
we define the few symbols used from the readline library as macros, so that
the code using them can be compiled without needing additional \&{\#ifdef}
lines. It turns out that |getline| and |add_history| are not used in calls,
but rather passed as function pointers to the |BufferedInput| constructor, so
the appropriate expansion for these macros is the null pointer.

@h <iostream>
@h <fstream>

@h "buffer.h"
@h "lexer.h"
@h "version.h"

@< Conditionally include the header files for the readline library @>=
#ifdef NREADLINE
#define readline nullptr
#define add_history nullptr
#define clear_history()
#else
#include <readline/readline.h>
#include <readline/history.h>
#endif

@ After a basic initialisation, our main program constructs unique instances
for various classes of the interpreter, and sets pointers to them so that
various compilation units can access them. Then in a loop it calls the parser
until it sets |verbosity<0|, which is done upon seeing the \.{quit} command.
We call the |reset| method of the lexical scanner before calling the parser,
which will discard any input that is left by a possible previous erroneous
input. This also already fetches a new line of input, or abandons the program
in case none can be obtained.

@h <unistd.h> // for |isatty|
@< Main program @>=

int main(int argc, char** argv)
{ using namespace std; using namespace atlas::interpreter;
@)
  @< Handle command line arguments @>

@/BufferedInput input_buffer(isatty(STDIN_FILENO) ? "expr> " : nullptr
                            ,use_readline ? readline : nullptr
			    ,use_readline ? add_history : nullptr);
  main_input_buffer= &input_buffer;
@/Hash_table hash; main_hash_table= &hash;
@/Lexical_analyser ana(input_buffer,hash,keywords,prim_names); lex=&ana;
@/Id_table main_table; @+ global_id_table=&main_table;
@/overload_table main_overload_table;
 @+ global_overload_table=&main_overload_table;
@)
  @< Initialise various parts of the program @>
@)
  cout << "This is 'realex', version " realex_version " (compiled on " @|
       << atlas::version::COMPILEDATE @| <<
").\nIt is the programmable interpreter interface to the library (version " @|
       << atlas::version::VERSION @| << ") of\n"
       << atlas::version::NAME << @| ". http://www.liegroups.org/\n";
@)
  last_value = shared_value (new tuple_value(0));
  last_type = void_type.copy();
   // |last_type| is a |type_ptr| defined in \.{evaluator.w}
  while (ana.reset()) // get a fresh line for lexical analyser, or quit
  { expr_p parse_tree;
    int old_verbosity=verbosity;
    ofstream redirect; // if opened, this will be closed at end of loop
    if (yyparse(&parse_tree,&verbosity)!=0)
      continue; // syntax error or non-expression
    if (verbosity!=0) // then some special action was requested
    { if (verbosity<0)
        break; // \.{quit} command
      if (verbosity==2 or verbosity==3)
        // indicates output redirection was requested
      { @< Open |redirect| to specified file, and if successful make
        |output_stream| point to it; otherwise |continue| @>
        verbosity=old_verbosity; // verbosity change was temporary
      }
      if (verbosity==1) //
        cout << "Expression before type analysis: " << *parse_tree << endl;
    }
    @< Analyse types and then evaluate and print, or catch runtime or other
       errors @>
    output_stream= &cout; // reset output stream if it was changed
  }
  clear_history();
  // clean up (presumably disposes of the lines stored in history)
  cout << "Bye.\n";
  return 0;
}

@ Here are several calls necessary to get various parts of this program off to
a good start, starting with the history and readline libraries, and setting a
comment convention. Initialising the constants in the Atlas library is no
longer necessary, as it is done automatically before |main| is called. Our own
compilation units do require explicit calling of their initialisation
functions.

@h "built-in-types.h"
@h "constants.h"
@< Initialise various parts of the program @>=
@< Initialise the \.{readline} library interface @>

@)ana.set_comment_delims('{','}');
@)initialise_evaluator(); initialise_builtin_types();

@ The function |id_completion_func| define in the \.{lexer} module will not be
plugged directly into the readline completion mechanism, but instead we
provide an alternative function for generating matches, which may pass the
above function to |rl_completion_matches| when it deems the situation
appropriate, or else returns |nullptr| to indicate that the default function,
completing on file names, should be used instead.

@< Definitions of local functions @>=
#ifndef NREADLINE
extern "C" char** do_completion(const char* text, int start, int end)
{
  if (start>0)
  { int i; char c; bool need_file=false;
    for (i=0; i<start; ++i)
      if (std::isspace(c=rl_line_buffer[i]))
        continue; // ignore space characters
      else if (c=='<' or c=='>')
        need_file=true;
      else
        break;

    if (need_file and i==start)
       // the text is preceded by one or more copies of \.<, \.>
      return nullptr; // signal that file name completion should be used
  }
  rl_attempted_completion_over = true;
    // don't try file name completion if we get here
  return rl_completion_matches(text,atlas::interpreter::id_completion_func);
}
#endif

@ The code concerning the \.{readline} library is excluded in cas
the \.{NREADLINE} flag is set.

@< Initialise the \.{readline} library interface @>=
#ifndef NREADLINE
  using_history();
  rl_completer_word_break_characters = lexical_break_chars;
  rl_attempted_completion_function = do_completion; // set up input completion

#endif

@ If a type error is detected by |analyse_types|, then it will have signalled
it and thrown a |runtime_error|; if that happens |type_OK| will remain |false|
and the runtime error is silently caught. If the result is an empty tuple, we
suppress printing of the uninteresting value.

@h <stdexcept>
@h "evaluator.h"

@< Analyse types and then evaluate and print... @>=
{ bool type_OK=false;
  try
  { expression_ptr e;
    last_type=analyse_types(*parse_tree,e); // move assignment of a |type_expr|
    type_OK=true;
    if (verbosity>0)
      cout << "Type found: " << last_type << endl @|
	   << "Converted expression: " << *e << endl;
    e->evaluate(expression_base::single_value);
    last_value=pop_value();
    static type_expr empty(empty_tuple());
    if (last_type!=empty)
      *output_stream << "Value: " << *last_value << endl;
    destroy_expr(parse_tree);
  }
  catch (runtime_error& err)
  { if (type_OK)
      cerr << "Runtime error:\n  " << err.what() << "\nEvaluation aborted.";
    else cerr << err.what();
    cerr << std::endl;
    reset_evaluator(); main_input_buffer->close_includes();
  }
  catch (logic_error& err)
  { cerr << "Internal error: " << err.what() << ", evaluation aborted.\n";
    reset_evaluator(); main_input_buffer->close_includes();
  }
  catch (exception& err)
  { cerr << err.what() << ", evaluation aborted.\n";
    reset_evaluator(); main_input_buffer->close_includes();
  }
}

@ For the moment the only command line argument accepted is \.{-nr}, which
indicates to not use the readline and history library in the input buffer.

@h <cstring>
@< Handle command line arguments @>=
bool use_readline = argc<2 or std::strcmp(argv[1],"-nr")!=0;


@ The |std::ofstream| object was already created earlier in the main loop,
but it will only be opened if we come here. If this fails then we report it
directly and |continue| to the next iteration of the main loop, which is more
practical at this point than throwing and catching an error.

@< Open |redirect| to specified file... @>=
{ redirect.open(ana.scanned_file_name() ,ios_base::out |
     (verbosity==2 ? ios_base::trunc : ios_base::@;app));
  if (redirect.is_open())
    output_stream = &redirect;
  else
  {@; cerr << "Failed to open " << ana.scanned_file_name() << endl;
    continue;
  }
}

@* Index.

% Local IspellDict: british
