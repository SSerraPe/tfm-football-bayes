% ============================================================
% fmetfm.cls  (English thesis/TFM class, UPC/FME inspired)
% Collision-proof with AMS classes
% Version: 2026-03-03
% ============================================================
\NeedsTeXFormat{LaTeX2e}[1995/06/01]
\ProvidesClass{fmetfm}[2026/03/03 v2.2 English thesis class (UPC/FME inspired)]

% ------------------------------------------------------------
% Base class
% ------------------------------------------------------------
\LoadClass[11pt,oneside,onecolumn]{amsbook}

% ------------------------------------------------------------
% Core packages (safe + minimal)
% ------------------------------------------------------------
\RequirePackage{graphicx}
\RequirePackage{amsmath,amssymb,amsthm}
\RequirePackage{hyperref}

% NOTE:
% - Do NOT load makeidx here (AMS classes define index helpers and may
%   load/expect index macros; makeidx can clash).
% - Let the user decide in main.tex whether they want an index.
% ------------------------------------------------------------

% ------------------------------------------------------------
% Metadata fields
% ------------------------------------------------------------
\newcommand{\autor}[1]{\gdef\@fme@author{#1}}
\newcommand{\treball}[1]{\gdef\@fme@worktype{#1}}
\newcommand{\advisor}[1]{\gdef\@fme@advisor{#1}}
\newcommand{\departament}[1]{\gdef\@fme@dept{#1}}

\gdef\@fme@author{Author name}
\gdef\@fme@worktype{Type of work}
\gdef\@fme@advisor{Advisor name}
\gdef\@fme@dept{Department name}

\providecommand{\advisorname}{Advisor:}
\providecommand{\keywordname}{\textbf{Key words:}}
\providecommand{\mscname}{\textbf{MSC (2020):}}

% ------------------------------------------------------------
% Title page
% ------------------------------------------------------------
\renewcommand{\maketitle}{%
  \cleardoublepage
  \thispagestyle{empty}%
  \begingroup
  \centering
  {\Large \sffamily Universitat Polit\`ecnica de Catalunya\par}%
  \vspace{0.35em}%
  {\Large \sffamily Facultat de Matem\`atiques i Estad\'istica\par}%
  \vspace{4.5em}%
  {\Large \@fme@worktype\par}%
  \vspace{4.5em}%
  {\huge\bfseries \@title\par}%
  \vspace{1.5em}%
  {\Large \@fme@author\par}%
  \vspace{6em}%
  {\large \advisorname\ \@fme@advisor\par}%
  \vfill
  {\large \@fme@dept\par}%
  \endgroup
  \cleardoublepage
}

% ------------------------------------------------------------
% Optional: Dedication page
% ------------------------------------------------------------
\newcommand{\dedicatoria}[1]{%
  \cleardoublepage
  \thispagestyle{empty}%
  \null\vspace*{2cm}%
  \hfill\begin{minipage}{0.65\textwidth}\raggedleft #1\end{minipage}%
  \cleardoublepage
}

% ------------------------------------------------------------
% Abstract pages + keyword/MSC helpers
% Use \providecommand to avoid redefinition clashes.
% ------------------------------------------------------------
\providecommand{\keywords}[1]{\par\medskip\noindent\keywordname\enspace #1\par}
\providecommand{\msc}[1]{\par\medskip\noindent\mscname\enspace #1\par}

% Legacy-friendly aliases
\providecommand{\keywordseng}[1]{\keywords{#1}}

\newenvironment{abstractpage}[1]{%
  \cleardoublepage
  \thispagestyle{empty}%
  \null
  \vspace*{1.5cm}%
  \begin{center}\Large\bfseries #1\end{center}%
  \vspace{1cm}%
}{%
  \cleardoublepage
}

\renewenvironment{abstract}{\begin{abstractpage}{Abstract}}{\end{abstractpage}}
\newenvironment{abstracteng}{\begin{abstractpage}{Abstract}}{\end{abstractpage}}

% ------------------------------------------------------------
% Optional pages: Preface + Notation
% ------------------------------------------------------------
\newenvironment{preface}{%
  \cleardoublepage
  \thispagestyle{empty}%
  \null\vspace*{1.5cm}%
  \begin{center}\Large\bfseries Preface\end{center}%
  \vspace{1cm}%
}{%
  \cleardoublepage
}

\newenvironment{notation}{%
  \cleardoublepage
  \thispagestyle{empty}%
  \null\vspace*{1.5cm}%
  \begin{center}\Large\bfseries Notation\end{center}%
  \vspace{1cm}%
}{%
  \cleardoublepage
}

% ------------------------------------------------------------
% Theorem environments
% ------------------------------------------------------------
\theoremstyle{plain}
\newtheorem{theorem}{Theorem}[chapter]
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{conjecture}[theorem]{Conjecture}

\theoremstyle{definition}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{example}[theorem]{Example}
\newtheorem{exercise}[theorem]{Exercise}

\theoremstyle{remark}
\newtheorem{remark}[theorem]{Remark}
\newtheorem{note}[theorem]{Note}

% ------------------------------------------------------------
% Page style
% ------------------------------------------------------------
\pagestyle{headings}

\endinput
