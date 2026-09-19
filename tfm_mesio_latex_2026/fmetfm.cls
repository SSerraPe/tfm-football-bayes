% ============================================================
% fmetfm.cls  (English thesis/TFM class, UPC/FME/MESIO inspired)
% Collision-proof with AMS classes
% Version: 2026-09 -- migrated to adopt the MESIO templateTFE2026
% (MESIOTFM.cls) page geometry, sans-serif fonts, colour scheme,
% coloured hyperref links, and cover-page field set (keywords, MSC,
% month/year, acknowledgements), while keeping the chapter-based
% amsbook structure the thesis actually needs (MESIOTFM.cls itself
% is a flat, single-file `article`-based template with no \chapter
% support -- not reused verbatim; see docs/... migration notes).
% ============================================================
\NeedsTeXFormat{LaTeX2e}[1995/06/01]
\ProvidesClass{fmetfm}[2026/09/18 v3.0 English thesis class (UPC/FME/MESIO inspired)]

% ------------------------------------------------------------
% Base class
% ------------------------------------------------------------
\LoadClass[11pt,twoside]{amsbook}

% ------------------------------------------------------------
% Core packages (safe + minimal)
% ------------------------------------------------------------
\RequirePackage{graphicx}
\RequirePackage{amsmath,amssymb,amsthm}
\RequirePackage{xcolor}
\RequirePackage{hyperref}

% ------------------------------------------------------------
% MESIO template: page geometry (adopted verbatim from MESIOTFM.cls)
% ------------------------------------------------------------
\RequirePackage[a4paper,left=20mm,top=18mm,right=20mm,bottom=18mm,%
  includeheadfoot,headheight=24mm,headsep=0mm,footskip=15mm]{geometry}

% ------------------------------------------------------------
% MESIO template: sans-serif document + math fonts
% ------------------------------------------------------------
\RequirePackage{sfmath}
\renewcommand*{\familydefault}{\sfdefault}

% ------------------------------------------------------------
% MESIO template: colour scheme + coloured hyperlinks
% ------------------------------------------------------------
\definecolor{blau}{RGB}{43,129,173}
\definecolor{gris}{RGB}{100,100,100}
\hypersetup{
  colorlinks=true,
  linkcolor=blau,
  citecolor=blau,
  urlcolor=blau
}

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

% New fields introduced by the MESIO templateTFE2026 that this class
% did not previously have: month/year (shown on the cover page) and
% acknowledgements (its own page, right after the cover page).
\newcommand{\monthYear}[1]{\gdef\@fme@monthyear{#1}}
\gdef\@fme@monthyear{Month, year}

\newcommand{\agraiments}[1]{\gdef\@fme@ack{#1}}
\gdef\@fme@ack{}

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
  {\huge\bfseries\color{blau} \@title\par}%
  \vspace{1.5em}%
  {\Large \@fme@author\par}%
  \vspace{6em}%
  {\large \advisorname\ \@fme@advisor\par}%
  \vspace{1.5em}%
  {\large \@fme@monthyear\par}%
  \vfill
  {\large \@fme@dept\par}%
  \endgroup
  \cleardoublepage
  \ifx\@fme@ack\@empty\else
    \thispagestyle{empty}%
    \null\vspace*{2cm}%
    \begin{center}\Large\bfseries\color{blau} Acknowledgements\end{center}%
    \vspace{1cm}%
    \@fme@ack
    \cleardoublepage
  \fi
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
% NOTE: amsbook.cls already defines \keywords (for its own \maketitle-driven
% front matter, which we do not use here -- fmetfm.cls installs its own
% \maketitle above). \providecommand would therefore silently keep amsbook's
% version and swallow the keywords line, so this one deliberately overrides
% it with \renewcommand instead.
\renewcommand{\keywords}[1]{\par\medskip\noindent\keywordname\enspace #1\par}
\providecommand{\msc}[1]{\par\medskip\noindent\mscname\enspace #1\par}

% Legacy-friendly aliases
\providecommand{\keywordseng}[1]{\keywords{#1}}
% MESIO template's field name for keywords (\paraulesclau) -- alias so
% either name works.
\providecommand{\paraulesclau}[1]{\keywords{#1}}

\newenvironment{abstractpage}[1]{%
  \cleardoublepage
  \thispagestyle{empty}%
  \null
  \vspace*{1.5cm}%
  \begin{center}\Large\bfseries\color{blau} #1\end{center}%
  \vspace{1cm}%
}{%
  \cleardoublepage
}

\renewenvironment{abstract}{\begin{abstractpage}{Abstract}}{\end{abstractpage}}
\newenvironment{abstracteng}{\begin{abstractpage}{Abstract}}{\end{abstractpage}}

% ------------------------------------------------------------
% Theorem environments (names/styles match MESIOTFM.cls exactly;
% counter shares [chapter] rather than the template's [section],
% since chapters -- not sections -- are this thesis's top-level
% division, matching the existing "Remark 4.13"-style numbering
% used throughout the document).
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
\newtheorem{observation}[theorem]{Observation}

% ------------------------------------------------------------
% Page style
% ------------------------------------------------------------
\pagestyle{headings}

\endinput
